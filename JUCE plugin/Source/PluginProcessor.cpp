/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
NewProjectAudioProcessor::NewProjectAudioProcessor()
#ifndef JucePlugin_PreferredChannelConfigurations
     : AudioProcessor (BusesProperties()
                     #if ! JucePlugin_IsMidiEffect
                      #if ! JucePlugin_IsSynth
                       .withInput  ("Input",  juce::AudioChannelSet::stereo(), true)
                      #endif
                       .withOutput ("Output", juce::AudioChannelSet::stereo(), true)
                     #endif
                       )
#endif
{
}

NewProjectAudioProcessor::~NewProjectAudioProcessor()
{
}

//==============================================================================
const juce::String NewProjectAudioProcessor::getName() const
{
    return JucePlugin_Name;
}

bool NewProjectAudioProcessor::acceptsMidi() const
{
   #if JucePlugin_WantsMidiInput
    return true;
   #else
    return false;
   #endif
}

bool NewProjectAudioProcessor::producesMidi() const
{
   #if JucePlugin_ProducesMidiOutput
    return true;
   #else
    return false;
   #endif
}

bool NewProjectAudioProcessor::isMidiEffect() const
{
   #if JucePlugin_IsMidiEffect
    return true;
   #else
    return false;
   #endif
}

double NewProjectAudioProcessor::getTailLengthSeconds() const
{
    return 0.0;
}

int NewProjectAudioProcessor::getNumPrograms()
{
    return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
                // so this should be at least 1, even if you're not really implementing programs.
}

int NewProjectAudioProcessor::getCurrentProgram()
{
    return 0;
}

void NewProjectAudioProcessor::setCurrentProgram (int index)
{
}

const juce::String NewProjectAudioProcessor::getProgramName(int index)
{
    return {};
}

void NewProjectAudioProcessor::changeProgramName(int index, const juce::String& newName)
{
}

//==============================================================================
void NewProjectAudioProcessor::prepareToPlay(double sampleRate, int samplesPerBlock)
{
    fftWindowSize = nextPowerOfTwo(samplesPerBlock+1);
    delayBuffer = AudioBuffer<float>(2, maxBufferSize);
    delayBuffer.clear();
    overlapBuffer = AudioBuffer<float>(2, samplesPerBlock);
    overlapBuffer.clear();
    internalOverlapBuffer = AudioBuffer<float>(1, samplesPerBlock);
    internalOverlapBuffer.clear();
    outputBuffer = AudioBuffer<float>(2, samplesPerBlock);
    outputBuffer.clear();
    IRAudioLeft.set_size(20000);
    IRAudioRight.set_size(20000);
    internalBufferX.set_size(samplesPerBlock);
    internalBufferIRReal.set_size(fftWindowSize);
    internalBufferIRImag.set_size(fftWindowSize);
    for (int i = 0; i < 20000; ++i) 
    {
        IRAudioLeft[i] = 1.0f;
        IRAudioRight[i] = 1.0f;
    }
    RemoveTailBelowThreshold(IRAudioLeft, -60);
    RemoveTailBelowThreshold(IRAudioRight, -60);
    GetUnisonPartitionedIRFrames(IRAudioLeft, fftWindowSize, samplesPerBlock, IRFramesImagLeft, IRFramesRealLeft, &nIRFrames, &outputSize);
    GetUnisonPartitionedIRFrames(IRAudioRight, fftWindowSize, samplesPerBlock, IRFramesImagRight, IRFramesRealRight, &nIRFrames, &outputSize);
}

void NewProjectAudioProcessor::releaseResources()
{
    // When playback stops, you can use this as an opportunity to free up any
    // spare memory, etc.
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool NewProjectAudioProcessor::isBusesLayoutSupported(const BusesLayout& layouts) const
{
#if JucePlugin_IsMidiEffect
    juce::ignoreUnused(layouts);
    return true;
#else
    // This is the place where you check if the layout is supported.
    // In this template code we only support mono or stereo.
    if (layouts.getMainOutputChannelSet() != juce::AudioChannelSet::mono()
        && layouts.getMainOutputChannelSet() != juce::AudioChannelSet::stereo())
        return false;

    // This checks if the input layout matches the output layout
#if ! JucePlugin_IsSynth
    if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
        return false;
#endif

    return true;
#endif
}
#endif

void NewProjectAudioProcessor::processBlock(juce::AudioBuffer<float>& buffer, juce::MidiBuffer& midiMessages)
{
    juce::ScopedNoDenormals noDenormals;
    auto totalNumInputChannels = getTotalNumInputChannels();
    auto totalNumOutputChannels = getTotalNumOutputChannels();
    int numberOfSamples = buffer.getNumSamples();

    for (auto i = totalNumInputChannels; i < totalNumOutputChannels; ++i)
        buffer.clear(i, 0, buffer.getNumSamples());

    for (int channel = 0; channel < totalNumInputChannels; ++channel)
    {
        auto* channelData = buffer.getWritePointer(channel);
        for (int n = 0; n < numberOfSamples; ++n)
        {
            delayBuffer.addSample(channel, delayBufferWriteHead[channel], channelData[n]);
            delayBufferWriteHead[channel]++;
            if (delayBufferWriteHead[channel] >= maxBufferSize)
                delayBufferWriteHead[channel] = 0;
        }
        if (delayBufferWriteHead[channel] < numberOfSamples)
        {
            delayBufferReadHead[channel] = delayBufferWriteHead[channel] - numberOfSamples + maxBufferSize;
        }
        else
        {
            delayBufferReadHead[channel] = delayBufferWriteHead[channel] - numberOfSamples;
        }
        for (int n = 0; n < numberOfSamples; ++n)
        {
            outputBuffer.setSample(channel, n, overlapBuffer.getSample(channel, n));
        }
        overlapBuffer.clear(channel, 0, numberOfSamples);
        int maxFrameIters = jmin(maxFrames, (int)nIRFrames);
        for (int frame = 0; frame < maxFrameIters; ++frame)
        {
            for (int n = 0; n < numberOfSamples; ++n)
            {
                internalBufferX[n] = delayBuffer.getSample(channel, delayBufferReadHead[channel]);
                delayBufferReadHead[channel]++;
                if (delayBufferReadHead[channel] >= maxBufferSize)
                    delayBufferReadHead[channel] = 0;
            }
            for (int n = 0; n < fftWindowSize; ++n)
            {
                if (channel == 0)
                {
                    internalBufferIRReal[n] = IRFramesRealLeft[frame * fftWindowSize + n];
                    internalBufferIRImag[n] = IRFramesImagLeft[frame * fftWindowSize + n];
                }
                else if (channel == 1)
                {
                    internalBufferIRReal[n] = IRFramesRealLeft[frame * fftWindowSize + n];
                    internalBufferIRImag[n] = IRFramesImagLeft[frame * fftWindowSize + n];
                }
            }
            FreqConvolute(internalBufferX, internalBufferIRReal, internalBufferIRImag, fftWindowSize, convOutBuffer, &outputSize);
            for (int n = 0; n < numberOfSamples; ++n)
            {
                float newOverlapVal = overlapBuffer.getSample(channel, n) + convOutBuffer[numberOfSamples + n]/maxFrameIters;
                overlapBuffer.addSample(channel, n, newOverlapVal);
                float newOutputVal = outputBuffer.getSample(channel, n) + convOutBuffer[n] / maxFrameIters;
                outputBuffer.addSample(channel, n, newOutputVal);
            }
            if (delayBufferReadHead[channel] < numberOfSamples)
                delayBufferReadHead[channel] = maxBufferSize - numberOfSamples;
            else
                delayBufferReadHead[channel] = delayBufferReadHead[channel] - numberOfSamples;
        }

        for (int n = 0; n < numberOfSamples; ++n)
        {
            channelData[n] = (1.0f - mix / 100.0f) * channelData[n] + (mix / 100.0f) * outputBuffer.getSample(channel, n);
        }

    }
}

//==============================================================================
bool NewProjectAudioProcessor::hasEditor() const
{
    return true; // (change this to false if you choose to not supply an editor)
}

juce::AudioProcessorEditor* NewProjectAudioProcessor::createEditor()
{
    return new NewProjectAudioProcessorEditor (*this);
}

//==============================================================================
void NewProjectAudioProcessor::getStateInformation (juce::MemoryBlock& destData)
{
    // You should use this method to store your parameters in the memory block.
    // You could do that either as raw data, or use the XML or ValueTree classes
    // as intermediaries to make it easy to save and load complex data.
}

void NewProjectAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    // You should use this method to restore your parameters from this memory block,
    // whose contents will have been created by the getStateInformation() call.
}

//==============================================================================
// This creates new instances of the plugin..
juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new NewProjectAudioProcessor();
}
