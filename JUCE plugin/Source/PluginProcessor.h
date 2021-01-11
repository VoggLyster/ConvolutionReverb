/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "FreqConvolute.h"
#include "coder_array.h"

//==============================================================================
/**
*/
class ConvolutionReverbAudioProcessor  : public juce::AudioProcessor
{
public:
    //==============================================================================
    ConvolutionReverbAudioProcessor();
    ~ConvolutionReverbAudioProcessor() override;

    //==============================================================================
    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;

   #ifndef JucePlugin_PreferredChannelConfigurations
    bool isBusesLayoutSupported (const BusesLayout& layouts) const override;
   #endif

    void processBlock (juce::AudioBuffer<float>&, juce::MidiBuffer&) override;

    //==============================================================================
    juce::AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;

    //==============================================================================
    const juce::String getName() const override;

    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool isMidiEffect() const override;
    double getTailLengthSeconds() const override;

    //==============================================================================
    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram (int index) override;
    const juce::String getProgramName (int index) override;
    void changeProgramName (int index, const juce::String& newName) override;

    //==============================================================================
    void getStateInformation (juce::MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;

    //==================
    void IRFileReady(AudioBuffer<float>* IRFileAudioBuffer);
    void PrepareIRFromBuffers();

    // Inherited via Listener
    //void parameterChanged(const String& parameterID, float newValue) override;

private:
    // Exposed params
    AudioProcessorValueTreeState parameters;

    std::atomic<float>* mixParameter = nullptr;
    std::atomic<float>* maxFramesParameter = nullptr;

    // Internal stuff
    int maxBufferSize = 1048576;
    AudioBuffer<float> delayBuffer;
    int delayBufferWriteHead[2] = { 0, 0 };
    int delayBufferReadHead[2] = { 0, 0 };
    AudioBuffer<float> overlapBuffer;
    AudioBuffer<float> internalOverlapBuffer;
    AudioBuffer<float> outputBuffer;
    coder::array<double, 1> IRFramesRealLeft;
    coder::array<double, 1> IRFramesRealRight;
    coder::array<double, 1> IRFramesImagLeft;
    coder::array<double, 1> IRFramesImagRight;
    coder::array<double, 1> IRAudioLeft;
    coder::array<double, 1> IRAudioRight;
    coder::array<double, 1> internalBufferX;
    coder::array<double, 1> internalBufferIRReal;
    coder::array<double, 1> internalBufferIRImag;
    coder::array<double, 1> convOutBuffer;
    double nIRFrames = 0;
    double outputSize = 0;
    int writeHead = 0;
    int fftWindowSize = 1024;
    int blockSize = 0;
    bool IRLoaded = false;

    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (ConvolutionReverbAudioProcessor)
};
