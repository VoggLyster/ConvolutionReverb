/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
ConvolutionReverbAudioProcessorEditor::ConvolutionReverbAudioProcessorEditor (ConvolutionReverbAudioProcessor& p, AudioProcessorValueTreeState& vts)
    : AudioProcessorEditor (&p), audioProcessor (p), valueTreeState(vts)
{
    // Make sure that before the constructor has finished, you've set the
    // editor's size to whatever you need it to be.
    mixLabel.setText("Mix", dontSendNotification);
    addAndMakeVisible(mixLabel);
    addAndMakeVisible(mixSlider);
    mixAttachment.reset(new AudioProcessorValueTreeState::SliderAttachment(valueTreeState, "mix", mixSlider));
    maxFramesLabel.setText("Max frames", dontSendNotification);
    addAndMakeVisible(maxFramesLabel);
    addAndMakeVisible(maxFramesSlider);
    maxFramesAttachment.reset(new AudioProcessorValueTreeState::SliderAttachment(valueTreeState, "maxFrames", maxFramesSlider));
    
    IRNameLabel.setText("No IR loaded", dontSendNotification);
    addAndMakeVisible(IRNameLabel);
    addAndMakeVisible(loadIRButton);
    loadIRButton.setButtonText("Load IR");
    loadIRButton.addListener(this);

    setSize(400, 300);

    formatManager.registerBasicFormats();
}

ConvolutionReverbAudioProcessorEditor::~ConvolutionReverbAudioProcessorEditor()
{
    
}

//==============================================================================
void ConvolutionReverbAudioProcessorEditor::paint (juce::Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll (getLookAndFeel().findColour (juce::ResizableWindow::backgroundColourId));

    g.setColour (juce::Colours::white);
    g.setFont (15.0f);
    //g.drawFittedText ("Hello World!", getLocalBounds(), juce::Justification::centred, 1);
}

void ConvolutionReverbAudioProcessorEditor::resized()
{
    int border = 20;
    int w = getWidth() - 2 * border;
    mixLabel.setBounds(border, border, w, labelHeight);
    mixSlider.setBounds(border, border + labelHeight, w, sliderHeight);
    maxFramesLabel.setBounds(border, border * 2 + labelHeight + sliderHeight, w, labelHeight);
    maxFramesSlider.setBounds(border, border * 2 + labelHeight * 2 + sliderHeight, w, sliderHeight);
    IRNameLabel.setBounds(border, border * 3 + labelHeight * 2 + sliderHeight * 2, w, labelHeight);
    loadIRButton.setBounds(border, border * 3 + labelHeight * 3 + sliderHeight * 2, w, sliderHeight);
}

void ConvolutionReverbAudioProcessorEditor::buttonClicked(Button*)
{
    FileChooser chooser("Select an impulse response to load",
        {},
        "*wav");
    if (chooser.browseForFileToOpen())
    {
        auto file = chooser.getResult();
        auto* reader = formatManager.createReaderFor(file);
        IRNameLabel.setText(file.getFileName(), dontSendNotification);

        if (reader != nullptr)
        {
            IRFileAudioBuffer = AudioBuffer<float>(reader->numChannels, reader->lengthInSamples);
            int numSamples = reader->lengthInSamples;
            reader->read(&IRFileAudioBuffer, 0, IRFileAudioBuffer.getNumSamples(), 0, true, true);
            audioProcessor.IRFileReady(&IRFileAudioBuffer);
        }
        delete(reader);
    }
}
