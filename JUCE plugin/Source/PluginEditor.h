/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
#include "PluginProcessor.h"

//==============================================================================
/**
*/
class ConvolutionReverbAudioProcessorEditor  : public juce::AudioProcessorEditor, public Button::Listener
{
public:
    ConvolutionReverbAudioProcessorEditor (ConvolutionReverbAudioProcessor&, AudioProcessorValueTreeState&);
    ~ConvolutionReverbAudioProcessorEditor() override;

    //==============================================================================
    void paint (juce::Graphics&) override;
    void resized() override;

private:
    AudioProcessorValueTreeState& valueTreeState;

    AudioFormatManager formatManager;
    std::unique_ptr<AudioFormatReaderSource> readerSource;

    int labelHeight = 20;
    int sliderHeight = 30;

    Label mixLabel;
    Slider mixSlider;
    std::unique_ptr<AudioProcessorValueTreeState::SliderAttachment> mixAttachment;

    Label maxFramesLabel;
    Slider maxFramesSlider;
    std::unique_ptr<AudioProcessorValueTreeState::SliderAttachment> maxFramesAttachment;

    Label IRNameLabel; 
    TextButton loadIRButton;

    // This reference is provided as a quick way for your editor to
    // access the processor object that created it.
    ConvolutionReverbAudioProcessor& audioProcessor;

    AudioBuffer<float> IRFileAudioBuffer;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (ConvolutionReverbAudioProcessorEditor)

        // Inherited via Listener
        virtual void buttonClicked(Button*) override;
};
