//
//  Academic License - for use in teaching, academic research, and meeting
//  course requirements at degree granting institutions only.  Not for
//  government, commercial, or other organizational use.
//
//  FreqConvolute.h
//
//  Code generation for function 'FreqConvolute'
//


#ifndef FREQCONVOLUTE_H
#define FREQCONVOLUTE_H

// Include files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
extern void FreqConvolute(coder::array<double, 1U> &input, const coder::array<
  double, 1U> &ir_frame_real, const coder::array<double, 1U> &ir_frame_imag,
  double fft_frame_size, coder::array<double, 1U> &out, double *out_size);
extern void FreqConvolute_initialize();
extern void FreqConvolute_terminate();
extern void GetUnisonPartitionedIRFrames(const coder::array<double, 1U> &IR,
  double frame_size, double buffer_size, coder::array<double, 1U>
  &ir_frames_real, coder::array<double, 1U> &ir_frames_imag, double *n_frames,
  double *output_size);
extern void RemoveTailBelowThreshold(coder::array<double, 1U> &IR, double
  threshold_dB);

#endif

// End of code generation (FreqConvolute.h)
