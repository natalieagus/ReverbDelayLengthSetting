//
//  SFM.hpp
//  offlineReverb
//
//  Created by Natalie Agus on 13/3/18.
//  Copyright © 2018 Hans. All rights reserved.
//

#ifndef SFM_hpp
#define SFM_hpp
#include <stdio.h>
#include <Accelerate/Accelerate.h>


class SFM
{
public:
    
    //constructor
    SFM(size_t N);
    
    // computes SFM value for entire sample
    // x        : time samples
    // returns  : SFM value
    float spectral_flatness_value(float* x);
    
    // get an array of SFM values from samples with size 'n' each on input x
    void spectral_flatness_value_array(float* x, float* SFM_array, int n);
    
    
private:
    
    float* input;
    // pointer to complex form array {{re, im}, {re, im}, ... }
    float* x_ptr_complex;
    
    float* power_spectra;
    size_t sequence_length;
    
    // Helpers for FFT operation
    float* inputMemory;
    float* outputMemory;
    DSPSplitComplex inputSplit;
    DSPSplitComplex outputSplit;
//    vDSP_DFT_Setup setup;
    FFTSetup setup_fft;
    
    // Compute geometric mean
    // data     : real value array
    float geometric_mean(float* data, size_t N);
    
    // Convert real input sequence x into complex form, stored in x_ptr_complex
    // Imaginary part is always zero
    // e.g: x = {1,2,3,4}, x_ptr_complex = {{1,0}, {2,0}, {3,0}, {4,0}}
    // This is required only if vDSP_DFT is used in fft()
    void convert_real_to_complex(float* x);
    
    // perform FFT on input
    void fft(float* input, size_t* outputLength, size_t inputLength);
    
    // computes SFM value
    float compute_spectral_flatness_value(float* x, size_t samples_length);
    
    // clear buffers
    void reset();
    
};


#endif /* SFM_hpp */
