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
    
    // x        : time samples
    // returns  : SFM value
    float spectral_flatness_value(float* x);
    
    
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
    vDSP_DFT_Setup setup;
    
    // Compute geometric mean
    // data     : real value array
    float geometric_mean(float* data, size_t N);
    
    // Convert real input sequence x into complex form, stored in x_ptr_complex
    // Imaginary part is always zero
    // e.g: x = {1,2,3,4}, x_ptr_complex = {{1,0}, {2,0}, {3,0}, {4,0}}
    void convert_real_to_complex(float* x);
    
    // perform FFT on input in DSP complex form : (DSPComplex*) x_ptr_complex
    void fft(DSPComplex input[]);
    
    void reset();
    
};


#endif /* SFM_hpp */
