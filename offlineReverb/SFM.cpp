//
//  SFM.cpp
//  offlineReverb
//
//  Created by Natalie Agus on 13/3/18.
//  Copyright © 2018 Hans. All rights reserved.
//

#include "SFM.hpp"



SFM::SFM(size_t N){

    printf("\nSFM class instantiated \n\n");
    
    // sequence length has to be power of two
    bool powerOfTwo = !(N==0) && !(N & (N-1));
    assert(powerOfTwo == true);
    
    this->sequence_length = N;
    
    // initialize the setup struct
//    this->setup = vDSP_DFT_zop_CreateSetup(NULL, this->sequence_length, vDSP_DFT_FORWARD);
    this->setup_fft = vDSP_create_fftsetup(log2(this->sequence_length), kFFTRadix2);

    
    // create complex array with size of 2N (for Re and Im each)
    this->x_ptr_complex = (float*) malloc(N*sizeof(float)*2);
    // create power spectra array
    this->power_spectra = (float*) malloc(N*sizeof(float));
    // create input copy array
    this->input = (float*) malloc(N*sizeof(float));
    
    
    // create input and output complex array for fft with size 2N (for Re and Im each)
    this->inputMemory = (float*) malloc(N*sizeof(float)*2);
    this->outputMemory = (float*) malloc(N*sizeof(float)*2);
    
    // link the input and output complex array for fft to object DSPSplitComplex for vDSP operations
    this->inputSplit.realp = inputMemory;
    this->inputSplit.imagp = inputMemory+N;
    this->outputSplit.realp = outputMemory;
    this->outputSplit.imagp = outputMemory+N;

    
    reset();

    
}

void SFM::reset(){
    memset(x_ptr_complex, 0, this->sequence_length*sizeof(float)*2);
    memset(power_spectra, 0, this->sequence_length*sizeof(float));
    memset(inputMemory, 0, this->sequence_length*sizeof(float)*2);
    memset(outputMemory, 0, this->sequence_length*sizeof(float)*2);
}

void SFM::convert_real_to_complex(float* x){
    memcpy(input, x, sequence_length*sizeof(float));
    
//    printf("\n sequence Length : %zu", sequence_length);
//
//    printf("\n Signal: \n");
//    for (int i = 0; i<sequence_length; i++){
//        printf(" %f " , x[i]);
//    }
//
//    printf("\n Signal copied: \n");
//    for (int i = 0; i<sequence_length; i++){
//        printf(" %f " , input[i]);
//    }
//
    
    vDSP_vswap(input, 1, x_ptr_complex, 2, sequence_length);
    
//    printf("\n Signal: \n");
//    for (int i = 0; i<sequence_length*2; i++){
//        printf("{Re: %f, Im : %f}, \n ", x_ptr_complex[i], x_ptr_complex[i+1]);
//        i ++;
//    }
//
    
    
}


inline void SFM::fft(float* input, size_t* outputLength, size_t inputLength) {

//    for (int i = 0; i<sequence_length; i++) printf("%f ,", input[i]);

    
    vDSP_ctoz((DSPComplex*) input, 2, &inputSplit, 1, inputLength/2);
    
//    for (int i = 0; i<sequence_length; i++){
//        printf("{ %f , %f } , ", inputSplit.realp[i], inputSplit.imagp[i]);
//    }
    
    //
    // this seems to do the same thing as the line above it. I rewrote it
    // because it seems likely to be faster.
    //memcpy(inputSplit.realp, input, sizeof(float)*this->sequence_length);
    //memset(inputSplit.imagp, 0,     sizeof(float)*this->sequence_length);
    
    //using DFT
//    vDSP_DFT_Execute(this->setup,
//                     inputSplit.realp, inputSplit.imagp,
//                     outputSplit.realp, outputSplit.imagp);
//
    
    // using out of place fft
    // vDSP_fft_zop(setup_fft, &inputSplit, 1, &outputSplit, 1, log2(this->sequence_length), FFT_FORWARD);
    // *outputLength = (sequence_length/2) + 1;
    
    // in place fft for real valued input
    vDSP_fft_zrip(setup_fft, &inputSplit, 1, log2(inputLength), FFT_FORWARD);
    *outputLength = inputLength/2;
}

// Geometric mean computation using split exponent and mantissa
float SFM::geometric_mean(float* data, size_t N)
{
    long long ex = 0;
    auto do_bucket = [data,&ex](size_t first,size_t last) -> float
    {
        float ans = 1.0;
        for ( ;first != last;++first)
        {
            int i;
//            printf("%f ", data[first]);
            ans *= std::frexp(data[first],&i);
            ex+=i;
        }
//        printf("ans : %f ", ans);
        return ans;
    };
    
    const size_t bucket_size = -std::log2( std::numeric_limits<float>::min() );
    std::size_t buckets = N / bucket_size;
    
    float invN = 1.0 / N;
    float m = 1.0;
    
    for (std::size_t i = 0;i < buckets;++i)
        m *= std::pow( do_bucket(i * bucket_size,(i+1) * bucket_size),invN );
    
    m*= std::pow( do_bucket( buckets * bucket_size, N ),invN );
    
    return std::pow( std::numeric_limits<float>::radix,ex * invN ) * m;
}

// correct the DC and Nyquist terms to follow the same probability
// density function as the complex terms, for the vDSP_fft_zrip
// packing that packs DC and Nyquist together in the zero position
// of the array
//
// see our paper titled
// The Probability Density Function of the Spectral Flatness Measure
// for an explanation.
void correctDCNyquist_zripPacked(DSPSplitComplex* fftZripOutput){
    // divide the real and imaginary parts of the zero element by sqrt(2)
    fftZripOutput->realp[0] *= M_SQRT1_2;
    fftZripOutput->imagp[0] *= M_SQRT1_2;
}


float SFM::compute_spectral_flatness_value(float* x, size_t length){
    
    
////         print the input signal
//        printf("\n Input Signal: \n");
//        for (int i = 0; i<length; i++){
//            printf("%f, ", *(x + i));
//        }
    
    //convert_real_to_complex(x);
    size_t outputLength;
    
    fft(x, &outputLength, length);
    
    // create a pointer as an alias to the output of the fft (for readability)
    DSPSplitComplex* fftResult = &inputSplit;
    
//     print the FFT output
//    printf("\n FFT result: \n");
//    for (int i = 0; i<outputLength; i++){
//        printf("{Re: %f, Im : %f}, \n ", fftResult->realp[i]/2, fftResult->imagp[i]/2);
//    }
    
    correctDCNyquist_zripPacked(fftResult);
    
    // print the FFT output after correcting the DC / Nyquist PDF
//    printf("\n FFT result: \n");
//    for (int i = 0; i<outputLength; i++){
//        printf("{Re: %f, Im : %f}, \n ", fftResult->realp[i], fftResult->imagp[i]);
//    }
    
//    printf("\n FFT result: \n");
//    for (int i = 0; i<outputLength; i++){
//        printf("{Re: %f, Im : %f}, \n ", fftResult->realp[i], fftResult->imagp[i]);
//    }
//
    assert(sequence_length%2 == 0);
    
    vDSP_zvabs(fftResult, 1, power_spectra, 1, outputLength);
    
//    printf("\n Magnitude spectra: ");
//    for (int i = 0; i<outputLength; i++){
//        printf("%f,  ", power_spectra[i]);
//    }

    vDSP_vsq(power_spectra, 1, power_spectra, 1, outputLength);
    
//    printf("\n Power spectra: ");
//    for (int i = 0; i<outputLength; i++){
//        printf("%f,  ", power_spectra[i]);
//    }
//
//    printf("\n");
    
    float SFM_numerator = 0.0f;
    float SFM_denominator = 0.0f;
    
    SFM_numerator = geometric_mean(power_spectra, outputLength);
    vDSP_meanv(power_spectra, 1, &SFM_denominator, outputLength);

//    printf("\n Geometric mean : %f ", SFM_numerator);
//    printf("\n Arithmetic mean: %f \n", SFM_denominator);
    
    reset();
    
    return SFM_numerator/SFM_denominator;
}


float SFM::spectral_flatness_value(float* x){
    return compute_spectral_flatness_value(x, this->sequence_length);
}


void SFM::spectral_flatness_value_array(float *x, float *SFM_array, int n){
    
    // n has to be power of two
    bool powerOfTwo = !(n==0) && !(n & (n-1));
    assert(powerOfTwo == true);
    
    for (int i = 0; i<(this->sequence_length/n); i++){
        SFM_array[i] = compute_spectral_flatness_value((x+(i*n)), n);
    }
    
}
