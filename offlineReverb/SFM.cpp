//
//  SFM.cpp
//  offlineReverb
//
//  Created by Natalie Agus on 13/3/18.
//  Copyright © 2018 Hans. All rights reserved.
//

#include "SFM.hpp"

SFM::SFM(size_t N){

    printf("\n \n SFM class instantiated\n");
    
    bool powerOfTwo = !(N==0) && !(N & (N-1));
    assert(powerOfTwo == true);
    
    this->sequence_length = N;
    
    // initialize the setup struct
    this->setup = vDSP_DFT_zop_CreateSetup(NULL, this->sequence_length, vDSP_DFT_FORWARD);
    
    // create complex array with size of 2N (for Re and Im each)
    this->x_ptr_complex = (float*) malloc(N*sizeof(float)*2);
    // create power spectra array
    this->power_spectra = (float*) malloc(N*sizeof(float));
    
    // create input and output complex array for fft with size 2N (for Re and Im each)
    this->inputMemory = (float*) malloc(N*sizeof(float)*2);
    this->outputMemory = (float*) malloc(N*sizeof(float)*2);
    
    // link the input and output complex array for fft to object DSPSplitComplex for vDSP operations
    this->inputSplit.realp = inputMemory;
    this->inputSplit.imagp = inputMemory+N;
    this->outputSplit.realp = outputMemory;
    this->outputSplit.imagp = outputMemory+N;

    
    clear_all();

    
}

void SFM::clear_all(){
    memset(x_ptr_complex, 0, this->sequence_length*sizeof(float)*2);
    memset(power_spectra, 0, this->sequence_length*sizeof(float));
    memset(inputMemory, 0, this->sequence_length*sizeof(float)*2);
    memset(outputMemory, 0, this->sequence_length*sizeof(float)*2);
}

void SFM::convert_real_to_complex(float* x){
    vDSP_vswap(x, 1, x_ptr_complex, 2, sequence_length);
}


inline void SFM::fft(DSPComplex input[]) {

    vDSP_ctoz(input, 2, &inputSplit, 1, this->sequence_length);
    
    //vDSP_DFT_Setup setup = vDSP_DFT_zop_CreateSetup(NULL, this->sequence_length, vDSP_DFT_FORWARD);
    
    vDSP_DFT_Execute(this->setup,
                     inputSplit.realp, inputSplit.imagp,
                     outputSplit.realp, outputSplit.imagp);
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


float SFM::spectral_flatness_value(float* x){
    convert_real_to_complex(x);
    fft((DSPComplex*)x_ptr_complex);
    
//    printf("\n FFT result: \n");
//    for (int i = 0; i<sequence_length; i++){
//        printf("{Re: %f, Im : %f}, \n ", outputSplit.realp[i], outputSplit.imagp[i]);
//    }
    
    assert(sequence_length%2 == 0);
    
    vDSP_zvabs(&this->outputSplit, 1, power_spectra, 1, sequence_length/2+1);
    
//    printf("\n Magnitude spectra: ");
//    for (int i = 0; i<sequence_length/2+1; i++){
//        printf("%f,  ", power_spectra[i]);
//    }

    vDSP_vsq(power_spectra, 1, power_spectra, 1, sequence_length/2+1);
    
//    printf("\n Power spectra: ");
//    for (int i = 0; i<sequence_length/2+1; i++){
//        printf("%f,  ", power_spectra[i]);
//    }
//
//    printf("\n");
    
    float SFM_numerator = 0.0f;
    float SFM_denominator = 0.0f;
    
    SFM_numerator = geometric_mean(power_spectra, sequence_length/2+1);
    vDSP_meanv(power_spectra, 1, &SFM_denominator, sequence_length/2+1);

    printf("\n Geometric mean : %f ", SFM_numerator);
    printf("\n Arithmetic mean: %f \n", SFM_denominator);
    
    clear_all();
    
    return SFM_numerator/SFM_denominator;
    
}
