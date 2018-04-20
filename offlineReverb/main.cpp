// offlineReverb.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include "FDN.h"
#include <ctime>
#include <string>
#include <fstream>
#include <cstdint>
#include "assert.h"
#include <Accelerate/Accelerate.h>
#include "SFM.hpp"
#include "LBQ.hpp"
#include "primes.h"

#define FS 44100

using std::cin;
using std::cout;
using std::endl;
using std::fstream;
using std::string;
using namespace std;

// type         : FDN size
// samples      : length of IR in samples
// delaySetting : delay setting type (see setDelayTime method in FDN.cpp), e.g: 1 is velvet delay line method
// output       : IR output
void impulseResponse(int type, int samples, float* output, DelayTimeAlgorithm algo){
    FDN reverb = FDN(type,algo);
    reverb.impulseResponse_output(samples, output);
    
}

// Writes multiple impulse response .csv files
// times    : the number of IR .csv files generated
void saveImpulseMultiple(int type, int samples, int times){
    srand((int)time(0));

    std::ofstream ofstream;
    std::ofstream* ofL = &ofstream;

    for (int i = 0; i < times; i++){


        //    clock_t begin = clock();
        std::string filenameL = "impulse";
        if (type > 0) filenameL += "H";
        else filenameL += "C";
        filenameL += std::to_string(abs(type));
        filenameL += "_";
        filenameL += std::to_string(abs(i));
        filenameL += "_";
        filenameL += ".csv";
        
        
        ofL->open(filenameL);
        
        FDN reverb = FDN(type,DelayTimeAlgorithm::velvetNoise);
        reverb.impulseResponse_write(samples, ofL);
        
        std::cout << "impulse saved for type " << type << ".\n";
        ofL->close();
    }
}

// Writes one impulse response .csv file
void saveImpulse(int type, int samples){
    
    std::ofstream ofstream;
    std::ofstream* ofL = &ofstream;
    
    //    clock_t begin = clock();
    std::string filenameL = "impulse";
    if (type > 0) filenameL += "H";
    else filenameL += "C";
    filenameL += std::to_string(abs(type));
    filenameL += ".csv";
    
    
    ofL->open(filenameL);
    
    
    
    FDN reverb = FDN(type,DelayTimeAlgorithm::velvetNoise);
    reverb.impulseResponse_write(samples, ofL);
    
    std::cout << "impulse saved for type " << type << ".\n";
    ofL->close();

}

/*
 *For array with single dimension, set either outer_loop or inner_loop to 1
 * Write double arrays to .csv file, with line break on each inner_loop
 */
void writeToFile(float* array, int outer_loop, int inner_loop, std::string filename = "SFM_"){
    
    std::ofstream ofstream;
    std::ofstream* of = &ofstream;
    filename += std::to_string(inner_loop*outer_loop);
    filename += ".csv";
    
    of->open(filename);
    
    
    for (int i = 0; i<outer_loop; i++){
        for (int j = 0; j<inner_loop; j++){
            *of << array[i*inner_loop + j] << ",";
        }
        *of << "\n";
    }
    
    std::cout << "array saved for filename: " << filename << "\n";
    
    of->close();
    
}
//
//void writeToFileDoubleArray(float** array, int outer_loop, int inner_loop, std::string filename = "SFM_"){
//
//    std::ofstream ofstream;
//    std::ofstream* of = &ofstream;
//    filename += std::to_string(outer_loop*inner_loop);
//    filename += ".csv";
//
//    of->open(filename);
//
//
//    for (int i = 0; i<outer_loop; i++){
//        for (int j = 0; j<inner_loop; j++){
//            *of << array[i][j] << ",";
//        }
//        *of << "\n";
//    }
//
//    std::cout << "array saved for filename: " << filename << "\n";
//
//    of->close();
//
//}

int main(int argc, char* argv[])
{
    //number of repetition to generate IR
    int iteration = 3;
    int FDN_Size = 16;
    
    // Take 3 seconds of lossless FDN output:  44100 * 3 = 132300
    // Length of signal must be power of 2 for FFT in SFM.cpp
    // Nearest power of 2 to 132300 is 2^17 = 131072
    // So set impulseLength (in samples) to 131072

    int impulseLength = 131072; // changed to smaller value for testing
    int windowLength = 1024;

    int lags = 2;

    bool powerOfTwo = !(impulseLength==0) && !(impulseLength & (impulseLength-1));
    assert(powerOfTwo == true);

    float* output = (float*) malloc(impulseLength*sizeof(float));
    SFM sfm = SFM(impulseLength);
    LBQ lbq_test = LBQ(impulseLength, lags);

    float* SFM_output = (float*) malloc(iteration*sizeof(float));
    float* LBQ_output = (float*) malloc(iteration*sizeof(float));
    float* mean_SFM = (float*) malloc(iteration*sizeof(float));
    float* stdev_SFM = (float*) malloc(iteration*sizeof(float));
    
    //separated into 2 values per iteration: SFM for early part and SFM for late part
    float* SFM_early_late_output = (float*) malloc(iteration*2*sizeof(float));
 
    float* SFM_output_window_array = (float*) malloc(iteration * impulseLength/windowLength*sizeof(float));
//    for (int i=0; i<iteration; i++)
//        SFM_output_window_array[i] = (float *)malloc(impulseLength/windowLength * sizeof(float));

    for (int i = 0 ; i<iteration ; i++){

        memset(output, 0, impulseLength*sizeof(float));
        impulseResponse(FDN_Size, impulseLength, output,DelayTimeAlgorithm::velvetNoise);
        
//        std::cout << "Signal: " ;
//        printf("{");
//        for (int i = 0; i<impulseLength-1; i++) printf("%f ,", output[i]);
//        printf("%f}", output[impulseLength-1]);

        SFM_output[i] = sfm.spectral_flatness_value(output);
        LBQ_output[i] = lbq_test.LBQtest(output);
        
        sfm.spectral_flatness_value_array(output, SFM_output_window_array+(impulseLength/windowLength * i), windowLength);
        
        // 2205 is 50 ms, the nearest power of two to this is 2048
        sfm.spectral_flatness_value_early_late(output, SFM_early_late_output+(2*i), SFM_early_late_output+(2*i + 1), windowLength);
        
//        printf("\n");
//        for (int k=0; k<impulseLength/windowLength; k++) printf("%f , ", *(SFM_output_window_array+(impulseLength/windowLength * i + k)));
//        printf("\n");
        
        vDSP_normalize(SFM_output_window_array+(impulseLength/windowLength * i), 1, NULL, 1, &mean_SFM[i], &stdev_SFM[i], impulseLength/windowLength);
    
        
        printf("Iteration %i : ", i );
        printf("SFM : %f LBQ : %f mean_SFM: %f stdev_SFM: %f \n \n", SFM_output[i], LBQ_output[i], mean_SFM[i], stdev_SFM[i]);
    }
    
    writeToFile(SFM_output, iteration, 1, "SFM_all");
    writeToFile(SFM_output_window_array, iteration, impulseLength/windowLength, "SFM_array");
    writeToFile(SFM_early_late_output, iteration, 2, "SFM_early_late");
    
}



