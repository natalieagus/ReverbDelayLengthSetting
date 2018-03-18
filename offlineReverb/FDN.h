//
//  FDN.h
//
//   A Feedback Delay Network Reverberator
//
//  offlineReverb
//
//  Created by Hans on 19/5/15.
//  Copyright (c) 2015 Hans. All rights reserved.
//

#ifndef __offlineReverb__FDN__
#define __offlineReverb__FDN__

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <Accelerate/Accelerate.h>

#endif /* defined(__offlineReverb__FDN__) */


#define NUMDELAYS 16
#define NUMDELAYSEXT 1024
#define AUDIOCHANNELS 1

enum DelayTimeAlgorithm {velvetNoise, velvetPrime1, velvetPrime2, primePrime, randomBasic, randomPrime, roomDimension};

class FDN
{
public:
    // constructor
    FDN();
    FDN(int type, DelayTimeAlgorithm algo);
    void setDecayTime(double rt60);
    ~FDN();
    
protected:
    float inputAttenuation;
    float matrixAttenuation;
    float wetPct,dryPct;
    float un_matrixAttenuation;
    float tapGains[NUMDELAYSEXT];
    float tapSigns[NUMDELAYSEXT];
    float temp1[NUMDELAYSEXT];
    float temp2[NUMDELAYSEXT];
    float temp3[NUMDELAYSEXT];
    float temp16[16];
    float ppxx_xn, xxpp, pnxx_xn, xxpn;
    
    float* delayBuffers;
    
    // read / write indices
    int rwIndices[NUMDELAYSEXT];
    int endIndices[NUMDELAYSEXT];
    int startIndices[NUMDELAYSEXT];
    int totalDelayTime;
    
    // delay times
    int delayTimes[NUMDELAYSEXT];
    
    // The Feedback matrix
    float feedbackMatrix[NUMDELAYS][NUMDELAYS];
    
    int random_val, randomSeed;
    float maxRand_f;
    
    float outputs[NUMDELAYSEXT];
    
    void updateRand();
    void resetReadIndices();
    void incrementIndices();
    void randomPermutation(int* a, int length, int audioChannels);
    void randomPermutation1Channel(int* a, int length, int audioChannels);
    void hadamardTransform(float* input, float* output, size_t length);
    void hadamard16(float* input, float* output);
    
    void setDelayTime(int delayType);
    
    // Delay time setting methods
    void setDelayTimesVelvetNoise();
    void setDelayTimesVelvetPrime1();
    void setDelayTimesVelvetPrime2();
    void setDelayTimesRandom();
    void setDelayTimesRandomPrime();
    void setDelayTimesRoomDimension();
    
    
    int rvType, numDelays;
    
public:
    void resetDelay(float decayTime);
    void processAudio(float* pInput, float* pOutputL, float* pOutputR);
    void impulseResponse_write(long numSamples, std::ofstream* outputFile);
    void impulseResponse_output(long numSamples, float* output);
};

