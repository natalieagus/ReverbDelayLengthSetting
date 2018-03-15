//
//  FDN.cpp
//  offlineReverb
//
//  Created by Hans on 19/5/15.
//  Copyright (c) 2015 Hans. All rights reserved.
//

#include "FDN.h"
#include <cstdlib>
#include <cstring>
#include <math.h>
#include <assert.h>
#include <cstdlib>
#include <ctime>
#include "primes.h"
//#define M_E 2.71828182845904523536
#define DENSITY_WINDOW_SIZE 882 // 20.0 * (44100.0 / 1000.0); (20 ms)
#define RV_MIN_DELAY_TIME 40
#define RV_MAX_DELAY_TIME 200

//#include <iostream>


/*
 Helper method to set which delay time setting to use
 Write new methods at the end of this file
 */
void FDN::setDelayTime(int delayType){
    switch(delayType){
        
            //continue adding cases here
        case DelayTimeAlgorithm::velvetNoise:  setDelayTimesVelvetNoise();
        case DelayTimeAlgorithm::velvetPrime1: setDelayTimesVelvetPrime1();
        case DelayTimeAlgorithm::velvetPrime2: setDelayTimesVelvetPrime2();
        case DelayTimeAlgorithm::randomBasic:  setDelayTimesRandom();
        case DelayTimeAlgorithm::randomPrime:  setDelayTimesRandomPrime();
            
        default:setDelayTimesVelvetNoise();
    }
}


FDN::~FDN(){
    if (delayBuffers) delete[] delayBuffers;
}

FDN::FDN(int rvType, DelayTimeAlgorithm algo)
{
    this->rvType = rvType;
    
    numDelays = abs(rvType);
    
    setDelayTime(algo);
    
    
    vDSP_sve((float*) delayTimes, 1, (float*) &totalDelayTime, numDelays);
//    printf("\n Total Delay length %i  ", totalDelayTime);
//    printf("\n No. of FDN numDelays %i \n", numDelays);
    
//    for (int i = 0; i<numDelays; i ++){
//        printf("Delay line %i length : %i \n", i, delayTimes[i]);
//    }
//
    //Settng delay buffers
    // declare and initialise the delay buffers
    delayBuffers = new float[totalDelayTime];
    for (int i = 0; i < totalDelayTime; i++) delayBuffers[i] = 0;
    
    // calculate start and end indices for each delay line
    startIndices[0] = 0;
    endIndices[0] = delayTimes[0];
    for (int i = 1; i < numDelays; i++){
        startIndices[i] = startIndices[i - 1] + delayTimes[i - 1];
        endIndices[i] = startIndices[i] + delayTimes[i];
    }
    
    // since we copy one input channel into many delay lines, we need to attenuate the input to avoid clipping.
    inputAttenuation = 1.0f / sqrt((float)numDelays);
    

    // we specify our feedback matrix with elements of 0, 1, and -1.  In order to keep it unitary,
    // we need to scale it down.
    // here we calculate the dividend needed to scale the feedback matrix in order to keep it unitary.
    matrixAttenuation = 1.0f / sqrt(numDelays);
    
    // randomise output tap signs using hadamard Transform to even out output
//    printf("\n tapsigns before: \n");
    for (int i = 0; i < numDelays; i++){
        tapSigns[i] = pow(-1.0, i * 31 % 71);
//        printf(" tap %d %f, \n ", i,tapSigns[i]);
    }

    hadamardTransform(tapSigns, temp3, numDelays);
    
    memcpy(tapSigns, temp3, sizeof(float) * numDelays);
    
//    printf("\n tapsigns after: \n");
//    for (int i = 0; i < numDelays; i++){
//        printf("tap %d %f, \n", i, tapSigns[i]);
//    }
    
    // compute decay time gain. Set to -ve value for lossless FDN
    resetDelay(-1.0f);
    
    
}



void FDN::impulseResponse_write(long numSamples, std::ofstream* outputFile){
    float left, right;
    float zero = 0.0f;
    float one = 1.0f;
    
    processAudio(&one, &left, &right);
    *outputFile << left << ",";
    for (int i = 1; i < numSamples; i++){
        processAudio(&zero, &left, &right);
        *outputFile << left << ",";
    }
    
    *outputFile << "\n";
}


void FDN::impulseResponse_output(long numSamples, float* output){
//    float left, right;
    float left = 0.0f;
    float right = 0.0f;
    float zero = 0.0f;
    float one = 1.0f;
    
    processAudio(&one, &left, &right);
//    printf("left %f \n", left);
    output[0] = left;
    
    for (int i = 1; i < numSamples; i++){
        processAudio(&zero, &left, &right);
//        printf("left %f \n", left);
        output[i] = left;
    }
    
}


inline void FDN::processAudio(float* pInput, float* pOutputL, float* pOutputR)
{
    // Read the Input
    float xn = *pInput;
    *pOutputL = 0;
//    *pOutputR = 0;
    int i;
    
    // attenuate input volume
    xn = (inputAttenuation * xn);
    
    // mono output copied to both channels
    for (i = 0; i < numDelays; i++) {
        // copy delay line signals to output buffer and attenuate
        outputs[i] = delayBuffers[rwIndices[i]] * tapGains[i];
//        printf("outputs[%i] %f tapsigns[%i] %f", i, outputs[i], i, tapSigns[i]);
        *pOutputL += outputs[i] * tapSigns[i];
    }
    
    // process the feedback matrix
    hadamardTransform(outputs, temp3, numDelays);
    for (int i = 0; i < numDelays; i++) delayBuffers[rwIndices[i]] = temp3[i] + xn;

    incrementIndices();
}

void swap(int* a, int b, int c){
    int tmp = a[b];
    a[b] = a[c];
    a[c] = tmp;
}

// randomly permute the elements of a without swapping any delay times between channels
void FDN::randomPermutation(int* a, int length, int audioChannels){
    for (int j = 0; j < audioChannels; j++){
        for (int i = j; i < length; i += audioChannels){
            updateRand();
            swap(a, i, j + (random_val % (length / audioChannels))*audioChannels);
        }
    }
}

void FDN::randomPermutation1Channel(int* a, int length, int audioChannels){
    int j = 0;
    for (int i = j; i < length; i += audioChannels){
        updateRand();
        swap(a, i, j + (random_val % (length / audioChannels))*audioChannels);
    }
}



void FDN::resetDelay(float decayTime)
{
    // clear the buffers
    memset(delayBuffers, 0, sizeof(float)* totalDelayTime);
    memset(outputs, 0, sizeof(float)* NUMDELAYSEXT);
    
    // init read indices
    resetReadIndices();
    
    setDecayTime(decayTime);
    randomSeed = 0;
}

inline void FDN::incrementIndices(){
    for (int i = 0; i < numDelays; i++){
        rwIndices[i]++;
        if (rwIndices[i] >= endIndices[i])
            rwIndices[i] = startIndices[i];
    }
}

//
//// see: http://software.intel.com/en-us/articles/fast-random-number-generator-on-the-intel-pentiumr-4-processor/
//// rather than returning an output, this function updates a class variable so that we only have to generate 1 random number for each sample.
inline void FDN::updateRand()
{
    randomSeed = (214013 * randomSeed + 2531011);
    random_val = (randomSeed >> 16) & 0x7FFF;
}



float max(float* window, int length){
    float max = 0.0;
    for (int i = 0; i < length; i++){
        if (abs(window[i]) > max) max = abs(window[i]);
    }
    return max;
}

int thresholdCount(float* window, int length, float threshold){
    int count = 0;
    for (int i = 0; i < length; i++){
        if (abs(window[i]) > threshold) count++;
    }
    return count;
}



// Helper for FHT
inline void FDN::hadamard16(float* input, float* output){
    // level 1
    // +
    temp16[0] = input[0] + input[8];
    temp16[1] = input[1] + input[9];
    temp16[2] = input[2] + input[10];
    temp16[3] = input[3] + input[11];
    temp16[4] = input[4] + input[12];
    temp16[5] = input[5] + input[13];
    temp16[6] = input[6] + input[14];
    temp16[7] = input[7] + input[15];
    // -
    temp16[8] = input[0] - input[8];
    temp16[9] = input[1] - input[9];
    temp16[10] = input[2] - input[10];
    temp16[11] = input[3] - input[11];
    temp16[12] = input[4] - input[12];
    temp16[13] = input[5] - input[13];
    temp16[14] = input[6] - input[14];
    temp16[15] = input[7] - input[15];
    
    // level 2.1
    //+
    output[0] = temp16[0] + temp16[4];
    output[1] = temp16[1] + temp16[5];
    output[2] = temp16[2] + temp16[6];
    output[3] = temp16[3] + temp16[7];
    //-
    output[4] = temp16[0] - temp16[4];
    output[5] = temp16[1] - temp16[5];
    output[6] = temp16[2] - temp16[6];
    output[7] = temp16[3] - temp16[7];
    
    //+
    output[8] = temp16[8] + temp16[12];
    output[9] = temp16[9] + temp16[13];
    output[10] = temp16[10] + temp16[14];
    output[11] = temp16[11] + temp16[15];
    //-
    output[12] = temp16[8] - temp16[12];
    output[13] = temp16[9] - temp16[13];
    output[14] = temp16[10] - temp16[14];
    output[15] = temp16[11] - temp16[15];
    
    // level 3
    // +
    temp16[0] = output[0] + output[2];
    temp16[1] = output[1] + output[3];
    // -
    temp16[2] = output[0] - output[2];
    temp16[3] = output[1] - output[3];
    
    // +
    temp16[4] = output[4] + output[6];
    temp16[5] = output[5] + output[7];
    // -
    temp16[6] = output[4] - output[6];
    temp16[7] = output[5] - output[7];
    
    // +
    temp16[8] = output[8] + output[10];
    temp16[9] = output[9] + output[11];
    // -
    temp16[10] = output[8] - output[10];
    temp16[11] = output[9] - output[11];
    
    // +
    temp16[12] = output[12] + output[14];
    temp16[13] = output[13] + output[15];
    // -
    temp16[14] = output[12] - output[14];
    temp16[15] = output[13] - output[15];
    
    // level 4
    output[0] = temp16[0] + temp16[1];
    output[1] = temp16[0] - temp16[1];
    output[2] = temp16[2] + temp16[3];
    output[3] = temp16[2] - temp16[3];
    
    output[4] = temp16[4] + temp16[5];
    output[5] = temp16[4] - temp16[5];
    output[6] = temp16[6] + temp16[7];
    output[7] = temp16[6] - temp16[7];
    
    output[8] = temp16[8] + temp16[9];
    output[9] = temp16[8] - temp16[9];
    output[10] = temp16[10] + temp16[11];
    output[11] = temp16[10] - temp16[11];
    
    output[12] = temp16[12] + temp16[13];
    output[13] = temp16[12] - temp16[13];
    output[14] = temp16[14] + temp16[15];
    output[15] = temp16[14] - temp16[15];
}

// input and output must not be the same array
inline void FDN::hadamardTransform(float* input, float* output, size_t length){
    size_t partitionSize = length;
    float* tmpIn = input;
    float* tmpOut = temp1;
    bool tmpState = true;
    
    // iteratively calculated recursion
    while (partitionSize > 16) {
        size_t halfPartitionSize = partitionSize >> 1;
        for (size_t i = 0; i < length; i += partitionSize){
            // copy all the lower terms into place
            memcpy(tmpOut+i,tmpIn+i,sizeof(float)*halfPartitionSize);
            memcpy(tmpOut+i+halfPartitionSize,tmpIn+i,sizeof(float)*halfPartitionSize);
            // sum all the higher terms into place
            for (size_t j=i; j<halfPartitionSize+i; j++) {
                size_t idx2 = j+halfPartitionSize;
                tmpOut[j] += tmpIn[idx2];
                tmpOut[idx2] -= tmpIn[idx2];
            }
        }
        
        // swap temp buffers to avoid using the same memory for reading and writing.
        tmpIn = tmpOut;
        tmpState = !tmpState;
        if (tmpState) tmpOut = temp1;
        else tmpOut = temp2;
        partitionSize >>= 1;
    }
    
    // base case
    for (int i = 0; i < length; i += 16) hadamard16(tmpIn + i, output + i);
    
    // scale the output to make it unitary
    vDSP_vsmul(output,1,&matrixAttenuation,output,1,(size_t) length);
}

void FDN::resetReadIndices(){
    for (int i = 0; i < numDelays; i++){
        rwIndices[i] = startIndices[i];
    }
}


// computes the appropriate feedback gain attenuation
// to get a decay envelope with the specified RT60 time (in seconds)
// from a delay line of the specified length.
//
// This formula comes from solving EQ 11.33 in DESIGNING AUDIO EFFECT PLUG-INS IN C++ by Will Pirkle
// which is attributed to Jot, originally.
double gain(double rt60, double delayLengthInSamples) {
    return pow(M_E, (-3.0 * delayLengthInSamples) / (rt60 * 44100.0));
}

// larger values decay more slowly
void FDN::setDecayTime(double rt60){
    if (rt60 < 0) {
        for (int i = 0; i < numDelays; i++){
            tapGains[i] = 1.0f;
        }
    }
    else if (rt60 == 0) {
        for (int i = 0; i < numDelays; i++){
            tapGains[i] = 0.0f;
        }
    }
    else {
        for (int i = 0; i < numDelays; i++) {
            tapGains[i] = gain(rt60, delayTimes[i]);
            //            tapGains[i] = 1.0f;
        }
    }
}


// Methods for setting delay times

void FDN::setDelayTimesVelvetNoise(){
    
    // generate randomised delay tap outputs. See (http://users.spa.aalto.fi/mak/PUB/AES_Jarvelainen_velvet.pdf)
    //    float maxDelayTime = 0.100f * 44100.0f;
    //    float minDelayTime = 0.007f * 44100.0f;
    float minDelayTime = RV_MIN_DELAY_TIME;
    float maxDelayTime = RV_MAX_DELAY_TIME;
    
    float outTapSpacing = (float)(maxDelayTime - minDelayTime) / (float)numDelays;
    randomSeed = std::rand() ;
    updateRand();
    // Set output tap times
    totalDelayTime = 0;
    
    
    
    float scale = 1.0;
    for (int i = 0; i < numDelays; i++){
        delayTimes[i] = minDelayTime + outTapSpacing*((float)i + 0.5f);
        float jitter = ((float)randomSeed / (float)RAND_MAX) * (outTapSpacing);
        delayTimes[i] += jitter;
        delayTimes[i] *= scale;
        totalDelayTime += delayTimes[i];
    }
    
    
    
    // randomly shuffle the order of delay times in the array
    randomPermutation1Channel(delayTimes, numDelays, 1);
    
}



// set even spacing and then find the nearest prime numbers to that.
void FDN::setDelayTimesVelvetPrime1(){
    
    float minDelayTime = RV_MIN_DELAY_TIME;
    float maxDelayTime = RV_MAX_DELAY_TIME;
    
    // find the appropriate spacing to make evenly spaced delay times
    float outTapSpacing = (float)(maxDelayTime - minDelayTime) / (float)numDelays;

    // Set output tap times
    totalDelayTime = 0;
    
    
    // set delay times to the nearest prime number to even spacing
    for (int i = 0; i < numDelays; i++){
        delayTimes[i] = minDelayTime + outTapSpacing*(float)i;
        delayTimes[i] = RV_nearestPrime(delayTimes[i]);
        totalDelayTime += delayTimes[i];
    }
    
    
    
    // randomly shuffle the order of delay times in the array
    randomPermutation1Channel(delayTimes, numDelays, 1);
    
}


/*
 * Do the complete velvet noise algorithm but then move the times to nearest
 * prime number
 */
void FDN::setDelayTimesVelvetPrime2(){
    
    // generate randomised delay tap outputs. See (http://users.spa.aalto.fi/mak/PUB/AES_Jarvelainen_velvet.pdf)
    //    float maxDelayTime = 0.100f * 44100.0f;
    //    float minDelayTime = 0.007f * 44100.0f;
    float minDelayTime = RV_MIN_DELAY_TIME;
    float maxDelayTime = RV_MAX_DELAY_TIME;
    
    float outTapSpacing = (float)(maxDelayTime - minDelayTime) / (float)numDelays;
    randomSeed = std::rand() ;
    updateRand();
    // Set output tap times
    totalDelayTime = 0;
    
    
    

    for (int i = 0; i < numDelays; i++){
        delayTimes[i] = minDelayTime + outTapSpacing*((float)i + 0.5f);
        float jitter = ((float)randomSeed / (float)RAND_MAX) * (outTapSpacing);
        delayTimes[i] += jitter;
        delayTimes[i] = RV_nearestPrime(delayTimes[i]);
        totalDelayTime += delayTimes[i];
    }
    
    
    
    // randomly shuffle the order of delay times in the array
    randomPermutation1Channel(delayTimes, numDelays, 1);
    
}

/*
 * returns true if x is in array
 *
 * @param array    an array of length "length"
 * @param x        the number to search for
 * @param length   length of array
 */
bool arrayContains(int* array, int x, size_t length){
    for(size_t i=0; i<length; i++)
        if(array[i]==x) return true;
    
    return false;
}


/*
 * return a random number in the range [lowerBound,upperBound]
 */
int randInRange(int lowerBound, int upperBound){
    assert(upperBound >= lowerBound);
    
    uint32_t newRand = arc4random();
    
    uint32_t rangeWidth = upperBound - lowerBound + 1;
    
    return lowerBound + (int)(newRand % rangeWidth);
}



/*
 * return a random prime number in the range [lowerBound,upperBound]
 */
int randPrimeInRange(int lowerBound, int upperBound){
    assert(upperBound >= lowerBound);
    
    uint32_t newRand = arc4random();
    
    uint32_t rangeWidth = upperBound - lowerBound + 1;
    
    int randomNumber = lowerBound + (int)(newRand % rangeWidth);
    
    return RV_nearestPrime(randomNumber);
}



/*
 * return the number of primes between lowerBound and upperBound
 */
int countPrimesInRange(int lowerBound, int upperBound){
    assert(upperBound >= lowerBound);
    assert(lowerBound >= 0);
    
    // reset lowerBound to the smallest prime >= lowerBound
    lowerBound = RV_minPrime(lowerBound);
    
    // reset upperBound to the largest prime <= upperBound
    upperBound = RV_maxPrime(upperBound);
    
    int primeCount = 1;
    int currentPrime = lowerBound;
    for(int i = lowerBound+1; i<=upperBound; i++){
        // if there is a new prime <= i
        int maxPrimeBoundedByI = RV_maxPrime(i);
        if (maxPrimeBoundedByI != currentPrime){
            currentPrime = maxPrimeBoundedByI;
            primeCount++;
        }
    }
    
    return primeCount;
}



/*
 * Set delay times to random numbers within bounds
 */
void FDN::setDelayTimesRandom(){
    
    // Set output tap times
    totalDelayTime = 0;
    
    // clear the delay times array
    memset(delayTimes, 0, sizeof(int)*numDelays);
    
    for (int i = 0; i < numDelays; i++){
        // generate a random delay time in the specified range
        int nextDelayTime = randInRange(RV_MIN_DELAY_TIME, RV_MAX_DELAY_TIME);
        
        // because we are not allowing duplicate delay times, we need to
        // ensure that the range is sufficiently wide so that it is possible
        // to fill the array without duplicating any values
        assert(numDelays <= RV_MAX_DELAY_TIME - RV_MIN_DELAY_TIME + 1);
        
        // if the time we got is already in the array, choose another one
        // until we find a unique one
        while (arrayContains(delayTimes, nextDelayTime, numDelays))
            nextDelayTime = randInRange(RV_MIN_DELAY_TIME, RV_MAX_DELAY_TIME);
        
        // store the unique random result to the ith position in the array
        delayTimes[i] = nextDelayTime;
    }
    
    
    
    // randomly shuffle the order of delay times in the array
    randomPermutation1Channel(delayTimes, numDelays, 1);
    
}



/*
 * Set delay times to random prime numbers within bounds
 */
void FDN::setDelayTimesRandomPrime(){
    
    // Set output tap times
    totalDelayTime = 0;
    
    // clear the delay times array
    memset(delayTimes, 0, sizeof(int)*numDelays);
    
    for (int i = 0; i < numDelays; i++){
        // generate a random delay time in the specified range
        int nextDelayTime = randInRange(RV_MIN_DELAY_TIME, RV_MAX_DELAY_TIME);
        
        // because we are not allowing duplicate delay times, we need to
        // ensure that the range is sufficiently wide so that it is possible
        // to fill the array without duplicating any values
        assert(numDelays <= countPrimesInRange(RV_MIN_DELAY_TIME, RV_MAX_DELAY_TIME));
        
        // if the time we got is already in the array, choose another one
        // until we find a unique one
        while (arrayContains(delayTimes, nextDelayTime, numDelays))
            nextDelayTime = randPrimeInRange(RV_MIN_DELAY_TIME, RV_MAX_DELAY_TIME);
        
        // store the unique random prime result to the ith position in the array
        delayTimes[i] = nextDelayTime;
    }
    
    
    
    // randomly shuffle the order of delay times in the array
    randomPermutation1Channel(delayTimes, numDelays, 1);
    
}
