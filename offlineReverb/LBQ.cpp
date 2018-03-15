//
//  LBQ.cpp
//  offlineReverb
//
//  Created by Natalie Agus on 13/3/18.
//  Copyright © 2018 Hans. All rights reserved.
//

#include "LBQ.hpp"

LBQ::LBQ(int N, int m){
    
    this->N = N;
    this->m = m;
    this->Q_coefficient = N * (N + 2);
    
    this->zero = 0.0f;
    this->one = 1.0f;
    
    this->padded_x = (float*) malloc((N+m)*sizeof(float));
    this->r_k = (float*) malloc(N*sizeof(float));
    this->temp = (float*) malloc(N*sizeof(float));
    this->LBQ_denominator = (float*) malloc(m*sizeof(float));
    
    float initVal = this->N-1; // N - 1
    float rampVal = -1.0f;
    vDSP_vramp(&initVal, &rampVal, LBQ_denominator, 1, m);
    
}



float LBQ::LBQtest(float* x){
    
    //copy to padded_x
    memcpy(padded_x, x, N * sizeof(float));
    
    //set the trailing to zero
    //memcpy will screw up unused padded_x so there's no point setting padded_x to zero for N+m before memcpy
    memset(padded_x+N, 0, m * sizeof(float));
    
//    for(int i =0; i<N+m; i++){
//        printf("Signal %f \n", padded_x[i]);
//    }

    //get sum squared of the signal
    vDSP_vsq(padded_x, 1, temp, 1, N);
    vDSP_sve(temp, 1, &sum_sq, N);
    
//    printf("Sum squared : %f \n", sum_sq);
    
    //autocorrelation to obtain r_k
    vDSP_conv(padded_x, 1, padded_x, 1, r_k, 1, m+1, N);
    vDSP_vsdiv(r_k, 1, &sum_sq, r_k, 1, N);
    
//    for(int i =0; i<m; i++){
//        printf("LBQ DENOM %i : %f \n", i, *(LBQ_denominator+i));
//        //        std::cout << output[i] << "\n";
//    }
    
//    for(int i =0; i<m; i++){
//        printf("r_k %f \n", *(r_k+1+i));
//        //        std::cout << output[i] << "\n";
//    }
    
    //squared autocorrelation
    vDSP_vsq(r_k, 1, r_k, 1, m+1);
    
//    for(int i =0; i<m; i++){
//        printf("squared r_k %f \n", *(r_k+1+i));
//        //        std::cout << output[i] << "\n";
//    }
    
    //divide by the denominator
    vDSP_vdiv(LBQ_denominator, 1, r_k+1, 1, temp, 1, m);
    
//    for(int i =0; i<m; i++){
//        printf("Q %i : %f \n", i, *(temp+i));
//        //        std::cout << output[i] << "\n";
//    }
//    
    //get Q-value
    vDSP_sve(temp, 1, &Q_val, m);
    
    return Q_val * Q_coefficient;
    
}
