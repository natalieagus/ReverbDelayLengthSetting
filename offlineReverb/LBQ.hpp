//
//  LBQ.hpp
//  offlineReverb
//
//  Created by Natalie Agus on 13/3/18.
//  Copyright © 2018 Hans. All rights reserved.
//

#ifndef LBQ_hpp
#define LBQ_hpp

#include <stdio.h>
#include <Accelerate/Accelerate.h>

class LBQ
{
public:
    LBQ(int N, int m);
    float LBQtest(float* x);
    
private:
    int N, m;
    float zero, one;
    float sum_sq;
    float Q_val;
    float Q_coefficient;
    
    float* padded_x;
    float* r_k;
    float* temp;
    float* LBQ_denominator;

};

#endif /* LBQ_hpp */
