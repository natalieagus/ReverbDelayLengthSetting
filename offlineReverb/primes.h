//
//  primes.h
//  offlineReverb
//
//  Created by hans anderson on 3/14/18.
//  anyone may use this file without restrictions of any kind.
//

#ifndef primes_h
#define primes_h

#ifdef __cplusplus
extern "C" {
#endif


/*
 * recursive binary search for x in RV_primes, starting at idx
 * returns the prime nearest to x
 */
int RV_binarySearch(int minIdx, int maxIdx, int x);




/*
 * returns the prime number nearest to x
 */
int RV_nearestPrime(int x);






#ifdef __cplusplus
}
#endif
    
#endif /* primes_h */
