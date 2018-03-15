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
     * recursive binary search for x in RV_primes
     * returns the prime nearest to x
     *
     * @param minIdx  lower bound of the search range in RV_primes
     * @param maxIdx  upper bound of the serach range in RV_primes
     * @param x       value to search for nearest primes
     */
    int RV_binarySearchGetValue(int minIdx, int maxIdx, int x);
    


    
    
    /*
     * recursive binary search for x in RV_primes, starting at idx
     * returns the index in RV_primes of the prime nearest to x
     */
    int RV_binarySearchGetIndex(int minIdx, int maxIdx, int x);
    
    
    

    /*
     * returns the prime number nearest to x
     */
    int RV_nearestPrime(int x);
    
    
    
    /*
     * returns the maximum prime number p such that p <= x
     */
    int RV_maxPrime(int x);



    /*
     * returns the minimum prime number p such that p >= x
     */
    int RV_minPrime(int x);



#ifdef __cplusplus
}
#endif
    
#endif /* primes_h */
