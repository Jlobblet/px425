/* Random number generator code for PX425 assignments.
 * Uses the Mersenne Twister algorithm.
 *
 * University of Warwick. Nov 2013
 *
 *  Original code can be found at the following link:
 *  http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html
 *
 * A C-program for MT19937, with initialization improved 2002/1/26.
 * Coded by Takuji Nishimura and Makoto Matsumoto.
 *
 * Before using, initialize the state by using init_genrand(seed)
 *
 * Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
 * All rights reserved. */

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* Constant vector a */
#define UPPER_MASK 0x80000000UL /* Most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* Least significant r bits */

unsigned long mt[N];   /* The array for the state vector  */
int mti = N + 1;           /* mti==N+1 means mt[N] is not initialized */

/* Initializes mt[N] with a seed */
void init_genrand(unsigned long s) {

    mti = N + 1;
    mt[0] = s & 0xffffffffUL;
    for (mti = 1; mti < N; mti++) {
        mt[mti] =
                (1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier.
         * In the previous versions, MSBs of the seed affect
         * only MSBs of the array mt[].
         * 2002/01/09 modified by Makoto Matsumoto */
        mt[mti] &= 0xffffffffUL;
        /* For >32 bit machines */
    }
}


/* Generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void) {
    unsigned long y;
    unsigned long mag01[2] = {0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) {
        /* Generate N words at one time */
        int kk;

        /* If init_genrand() has not been called,
         * a default initial seed is used */
        if (mti == N + 1) {
            init_genrand(5489UL);
        }

        for (kk = 0; kk < N - M; kk++) {
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (; kk < N - 1; kk++) {
            y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
            mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
        mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* Generates a random number on [0,1)-real-interval */
double genrand(void) {
    return genrand_int32() * (1.0 / 4294967296.0);
    /* divided by 2^32 */
}

