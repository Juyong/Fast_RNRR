/*

BSD license:

Copyright (c) 2008 Paul Hsieh
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.

    Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

    Neither the name of randutil nor the names of its contributors may be used
    to endorse or promote products derived from this software without
    specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

*/

/*

    Randutil is meant for creating very accurate pseudo random numbers.  The
    problem with the typical methods is that they are either approximations
    or limited to small ranges.  The functions given below use double
    precision floating point for all internal calculations without introducing
    any runaway accuracy problems.  Hence all output should considered to be
    accurate within a few ulps (you can assume a handful of round off errors.)

    Where size_t is used, only values which are round trip convertable through
    a double precision floating point number are valid.  On most systems this
    is either 32 or 53 bits worth.


        int randbiased (double x);

        Returns 1 with a probability of x otherwise a 0.  It is assumed that
        x is between 0.0 and 1.0 inclusive, otherwise x will behave as if
        clamped to that range.


        size_t randslot (const double slots[n-1], size_t n);

        Picks an index from 0 to n-1 inclusive, where the array slots[] is
        at least n-1 elements, and each entry must be a monotonically non-
        decreasing sequence of probabilities.  If the slots array does not
        follow the conditions required, this function may be very slow and
        will output meaningless results.


        size_t randrange (size_t n);

        Picks a value from 0 to n-1 inclusive with equal expected
        likelihood.

        This is commonly implemented as (PRNG_FUNCTION() % n) which is
        inescapably inaccurate if n does not divide evenly into
        (1+PRNG_MAXIMUM) let alone if n is larger than 1+PRNG_MAXIMUM.

        TODO: switch to iterated method if n is too large for double
        precision accuracy.


        double drand (void);

        Picks a uniformly value between 0.0 and 1.0 .


        unsigned int randreseed (void);

        Reseed the random number generator and return the seed used.

        TODO: Generalize by allowing for other sources of entropy via callbacks
        and try to comply with the Fortuna mechanism.

 */

#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include "randutil.h"
#include <iostream>

/* =============================================================================================== */

/*
 *	Install your PRNG so that PRNG_SEED(src, len), PRNG_FUNCTION() and PRNG_MAXIMUM are defined.
 *
 *  PRNG_FUNCTION() outputs an integer type whose range is covered by the double type.  Its range
 *  is between 0 and PRNG_MAXIMUM inclusive.
 *
 *  Like ANSI C's RAND_MAX it is assumed that PRNG_MAXIMUM >= 32767 otherwise drand() will not
 *  work as intended.
 */

static void xseed (const void * src, size_t len) {
    unsigned int seed = 0;
    if (len > sizeof (seed)) len = sizeof (seed);
    if (len > 0) memcpy (&seed, src, len);
    srand (seed);
}

#define PRNG_SEED(src, len) xseed (src, len)
#define PRNG_FUNCTION() rand()
#define PRNG_MAXIMUM RAND_MAX

/* =============================================================================================== */

#define RS_SCALE (1.0 / (1.0 + PRNG_MAXIMUM))

/*
    High quality random number generation functions.  Each function introduces no accuracy issues
    beyond inherent double precision accuracy limits.  Hence the probability distributions should
    alway be considered accurate to within a few ulp.
 */

/*  Returns 0 or 1 with the probability x for picking 0. */

int randbiased (double x) {
    for (;;) {
        double p = PRNG_FUNCTION () * RS_SCALE;
        if (p >= x) return 0;
        if (p + RS_SCALE <= x) return 1;
        /* p < x < p+RS_SCALE */
        x = (x - p) * (1.0 + PRNG_MAXIMUM);
    }
}

/* A non-overflowing average function */
#define average2scomplement(x,y) ((x) & (y)) + (((x) ^ (y))/2)

/* Returns a random number between 0 and n-1 inclusive according to the distribution given in slots[]  */

size_t randslot (const double slots[/* n-1 */], size_t n) {
    double xhi;

    /* Select a random range [x,x+RS_SCALE) */
    double x = PRNG_FUNCTION () * RS_SCALE;

    /* Perform binary search to find the intersecting slot */
    size_t hi = n-2, lo = 0, mi;
    while (hi > lo) {
        mi = average2scomplement (lo, hi);
        if (x >= slots[mi]) lo = mi+1; else hi = mi;
    }

    /* Taking slots[-1]=0.0, this is now true: slots[lo-1] <= x < slots[lo] */

    /* If slots[lo-1] <= x < x+RS_SCALE <= slots[lo] then
       any point in [x,x+RS_SCALE) is in [slots[lo-1],slots[lo]) */

    if ((xhi = x+RS_SCALE) <= slots[lo]) return lo;

    /* Otherwise x < slots[lo] < x+RS_SCALE */

    for (;;) {
        /* x < slots[lo] < xhi */
        if (randbiased ((slots[lo] - x) / (xhi - x))) return lo;
        x = slots[lo];
        lo++;
        if (lo >= n-1) return n-1;

        if (xhi <= slots[lo]) {
            /* slots[lo-1] = x <= xhi <= slots[lo] */
            return lo;
        }
    }
}

/* Returns a (uniformly) randomly chosen number between 0 and n-1 inclusive. */

size_t randrange (size_t n) {
    double resolution = n * RS_SCALE;
    double x = resolution * PRNG_FUNCTION (); /* x in [0,n) */
    double xhi = x + resolution;
    size_t lo = (size_t) floor (x);

    for (;;) {
        lo++;
        if (lo >= xhi || randbiased ((lo - x) / (xhi - x))) return lo-1;
        x = lo;
    }
}

/* Returned a random fraction between 0.0 inclusive and 1.0 exclusive. */

#if (PRNG_MAXIMUM >= 0x7FFFFFFF && DBL_MANT_DIG <= 62) || (PRNG_MAXIMUM >= 0x7FFF && DBL_MANT_DIG <= 30)
/* GCC */
double drand (void) {
    double d;
    do {
        d = ((PRNG_FUNCTION () * RS_SCALE) + PRNG_FUNCTION ()) * RS_SCALE;
    } while (d >= 1.0); /* In case of round off error. */
    return d;
}
#elif PRNG_MAXIMUM < 208064 && DBL_MANT_DIG <= 53
/* MSVC, WATCOM, Turbo C, etc */
double drand (void) {
    double d;
    do {
        d = ((((PRNG_FUNCTION () * RS_SCALE) + PRNG_FUNCTION ()) * RS_SCALE + PRNG_FUNCTION ()) * RS_SCALE + PRNG_FUNCTION ()) * RS_SCALE;
    } while (d >= 1.0); /* In case of round off error. */
    return d;
}
#else
static int drandIterationsRequired (void) {
    double m, n = 0;
    int c = -1;
    do {
        m = n;
        n = (n + PRNG_MAXIMUM) * (1.0 / (1.0 + PRNG_MAXIMUM));
        c++;
    } while (n > m); /* Accuracy limitations will eventually force n <= m */
    return c;
}

static int drandIterationsCount = 0;

double drand (void) {
    double d;

    if (0 == drandIterationsCount)
        drandIterationsCount = drandIterationsRequired ();

    do {
        switch (drandIterationsCount) {
            case 1:	d = PRNG_FUNCTION () * RS_SCALE;
                    break;
            case 2: d = (PRNG_FUNCTION () * RS_SCALE + PRNG_FUNCTION ()) * RS_SCALE;
                    break;
            case 3: d = ((PRNG_FUNCTION () * RS_SCALE + PRNG_FUNCTION ()) * RS_SCALE + PRNG_FUNCTION ()) * RS_SCALE;
                    break;
            case 4: d = (((PRNG_FUNCTION () * RS_SCALE + PRNG_FUNCTION ()) * RS_SCALE + PRNG_FUNCTION ()) * RS_SCALE + PRNG_FUNCTION ()) * RS_SCALE;
                    break;
            case 5: d = ((((PRNG_FUNCTION () * RS_SCALE + PRNG_FUNCTION ()) * RS_SCALE + PRNG_FUNCTION ()) * RS_SCALE + PRNG_FUNCTION ()) * RS_SCALE + PRNG_FUNCTION ()) * RS_SCALE;
                    break;
            default: {
                int i;
                for (d = 0.0, i = 0; i < drandIterationsCount; i++)
                    d = (d + PRNG_FUNCTION ()) * RS_SCALE;
                break;
            }
        }
    } while (d >= 1.0); /* In case of round off error. */
    return d;
}
#endif

/*
 *  Static entropy state.
 */

static struct {
    int which;
    union {
        time_t v;
        unsigned int seed;
    } t;
    union {
        clock_t v;
        unsigned int seed;
    } c;
    union {
        unsigned int v;
        unsigned int seed;
    } counter;
} entropy = {0};

static unsigned char * p = (unsigned char *) (&entropy + 1);
static unsigned int accSeed = 0;

/* Set the PRNG seed from some sources of entropy. */

unsigned int randreseed (void) {
    if (p == ((unsigned char *) (&entropy + 1))) {
        switch (entropy.which) {
            case 0:	entropy.t.v += time (NULL);
                    accSeed ^= entropy.t.v;
                    break;
            case 1:	entropy.c.v += clock();
                    break;
            case 2:	entropy.counter.v++;
                    break;
        }
        entropy.which = (entropy.which + 1) % 3;
        p = (unsigned char *) &entropy.t;
    }
    accSeed = ((accSeed * (UCHAR_MAX+2U)) | 1) + (int) *p;
    p++;
    PRNG_SEED (&accSeed, sizeof (accSeed));
    return accSeed;
}
