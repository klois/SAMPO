/*
 * @author      Klaus Kofler
 * @date		07/01/2013
 */
#pragma once

#include "device_types.h"
#include "agent.h"

// taken from GPU Gems3

// S1, S2, S3, and M are all constants, and z is part of the
// private per-thread generator state.
unsigned TausStep(unsigned z, int S1, int S2, int S3, unsigned M)
{
	unsigned b=(((z << S1) ^ z) >> S2);
	return (((z & M) << S3) ^ b);
}

// A and C are constants
unsigned LCGStep(unsigned z, unsigned A, unsigned C)
{
	return (A*z+C);
}

/*
 * uniform random number generator
 * @param z1 random integer
 * @param z2 random integer
 * @param z3 random integer
 * @param z4 random integer
 */
REAL HybridTaus(__constant struct Seeds* seeds)
{
	UINT gid = get_global_id(0);
	// Tausworthe state components should be greater than 128
	UINT z1 = seeds->z1 + 128u;
	UINT z2 = seeds->z2 + 128u + gid;
	UINT z3 = seeds->z3 + 128u;
	// multiplying by gid in order to have different seed for each thread
	UINT z4 = seeds->z4 * gid;

	// Combined period is lcm(p1,p2,p3,p4)~ 2^121
	return 2.3283064365387e-10f * (              // Periods
		TausStep(z1, 13, 19, 12, 4294967294UL) ^  // p1=2^31-1
		TausStep(z2, 2, 25, 4, 4294967288UL) ^    // p2=2^30-1
		TausStep(z3, 3, 11, 17, 4294967280UL) ^   // p3=2^28-1
		LCGStep(z4, 1664525, 1013904223UL)       // p4=2^32
	);
}

/*
 * BoxMuller algorithm to generate a normal distributed random number out of two uniform distributed random numbers
 * @param seeds an array with at least two seeds
 */
REAL BoxMuller(__constant struct Seeds* seeds, REAL mean, REAL stdDev)
{
	UINT offset = get_global_size(0);
	REAL u0=HybridTaus(&seeds[0]), u1=HybridTaus (&seeds[1]);
	REAL r=sqrt(-2*log(u0));
	REAL theta=2*M_PI_F*u1;
	return (r*sin(theta)) * stdDev + mean;
}
