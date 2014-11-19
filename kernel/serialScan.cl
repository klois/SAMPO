/*
 * @author      Klaus Kofler
 * @date		07/01/2013
 */

#include "device_types.h"
#include "agent.h"

/*
 * calculates a prefix sum, using only a single thread. The prefix sum is stored at the beginning of prefixSum
 * @param agentStates The array holding all agents' state information
 * @param prefixSum Oversized array, in which's beginning the prefix sum will be stored
 * @param pop the properties (number of agents in certain state) of the current population
 * @param offset the offset relative to the thread index, where the prefix sum should be stored in prefixSum
 * @param start the starting index for the prefix sum in agents
 * @param upperBound the end of the prefix sum in agents
 * @param match the state which will be treated as a 1 for the prefix sum. All other elements will be threaded as 0
 */
__kernel void serialScan(__global struct AgentState* agentStates, __global INT* prefixSum, __constant struct Population* pop,
		UINT offset, UINT start, UINT upperBound, enum State match) {
	INT curr = -1;

	for(UINT idx = start; idx <= upperBound; ++idx) {
		bool valid = !((idx >= pop->eggs.end && idx < pop->larvae.start) ||
					   (idx >= pop->larvae.end && idx < pop->pupae.start) ||
					   (idx >= pop->pupae.end && idx < pop->immatures.start) ||
					   (idx >= pop->immatures.end && idx < pop->mateSeekings.start) ||
					   (idx >= pop->mateSeekings.end && idx < pop->bmSeekings.start) ||
					   (idx >= pop->bmSeekings.end && idx < pop->bmDigestings.start) ||
					   (idx >= pop->bmDigestings.end && idx < pop->gravids.start) ||
					   (idx >= pop->gravids.end));

//		if(match == EGG)
//			printf("%d %s %d %d\n", idx, (valid ? "VALID" : "NOT valid"), curr, valid ? agentStates[idx].state : -1);

		if(valid && !agentStates[idx].dead && (agentStates[idx].state == match)) ++curr;

		prefixSum[idx - offset] = curr;

//		if(idx == 100) printf("screw you %d %d %d\n", match, idx, curr);
	}
}


/*
 * version that can calculate the indices of all states in only a single call. Incompatible interface with parallel prefix sum and therefore unused
 * gives an overal performance improvement of approximately 5% on CPUs
 * @param agentStates The array holding all agents' state information
 * @param prefixSum Oversized array, in which's beginning the prefix sum will be stored
 * @param pop the properties (number of agents in certain state) of the current population
 * @param offset the offset relative to the thread index, where the prefix sum should be stored in prefixSum
 * @param start the starting index for the prefix sum in agents
 * @param upperBound the end of the prefix sum in agents
 * @param match the state which will be treated as a 1 for the prefix sum. All other elements will be threaded as 0
 */
__kernel void serialScan0(__global struct AgentState* agentStates, __global INT* prefixSum, __constant struct Population* pop,
		UINT offset, UINT start, UINT upperBound, enum State match) {
	INT egg = -1;
	INT larva = -1;
	INT pupa = -1;
	INT immature = -1;
	INT mateSeeking = -1;
	INT bloodMealSeeking = -1;
	INT bloodMealDigesting = -1;
	INT gravid = -1;

	for(UINT idx = start; idx <= upperBound+1; ++idx) {
		bool valid = !((idx >= pop->eggs.end && idx < pop->larvae.start) ||
					   (idx >= pop->larvae.end && idx < pop->pupae.start) ||
					   (idx >= pop->pupae.end && idx < pop->immatures.start) ||
					   (idx >= pop->immatures.end && idx < pop->mateSeekings.start) ||
					   (idx >= pop->mateSeekings.end && idx < pop->bmSeekings.start) ||
					   (idx >= pop->bmSeekings.end && idx < pop->bmDigestings.start) ||
					   (idx >= pop->bmDigestings.end && idx < pop->gravids.start) ||
					   (idx >= pop->gravids.end));

//		if(match == EGG)
//			printf("%d %s %d %d\n", idx, (valid ? "VALID" : "NOT valid"), curr, valid ? agentStates[idx].state : -1);

		if(valid && !agentStates[idx].dead) {
			if(agentStates[idx].state == EGG) {
				++egg;
				prefixSum[idx - offset] = egg;
			}
			if(agentStates[idx].state == LARVA) {
				++larva;
				prefixSum[idx - offset] = larva;
			}
			if(agentStates[idx].state == PUPA) {
				++pupa;
				prefixSum[idx - offset] = pupa;
			}
			if(agentStates[idx].state == IMMATURE) {
				++immature;
				prefixSum[idx - offset] = immature;
			}
			if(agentStates[idx].state == MATESEEKING) {
				++mateSeeking;
				prefixSum[idx - offset] = mateSeeking;
			}
			if(agentStates[idx].state == BMS) {
				++bloodMealSeeking;
				prefixSum[idx - offset] = bloodMealSeeking;
			}
			if(agentStates[idx].state == BMD) {
				++bloodMealDigesting;
				prefixSum[idx - offset] = bloodMealDigesting;
			}
			if(agentStates[idx].state == GRAVID) {
				++gravid;
				prefixSum[idx - offset] = gravid;
			}

		} else {
			if(idx == pop->eggs.end)
				prefixSum[idx - offset] = egg;
			if(idx == pop->larvae.end)
				prefixSum[idx - offset] = larva;
			if(idx == pop->pupae.end) {
				prefixSum[idx - offset] = pupa;
//printf("writing %d to %d\n", pupa, idx - offset);
			}
			if(idx == pop->immatures.end)
				prefixSum[idx - offset] = immature;
			if(idx == pop->mateSeekings.end)
				prefixSum[idx - offset] = mateSeeking;
			if(idx == pop->gravids.end+1)
				prefixSum[idx - offset] = bloodMealSeeking;
			if(idx == pop->bmDigestings.end)
				prefixSum[idx - offset] = bloodMealDigesting;
			if(idx == pop->gravids.end)
				prefixSum[idx - offset] = gravid;
		}


//		if(idx == 100) printf("screw you %d %d %d\n", match, idx, curr);
	}
}
