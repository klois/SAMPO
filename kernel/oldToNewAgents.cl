/*
 * @author      Klaus Kofler
 * @date		07/01/2013
 */

#include "device_types.h"
#include "agent.h"

/*
void printRange(struct AgentRange* range) {
	printf("[%d\t%d]\t%d\r\n", range->start, range->end, range->end - range->start);
}

void printPopulation(struct Population* pop) {
	printf("Eggs \t\t"); printRange(&pop->eggs);
	printf("Larvae \t\t"); printRange(&pop->larvae);
	printf("Pupae \t\t"); printRange(&pop->pupae);
	printf("Immatures \t"); printRange(&pop->immatures);
	printf("MateSeeking \t"); printRange(&pop->mateSeekings);
	printf("BMS \t\t"); printRange(&pop->bmSeekings);
	printf("BMD \t\t"); printRange(&pop->bmDigestings);
	printf("Gravids \t"); printRange(&pop->gravids);
}
*/

void calcNewRange(struct AgentRange* newRange, UINT startingPoint, UINT size) {
	newRange->start = ((startingPoint + LOCAL_SIZE) / LOCAL_SIZE) * LOCAL_SIZE;
	newRange->end = newRange->start + size;
}

/*
 * Kernel to copy the old agents to the new position in newAgents using the various prefixSum to determine the index in the new array. It also writes the
 * properties of the next generation's population to position 1 in population
 * @param oldAgents array holding all agents of the current iteration
 * @param oldAgendAges array holding all agents' age information of the current iteration
 * @param oldAgendStates array holding all agents' state information of the current iteration
 * @param newAgents array that holds the eggs generated in the current iteration and to which the still living agents of the current iteration will be added
 * @param newAgentAges array that holds the age of the eggs generated in the current iteration and to which the age information of the still living
 * 			agents of the current iteration will be added
 * @param newAgentStates array that holds the state of the eggs generated in the current iteration and to which the state information of the still living
 * 			agents of the current iteration will be added
 * @param prefixSum1 array with a prefix sum for all eggs, immatures and blood meal digestings of the current iteration
 * @param prefixSum2 array with a prefix sum for all larvae, mate seekings and gravids of the current iteration
 * @param prefixSum3 array with a prefix sum for all pupae and blood meal seekings of the current iteration
 * @param population the properties (number of agents in certain state) of the current population at positon 0. The properties of the next iteraiont's
 * 			population will be written to position 1.
 * @param enb eggs and biomass struct to get the number of newly generated eggs
 */
__kernel void oldToNewAgents(__global struct Agent* oldAgents, __global struct AgentAge* oldAgentAges, __global struct AgentState* oldAgentStates,
		__global struct Agent* newAgents, __global struct AgentAge* newAgentAges, __global struct AgentState* newAgentStates,
		__global INT* prefixSum1, __global INT* prefixSum2, __global INT* prefixSum3,
		__global struct Population* population,	__constant struct EggsNbiomass* enb) {

	UINT gid = get_global_id(0);
	struct Population pop = population[0];

	if(gid == pop.gravids.end) { // write new properties array
#define USE_PRIVATE 0
#if USE_PRIVATE
		// NVIDIA cannot handle this version
		struct Population privateCopy;
		struct Population* offset = &privateCopy;
#else
		// good for all, max 1% overall slowdown
		__global struct Population* offset = &population[1];
#endif
		offset->eggs.start = 0;
		offset->eggs.end = prefixSum1[pop.eggs.end] + 1 + enb->newEggs;
//		printf("old %d, new %d + %d\n", pop.eggs.end, prefixSum1[pop.eggs.end-1], environment->newEggs);
/*
for(INT i = 0; i < environment->newEggs; ++i)
	if(newAgents[i].state != EGG)
		printf("\t%d is not an egg %d\n", i, newAgents[i].state);
*/
		// assure a space between states of at leas one in order to avoid overwriting of last element of prefix sum
		offset->larvae.start = ((offset->eggs.end + LOCAL_SIZE) / LOCAL_SIZE) * LOCAL_SIZE;
		offset->larvae.end = offset->larvae.start + prefixSum2[pop.larvae.end] + 1;
//		calcNewRange(&offset->larvae, offset->eggs.end, prefixSum2[pop.larvae.end]);
		offset->pupae.start = ((offset->larvae.end + LOCAL_SIZE) / LOCAL_SIZE) * LOCAL_SIZE;
		offset->pupae.end = offset->pupae.start + prefixSum3[pop.pupae.end-pop.larvae.start] + 1;
//		calcNewRange(&offset->pupae, offset->larvae.end, prefixSum3[pop.pupae.end]);
		offset->immatures.start = ((offset->pupae.end + LOCAL_SIZE) / LOCAL_SIZE) * LOCAL_SIZE;
		offset->immatures.end = offset->immatures.start + prefixSum1[pop.immatures.end] + 1;
//		calcNewRange(&offset->immatures, offset->pupae.end, prefixSum1[pop.immatures.end]);
		offset->mateSeekings.start = ((offset->immatures.end + LOCAL_SIZE) / LOCAL_SIZE) * LOCAL_SIZE;
		offset->mateSeekings.end = offset->mateSeekings.start + prefixSum2[pop.mateSeekings.end] + 1;
//		calcNewRange(&offset->mateSeekings, offset->immatures.end, prefixSum2[pop.mateSeekings.end]);
		offset->bmSeekings.start = ((offset->mateSeekings.end + LOCAL_SIZE) / LOCAL_SIZE) * LOCAL_SIZE;
		offset->bmSeekings.end = offset->bmSeekings.start + prefixSum3[pop.gravids.end-pop.larvae.start] + 1;
//		calcNewRange(&offset->bmSeekings, offset->mateSeekings.end, prefixSum3[pop.gravids.end]);
		offset->bmDigestings.start = ((offset->bmSeekings.end + LOCAL_SIZE) / LOCAL_SIZE) * LOCAL_SIZE;
		offset->bmDigestings.end = offset->bmDigestings.start + prefixSum1[pop.bmDigestings.end] + 1;
//		calcNewRange(&offset->bmDigestings, offset->bmSeekings.end, prefixSum1[pop.bmDigestings.end]);
		offset->gravids.start = ((offset->bmDigestings.end + LOCAL_SIZE) / LOCAL_SIZE) * LOCAL_SIZE;
		offset->gravids.end = offset->gravids.start + prefixSum2[pop.gravids.end] + 1;
//		calcNewRange(&offset.gravids, offset.bmDigestings.end, prefixSum2[pop.gravids.end]);
#if USE_PRIVATE
		population[1] = privateCopy;
#endif
	}

	bool valid = !((gid >= pop.eggs.end && gid < pop.larvae.start) ||
				   (gid >= pop.larvae.end && gid < pop.pupae.start) ||
				   (gid >= pop.pupae.end && gid < pop.immatures.start) ||
				   (gid >= pop.immatures.end && gid < pop.mateSeekings.start) ||
				   (gid >= pop.mateSeekings.end && gid < pop.bmSeekings.start) ||
				   (gid >= pop.bmSeekings.end && gid < pop.bmDigestings.start) ||
				   (gid >= pop.bmDigestings.end && gid < pop.gravids.start) ||
				   (gid >= pop.gravids.end));

	if(!valid) return;

	struct AgentState agentState = oldAgentStates[gid];

	if(agentState.dead) {
//printf("%d is dead\n", gid);
		return;
	}

	struct Population offset;
	offset.eggs.start = 0;
	offset.eggs.end = enb->newEggs;

	enum State state = agentState.state;

	if(state == EGG) {
		// start copying in old agents after the newly generated eggs
		UINT newIdx = offset.eggs.end + prefixSum1[gid];
		newAgents[newIdx] = oldAgents[gid];
		newAgentAges[newIdx] = oldAgentAges[gid];
		newAgentStates[newIdx] = oldAgentStates[gid];
		return;
	}

	// add all moved eggs to offset
	offset.eggs.end += (pop.eggs.end > 0) ? prefixSum1[pop.eggs.end] + 1 : 0; // TODO check if

	// round offset to the next multiple of LOCAL_SIZE
	offset.larvae.start = ((offset.eggs.end + LOCAL_SIZE) / LOCAL_SIZE) * LOCAL_SIZE;

	if(state == LARVA) {
		UINT newIdx = offset.larvae.start + prefixSum2[gid];
		newAgents[newIdx] = oldAgents[gid];
		newAgentAges[newIdx] = oldAgentAges[gid];
		newAgentStates[newIdx] = oldAgentStates[gid];
		return;
	}

	// add all moved larvae to offset
	offset.larvae.end = offset.larvae.start + prefixSum2[pop.larvae.end] + 1;
	// round offset to the next multiple of LOCAL_SIZE
	offset.pupae.start = ((offset.larvae.end + LOCAL_SIZE) / LOCAL_SIZE) * LOCAL_SIZE;

	if(state == PUPA) {
//printf("bmd %d %d\n", gid, offset.larvae.start);
		UINT newIdx = offset.pupae.start + prefixSum3[gid-pop.larvae.start];
		newAgents[newIdx] = oldAgents[gid];
		newAgentAges[newIdx] = oldAgentAges[gid];
		newAgentStates[newIdx] = oldAgentStates[gid];
		return;
	}

	// add all moved pupae to offset
	offset.pupae.end = offset.pupae.start + prefixSum3[pop.pupae.end-pop.larvae.start] + 1;
	// round offset to the next multiple of LOCAL_SIZE
	offset.immatures.start = ((offset.pupae.end + LOCAL_SIZE) / LOCAL_SIZE) * LOCAL_SIZE;

	if(state == IMMATURE) {
		UINT newIdx = offset.immatures.start + prefixSum1[gid];
		newAgents[newIdx] = oldAgents[gid];
		newAgentAges[newIdx] = oldAgentAges[gid];
		newAgentStates[newIdx] = oldAgentStates[gid];
		return;
	}

	// add all moved immatures to offset
	offset.immatures.end = offset.immatures.start + prefixSum1[pop.immatures.end] + 1;
	// round offset to the next multiple of LOCAL_SIZE
	offset.mateSeekings.start = ((offset.immatures.end + LOCAL_SIZE) / LOCAL_SIZE) * LOCAL_SIZE;

	if(state == MATESEEKING) {
//printf("%d : %d \n", prefixSum1[gid], agent.state);
		UINT newIdx = offset.mateSeekings.start + prefixSum2[gid];
		newAgents[newIdx] = oldAgents[gid];
		newAgentAges[newIdx] = oldAgentAges[gid];
		newAgentStates[newIdx] = oldAgentStates[gid];
		return;
	}

	// add all moved mate seekings to offset
	offset.mateSeekings.end = offset.mateSeekings.start + prefixSum2[pop.mateSeekings.end] + 1;
	// round offset to the next multiple of LOCAL_SIZE
	offset.bmSeekings.start = ((offset.mateSeekings.end + LOCAL_SIZE) / LOCAL_SIZE) * LOCAL_SIZE;

	if(state == BMS) {
		UINT newIdx = offset.bmSeekings.start + prefixSum3[gid-pop.larvae.start];
		newAgents[newIdx] = oldAgents[gid];
		newAgentAges[newIdx] = oldAgentAges[gid];
		newAgentStates[newIdx] = oldAgentStates[gid];
		return;
	}

	// add all moved blood meal seekings to offset
	offset.bmSeekings.end = offset.bmSeekings.start + prefixSum3[pop.gravids.end-pop.larvae.start] + 1;
	// round offset to the next multiple of LOCAL_SIZE
	offset.bmDigestings.start = ((offset.bmSeekings.end + LOCAL_SIZE) / LOCAL_SIZE) * LOCAL_SIZE;

	if(state == BMD) {
		UINT newIdx = offset.bmDigestings.start + prefixSum1[gid];
		newAgents[newIdx] = oldAgents[gid];
		newAgentAges[newIdx] = oldAgentAges[gid];
		newAgentStates[newIdx] = oldAgentStates[gid];
		return;
	}

	// add all moved blood meal seekings to offset
	offset.bmDigestings.end = offset.bmDigestings.start + prefixSum1[pop.bmDigestings.end] + 1;
	// round offset to the next multiple of LOCAL_SIZE
	offset.gravids.start = ((offset.bmDigestings.end + LOCAL_SIZE) / LOCAL_SIZE) * LOCAL_SIZE;

	if(state == GRAVID) {
		UINT newIdx = offset.gravids.start + prefixSum2[gid];
		newAgents[newIdx] = oldAgents[gid];
		newAgentAges[newIdx] = oldAgentAges[gid];
		newAgentStates[newIdx] = oldAgentStates[gid];
	}

}
