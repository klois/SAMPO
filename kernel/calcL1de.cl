/*
 * @author      Klaus Kofler
 * @date		07/01/2013
 */

#include "device_types.h"
#include "agent.h"

/*
 * Calculates partial sums of the larvae one day equivalents in agents, i.e. the total age of all larvae in days
 * One result per group is generated and stored in intermediateL1de
 * @param agents The array holding the required information of all agents
 * @param pop A single structure holding the needed information about the agent populations at position 1. Will be copied to position 0. This is needed to
 * 				avoid synchronization in the update kernel
 * @param intermediateL1de Array will be filled with partial sum of days for each group.
 */
__kernel void calcL1dePerGroup(__global struct AgentAge* agents, __constant struct Population* pop, __global UINT* intermediateL1de) {
	UINT gid = get_global_id(0);
	UINT lid = get_local_id(0);
	UINT offset = pop[1].larvae.start;
	UINT group = get_group_id(0);

	if(pop[1].larvae.end > offset) {
		__local UINT larvae1Dequiv[LOCAL_SIZE];
		larvae1Dequiv[lid] = 0; //initialize with 0, not all platforms do this by default

		for(UINT i = offset + gid; i < pop[1].larvae.end; i += get_global_size(0)) {
			larvae1Dequiv[lid] += (UINT)(agents[i].ageInHours/24.0f);
		}

		// reduction in local memory
		for(UINT currsize = LOCAL_SIZE/2; currsize > 0; currsize /= 2)
		{
			barrier(CLK_LOCAL_MEM_FENCE);
			if(lid < currsize) {
				larvae1Dequiv[lid] += larvae1Dequiv[lid + currsize];
			}
		}

		if(lid != 0) return;

		intermediateL1de[group] = larvae1Dequiv[0];
	}
}


/*
 * Calculates the tolal sums of the todal larva one day equivalents using the patial sum generated by calcFemalesPerGroup
 * @param gPop A single structure holding the needed information about the agent populations
 * @param intermediateL1de Array holding the partial sums as generated by calcL1dePerGroup.
 * @param numberOfSlots The number of partial sums that are provided by calcL1dePerGroup
 */
__kernel void calcL1deTotal(__global struct Population* gPop, __global UINT* intermediateL1de, UINT numberOfSlots) {
	// kernel is inteded to be executed by only one group in order to get the final result
	UINT id = get_local_id(0);
	UINT larvae1DequivPrivate = 0;
	bool anyLarvae = (gPop[1].larvae.end - gPop[1].larvae.start) != 0;

	if(anyLarvae) {
		__local UINT larvae1Dequiv[LOCAL_SIZE];
		larvae1Dequiv[id] = id < numberOfSlots ? intermediateL1de[id] : 0;

		// reduction in local memory
		for(UINT currsize = LOCAL_SIZE/2; currsize > 0; currsize /= 2)
		{
			barrier(CLK_LOCAL_MEM_FENCE);
			if(id < currsize) {
				larvae1Dequiv[id] += larvae1Dequiv[id + currsize];
			}
		}
		larvae1DequivPrivate = larvae1Dequiv[0];
	}

	if(id != 0) return;
//printf("population size on device: %d\n", sizeof(struct Population));
	struct Population pop = gPop[1];

	pop.numLarvae1DayEquiv = larvae1DequivPrivate;

/*
	pop.gravids.end = sizeof(struct Population);
	pop.gravids.start = 0;

	pop.bmDigestings.end = sizeof(struct Environment);
	pop.bmDigestings.start = 0;

	pop.bmSeekings.end = sizeof(struct Agent);
	pop.bmSeekings.start = 0;
*/

	// copy new population to place where everybody reads it
	gPop[0] = pop;
}
