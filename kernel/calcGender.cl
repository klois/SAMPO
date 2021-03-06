/*
 * @author      Klaus Kofler
 * @date		07/01/2013
 */

#include "device_types.h"
#include "agent.h"

/*
 * Calculates partial sums of the number of females in agents (considering only adult agents) as well as the number of potentially infected females
 * One result per group is generated and stored in intermediateFemales
 * @param agents The array holding the needed infromation of all agents
 * @param pop A single structure holding the needed information about the agent populations
 * @param intermediateFemales Array will be filled with partial sums for each group. In the range [0, #groups] the number of females will be stored,
 * 							  in the range [#groups, 2*#groups] each group stores the partial sum of the potential infected females
 */
__kernel void calcFemalesPerGroup(__global struct AgentAge* agents, __constant struct Population* pop, __global UINT* intermediateFemales) {
	UINT gid = get_global_id(0);
	UINT lid = get_local_id(0);
	UINT offset = pop->immatures.start;
	UINT group = get_group_id(0);

	__local UINT females[LOCAL_SIZE];
	__local UINT potentiallyIfective[LOCAL_SIZE];
	females[lid] = 0; 				//initialize with 0, not all platforms do this by default
	potentiallyIfective[lid] = 0;

	for(UINT i = offset + gid; i < pop->gravids.end; i += get_global_size(0)) {
		// do not consider padding area between states
		bool count = !((i >= pop->immatures.end && i < pop->mateSeekings.start) ||
					   (i >= pop->mateSeekings.end && i < pop->bmSeekings.start) ||
					   (i >= pop->bmSeekings.end && i < pop->bmDigestings.start) ||
					   (i >= pop->bmDigestings.end && i < pop->gravids.start));

		if(count) {
			if(agents[i].isFemale) {
				++females[lid];
				potentiallyIfective[lid] += agents[i].cumulativeSporogonicDevelopment >= 1.0f;
			}
		}
	}

	// reduction in local memory
	for(UINT currsize = LOCAL_SIZE/2; currsize > 0; currsize /= 2)
	{
		barrier(CLK_LOCAL_MEM_FENCE);
		if(lid < currsize) {
			females[lid] += females[lid + currsize];
			potentiallyIfective[lid] += potentiallyIfective[lid + currsize];
		}

	}

	if(lid != 0) return;

	intermediateFemales[group] = females[0];
	intermediateFemales[group + get_num_groups(0)] = potentiallyIfective[0];
}

/*
 * Calculates the tolal sums of the number of females in agents (considering only adult agents) as well as the number of potentially infected females using the
 * patial sum generated by calcFemalesPerGroup
 * @param pop A single structure holding the needed information about the agent populations
 * @param intermediateFemales Array holding the partial sums as generated by calcFemalesPerGroup.
 */
__kernel void calcFemalesTotal(__global struct Population* pop, __global UINT* intermediateFemales) {
	// kernel is inteded to be executed by only one group in order to get the final result
	UINT id = get_local_id(0);

	__local UINT females[LOCAL_SIZE];
	__local UINT potentiallyInfective[LOCAL_SIZE];
	females[id] = intermediateFemales[id];
	potentiallyInfective[id] = intermediateFemales[id + LOCAL_SIZE];

	// reduction in local memory
	for(UINT currsize = LOCAL_SIZE/2; currsize > 0; currsize /= 2)
	{
		barrier(CLK_LOCAL_MEM_FENCE);
		if(id < currsize) {
			females[id] += females[id + currsize];
			potentiallyInfective[id] += potentiallyInfective[id + currsize];
		}
	}

	if(id != 0) return;

	pop->numFemale = females[0];
	pop->numPotentiallyInfective = potentiallyInfective[0];
}
