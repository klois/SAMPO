#include "device_types.h"
#include "agent.h"

// optimizations using local memory in a two step approach abandoned since it is called only one time after the simulation

/*
 * Calculates a histogram of the age of all adult mosquitoes. One bin per day, cropped at 99 days.
 * @param agents Array holding all agents' age
 * @param pop the properties (number of agents in certain state) of the current population. The correct information is found at position 1 of this buffer since
 * 				it is calculated after the oldToNewAgents and before the calcL1de kernel
 * @param histogram buffer where the histogram will be stored
 */
__kernel void ageHistogram(__global struct AgentAge* agents, __constant struct Population* pop, __global UINT* histogram) {
	UINT gid = get_global_id(0);

	UINT offset = pop[1].immatures.start;
	UINT group = get_group_id(0);

	if(gid < offset || gid >= pop[1].gravids.end)
		return;

	bool count = !((gid >= pop[1].immatures.end && gid < pop[1].mateSeekings.start) ||
				   (gid >= pop[1].mateSeekings.end && gid < pop[1].bmSeekings.start) ||
				   (gid >= pop[1].bmSeekings.end && gid < pop[1].bmDigestings.start) ||
				   (gid >= pop[1].bmDigestings.end && gid < pop[1].gravids.start));

	if(!count) return;

	struct AgentAge agent = agents[gid ];

	UINT bin = min((UINT)(agent.ageInHours/24.0f), 99u);

	atomic_inc(&(histogram[bin]));

}
