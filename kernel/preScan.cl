/*
 * @author      Klaus Kofler
 * @date		07/01/2013
 */

#include "device_types.h"
#include "agent.h"

/*
 * generates a flag array with 1 at all places where the corresponding agent is alive and its state matches one of the states in match\
 * @param agentStates the agents' state array
 * @param pop the properties (number of agents in certain state) of the current population
 * @param data a array that will be filled with 0 and 1, depending on the value of agents.state and match
 * @param flags an array that will be 0 in the valid areas and 1 in the others
 * @param start the index in agentStates from where the prefix sum should start counting
 * @param break1 the end of the 1st area
 * @param break2 the end of the 2nd area
 * @param end the number of elements in agents/flags that are processed, at the same time, the end of the last area
 * @param match an integer that has a 1 at each bit that corresponds to a state that should be matched
 */
void __kernel preScan(__global struct AgentState* agentStates, __constant struct Population* pop, __global INT* data, __global UINT* flags,
		UINT start, UINT break1, UINT break2, UINT end, INT match) {
	UINT gid = get_global_id(0);
	
	// increment bounds by one because this field is used to read the number of elements in the corresponding bin
	++end;
	++break1;
	++break2;
//if(gid == 0) printf("device b1 %d, b3 %d, end %d\n", break1, break2, end);
	if(gid > end) return;
	if(gid < start) return;
/*
	if(gid == end) {
		flags[gid] = 1; // TODO check if this is really needed. In this way, 2 additional elements are accessed
//		return;
	}
*/
//	flags[gid] = (gid == break1) || (gid == break2) || (gid == end);

	bool valid = !((gid >= pop->eggs.end && gid < pop->larvae.start) ||
				   (gid >= pop->larvae.end && gid < pop->pupae.start) ||
				   (gid >= pop->pupae.end && gid < pop->immatures.start) ||
				   (gid >= pop->immatures.end && gid < pop->mateSeekings.start) ||
				   (gid >= pop->mateSeekings.end && gid < pop->bmSeekings.start) ||
				   (gid >= pop->bmSeekings.end && gid < pop->bmDigestings.start) ||
				   (gid >= pop->bmDigestings.end && gid < pop->gravids.start) ||
				   (gid >= pop->gravids.end));


	UINT idx = gid - start;
	flags[idx] = (gid >= break1) + (gid >= break2) + 1;
	if(!valid) {
		data[idx] = 0;
		return;
	}

	struct AgentState agentState = agentStates[gid];
	if(!agentState.dead && ((agentState.state & match) != 0))
		data[idx] = 1;
	else
		data[idx] = 0;

//	printf("%d gid %d, match %d, state %d\n", (agents[gid].state & match) != 0, gid, match, agents[gid].state);
}
