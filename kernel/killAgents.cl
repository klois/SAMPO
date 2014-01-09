/*
 * @author      Klaus Kofler
 * @date		07/01/2013
 */

#include "device_types.h"
#include "agent.h"
#include "gpuRand.h"

UINT calcDiff(struct AgentRange range) {
	return range.end - range.start;
}

/*
 * initializes the environment with several eggs of agents
 * @param agents the agents age array to read the age
 * @param agentStates array holding the sate information of all agents
 * @param pop the properties (number of agents in certain state) of the current population
 * @param enb this kernel resets the newEggs counter and calculates the total biomass
 * @param bnc structure to store the informations about bites and cycles. Will be nulled in this kernel
 * @param seeds Seeds to be used for the random number generator on the device
 * @param carryingCapacity the environment's carrying capacity
 */
__kernel void killAgents(__global struct AgentAge* agents, __global struct AgentState* agentStates, __constant struct Population* pop,
		__global struct EggsNbiomass* enb, __global struct BitesNcycles* bnc, __constant struct Seeds* seeds, REAL carryingCapacity) {
	UINT gid = get_global_id(0);

	bool valid = !((gid >= pop->eggs.end && gid < pop->larvae.start) ||
				   (gid >= pop->larvae.end && gid < pop->pupae.start) ||
				   (gid >= pop->pupae.end && gid < pop->immatures.start) ||
				   (gid >= pop->immatures.end && gid < pop->mateSeekings.start) ||
				   (gid >= pop->mateSeekings.end && gid < pop->bmSeekings.start) ||
				   (gid >= pop->bmSeekings.end && gid < pop->bmDigestings.start) ||
				   (gid >= pop->bmDigestings.end && gid < pop->gravids.start) ||
				   (gid >= pop->gravids.end));

	if(gid == 0) {
		// reset newEgg counter for next iteration of update
		enb->newEggs = 0;

		// calculate the total biomass
		UINT totalBiomass = pop->numLarvae1DayEquiv + (pop->eggs.end) + (pop->pupae.end - pop->pupae.start);
		enb->totalBiomass = totalBiomass;

		// null BitesNcycles
		bnc->numBitesReported = 0;
		bnc->numInfectBitesReported = 0;
		bnc->numCyclesReported = 0;
		bnc->sumCyclesReported = 0;
	}

	if(!valid) return;

	struct AgentAge agent = agents[gid];
	UINT ageInDays = (UINT)min(agent.ageInHours / 24.0f, 99.0f);

	REAL probability = HybridTaus(&seeds[0]); // TODO random
/*
if(gid >= pop->eggs.start && gid < pop->eggs.end)
	if(agents[gid].state != EGG) printf("%d Not an egg %d [%d %d]\n", gid, agents[gid].state, pop->eggs.start, pop->eggs.end);

if(gid >= pop->larvae.start && gid < pop->larvae.end)
	if(agents[gid].state != LARVA) printf("%d Not a larva %d [%d %d]\n", gid, agents[gid].state, pop->larvae.start, pop->larvae.end);

if(gid >= pop->pupae.start && gid < pop->pupae.end)
	if(agents[gid].state != PUPA) printf("%d Not a pupa %d [%d %d]\n", gid, agents[gid].state, pop->pupae.start, pop->pupae.end);

if(gid >= pop->immatures.start && gid < pop->immatures.end)
	if(agents[gid].state != IMMATURE) printf("%d Not immatue %d [%d %d]\n", gid, agents[gid].state, pop->immatures.start, pop->immatures.end);

if(gid >= pop->mateSeekings.start && gid < pop->mateSeekings.end)
	if(agents[gid].state != MATESEEKING) printf("%d Not mate seeking %d [%d %d]\n", gid, agents[gid].state, pop->mateSeekings.start, pop->mateSeekings.end);

if(gid >= pop->bmSeekings.start && gid < pop->bmSeekings.end)
	if(agents[gid].state != BMS) printf("%d Not blood meal seeking %d [%d %d]\n", gid, agents[gid].state, pop->bmSeekings.start, pop->bmSeekings.end);

if(gid >= pop->bmDigestings.start && gid < pop->bmDigestings.end)
	if(agents[gid].state != BMD) printf("%d Not blood meal digesting %d [%d %d]\n", gid, agents[gid].state, pop->bmDigestings.start, pop->bmDigestings.end);

if(gid >= pop->gravids.start && gid < pop->gravids.end)
	if(agents[gid].state != GRAVID) printf("%d Not gravid %d [%d %d]\n", gid, agents[gid].state, pop->gravids.start, pop->gravids.end);
*/

	// killing adults
	if(gid >= pop->immatures.start) {
		REAL a, B, s;
		a = 0.1f;
		B = 25.0f;
		s = 0.1f;
		REAL DMR = (a * exp(ageInDays/B)) / (1.0f + (a * B * s * (exp(ageInDays/B) - 1.0f)) );
		REAL DSR = 1.0f - DMR;	// DSR = Daily Survival Rate
		REAL HSR = pow(DSR, 1.0f/24.0f);	// HSR = Hourly Survival Rate
		REAL HMR = 1.0f - HSR;	// HMR = Hourly Mortality Rate

		if (probability <= HMR) {	// Kill this agent
			agentStates[gid].dead = true;
		}

		//agents[gid] = agent;
		return;
	}

	if(gid < pop->eggs.end) {
		REAL DMR = 0.1f;	// DMR = Daily Mortality Rate
		REAL DSR = 1.0f - DMR;	// DSR = Daily Survival Rate
		REAL HSR = pow(DSR, 1.0f/24.0f);	// HSR = Hourly Survival Rate
		REAL HMR = 1.0f - HSR;	// HMR = Hourly Mortality Rate

		if (probability <= HMR) {	// Kill this agent
			agentStates[gid].dead = true;
		}

		//agents[gid] = agent;
		return;
	}

	if(gid < pop->larvae.end) {
		REAL rainfallCoefficient = 1.0f;
		REAL a = 0.1f;
		REAL DMR = a * exp( (REAL)pop->numLarvae1DayEquiv / (ageInDays * carryingCapacity * rainfallCoefficient) );
		DMR = min(DMR, 0.8f);	// Clip if larger than 1.0

		REAL DSR = 1.0f - DMR;	// DSR = Daily Survival Rate
		REAL HSR = pow(DSR, 1.0f/24.0f);	// HSR = Hourly Survival Rate
		REAL HMR = 1.0f - HSR;	// HMR = Hourly Mortality Rate

		if (probability <= HMR) {	// Kill this agent
			agentStates[gid].dead = true;
		}

		//agents[gid] = agent;
		return;
	}

	// pupae state
	REAL DMR = 0.1f;	// DMR = Daily Mortality Rate
	REAL DSR = 1.0f - DMR;	// DSR = Daily Survival Rate
	REAL HSR = pow(DSR, 1.0f/24.0f);	// HSR = Hourly Survival Rate
	REAL HMR = 1.0f - HSR;	// HMR = Hourly Mortality Rate

	if (probability <= HMR) {	// Kill this agent
		agentStates[gid].dead = true;
	}

	//agents[gid] = agent;

	return;
}
