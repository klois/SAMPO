/*
 * @author      Klaus Kofler
 * @date		07/01/2013
 */

#define STR_VALUE(arg) #arg
#define INC(name) STR_VALUE(name)
#define SPECIES_INC INC(SPECIES)
#include SPECIES_INC
#include "device_types.h"
#include "agent.h"
#include "gpuRand.h"

// default values for probabilities if not defined in agents specific header file
#ifndef BITE_HUMAN_PROBABILITY
#define BITE_HUMAN_PROBABILITY 1.0f
#endif
#ifndef BITE_INDOOR_PROBABILITY
#define BITE_INDOOR_PROBABILITY 1.0f
#endif
#ifndef REST_INDOOR_PROBABILITY
#define REST_INDOOR_PROBABILITY 1.0f
#endif

/*
 * initializes the environment with several eggs of agents
 * @param agents the agents that will be updated
 * @param agentAges age information of the agents to update
 * @param agentStates state information of the agents to update
 * @param pop the properties (number of agents in certain state) of the current population
 * @param environment the environment, used to get the carrying capacity and properties of interventions
 * @param temperature the temperature in this timestep as a floating point number
 * @param bnc single struct that will be filled with the information about the bites and performed cycles using atomic operations
 * @param enb structure to store the number of laid eggs and to read the total biomass
 * @param seeds seeds for the random number generator
 * @param elapsedTimeInHours the time in hours since the last update
 * @param worldTime the current world time ranging from 0 to 23
 */
__kernel void updateAgents(__global struct Agent* agents, __global struct AgentAge* agentAges, __global struct AgentState* agentStates,
		__constant struct Population* pop, struct Environment environment, Temperature temperature, __global struct BitesNcycles* bnc,
		 __global struct EggsNbiomass* enb, __constant struct Seeds* seeds,	REAL elapsedTimeInHours, UINT worldTime) {

	UINT gid = get_global_id(0);

	bool valid = !((gid >= pop->eggs.end && gid < pop->larvae.start) ||
				   (gid >= pop->larvae.end && gid < pop->pupae.start) ||
				   (gid >= pop->pupae.end && gid < pop->immatures.start) ||
				   (gid >= pop->immatures.end && gid < pop->mateSeekings.start) ||
				   (gid >= pop->mateSeekings.end && gid < pop->bmSeekings.start) ||
				   (gid >= pop->bmSeekings.end && gid < pop->bmDigestings.start) ||
				   (gid >= pop->bmDigestings.end && gid < pop->gravids.start) ||
				   (gid >= pop->gravids.end));

	if(!valid) return;

	struct AgentState agentState = agentStates[gid];
	if(agentState.dead) return;

	struct Agent agent = agents[gid];
	struct AgentAge agentAge = agentAges[gid];

	agentAge.ageInHours += elapsedTimeInHours;
	agentAge.hoursInState += elapsedTimeInHours;

	//      if(environment.temperature > 16.0f && agent.bloodmealCount > 0) {
	//              agent.cumulativeSporogonicDevelopment += 1.0f/((111.0f/(environment.temperature-16.0f))*14.0);
	//      }

	agentAge.cumulativeSporogonicDevelopment += (temperature > 16.0f && agent.humanBloodmealCount > 0u) *
			(1.0f/((111.0f/(temperature-16.0f))*14.0)); // TODO how would biteHuman interact here

	if(gid < pop->eggs.end) {
		if(agentAge.hoursInState >= agent.delay && eggsTransitionTime(worldTime)) {
			// larva state enter
			agentAge.ageInHours = 24.0f;
			agentAge.hoursInState = 0.0f;
			// agent.timeEntered = environment.totalTime; not used
			// agent.cumulativeLarvalDleay = 0.0f; already set at egg creation
			agent.delay = larvaDelay(seeds+1);

			agentState.state = LARVA;
			agents[gid] = agent;
			agentStates[gid] = agentState;
		}
		agentAges[gid] = agentAge;
		return;
	}

	if(gid < pop->larvae.end) {
		if (environment.larvacideValue > 0.0f && agentAge.hoursInState == elapsedTimeInHours) { // TODO check if first time to update the larva only
			REAL larvacide = HybridTaus(&seeds[2]);
			if (larvacide <= environment.larvacideValue) {
				agentState.dead = true;	// Mark this agent as 'dead'
				agentStates[gid] = agentState;
				return;
			}
		}

		agent.cumulativeLarvalDelay += larvaDevelopment(temperature);

		if(agent.cumulativeLarvalDelay >= agent.delay && larvaeTransitionTime(worldTime)) {
			// pupa state enter
			agentAge.hoursInState = 0.0f;
			// agent.timeEntered = environment.totalTime; not used
			agent.delay = pupaDelay(temperature); //temp depended rate for pupa...

			agentState.state = PUPA;
			agentStates[gid] = agentState;
		}
		agents[gid] = agent;
		agentAges[gid] = agentAge;
		return;
	}

	if(gid < pop->pupae.end) {
		if((agentAge.hoursInState >= agent.delay) && pupaeTransitionTime(worldTime)) {
			// immature state enter
			agentAge.ageInHours = 0.0f;
			agentAge.hoursInState = 0.0f;
			// agent.timeEntered = environment.totalTime; not used
			agent.delay = immatureDelay(temperature); // temperature depended rate for ImmatureAdult...

			agentState.state = IMMATURE;
			agents[gid] = agent;
			agentStates[gid] = agentState;
		}
		agentAges[gid] = agentAge;
		return;
	}

	if(gid < pop->immatures.end) {
		if(agentAge.hoursInState >= agent.delay && immaturesTransitionTime(worldTime)) {
			agent.delay = mateSeekingDelay(temperature);
			// mate seeking state enter
			agentAge.hoursInState = 0.0f; // actually not done in java program and probably not used. But just to be save...

			agentState.state = MATESEEKING;
			agents[gid] = agent;
			agentStates[gid] = agentState;
		}
		agentAges[gid] = agentAge;
		return;
	}

	if(gid < pop->mateSeekings.end) {
		if((agentAge.hoursInState >= agent.delay) && agentAge.isFemale && mateSeekingsTransitionTime(worldTime)) {
			// blood meal seeking state enter
			agentAge.hoursInState = 0.0f; // actually not done in java program and probably not used. But just to be save...
			agent.delay = bloodMealSeekingDelay(temperature);

			agentState.state = BMS;
			agentStates[gid] = agentState;
			agents[gid] = agent;
		}
		agentAges[gid] = agentAge;
		return;
	}

	if(gid < pop->bmSeekings.end) {
		++agent.cycleLength;

		if((agentAge.hoursInState >= agent.delay) && bloodMealSeekingsTransitionTime(worldTime)) {
			REAL bloodmealSuccessProbability = HybridTaus(&seeds[1]); //TODO random
			if(bloodmealSuccessProbability <= environment.bloodmealSuccess) {

				if(environment.ITNValue * BITE_HUMAN_PROBABILITY * BITE_INDOOR_PROBABILITY > 0.0f) {
					REAL ITN = HybridTaus(&seeds[2]); //TODO random
					if(ITN <= environment.ITNValue * BITE_HUMAN_PROBABILITY * BITE_INDOOR_PROBABILITY) {
						agentState.dead = true;
						agentStates[gid] = agentState;
						return;
					}
				}

				int bitesAhuman = BITE_HUMAN_PROBABILITY < 1.0f ? (HybridTaus(&seeds[3]) <= BITE_HUMAN_PROBABILITY) : 1;

				agent.humanBloodmealCount += 1u * bitesAhuman;

				atomic_inc(&bnc->numBitesReported);

				// blood meal digesting state enter
				agentAge.hoursInState = 0.0f;
				// agent.timeEntered = environment.totalTime; not used
				agent.delay = bloodMealDigestingDelay(temperature);
				agentState.state = BMD;

				if((agentAge.cumulativeSporogonicDevelopment >= 1.0f) && bitesAhuman)
					atomic_inc(&(bnc->numInfectBitesReported));

				agentStates[gid] = agentState;
			}
		}

		agents[gid] = agent;
		agentAges[gid] = agentAge;
		return;
	}

	if(gid < pop->bmDigestings.end) {
		if(environment.IRSValue * REST_INDOOR_PROBABILITY > 0.0f) {
			REAL ITN = HybridTaus(&seeds[2]); //TODO random
			if(ITN <= environment.IRSValue * REST_INDOOR_PROBABILITY) {
				agentState.dead = true;
				agentStates[gid] = agentState;
				return;
			}
		}

		++agent.cycleLength;

		if(agentAge.hoursInState >= agent.delay && bloodMealDigestingsTransitionTime(worldTime)) {
			if(agent.availableEggs <= 0 && agentAge.isFemale) { // generate eggs
				agent.availableEggs = generateEggs(seeds+1, &agent.numEggBatches);
			}

			// gravid state enter
			agentAge.hoursInState = 0.0f;
			// agent.timeEntered = environment.totalTime; not used
			agent.ovipositionAttemps = 0u;

			agentState.state = GRAVID;
			agentStates[gid] = agentState;
		}

		agents[gid] = agent;
		agentAges[gid] = agentAge;
		return;
	}

	// gravid state
	++agent.cycleLength;
	if(gravidsEggLayTime(worldTime)) {
		bool shouldLayEggs = HybridTaus(&seeds[1]) < 0.25f; // TODO random

		if(shouldLayEggs) {
			if(environment.oviTrapValue > 0.0f) {
				REAL OVITrap = HybridTaus(&seeds[2]); // TODO random
				if(OVITrap <= environment.oviTrapValue) {
					agentState.dead = true;
					agent.availableEggs = 0;
				}
			}

			if(!agentState.dead) {
				//lay eggs
				agent.availableEggs = layEggs(&agent, &agentState, enb, seeds, environment.carryingCapacity);
			}

			++agent.ovipositionAttemps;
			//++agent.podsSampled; never read

			if(agent.availableEggs == 0) {
				// cycle occured
				atomic_inc(&(bnc->numCyclesReported));
				atomic_add(&(bnc->sumCyclesReported), agent.cycleLength);

				agent.cycleLength = 0;

				// blood meal seeking state enter
				agentAge.hoursInState = 0.0f; // actually not done in java program and probably not used. But just to be save...
				agent.delay = bloodMealSeekingDelay(temperature);
				agentState.state = BMS;
				agentStates[gid] = agentState;
			}
		}
	}
	agents[gid] = agent;
	agentAges[gid] = agentAge;

}
