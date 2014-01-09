/*
 * @author      Klaus Kofler
 * @date		07/01/2013
 */

// including species header as defined via command line argument
#define STR_VALUE(arg) #arg
#define INC(name) STR_VALUE(name)
#define SPECIES_INC INC(SPECIES)
#include SPECIES_INC


#include "device_types.h"
#include "agent.h"
#include "gpuRand.h"

/*
 * initializes the environment with several eggs of agents
 * @param agents the agents array where the new agents will be written to
 * @param agentAges the agentAges array where the new agents' age information will be written to
 * @param agentStates the agentStates array where the new agents' state information will be written to
 * @param enb EggsNbiomass structure containing also the number of agent eggs that will be put in agents in this step
 * @param temperature the temperature in this timestep as a floating point number
 * @param seeds seeds for the random number generator
 */
__kernel void initAgents(__global struct Agent* agents, __global struct AgentAge* agentAges, __global struct AgentState* agentStates,
		__constant struct EggsNbiomass* enb, Temperature temperature, __constant struct Seeds* seeds) {
	UINT initialAgentCount = enb->newEggs;
	UINT id = get_global_id(0);

	if(id >= initialAgentCount)
		return;

	//initialize the initial population
	struct AgentState as;
	as.dead = false;
	as.state = EGG;
	agentStates[id] = as;

	struct AgentAge agentAge;
	agentAge.ageInHours = 0.0f;
	agentAge.hoursInState = 0.0f;

	agentAge.isFemale = HybridTaus(&seeds[3]) <= 0.5; //TODO random;
	agentAges[id] = agentAge;

	struct Agent agent;
	agent.humanBloodmealCount = 0u;

	agent.availableEggs = 0u;
	agent.cycleLength = 0u;
	agent.cumulativeLarvalDelay = 0.0f; // set field needed in larva state already now

	REAL probability = HybridTaus(&seeds[2]); //TODO random defines "bin"

	agent.delay = eggsIncubation(temperature) + eggsHatching(probability); //set time for incubation + hatching

	// MosquitoAgent constructor
//	agent.availableEggs = 0; written before read
	agent.numEggBatches = 0u;

	// store generated agent in global memory
	agents[id] = agent;
}

/*
 * only to initialize a testing environment!!!
 * initializes the environment with several agents of all types
 * @param agents the agents array where the new agents will be written to
 * @param initialAgentcount the number of agents that will be put in agents in this step
 */
/*
__kernel void initTestAgents(__global struct Agent* agents, __global struct Environment* environment) {
	UINT initialAgentCount = environment->newEggs;
	UINT id = get_global_id(0);
	if(id >= initialAgentCount)
		return;

	//initialize the initial population
	struct Agent agent;
	agent.dead = false;
	agent.cycleLength = 0u;

	agent.hoursInState = 0.0f;

	if(id<64) {
		agent.ageInHours = 0.0f;
		agent.availableEggs = 0u;
		agent.bloodmealCount = 0.0f;
		agent.isFemale = (id % 2 == 0) ? true : false;

		agent.state = EGG;
		agent.delay = 24.0f; // TODO add randomness
		agents[id] = agent;
		return;
	}

	if(id < 128) {
		agent.ageInHours = 24.0f;
		agent.availableEggs = 0u;
		agent.cumulativeLarvalDelay = 0.0f;
		agent.bloodmealCount = 0u;
		agent.isFemale = (id % 2 == 0) ? true : false;

		agent.state = LARVA;
		agent.delay = 1.0f; // TODO add randomness
		agents[id] = agent;
		return;
	}

	if(id < 192) {
		agent.availableEggs = 0u;
		agent.bloodmealCount = 0.0f;
		agent.isFemale = (id % 2 == 0) ? true : false;

		agent.state = PUPA;
		agent.delay = 24.0f; // TODO add randomness
		agents[id] = agent;
		return;
	}

	if(id < 256) {
		agent.availableEggs = 0u;

		agent.cumulativeSporogonicDevelopment = 0.0f;
		agent.bloodmealCount = 0.0f;
		agent.isFemale = (id % 2 == 0) ? true : false;

		agent.state = IMMATURE;
		agent.delay = 24.0f; // TODO add randomness
		agents[id] = agent;
		return;
	}

	if(id < 320) {
		agent.availableEggs = 0u;

		agent.cumulativeSporogonicDevelopment = 0.0f;
		agent.bloodmealCount = 0.0f;
		agent.isFemale = (id % 2 == 0) ? true : false;

		agent.state = MATESEEKING;
		agent.delay = 1.0f;
		agents[id] = agent;
		return;
	}

	if(id < 384) {
		agent.availableEggs = 0u;

		agent.cumulativeSporogonicDevelopment = 0.0f;
		agent.bloodmealCount = 0.0f;
		agent.isFemale = (id % 2 == 0) ? true : false;

		agent.state = BMS;
		agent.delay = 1.0f;
		agents[id] = agent;
		return;
	}

	if(id < 448) {
		agent.availableEggs = 0u;

		agent.cumulativeSporogonicDevelopment = 0.0f;
		agent.bloodmealCount = 0.0f;
		agent.isFemale = (id % 2 == 0) ? true : false;

		agent.numEggBatches = 0u;

		agent.state = BMD;
		agent.delay = 48.0f; // TODO add randomness
		agents[id] = agent;
		return;
	}

//	if(id < 448) {
		agent.availableEggs = 0u;

		agent.bloodmealCount = 0.0f;
		agent.isFemale = (id % 2 == 0) ? true : false;
		agent.cumulativeSporogonicDevelopment = (agent.isFemale && id % 3 == 0) ? 1.0f : 0.0f;

		agent.ovipositionAttemps = 0u;

		agent.state = GRAVID;
		agents[id] = agent;
		return;
//	}

	// store generated agent in global memory
}
*/
