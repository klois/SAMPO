/*
 * @author      Klaus Kofler
 * @date		07/01/2013
 */
#pragma once

#include "device_types.h"
#include "agent.h"
#include "gpuRand.h"

// defines the properties for agents
#define BITE_HUMAN_PROBABILITY 1.0f
#define BITE_INDOOR_PROBABILITY 1.0f
#define REST_INDOOR_PROBABILITY 1.0f

///////////////////////////////////////////////////////////////////////////// EGG /////////////////////////////////////////////////////////////////////////////

/*
 * eggsIncubation and eggsHatching are summed together to form the egg's delay. An egg develops into a larva when it's age in our exceeds this delay and
 * world time is in the range specified by eggsTransitionTime
 */

/*
 * Incubation time for an egg. Depends on the environment's temperature. Transition of an egg to a larva occurs when the age in hours of an agent exceeds
 * the sum of eggsIncubation and eggsHatching
 * @param temperature The current temperature of the environment
 * @return the eggs incubation time, dependent on the temperature
 */
REAL eggsIncubation(Temperature temperature) {
	return (-0.923f*temperature+60.923f);
}

/*
 * Hatching time for an egg. Depends on a random variable with range [0, 1]. Transition of an egg to a larva occurs when the age in hours of an agent exceeds
 * the sum of eggsIncubation and eggsHatching
 * @param probability A random variable between 0 and 1
 * @return the eggs hatching value dependent on a random number
 */
REAL eggsHatching(REAL probability) {
	REAL m = 0.0f;
	REAL b = 0.0f;
	if (probability < 0.5f){
		m= 48.0f;
		b= 0.0f;
	} else if (probability < 0.85f){
		m= 68.5714f;
		b= -10.2857f;
	} else if (probability < 0.90f){
		m= 480.0f;
		b= -360.0f;
	} else if (probability < 0.94f){
		m= 600.00f;
		b= -468.0f;
	} else {
		m= 2400.0f;
		b= -2160.0f;
	}

	return (UINT)(m*probability+b+0.5f); // round to nearest integer
}

/*
 * determines if the current world time is inside the egg state transition time range. A transition from egg to larva can only occur if this function
 * evaluates to true
 * @param worldTime The current time in the environment, ranging form 0 to 23
 * @return true if the world time is in the egg state transition range, false otherwise
 */
bool eggsTransitionTime(UINT worldTime) {
	return true;
}

//////////////////////////////////////////////////////////////////////////// LARVA ////////////////////////////////////////////////////////////////////////////

/*
 * a larva develops to a pupa if the sum of devRates exceeds the larva delay. The delay is specified using the larvaDelay function. It is temperature
 * independent and based on a random number who's seeds are an argument to the function. In each iteration, a temperature dependent development rate is added
 * to the lara's development. This value is specified by the larvaDevRate function. Once the development exceeds the delay specified by the function larva
 * delay and world time is in the range specified by larvaeTransitionTime, the state transition into a pupa occurs
 */

/*
 * Specifies the larva Delay. Larva will transform into a pupa once the development exceeds the delay.
 * This version uses a normal distributed random number
 * @param seeds Seeds for the OpenCL device random number generator
 * @return The larva delay depending on a normal distribution
 */
REAL larvaDelay(__constant struct Seeds* seeds) {
	REAL mean = 1.0f;
	REAL stdDev = 0.1f;
	return BoxMuller(seeds, mean, stdDev); // resembles Simulation.randomNormalDouble(mean, stdDev)
}

// development rates to be used to calculate the larva development in each time step
__constant REAL devRates[] = { // from 12 to 43 degrees C if greater wrap to these limits
			0.000175044, 0.000549094,0.00109031,0.001494158,0.001752487,0.001962866,0.002171238,0.002392989,
			0.002633689,0.00289619,0.003182685,0.003495287,0.003836204,0.004207789,0.004612566,0.005053243,0.005532729,0.006054144,
			0.006620836,0.007236393,0.007904657,0.008629745,0.009416051,0.010267388,0.011113704,0.011960019,0.012806335,0.010401163,
			0.007995991,0.00033845,0.00000533086,0.0000000840986 };

/*
 * calculates the larva development for a single time step
 * @param temperature The current temperature of the environment
 * @return the temperature dependent development of a larva
 */
REAL larvaDevelopment(Temperature temperature) {
	Temperature clampedTemp = min(max(temperature, 12.0f), 48.0f); // clamp temperature between 12 and 48

	return devRates[((int)clampedTemp) - 12];
}

/*
 * determines if the current world time is inside the larva state transition time range. A transition from larva to pupa can only occur if this function
 * evaluates to true
 * @param worldTime The current time in the environment, ranging form 0 to 23
 * @return true if the world time is in the larva state transition range, false otherwise
 */
bool larvaeTransitionTime(UINT worldTime) {
	return (worldTime <= 6u || worldTime >= 18u); // TODO check
}

//////////////////////////////////////////////////////////////////////////// PUPA /////////////////////////////////////////////////////////////////////////////

/*
 * Pupae transform into immatures once their time in state (in hours) exceeds the delay defined by pupaDelay and world time is in the range specified by
 * pupaeTransitionTime
 */

/*
 * Defines the temperature dependent pupa delay. Transition to immatures happens once the time in state exceeds the delay
 * @param temperature The current temperature of the environment
 * @return the temperature dependent delay of the pupa
 */
REAL pupaDelay(Temperature temperature) {
	return -0.923*temperature+60.923;
}

/*
 * determines if the current world time is inside the pupa state transition time range. A transition from pupa to immature can only occur if this function
 * evaluates to true
 * @param worldTime The current time in the environment, ranging form 0 to 23
 * @return true if the world time is in the pupa state transition range, false otherwise
 */
bool pupaeTransitionTime(UINT worldTime) {
	return (worldTime <= 6u || worldTime >= 18u);
}

///////////////////////////////////////////////////////////////////////// IMMATURES ///////////////////////////////////////////////////////////////////////////

/*
 * Immatures develop into mate seeking agents only if the are female. To develop the time in this state (in hours) must exceed the delay specified by
 * immaturesDelay and the world time must be in the range specified by immaturesTransitionTime
 */

/*
 * Defines the temperature dependent immature delay. Transition to mate seeking agents happens once the time in state exceeds the delay
 * @param temperature The current temperature of the environment
 * @return the temperature dependent delay of the immature
 */
REAL immatureDelay(Temperature temperature) {
	return -2.667f * temperature + 120.0f;
}

/*
 * determines if the current world time is inside the immatures state transition time range. A transition from immature to mate seeking can only occur if this
 * function evaluates to true
 * @param worldTime The current time in the environment, ranging form 0 to 23
 * @return true if the world time is in the immature state transition range, false otherwise
 */
bool immaturesTransitionTime(UINT worldTime) {
	return true;
}

/////////////////////////////////////////////////////////////////////// MATE SEEKING //////////////////////////////////////////////////////////////////////////

/*
 * Mate seeking agents develop into blood meal seeking agents when the time in this state (in hours) exceeds the delay specified by bloodMealSeekingDelay
 * and the world time is in the range specified by mateSeekingsTransitionTime
 */

/*
 * Defines the temperature dependent mate seeking delay. Transition to blood meal seeking agents happens once the time in state exceeds the delay
 * @param temperature The current temperature of the environment
 * @return the temperature dependent delay of the mate seeking agent
 */
REAL mateSeekingDelay(Temperature temperature) {
	return 1.0f;
}

/*
 * determines if the current world time is inside the mate seekings state transition time range. A transition from immature to mate seeking can only occur if
 * this function evaluates to true
 * @param worldTime The current time in the environment, ranging form 0 to 23
 * @return true if the world time is in the mate seekings state transition range, false otherwise
 */
bool mateSeekingsTransitionTime(UINT worldTime) {
	return (worldTime == 18);
}


///////////////////////////////////////////////////////////////////// BLOOD MEAL SEEKING /////////////////////////////////////////////////////////////////////

/*
 * Blood meal seeking agents develop into blood meal digesting agents when the time in this state (in hours) exceeds the delay specified by
 * bloodMealSeekingDelay and the world time is in the range specified by bloodMealSeekingsTransitionTime
 */

/*
 * Defines the temperature dependent blood meal seeking delay. Transition to blood meal digesting agents happens once the time in state exceeds the delay
 * @param temperature The current temperature of the environment
 * @return the temperature dependent delay of the blood meal seeking agent
 */
REAL bloodMealSeekingDelay(Temperature temperature) {
	return 0.0f;
}

/*
 * determines if the current world time is inside the larva state transition time range. A transition from immature to mate seeking can only occur if this
 * function evaluates to true
 * @param worldTime The current time in the environment, ranging form 0 to 23
 * @return true if the world time is in the mate seekings state transition range, false otherwise
 */
bool bloodMealSeekingsTransitionTime(UINT worldTime) {
	return (worldTime <= 6u || worldTime >= 18u);
}

///////////////////////////////////////////////////////////////////// BLOOD MEAL DIGESTING ////////////////////////////////////////////////////////////////////

/*
 * Blood meal digesting agents develop into gravid agents when the time in this state (in hours) exceeds the delay specified by bloodMealDigestingDelay and
 * the world time is in the range specified by bloodMealSeekingsTransitionTime. Durint the state transition, the function generateEggs is called in order to
 * determine how many eggs the gravid agent will hold
 */

/*
 * Defines the temperature dependent blood meal digesting delay. Transition to gravid agents happens once the time in state exceeds the delay
 * @param temperature The current temperature of the environment
 * @return the temperature dependent delay of the blood meal seeking agent
 */
REAL bloodMealDigestingDelay(Temperature temperature) {
	return -1.231f * temperature + 77.231f;
}

/*
 * determines if the current world time is inside the larva state transition time range. A transition from immature to mate seeking can only occur if this
 * function evaluates to true
 * @param worldTime The current time in the environment, ranging form 0 to 23
 * @return true if the world time is in the mate seekings state transition range, false otherwise
 */
bool bloodMealDigestingsTransitionTime(UINT worldTime) {
	return true;
}

/*
 * Calculates the number of eggs that an agent hold when evolving from blood meal digesting to gravid. New eggs are only generated if the agent does not already
 * contain eggs
 */
UINT generateEggs(__constant struct Seeds* seeds, UINT* numEggBatches) {
	REAL eggBatchSize = 170.0f;
	REAL stdDev = 30.0f;
	UINT eggs = (UINT)(BoxMuller(seeds, eggBatchSize, stdDev) + 0.5f); //  (int)Math.round(Simulation.randomNormalDouble(170.0, 30.0));
	if (eggs < 0) eggs = 0;
	eggs = (UINT)(eggs * pow(.8f, (REAL)*numEggBatches) + 0.5f);//(int)Math.round((eggs * pow(.8f, (REAL)batch)));

	++(*numEggBatches);
	return eggs; // Normally distributed now.
}

//////////////////////////////////////////////////////////////////////////// GRAVID ///////////////////////////////////////////////////////////////////////////

/*
 * Gravid agents develop try to lay eggs at a certain time of day, defined by gravidsEggLayTime. The amount of eggs laid is determined by the function layEggs.
 * A gravid agent will transition into a blood meal seeking one after all its eggs are laid.
 */

/*
 * determines if the current world time is inside the egg laying time of the gravid agent
 * @param worldTime The current time in the environment, ranging form 0 to 23
 * @return true if the world time is in the mate seekings state transition range, false otherwise
 */
bool gravidsEggLayTime(UINT worldTime) {
	return (worldTime <= 6u || worldTime >= 18u);
}

/*
 * Function which determines how many eggs are laid by a gravid agent. It increases the global new eggs counter by the amount of laid eggs and returns the
 * number of eggs remaining in the agent.
 * @param agent The gravid agent laying the eggs. Holds the fields availableEggs, cycleLength, bloodmealCount, delay, numEggBatches and ovipositionAttemps.
 * 				The fields ovipositionAttemps and availableEggs will be read
 * @param agentState The state of the current agent. Holds the fields dead and state. The field dead will be updated in case of OVITrap capture
 * @param enb A struture holding global information about the new eggs and total biomass in the environment. It consists of totalBiomass and newEggs. The
 * 				field totalBiomass will be read, newEggs will be updated
 * @param seeds Seeds for the OpenCL device random number generator
 * @param caryingCapacity The current carryingCapacity of the environment
 * @return the eggs remaining in the agent
 */
INT layEggs(__private struct Agent* agent, __private struct AgentState* agentState, __global struct EggsNbiomass* enb, __constant struct Seeds* seeds,
		REAL carryingCapacity) {
	UINT totalBiomass = enb->totalBiomass;
	REAL normalizedBiomass = ((REAL)totalBiomass) / (((REAL)(agent->ovipositionAttemps) + 1.0f) * carryingCapacity);
	REAL potentialBiomass = max(0.0f, 1.0f - normalizedBiomass);
	UINT potentialEggs = (UINT)ceil((REAL)agent->availableEggs * potentialBiomass);

	//add new eggs, single counter for the entire environment, hence an atomic operation must be used
	UINT idx = atomic_add(&(enb->newEggs), potentialEggs);

	return agent->availableEggs - potentialEggs;
}
