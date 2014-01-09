/*
 * @author      Klaus Kofler
 * @date		07/01/2013
 */
#pragma once

#define DEVICE_TYPE ICL_GPU

#ifdef __APPLE__
#if DEVICE_TYPE == ICL_CPU
#define LOCAL_SIZE 1
#else
#define LOCAL_SIZE 64
#endif
#else
#define LOCAL_SIZE 64
#endif


typedef REAL Temperature;

struct Environment {
//	Temperature* temperatureProfile; separate array
	REAL bloodmealSuccess;

	REAL ITNValue;
//	REAL oviValue;  		unused
//	REAL immValue;	unused
	REAL IRSValue;

	// aquatic fields
	REAL carryingCapacity;
	REAL larvacideValue;
	REAL oviTrapValue;
};

struct EggsNbiomass {
	INT totalBiomass;
	UINT newEggs;
};

// fields moved over from environment since they have to be read in each iteration;
struct BitesNcycles {
	UINT numCyclesReported;
	UINT sumCyclesReported;

	UINT numBitesReported;
	UINT numInfectBitesReported;
};

typedef UINT EnvironmentId;

enum State {
	EGG = 1,
	LARVA = 2,
	PUPA = 4,
	IMMATURE = 8,
	MATESEEKING = 16,
	BMS = 32,
	BMD = 64,
	GRAVID = 128
};

#define MAX_OVIPOSITION_ATTEMPTS 3

struct Agent {
//	struct EnvironmentId environment; // can be dropped if there is only one environment
//	REAL ageInHours; // might better to use INT
//	bool dead;		 // moved to AgentState

	UINT availableEggs;
	UINT cycleLength;

	// AnophelesGambiaeAgent fields
	UINT humanBloodmealCount;

	// state fields
//	enum State state; // moved to AgentState
//	REAL hoursInState; moved to AgentAge
//	REAL timeEnteredInState; seems not to be needed
	REAL delay;

	UINT numEggBatches; // only needed in BMD, has to survive BMS and GRAVID
	union {
		UINT ovipositionAttemps; // gravid state field
		REAL cumulativeLarvalDelay; // larva state field
	};
};

// removed from agent struct for performance reasons
struct AgentState {
	bool dead;
	enum State state;
};

struct AgentAge {
	REAL ageInHours; // might better to use INT
	REAL hoursInState;

	bool isFemale;
	REAL cumulativeSporogonicDevelopment;
};

struct AgentRange {
	UINT start;
	UINT end;
};

struct Population {
	struct AgentRange eggs;
	struct AgentRange larvae;
	struct AgentRange pupae;
	struct AgentRange immatures;
	struct AgentRange mateSeekings;
	struct AgentRange bmSeekings;
	struct AgentRange bmDigestings;
	struct AgentRange gravids;

	// fields needed for stats
	UINT numLarvae1DayEquiv;
	UINT numFemale;
	UINT numPotentiallyInfective;
};

struct Stats {
	// aquatic
	UINT numEggs;
	UINT numLarvae;
	UINT numPupae;
	UINT numLarvae1DazEquiv;
	UINT numBiomass;

	// adult
	UINT numImmature;
	UINT numMating;
	UINT numBMS;
	UINT numBMD;
	UINT numOVI;
	UINT numFemales;
	UINT numPotentiallyInfective;
	UINT numMales;
};

//store random seeds generated on cpu to generate ranodm numbers on GPU
struct Seeds {
	UINT z1;
	UINT z2;
	UINT z3;
	UINT z4;
};
