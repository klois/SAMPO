/*
 * @author      Klaus Kofler
 * @date		07/01/2013
 */

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#include "abms.h"
#include "scan.h"

#define FACTOR 320 

REAL hoursInTimeStep = 1.0f;
long long maxSteps = 8760;
struct Environment* environment;
UINT initialAgentCount = 100 * FACTOR;
REAL carryingCapacity = 500.0f * (REAL)FACTOR;
Temperature* temperature;

UINT nSeeds = 4u; // number of seeds for random number generator

#define DEBUG 1

UINT buffSize = max(max(2*LOCAL_SIZE, 100), maxAgents);

//TODO: Fix this for MacOS X. MacOS X seemingly requires this to be an absolute path (relative path doesn't seem to work)
#define KENRNEL_INCLUDE_PATH "include"
char kernelBuildArgs[512];

// timing
icl_event* initEvent;
icl_event* killEvent;
icl_event* updateEvent;
icl_event* preScanEvent;
icl_event* oldToNewEvent;

#if TIMING
double initTime;
icl_timer* l1deTimer;
icl_timer* genderTimer;
icl_timer* randTimer;
double killTime;
double updateTime;
double preScanTime;
icl_timer* scanTimer;
double oldToNewTime;
#endif

// swap buffers
void swap(icl_buffer** a, icl_buffer** b)
{
	icl_buffer* tmp = *a;
	*a = *b;
	*b = tmp;
}

UINT calcNum(struct AgentRange range) {
	return range.end - range.start;
}

void printRange(struct AgentRange* range) {
	printf("[%d\t%d]\t%d\r\n", range->start, range->end, range->end - range->start);
}

void printStats(struct Stats* stats, struct BitesNcycles* bnc, UINT i) {
	UINT numMature = stats->numImmature + stats->numMating + stats->numBMS + stats->numBMD + stats->numOVI;

	printf("%0.1f, %0.1f", hoursInTimeStep * i, temperature[i]);

	printf(", %d, %d, %d, %d, %d, %d, %d, %d, %d", numMature, stats->numImmature, stats->numMating, stats->numBMS, stats->numBMD, stats->numOVI,
			stats->numFemales, stats->numPotentiallyInfective, stats->numMales);
	printf(", %d, %d, %d, %d, %d, ", stats->numEggs, stats->numLarvae, stats->numPupae, stats->numLarvae1DazEquiv, stats->numBiomass);

//	printf(", %d\r\n", stats->numEggs + stats->numLarvae + stats->numPupae + stats->numImmature + stats->numMating + stats->numBMS + stats->numBMD + stats->numOVI);
	printf("%d, %d, %d\r\n", bnc->numCyclesReported > 0 ? bnc->sumCyclesReported / bnc->numCyclesReported : 0,
			bnc->numBitesReported, bnc->numInfectBitesReported);
	// 290, 58, 54, 60, 56, 62, 145, 10, 145, 60, 56, 62, 56, 178

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

int readSingleValFile(const char* path, REAL** pointer) {
	*pointer = (REAL*)malloc(sizeof(REAL) * maxSteps);
	FILE* temp = fopen(path, "r");

	if(!temp) {
		printf("Cannot open %s. Check temperature file path\n", path);
		return -1;
	}

	for(UINT i = 0u; i < maxSteps; ++i)
		fscanf(temp, "%f", &(*pointer)[i]);

	fclose(temp);
	return 0;
}

int readEnvironment(const char* path) {
	environment = (struct Environment*)malloc(sizeof(struct Environment) * maxSteps);
	FILE* interventions = fopen(path, "r");

	if(!interventions) {
		printf("Cannot open %s. Check environment file path\n", path);
		return -1;
	}

	REAL coverage, effectiveness;

	// TODO check how arguments should be used
	for(UINT i = 0u; i < maxSteps; ++i) {
		fscanf(interventions, "%f, %f, ", &coverage, &effectiveness);
		environment[i].IRSValue = coverage * effectiveness;
		fscanf(interventions, "%f, %f, ", &coverage, &effectiveness);
		environment[i].ITNValue = coverage * effectiveness;
		fscanf(interventions, "%f, %f, ", &coverage, &effectiveness);
		environment[i].larvacideValue = coverage * effectiveness;
		fscanf(interventions, "%f, %f, ", &coverage, &effectiveness);
		environment[i].oviTrapValue = coverage * effectiveness;
		fscanf(interventions, "%f", &(environment[i].bloodmealSuccess));
		environment[i].carryingCapacity = carryingCapacity; // TODO change this to reading from file, once an appropriate file is available
	}
	fclose(interventions);

	return 0;
}


void generateSeeds(icl_buffer* seedsD, struct Seeds* seedsH, icl_device* dev) {
#if TIMING
	clFinish(dev->queue);
	icl_start_timer(randTimer);
#endif
	for(UINT i = 0; i < nSeeds; ++i) {
		seedsH[i].z1 = rand();
		seedsH[i].z2 = rand();
		seedsH[i].z3 = rand();
		seedsH[i].z4 = rand();
	}

	icl_write_buffer(seedsD, CL_TRUE, nSeeds * sizeof(struct Seeds), seedsH, NULL, NULL);

#if TIMING
	clFinish(dev->queue);
	icl_stop_timer(randTimer);
#endif
}

void storeHistogram(icl_buffer* agentAges, icl_buffer* pop, icl_buffer* hist, icl_device* dev) {
	size_t localWorkSize = LOCAL_SIZE;
	size_t globalWorkSize = ((maxAgents + LOCAL_SIZE) / LOCAL_SIZE) * LOCAL_SIZE; // one additional thread to write properties in oldToNew agents
	icl_kernel* ageHist = icl_create_kernel(dev, "kernel/ageHistogram.cl", "ageHistogram", kernelBuildArgs, ICL_SOURCE);

	UINT* histogram = (UINT*)calloc(100, sizeof(UINT));
	// set hist buffer to zero
	icl_write_buffer(hist, CL_TRUE, 100*sizeof(UINT), histogram, NULL, NULL);

	icl_run_kernel(ageHist, 1, &globalWorkSize, &localWorkSize, NULL, NULL, 3,
			(size_t)0, agentAges,
			(size_t)0, pop,
			(size_t)0, hist);

	icl_read_buffer(hist, CL_TRUE, 100*sizeof(UINT), histogram, NULL, NULL);

	FILE* hFile = fopen("histogram.txt", "w");
	for(UINT i = 0u; i < 100u; ++i)
		fprintf(hFile, "%d, ", histogram[i]);

	fclose(hFile);
	icl_release_kernel(ageHist);
	free(histogram);
}

icl_kernel* createInitialPopulation(icl_buffer* agents, icl_buffer* agentAges, icl_buffer* agentStates, struct Population* pop, icl_buffer* enbD,
		icl_buffer* seeds, icl_device* dev) {
	icl_kernel* init = icl_create_kernel(dev, "kernel/initAgents.cl", "initAgents", kernelBuildArgs, ICL_SOURCE);

	size_t localSize = LOCAL_SIZE;
	// overprovisioning, actual number known only on device
	size_t globalSize = ((maxAgents + LOCAL_SIZE - 1) / LOCAL_SIZE) * LOCAL_SIZE;

	icl_run_kernel(init, 1, &globalSize, &localSize, NULL, initEvent, 6,
			(size_t)0, (void *)agents,
			(size_t)0, (void *)agentAges,
			(size_t)0, (void *)agentStates,
			(size_t)0, (void *)enbD,
			sizeof(Temperature), &temperature[0],
			(size_t)0, (void *)seeds);

	pop->eggs.start = 0u;
	pop->eggs.end = initialAgentCount;

	UINT initVal = initialAgentCount + 1; //need offest of one in order to avoid overwriting of prefixSums
	pop->larvae.start = initVal;
	pop->pupae.start = initVal;
	pop->immatures.start = initVal;
	pop->mateSeekings.start = initVal;
	pop->bmSeekings.start = initVal;
	pop->bmDigestings.start = initVal;
	pop->gravids.start = initVal;

	pop->larvae.end = initVal;
	pop->pupae.end = initVal;
	pop->immatures.end = initVal;
	pop->mateSeekings.end = initVal;
	pop->bmSeekings.end = initVal;
	pop->bmDigestings.end = initVal;
	pop->gravids.end = initVal;

//	clFinish(dev->queue);
#if TIMING
	clWaitForEvents(1, initEvent->event);
	initTime += icl_profile_event(initEvent, MEASURE_START, ICL_FINISHED, ICL_MILLI);
#endif

	return init;
}

struct Stats calcStats(icl_buffer* agentAges, icl_buffer* popD, icl_buffer* buff, icl_device* dev,
		icl_kernel* calcL1dePerGroup, icl_kernel* calcL1deTotal,
		icl_kernel* calcFemalesPerGroup, icl_kernel* calcFemalesTotal,
		struct Population* popH) {
	size_t localSize = LOCAL_SIZE;
	UINT upperBound = max(popH->immatures.end, max(popH->bmSeekings.end, max(popH->bmDigestings.end, popH->gravids.end)));
	UINT adultsRange = upperBound - popH->immatures.start;

#if TIMING
	clFinish(dev->queue);
	icl_start_timer(l1deTimer);
#endif

	size_t globalSize = LOCAL_SIZE * LOCAL_SIZE;

	icl_run_kernel(calcL1dePerGroup, 1, &globalSize, &localSize, NULL, NULL, 3,
			(size_t)0, (void *)agentAges,
			(size_t)0, (void *)popD,
			(size_t)0, (void *)buff);


	// called always, since it also copies the properties array
	icl_run_kernel(calcL1deTotal, 1, &localSize, &localSize, NULL, NULL, 3,
			(size_t)0, (void *)popD,
			(size_t)0, (void *)buff,
			sizeof(UINT), &localSize);

#if TIMING
	clFinish(dev->queue);
	icl_stop_timer(l1deTimer);

	icl_start_timer(genderTimer);
#endif

	icl_run_kernel(calcFemalesPerGroup, 1, &globalSize, &localSize, NULL, NULL, 3,
			(size_t)0, (void *)agentAges,
			(size_t)0, (void *)popD,
			(size_t)0, (void *)buff);

	icl_run_kernel(calcFemalesTotal, 1, &localSize, &localSize, NULL, NULL, 2,
			(size_t)0, (void *)popD,
			(size_t)0, (void *)buff);

/*		globalSize = ((adultsRange + LOCAL_SIZE - 1) / LOCAL_SIZE) * LOCAL_SIZE;
	// calculate a histogram with the age of all adults
	icl_run_kernel(ageHistogram, 1, &localSize, &localSize, NULL, NULL, 3,
			(size_t)0, (void *)agents,
			(size_t)0, (void *)popD,
			(size_t)0, (void *)buff);*/

	icl_read_buffer(popD, CL_TRUE, sizeof(struct Population), popH, NULL, NULL);

#if TIMING
	clFinish(dev->queue);
	icl_stop_timer(genderTimer);
#endif

	struct Stats stats;
	stats.numEggs = calcNum(popH->eggs);
	stats.numLarvae = calcNum(popH->larvae);
	stats.numPupae = calcNum(popH->pupae);
	stats.numImmature = calcNum(popH->immatures);
	stats.numMating = calcNum(popH->mateSeekings);
	stats.numBMS = calcNum(popH->bmSeekings);
	stats.numBMD = calcNum(popH->bmDigestings);
	stats.numOVI = calcNum(popH->gravids);

	stats.numLarvae1DazEquiv = (stats.numLarvae > 0) ? popH->numLarvae1DayEquiv : 0;
	stats.numBiomass = stats.numLarvae1DazEquiv + stats.numEggs + stats.numPupae;

	stats.numFemales = adultsRange > 0 ? popH->numFemale : 0;
	stats.numPotentiallyInfective = adultsRange > 0 ? popH->numPotentiallyInfective : 0;
	stats.numMales = (stats.numImmature + stats.numMating + stats.numBMS + stats.numBMD + stats.numOVI) - stats.numFemales;

	UINT nAgents = stats.numEggs + stats.numLarvae + stats.numPupae + stats.numImmature + stats.numMating + stats.numBMS + stats.numBMD + stats.numOVI;
	assert(nAgents <= maxAgents * 0.95 && "not enough capacity");

	assert(nAgents > 0 && "No more agents left");

	return stats;
}

void update(icl_buffer* agents, icl_buffer* agentAges, icl_buffer* agentStates, icl_buffer* enbD, icl_buffer* bnc, icl_buffer* popD, icl_buffer* seeds,
		icl_kernel* killAgents, icl_kernel* updateAgents, UINT worldTime, struct Stats* stats, struct Population* popH, UINT i) {

	size_t localWorkSize = LOCAL_SIZE;
	size_t globalWorkSize = ((popH->gravids.end + LOCAL_SIZE - 1) / LOCAL_SIZE) * LOCAL_SIZE;

	// killing some agents
	icl_run_kernel(killAgents, 1, &globalWorkSize, &localWorkSize, NULL, killEvent, 7,
			(size_t)0, (void *)agentAges,
			(size_t)0, (void *)agentStates,
			(size_t)0, (void *)popD,
			(size_t)0, (void *)enbD,
			(size_t)0, (void *)bnc,
			(size_t)0, (void *)seeds,
			sizeof(REAL), &environment[i].carryingCapacity);

	// update the states of all agents
	icl_run_kernel(updateAgents, 1, &globalWorkSize, &localWorkSize, NULL, updateEvent, 11,
			(size_t)0, (void *)agents,
			(size_t)0, (void *)agentAges,
			(size_t)0, (void *)agentStates,
			(size_t)0, (void *)popD,
			sizeof(struct Environment), &environment[i],
			sizeof(Temperature), &temperature[i],
			(size_t)0, (void *)bnc,
			(size_t)0, (void *)enbD,
			(size_t)0, (void *)seeds,
			sizeof(REAL), &hoursInTimeStep,
			sizeof(UINT), &worldTime);
}

void sequentialPrefixScan(icl_buffer* agentStates, icl_buffer* pop, icl_buffer* prefixSum,
		icl_kernel* serialScan, UINT offset, UINT start, UINT end, INT match) {
	size_t scanLocalWorkSize = 1;
	size_t scanGlobalWorkSize = 1;

	icl_run_kernel(serialScan, 1, &scanGlobalWorkSize, &scanLocalWorkSize, NULL, NULL, 7,
				(size_t)0, (void *)agentStates,
				(size_t)0, (void *)prefixSum,
				(size_t)0, (void *)pop,
				sizeof(UINT), &(offset),
				sizeof(UINT), &(start),
				sizeof(UINT), &(end),
				sizeof(enum State), &match);
}

void parallelPrefixScan(icl_buffer* agentStates, icl_buffer* pop, icl_buffer* prefixSum, icl_buffer* flags,
		icl_kernel* preScan, UINT start, UINT break1, UINT break2, UINT end, INT match) {
	size_t localWorkSize = LOCAL_SIZE;
	size_t globalWorkSize = ((end + LOCAL_SIZE) / LOCAL_SIZE) * LOCAL_SIZE;
	icl_run_kernel(preScan, 1, &globalWorkSize, &localWorkSize, NULL, preScanEvent, 9,
			(size_t)0, (void *)agentStates,
			(size_t)0, (void *)pop,
			(size_t)0, (void *)prefixSum,
			(size_t)0, (void *)flags,
			sizeof(UINT), &(start),
			sizeof(UINT), &(break1),
			sizeof(UINT), &(break2),
			sizeof(UINT), &(end),
			sizeof(INT), &match);

	scan(prefixSum, flags, end+1-start);

#if TIMING
	clWaitForEvents(1u, preScanEvent->event);
	preScanTime += icl_profile_event(preScanEvent, MEASURE_START, ICL_FINISHED, ICL_MILLI);
#endif
}

void createNewAgents(icl_buffer* agents, icl_buffer* agentAges, icl_buffer* agentStates, icl_buffer* newAgents, icl_buffer* newAgentAges,
		icl_buffer* newAgentStates, icl_buffer* enbD, icl_buffer* popD,
		icl_buffer* seeds, icl_buffer* prefixSum1, icl_buffer* prefixSum2, icl_buffer* prefixSum3, icl_buffer* buff,
		icl_kernel* createEggs, icl_kernel* serialScan, icl_kernel* preScan, icl_kernel* oldToNewAgents, struct Population* popH, icl_device* dev, UINT i) {
	// TODO add check for over limit size
	assert(1);

	size_t localWorkSize = LOCAL_SIZE;
	size_t globalWorkSize = ((maxAgents + LOCAL_SIZE) / LOCAL_SIZE) * LOCAL_SIZE; // one additional thread to write properties in oldToNew agents

	icl_run_kernel(createEggs, 1, &globalWorkSize, &localWorkSize, NULL, initEvent, 6,
			(size_t)0, (void *)newAgents,
			(size_t)0, (void *)newAgentAges,
			(size_t)0, (void *)newAgentStates,
			(size_t)0, (void *)enbD,
			sizeof(Temperature), &temperature[i],
			(size_t)0, (void *)seeds);

#if TIMING
	clFinish(dev->queue);
	icl_start_timer(scanTimer);
#endif

#if DEVICE_TYPE != ICL_CPU
#define PARALLEL_SCAN
//	printf("host   b1 %d, b3 %d, end %d\n", popH->eggs.end, popH->immatures.end, popH->bmDigestings.end);
	parallelPrefixScan(agentStates, popD, prefixSum1, buff, preScan,
			0, popH->eggs.end, popH->immatures.end, popH->bmDigestings.end, EGG | IMMATURE | BMD);

	parallelPrefixScan(agentStates, popD, prefixSum2, buff, preScan,
			0, popH->larvae.end, popH->mateSeekings.end, popH->gravids.end, LARVA | MATESEEKING | GRAVID);

	parallelPrefixScan(agentStates, popD, prefixSum3, buff, preScan,
			popH->larvae.start, popH->pupae.end, popH->gravids.end, popH->gravids.end, PUPA | BMS);

#else
	sequentialPrefixScan(agentStates, popD, prefixSum1, serialScan,
			0, popH->eggs.start, popH->eggs.end, EGG);

	sequentialPrefixScan(agentStates, popD, prefixSum2, serialScan,
			0, popH->eggs.start, popH->larvae.end, LARVA);

	sequentialPrefixScan(agentStates, popD, prefixSum3, serialScan,
			popH->larvae.start, popH->larvae.start, popH->pupae.end, PUPA);

	sequentialPrefixScan(agentStates, popD, prefixSum1, serialScan,
			0, popH->pupae.start, popH->immatures.end, IMMATURE);

	sequentialPrefixScan(agentStates, popD, prefixSum2, serialScan,
			0, popH->immatures.start, popH->mateSeekings.end, MATESEEKING);

	sequentialPrefixScan(agentStates, popD, prefixSum3, serialScan,
			popH->larvae.start, popH->mateSeekings.start, popH->gravids.end, BMS);

	sequentialPrefixScan(agentStates, popD, prefixSum1, serialScan,
			0, popH->bmSeekings.start, popH->bmDigestings.end, BMD);

	sequentialPrefixScan(agentStates, popD, prefixSum2, serialScan,
			0, popH->bmDigestings.start, popH->gravids.end, GRAVID);
#endif

#if TIMING
	clFinish(dev->queue);
	icl_stop_timer(scanTimer);
#endif

	globalWorkSize = ((popH->gravids.end + LOCAL_SIZE) / LOCAL_SIZE) * LOCAL_SIZE; // one additional thread to write properties in oldToNew agents

	// sort agents form old array into new array
	icl_run_kernel(oldToNewAgents, 1, &globalWorkSize, &localWorkSize, NULL, oldToNewEvent, 11,
			(size_t)0, (void *)agents,
			(size_t)0, (void *)agentAges,
			(size_t)0, (void *)agentStates,
			(size_t)0, (void *)newAgents,
			(size_t)0, (void *)newAgentAges,
			(size_t)0, (void *)newAgentStates,
			(size_t)0, (void *)prefixSum1,
			(size_t)0, (void *)prefixSum2,
			(size_t)0, (void *)prefixSum3,
			(size_t)0, (void *)popD,
			(size_t)0, (void *)enbD);
}

void run(icl_buffer* enbD, icl_buffer* bnc, icl_buffer* agents, icl_buffer* agentAges, icl_buffer* agentStates,
		icl_buffer* seedsD, struct Population* popH, struct Seeds* seedsH, icl_device* dev, icl_kernel* createEggs) {
	icl_buffer* newAgents = icl_create_buffer(dev, CL_MEM_READ_WRITE, maxAgents * sizeof(struct Agent));
	icl_buffer* newAgentAges = icl_create_buffer(dev, CL_MEM_READ_WRITE, maxAgents * sizeof(struct AgentAge));
	icl_buffer* newAgentStates = icl_create_buffer(dev, CL_MEM_READ_WRITE, maxAgents * sizeof(struct AgentState));
	icl_buffer* popD = icl_create_buffer(dev, CL_MEM_READ_WRITE, sizeof(struct Population) * 2); // using double buffering
	icl_buffer* buff = icl_create_buffer(dev, CL_MEM_READ_WRITE, buffSize * sizeof(UINT));

	icl_buffer* prefixSum1 = icl_create_buffer(dev, CL_MEM_READ_WRITE, maxAgents * sizeof(INT));
	icl_buffer* prefixSum2 = icl_create_buffer(dev, CL_MEM_READ_WRITE, maxAgents * sizeof(INT));
	icl_buffer* prefixSum3 = icl_create_buffer(dev, CL_MEM_READ_WRITE, maxAgents * sizeof(INT));

	// write initial population to popD[1]. calcStats will read it from there and copy to popD[0]
	icl_write_buffer_offset(popD, CL_TRUE, sizeof(struct Population), sizeof(struct Population), popH, NULL, NULL);

	icl_kernel* calcL1dePerGroup = icl_create_kernel(dev, "kernel/calcL1de.cl", "calcL1dePerGroup", kernelBuildArgs, ICL_SOURCE);
	icl_kernel* calcL1deTotal = icl_create_kernel(dev, "kernel/calcL1de.cl", "calcL1deTotal", kernelBuildArgs, ICL_SOURCE);
	icl_kernel* calcFemalesPerGroup = icl_create_kernel(dev, "kernel/calcGender.cl", "calcFemalesPerGroup", kernelBuildArgs, ICL_SOURCE);
	icl_kernel* calcFemalesTotal = icl_create_kernel(dev, "kernel/calcGender.cl", "calcFemalesTotal", kernelBuildArgs, ICL_SOURCE);

//	icl_kernel* ageHistogram = icl_create_kernel(dev, "kernel/ageHistogram.cl", "ageHistogram", KERNEL_BUILD_MACRO, ICL_SOURCE);
	icl_kernel* killAgents = icl_create_kernel(dev, "kernel/killAgents.cl", "killAgents", kernelBuildArgs, ICL_SOURCE);
	icl_kernel* updateAgents = icl_create_kernel(dev, "kernel/updateAgents.cl", "updateAgents", kernelBuildArgs, ICL_SOURCE);
	icl_kernel* oldToNewAgents = icl_create_kernel(dev, "kernel/oldToNewAgents.cl", "oldToNewAgents", kernelBuildArgs, ICL_SOURCE);

	icl_kernel* serialScan = icl_create_kernel(dev, "kernel/serialScan.cl", "serialScan", kernelBuildArgs, ICL_SOURCE);
	icl_kernel* preScan =  icl_create_kernel(dev, "kernel/preScan.cl", "preScan", kernelBuildArgs, ICL_SOURCE);

	// create temporary buffers for prefix scan
#ifdef PARALLEL_SCAN
#if LOCAL_SIZE == 1
#define SCAN_SIZE 1
#else
#define SCAN_SIZE LOCAL_SIZE*2
#endif
	segmented_scan_init(SCAN_SIZE, maxAgents, dev, kernelBuildArgs, ICL_SOURCE);
#endif

	for(UINT currentStep = 0; currentStep < maxSteps; ++currentStep) {
		struct Stats stats = calcStats(agentAges, popD, buff, dev, calcL1dePerGroup, calcL1deTotal, calcFemalesPerGroup, calcFemalesTotal, popH);
#if DEBUG
		struct BitesNcycles bncH;
		icl_read_buffer(bnc, CL_TRUE, sizeof(struct BitesNcycles), &bncH, NULL, NULL);
		printStats(&stats, &bncH, currentStep);
#endif

		// generate random seedsD and copy them to device
		generateSeeds(seedsD, seedsH, dev);

		update(agents, agentAges, agentStates, enbD, bnc, popD, seedsD, killAgents, updateAgents,
				(UINT)(currentStep * hoursInTimeStep)%24u, &stats, popH, currentStep);

		createNewAgents(agents, agentAges, agentStates, newAgents, newAgentAges, newAgentStates, enbD, popD, seedsD,
				prefixSum1, prefixSum2, prefixSum3, buff, createEggs, serialScan, preScan, oldToNewAgents, popH, dev, currentStep);

		swap(&agents, &newAgents);
		swap(&agentAges, &newAgentAges);
		swap(&agentStates, &newAgentStates);
#if TIMING
		clFinish(dev->queue);
		initTime += icl_profile_event(initEvent, MEASURE_START, ICL_FINISHED, ICL_MILLI);
		killTime += icl_profile_event(killEvent, MEASURE_START, ICL_FINISHED, ICL_MILLI);
		updateTime += icl_profile_event(updateEvent, MEASURE_START, ICL_FINISHED, ICL_MILLI);
		oldToNewTime += icl_profile_event(oldToNewEvent, MEASURE_START, ICL_FINISHED, ICL_MILLI);
#endif
	}
	icl_read_buffer_offset(popD, CL_TRUE, sizeof(struct Population), sizeof(struct Population), popH, NULL, NULL);
	printPopulation(popH);

#if TIMING
	printf(" preScan\t\t%f\n", preScanTime);
#endif

#ifdef PARALLEL_SCAN
	segmented_scan_release();
#endif

	storeHistogram(agentAges, popD, buff, dev);

	icl_release_kernel(calcL1dePerGroup);
	icl_release_kernel(calcL1deTotal);
	icl_release_kernel(calcFemalesPerGroup);
	icl_release_kernel(calcFemalesTotal);
	icl_release_kernel(killAgents);
	icl_release_kernel(updateAgents);
	icl_release_kernel(createEggs);
	icl_release_kernel(oldToNewAgents);
	icl_release_buffers(5, newAgents, newAgentAges, newAgentStates, popD, buff, seedsD);
	icl_release_buffers(3, prefixSum1, prefixSum2, prefixSum3);

	icl_release_buffers(5, enbD, bnc, agents, agentAges, agentStates);
}

int readArguments(int argc, char **argv, char* temperaturePath, char* environmentPath) {
	if(argc > 1) {
		if(strcmp(argv[1], "-h") == 0) {
			printf("Usage: bin\\abms temperatueFileName environmentFileName speciesHeaderFileName\n");
			return -1;
		}
	}

	if(argc > 3)
		sprintf(kernelBuildArgs, "-I%s -DSPECIES=%s", KENRNEL_INCLUDE_PATH, argv[3]);
	else
		sprintf(kernelBuildArgs, "-I%s -DSPECIES=%s", KENRNEL_INCLUDE_PATH, "properties.h");

	if(argc > 2)
		sprintf(environmentPath, "%s",  argv[2]);
	else
		sprintf(environmentPath, "%s", "environment.txt");

	if(argc > 1)
		sprintf(temperaturePath, "%s", argv[1]);
	else
		sprintf(temperaturePath, "%s", "temp.txt");

	printf("Temperature:\t%s\nEnvironment:\t%s\nSpecies:\t%s\n", temperaturePath, environmentPath, argc > 3 ? argv[3] : "properties.h");

	return 0;
}


int main (int argc, char **argv)
{
	char temperaturePath[512],  environmentPath[512];

	if(readArguments(argc, argv, temperaturePath, environmentPath))
		return 0;
	// setup input
/*
//	environmentH.immValue = 0.0f;
//	environmentH.oviValue = 0.0f;
	environment[0].IRSValue = 0.0f;
	environment[0].ITNValue = 0.0f;
	environment[0].bloodmealSuccess = 0.25f;
//	environmentH.step = 0;

	//aquatic fields
	environment[0].oviTrapValue = 0.0f;
	environment[0].larvacideValue = 0.0f;
*/

	struct EggsNbiomass enbH;
	// number of initial eggs in environment
	enbH.newEggs = initialAgentCount;
	enbH.totalBiomass = initialAgentCount;

	srand(time(NULL));

	// init ocl
	icl_init_devices(DEVICE_TYPE);

	if (icl_get_num_devices() != 0)
	{
		icl_device* dev = icl_get_device(0);

		icl_print_device_short_info(dev);

		initEvent = icl_create_event();
		killEvent = icl_create_event();
		updateEvent = icl_create_event();
		oldToNewEvent = icl_create_event();
#ifdef PARALLEL_SCAN
		preScanEvent = icl_create_event();
#endif

#if TIMING
		initTime = 0.0;
		l1deTimer = icl_init_timer(ICL_MILLI);
		genderTimer = icl_init_timer(ICL_MILLI);
		randTimer = icl_init_timer(ICL_MILLI);
		killTime = 0.0;
		updateTime = 0.0;
		scanTimer = icl_init_timer(ICL_MILLI);
		oldToNewTime = 0.0;
#endif

		icl_timer* totalTime = icl_init_timer(ICL_MILLI);
		icl_start_timer(totalTime);

		// seeds for random number generation
		icl_buffer* seedsD = icl_create_buffer(dev, CL_MEM_READ_ONLY, nSeeds * sizeof(struct Seeds));
		struct Seeds* seedsH = (struct Seeds*)malloc(nSeeds * sizeof(struct Seeds));

		//create initial population
		struct Population population;
		icl_buffer* enbD = icl_create_buffer(dev, CL_MEM_READ_WRITE, sizeof(struct EggsNbiomass));
		icl_buffer* bnc = icl_create_buffer(dev, CL_MEM_READ_WRITE, sizeof(struct BitesNcycles));
		icl_buffer* agents = icl_create_buffer(dev, CL_MEM_READ_WRITE, maxAgents * sizeof(struct Agent));
		icl_buffer* agentAges = icl_create_buffer(dev, CL_MEM_READ_WRITE, maxAgents * sizeof(struct AgentAge));
		icl_buffer* agentStates = icl_create_buffer(dev, CL_MEM_READ_WRITE, maxAgents * sizeof(struct AgentState));

		icl_write_buffer(enbD, CL_TRUE, sizeof(struct EggsNbiomass), &enbH, NULL, NULL);

		struct BitesNcycles bncH;
		bncH.numCyclesReported = 0;
		bncH.sumCyclesReported = 0;
		bncH.numBitesReported = 0;
		bncH.numInfectBitesReported = 0;
		icl_write_buffer(bnc, CL_TRUE, sizeof(struct BitesNcycles), &bncH, NULL, NULL);
		icl_write_buffer(enbD, CL_TRUE, sizeof(struct EggsNbiomass), &enbH, NULL, NULL);

		if(readSingleValFile(temperaturePath, &temperature) < 0)
			return -1;
		if(readEnvironment(environmentPath) < 0)
			return -1;

		generateSeeds(seedsD, seedsH, dev);
		icl_kernel* createEggs = createInitialPopulation(agents, agentAges, agentStates, &population, enbD, seedsD, dev);

		run(enbD, bnc, agents, agentAges, agentStates, seedsD, &population, seedsH, dev, createEggs);

		clFinish(dev->queue);
		free(seedsH);

		icl_stop_timer(totalTime);
		icl_release_events(4, initEvent, killEvent, updateEvent, oldToNewEvent);
#ifdef PARALLEL_SCAN
		icl_release_event(preScanEvent);
#endif
#if TIMING
		printf("init\t\t%lf\nl1de\t\t%lf\ngender\t\t%lf\nrand\t\t%lf\nkill\t\t%lf\nupdate\t\t%lf\nscan\t\t%lf\noldToNew\t%lf\n",
				initTime, l1deTimer->current_time, genderTimer->current_time, randTimer->current_time,
				killTime, updateTime, scanTimer->current_time, oldToNewTime);

		icl_release_timer(l1deTimer);
		icl_release_timer(genderTimer);
		icl_release_timer(randTimer);
		icl_release_timer(scanTimer);
#endif
		icl_release_devices();
		free(temperature);
		free(environment);

		printf("Execution time: %f ms\n", totalTime->current_time);
		icl_release_timer(totalTime);
	} else {
		printf("Cannot find any OpenCL device of %s\n", DEVICE_TYPE == ICL_CPU ? "type CPU" : DEVICE_TYPE == ICL_GPU ? "type GPU" : "requested type");
	}

	return 0;
}

