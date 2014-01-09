/*
 * @author      Klaus Kofler
 * @date		07/01/2013
 */

#include <stdio.h>
#include "host_types.h"

#include "lib_icl.h"
#include "scan.h"
#include "math.h"

icl_device* dev;

icl_buffer* preSumArray;
icl_buffer* keySumArray;
icl_buffer* postSumArray;

icl_kernel* perBlockScanByKey;
icl_kernel* intraBlockInclusiveScanByKey;
icl_kernel* perBlockAdditionByKey;

size_t wx;

icl_event* perBlockScanEvent;
icl_event* intraBlockEvent;
icl_event* perBlockAdditionEvent;

#if TIMING
double perBlockScanTime;
double intraBlockTime;
double perBlockAdditionTime;
icl_timer* timer;
#endif

void scan(icl_buffer *data, icl_buffer *flag, UINT n) {
	// use size for actual n, not overapproximation as in allocation
	UINT numWorkGroups = ((n + wx - 1) / wx);
	UINT sizeScanBuff = ((numWorkGroups + wx -1) / wx) * wx;

	size_t gx = numWorkGroups * wx;
	INT init = -1;

#if TIMING
	clFinish(dev->queue);
	icl_start_timer(timer);
#endif

	icl_run_kernel(perBlockScanByKey, 1, &gx, &wx, NULL, perBlockScanEvent, 8,
			(size_t)0, (void *)flag,
			(size_t)0, (void *)data,
			sizeof(INT), &init,
			sizeof(UINT), &n,
			sizeof(UINT) * wx, NULL,
			sizeof(INT) * wx, NULL,
			(size_t)0, (void *)keySumArray,
			(size_t)0, (void *)preSumArray);
/*
    V_OPENCL( kernels[0].setArg( 0, firstKey.getBuffer()), "Error setArg kernels[ 0 ]" ); // Input keys
    V_OPENCL( kernels[0].setArg( 1, firstValue.getBuffer()),"Error setArg kernels[ 0 ]" ); // Input buffer
    V_OPENCL( kernels[0].setArg( 2, result.getBuffer( ) ), "Error setArg kernels[ 0 ]" ); // Output buffer
    V_OPENCL( kernels[0].setArg( 3, init ),                 "Error setArg kernels[ 0 ]" ); // Initial value exclusive
    V_OPENCL( kernels[0].setArg( 4, numElements ),          "Error setArg kernels[ 0 ]" ); // Size of scratch buffer
    V_OPENCL( kernels[0].setArg( 5, ldsKeySize, NULL ),     "Error setArg kernels[ 0 ]" ); // Scratch buffer
    V_OPENCL( kernels[0].setArg( 6, ldsValueSize, NULL ),   "Error setArg kernels[ 0 ]" ); // Scratch buffer
    V_OPENCL( kernels[0].setArg( 7, *binaryPredicateBuffer),"Error setArg kernels[ 0 ]" ); // User provided functor
    V_OPENCL( kernels[0].setArg( 8, *binaryFunctionBuffer ),"Error setArg kernels[ 0 ]" ); // User provided functor
    V_OPENCL( kernels[0].setArg( 9, *keySumArray ),         "Error setArg kernels[ 0 ]" ); // Output per block sum
    V_OPENCL( kernels[0].setArg(10, *preSumArray ),         "Error setArg kernels[ 0 ]" ); // Output per block sum
    V_OPENCL( kernels[0].setArg(11, doExclusiveScan ),      "Error setArg kernels[ 0 ]" ); // Exclusive scan?
*/

	UINT workPerThread = sizeScanBuff / wx;

	icl_run_kernel(intraBlockInclusiveScanByKey, 1, &wx, &wx, NULL, intraBlockEvent, 7,
			(size_t)0, (void *)keySumArray,
			(size_t)0, (void *)preSumArray,
			(size_t)0, (void *)postSumArray,
			sizeof(UINT), &numWorkGroups,
			sizeof(UINT) * wx, NULL,
			sizeof(INT) * wx, NULL,
			sizeof(UINT), &workPerThread);

/*
INT* output = (UINT*)malloc(n * sizeof(INT));
icl_read_buffer(keySumArray, CL_TRUE, n * sizeof(INT), output, NULL, NULL);
for(int i = 0; i < n; ++i)
	printf("%d ", output[i]);
printf("\n"); */
/*
    V_OPENCL( kernels[1].setArg( 0, *keySumArray ),         "Error setArg kernels[ 1 ]" ); // Input keys
    V_OPENCL( kernels[1].setArg( 1, *preSumArray ),         "Error setArg kernels[ 1 ]" ); // Input buffer
    V_OPENCL( kernels[1].setArg( 2, *postSumArray ),        "Error setArg kernels[ 1 ]" ); // Output buffer
    V_OPENCL( kernels[1].setArg( 3, numWorkGroupsK0 ),      "Error setArg kernels[ 1 ]" ); // Size of scratch buffer
    V_OPENCL( kernels[1].setArg( 4, ldsKeySize, NULL ),     "Error setArg kernels[ 1 ]" ); // Scratch buffer
    V_OPENCL( kernels[1].setArg( 5, ldsValueSize, NULL ),   "Error setArg kernels[ 1 ]" ); // Scratch buffer
    V_OPENCL( kernels[1].setArg( 6, workPerThread ),        "Error setArg kernels[ 1 ]" ); // User provided functor
    V_OPENCL( kernels[1].setArg( 7, *binaryPredicateBuffer ),"Error setArg kernels[ 1 ]" ); // User provided functor
    V_OPENCL( kernels[1].setArg( 8, *binaryFunctionBuffer ),"Error setArg kernels[ 1 ]" ); // User provided functor
*/

	icl_run_kernel(perBlockAdditionByKey, 1, &gx, &wx, NULL, perBlockAdditionEvent, 5,
			(size_t)0, (void *)keySumArray,
			(size_t)0, (void *)postSumArray,
			(size_t)0, (void *)flag,
			(size_t)0, (void *)data,
			sizeof(UINT), &n);
/*
    V_OPENCL( kernels[2].setArg( 0, *keySumArray ),         "Error setArg kernels[ 2 ]" ); // Input buffer
    V_OPENCL( kernels[2].setArg( 1, *postSumArray ),        "Error setArg kernels[ 2 ]" ); // Input buffer
    V_OPENCL( kernels[2].setArg( 2, firstKey.getBuffer()), "Error setArg kernels[ 2 ]" ); // Output buffer
    V_OPENCL( kernels[2].setArg( 3, result.getBuffer()),   "Error setArg kernels[ 2 ]" ); // Output buffer
    V_OPENCL( kernels[2].setArg( 4, numElements ),          "Error setArg kernels[ 2 ]" ); // Size of scratch buffer
    V_OPENCL( kernels[2].setArg( 5, *binaryPredicateBuffer ),"Error setArg kernels[ 2 ]" ); // User provided functor
    V_OPENCL( kernels[2].setArg( 6, *binaryFunctionBuffer ),"Error setArg kernels[ 2 ]" ); // User provided functor
*/


#if TIMING
	clFinish(dev->queue);
	icl_stop_timer(timer);
	perBlockScanTime += icl_profile_event(perBlockScanEvent, MEASURE_START, ICL_FINISHED, ICL_MILLI);
	intraBlockTime += icl_profile_event(intraBlockEvent, MEASURE_START, ICL_FINISHED, ICL_MILLI);
	perBlockAdditionTime += icl_profile_event(perBlockAdditionEvent, MEASURE_START, ICL_FINISHED, ICL_MILLI);
#endif
}

void segmented_scan_init(size_t _wx, UINT maxN, icl_device *_dev, const char* build_options, icl_create_kernel_flag flag) {
/*
	: clw(clw), wx(wx),
	m0(0), m1(0), m2(0), m3(0),
	k0(0), k1(0), k2(0) {
*/
	dev = _dev;

	wx = _wx;

	// overapproximation for allocation, using maximum allowed size for n
	UINT numWorkGroups = ((maxN + wx - 1) / wx);
	UINT sizeScanBuff = ((numWorkGroups + wx -1) / wx) * wx;

	const UINT buf_1 = sizeof(UINT)*sizeScanBuff;

	preSumArray = icl_create_buffer(dev, CL_MEM_READ_WRITE, buf_1);
	postSumArray = icl_create_buffer(dev, CL_MEM_READ_WRITE, buf_1);
	keySumArray = icl_create_buffer(dev, CL_MEM_READ_WRITE, buf_1);

	perBlockScanByKey = icl_create_kernel(dev, "kernel/boltScan.cl", "perBlockScanByKey", build_options, flag);
	intraBlockInclusiveScanByKey = icl_create_kernel(dev, "kernel/boltScan.cl", "intraBlockInclusiveScanByKey", build_options, flag);
	perBlockAdditionByKey = icl_create_kernel(dev, "kernel/boltScan.cl", "perBlockAdditionByKey", build_options, flag);

	perBlockScanEvent = icl_create_event();
	intraBlockEvent = icl_create_event();
	perBlockAdditionEvent = icl_create_event();
#if TIMING
	perBlockScanTime = 0;
	intraBlockTime = 0;
	perBlockAdditionTime = 0;
	timer = icl_init_timer(ICL_MILLI);
#endif
}

void segmented_scan_release() {
	icl_release_kernel(perBlockScanByKey);
	icl_release_kernel(intraBlockInclusiveScanByKey);
	icl_release_kernel(perBlockAdditionByKey);

	icl_release_buffers(3, preSumArray, postSumArray, keySumArray);

	icl_release_events(3, perBlockScanEvent, intraBlockEvent, perBlockAdditionEvent);

#if TIMING
	printf(" perBlockScanTime \t%f\n intraBlockTime \t%f\n perBlockAdditionTime \t%f\n", perBlockScanTime, intraBlockTime, perBlockAdditionTime);
	printf(" scan time %f\n", timer->current_time);
	icl_release_timer(timer);
#endif
}
