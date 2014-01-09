#include "device_types.h"

/***************************************************************************
*   Copyright 2012 - 2013 Advanced Micro Devices, Inc.                                     
*                                                                                    
*   Licensed under the Apache License, Version 2.0 (the "License");   
*   you may not use this file except in compliance with the License.                 
*   You may obtain a copy of the License at                                          
*                                                                                    
*       http://www.apache.org/licenses/LICENSE-2.0                      
*                                                                                    
*   Unless required by applicable law or agreed to in writing, software              
*   distributed under the License is distributed on an "AS IS" BASIS,              
*   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.         
*   See the License for the specific language governing permissions and              
*   limitations under the License.                                                   
***************************************************************************/

/******************************************************************************
 *  Kernel 0
 *****************************************************************************/
__kernel void perBlockScanByKey(
	__global UINT *keys,
	__global INT *vals,
	INT init,
    const uint vecSize,
    __local UINT *ldsKeys,
    __local INT *ldsVals,
    __global UINT *keyBuffer,
    __global INT *valBuffer)
{
    size_t gloId = get_global_id( 0 );
    size_t groId = get_group_id( 0 );
    size_t locId = get_local_id( 0 );
    size_t wgSize = get_local_size( 0 );
    //printf("gid=%i, lTid=%i, gTid=%i\n", groId, locId, gloId);

    // if exclusive, load gloId=0 w/ init, and all others shifted-1
    UINT key;
    INT val;
	  if (gloId < vecSize){
          //vType inVal = vals[gloId];
          //val = (INT) (*unaryOp)(inVal);
          key = keys[ gloId ];
          val = vals[ gloId ];
          ldsKeys[ locId ] = key;
          ldsVals[ locId ] = val;
     }
    // Computes a scan within a workgroup
    // updates vals in lds but not keys
    INT sum = val;
    for( size_t offset = 1; offset < wgSize; offset *= 2 )
    {
        barrier( CLK_LOCAL_MEM_FENCE );
        if (locId >= offset )
        {
            UINT key2 = ldsKeys[ locId - offset];
            if( ( key == key2 )  )
            {
                INT y = ldsVals[ locId - offset ];
                sum = ( sum + y );
            }
        }
        barrier( CLK_LOCAL_MEM_FENCE );
        ldsVals[ locId ] = sum;
    }
    
  	barrier( CLK_LOCAL_MEM_FENCE ); // needed for large data types
	  //  Abort threads that are passed the end of the input vector
    if (gloId >= vecSize) return; 

    // Each work item writes out its calculated scan result, relative to the beginning
    // of each work group
    UINT curkey, prekey;

    vals[ gloId ] = sum + init;

    
    if (locId == 0)
    {
        keyBuffer[ groId ] = ldsKeys[ wgSize-1 ];
        valBuffer[ groId ] = ldsVals[ wgSize-1 ];
    }
}


/******************************************************************************
 *  Kernel 1
 *****************************************************************************/
__kernel void intraBlockInclusiveScanByKey(
	__global UINT *keySumArray,
	__global INT *preSumArray,
	__global INT *postSumArray,
	const uint vecSize,
	__local UINT *ldsKeys,
	__local INT *ldsVals,
    const uint workPerThread)
{
    size_t groId = get_group_id( 0 );
    size_t gloId = get_global_id( 0 );
    size_t locId = get_local_id( 0 );
    size_t wgSize = get_local_size( 0 );
    uint mapId  = gloId * workPerThread;

    // do offset of zero manually
    uint offset;
    UINT key;
    INT workSum;
    if (mapId < vecSize)
    {
        UINT prevKey;

        // accumulate zeroth value manually
        offset = 0;
        key = keySumArray[ mapId+offset ];
        workSum = preSumArray[ mapId+offset ];
        postSumArray[ mapId+offset ] = workSum;

        //  Serial accumulation
        for( offset = offset+1; offset < workPerThread; offset += 1 )
        {
            prevKey = key;
            key = keySumArray[ mapId+offset ];
            if (mapId+offset<vecSize )
            {
                INT y = preSumArray[ mapId+offset ];
                if ( (key == prevKey ) )
                {
                    workSum = ( workSum + y );
                }
                else
                {
                    workSum = y;
                }
                postSumArray[ mapId+offset ] = workSum;
            }
        }
    }
    barrier( CLK_LOCAL_MEM_FENCE );
    INT scanSum = workSum;
    offset = 1;
    // load LDS with register sums
  	ldsVals[ locId ] = workSum;
    ldsKeys[ locId ] = key;
    // scan in lds
    for( offset = offset*1; offset < wgSize; offset *= 2 )
    {
        barrier( CLK_LOCAL_MEM_FENCE );
        if (mapId < vecSize)
        {
            if (locId >= offset  )
            {
                INT y = ldsVals[ locId - offset ];
                UINT key1 = ldsKeys[ locId ];
                UINT key2 = ldsKeys[ locId-offset ];
                if ( ( key1 == key2 ) )
                {
                    scanSum = ( scanSum + y );
                }
			        	else
				            scanSum = ldsVals[ locId ];
            }
			
        }
		    barrier( CLK_LOCAL_MEM_FENCE );
        ldsVals[ locId ] = scanSum;
    } // for offset
    barrier( CLK_LOCAL_MEM_FENCE );
    
    // write final scan from pre-scan and lds scan
    for( offset = 0; offset < workPerThread; offset += 1 )
    {
        barrier( CLK_GLOBAL_MEM_FENCE );

        if (mapId < vecSize && locId > 0)
        {
            INT y = postSumArray[ mapId+offset ];
            UINT key1 = keySumArray[ mapId+offset ]; // change me
            UINT key2 = ldsKeys[ locId-1 ];
            if ( ( key1 == key2 ) )
            {
                INT y2 = ldsVals[locId-1];
                y = ( y + y2 );
            }
            postSumArray[ mapId+offset ] = y;
        } // thread in bounds
    } // for 
} // end kernel


/******************************************************************************
 *  Kernel 2
 *****************************************************************************/
__kernel void perBlockAdditionByKey(
	__global UINT *keySumArray,
	__global INT *postSumArray,
	__global UINT *keys,
	__global INT *vals,
    const uint vecSize)
{
    size_t gloId = get_global_id( 0 );
    size_t groId = get_group_id( 0 );
    size_t locId = get_local_id( 0 );

    //  Abort threads that are passed the end of the input vector
    if( gloId >= vecSize )
        return;
        
    INT scanResult = vals[ gloId ];

    // accumulate prefix
    UINT key1 = keySumArray[ groId-1 ];
    UINT key2 = keys[ gloId ];
    if (groId > 0 && ( key1 == key2 ) )
    {
        INT postBlockSum = postSumArray[ groId-1 ];
        INT newResult = ( scanResult + postBlockSum );
        vals[ gloId ] = newResult;
    }
}
