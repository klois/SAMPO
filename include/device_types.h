/*
 * @author      Klaus Kofler
 * @date		07/01/2013
 */
#pragma once

#define REAL float
#define INT int
#define UINT uint
#define BOOL bool

// only for eclipse/VS for syntax highlighting and to avoid warnings
#ifndef __OPENCL_VERSION__
/* Define out keywords causing errors */
#define __kernel
#define __global
#define __constant
#define __local
#define __private
#define bool int
#define true 1
#define false 0
#define CLK_LOCAL_MEM_FENCE 32
#define CLK_GLOBAL_MEM_FENCE 64
#define align /*__declspec(align(8)) */
#else
#define align /*__attribute__ ((aligned))*/
#endif
