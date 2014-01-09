/*
 * @author      Biagio Cosenza
 * @date		04/01/2013
 */

#include "abms.h"

#include "scan.h"
#include "math.h"

icl_kernel *scan_pow2;
icl_kernel *scan_pad_to_pow2;
icl_kernel *upsweep_subarrays;
icl_kernel *downsweep_subarrays;
icl_device *dev;

icl_buffer* d_data;
icl_buffer* d_part;
icl_buffer* d_flag;

icl_buffer* d_data2;
icl_buffer* d_part2;
icl_buffer* d_flag2;

size_t wx; // workgroup size
size_t m;     // length of each subarray ( = wx*2 )


void scan_host(UINT *data, UINT *flag, UINT n) {
	/*
	icl_buffer d_data = clw.dev_malloc(sizeof(UINT)*k*m);
	icl_buffer d_part = clw.dev_malloc(sizeof(UINT)*k*m);
	icl_buffer d_flag = clw.dev_malloc(sizeof(UINT)*k*m);
	*/

	UINT k = (UINT) ceil((float)n/(float)m);
	const UINT buf_1 = sizeof(UINT)*k*m;
	icl_buffer* d_data = icl_create_buffer(dev, CL_MEM_READ_WRITE, buf_1);
	icl_buffer* d_part = icl_create_buffer(dev, CL_MEM_READ_WRITE, buf_1);
	icl_buffer* d_flag = icl_create_buffer(dev, CL_MEM_READ_WRITE, buf_1);
	
	/*
	m0 += clw.memcpy_to_dev(d_data, sizeof(UINT)*n, data);
	m1 += clw.memcpy_to_dev(d_part, sizeof(UINT)*n, flag);
	m2 += clw.memcpy_to_dev(d_flag, sizeof(UINT)*n, flag);
	*/
	const UINT buf_2 = sizeof(UINT)*n;
	icl_write_buffer(d_data, CL_FALSE, buf_2, data, NULL, NULL);
	icl_write_buffer(d_part, CL_FALSE, buf_2, flag, NULL, NULL);
	icl_write_buffer(d_flag, CL_FALSE, buf_2, data, NULL, NULL);


	recursive_scan(d_data, d_part, d_flag, n);
	
	/* m3 += clw.memcpy_from_dev(d_data, sizeof(UINT)*n, data); */
	icl_read_buffer(d_data, CL_FALSE, buf_2, data, NULL, NULL);

	/*
	clw.dev_free(d_data);
	clw.dev_free(d_part);
	clw.dev_free(d_flag);
	*/
	icl_release_buffers(3, d_data, d_part, d_flag);
}

void scan(icl_buffer *data, icl_buffer *flag, UINT n) {
	/*
	UINT k = (UINT) ceil((float)n/(float)m);
	cl_mem d_data = clw.dev_malloc(sizeof(UINT)*k*m);
	cl_mem d_part = clw.dev_malloc(sizeof(UINT)*k*m);
	cl_mem d_flag = clw.dev_malloc(sizeof(UINT)*k*m);
	*/
/*	UINT k = (UINT) ceil((float)n/(float)m);
	const UINT buf_1 = sizeof(UINT)*k*m;
	icl_buffer* d_data = icl_create_buffer(dev, CL_MEM_READ_WRITE, buf_1);
	icl_buffer* d_part = icl_create_buffer(dev, CL_MEM_READ_WRITE, buf_1);
	icl_buffer* d_flag = icl_create_buffer(dev, CL_MEM_READ_WRITE, buf_1);
*/
	/*
	clw.copy_buffer(data, d_data, sizeof(UINT)*n);
	clw.copy_buffer(flag, d_part, sizeof(UINT)*n);
	clw.copy_buffer(flag, d_flag, sizeof(UINT)*n);
	*/
	const UINT buf2 = sizeof(UINT)*n;
	icl_copy_buffer(data, d_data, buf2, NULL, NULL);
	icl_copy_buffer(flag, d_part, buf2, NULL, NULL);
	icl_copy_buffer(flag, d_flag, buf2, NULL, NULL);

	recursive_scan(d_data, d_part, d_flag, n);

/*	clw.copy_buffer(d_data, data, sizeof(UINT)*n); */
	icl_copy_buffer(d_data, data, buf2, NULL, NULL);

	/*
	clw.dev_free(d_data);
	clw.dev_free(d_part);
	clw.dev_free(d_flag);
	*/
//	icl_release_buffers(3, d_data, d_part, d_flag);
}

void recursive_scan(icl_buffer *d_data, icl_buffer *d_part, icl_buffer *d_flag, UINT n) {
	UINT k = (UINT) ceil((float)n/(float)m);
	//size of each subarray stored in local memory
	UINT bufsize = m;
	if (k == 1) {
		/*
		clw.kernel_arg(scan_pad_to_pow2,
			d_data,  d_part,  d_flag,
			bufsize, bufsize, bufsize,
			n);
		k0 += clw.run_kernel_with_timing(scan_pad_to_pow2, 1, &wx, &wx); // dim = 1
		*/
				
		icl_run_kernel(scan_pad_to_pow2, 1, &wx, &wx, NULL, NULL, 7,
			(size_t)0, (void *)d_data,
			(size_t)0, (void *)d_part,
			(size_t)0, (void *)d_flag,
			sizeof(UINT) * bufsize, NULL,
			sizeof(UINT) * bufsize, NULL,
			sizeof(UINT) * bufsize, NULL,
			sizeof(UINT), &n
		);

	} else {
		size_t gx = k * wx;

		/*
		cl_mem d_data2 = clw.dev_malloc(sizeof(UINT)*k);
		cl_mem d_part2 = clw.dev_malloc(sizeof(UINT)*k);
		cl_mem d_flag2 = clw.dev_malloc(sizeof(UINT)*k);
		*/
/*		UINT buf_size = sizeof(UINT)*k;
		icl_buffer* d_data2 = icl_create_buffer(dev, CL_MEM_READ_WRITE, buf_size);
		icl_buffer* d_part2 = icl_create_buffer(dev, CL_MEM_READ_WRITE, buf_size);
		icl_buffer* d_flag2 = icl_create_buffer(dev, CL_MEM_READ_WRITE, buf_size);
*/
		/*
		clw.kernel_arg(upsweep_subarrays,
			d_data,  d_part,  d_flag,
			d_data2, d_part2, d_flag2,
			bufsize, bufsize, bufsize,
			n);
		k1 += clw.run_kernel_with_timing(upsweep_subarrays, 1, &gx, &wx);
		*/

		icl_run_kernel(upsweep_subarrays, 1, &gx, &wx, NULL, NULL, 10,
			(size_t)0, (void *)d_data,
			(size_t)0, (void *)d_part,
			(size_t)0, (void *)d_flag,
			(size_t)0, (void *)d_data2,
			(size_t)0, (void *)d_part2,
			(size_t)0, (void *)d_flag2,
			sizeof(UINT) * bufsize, NULL,
			sizeof(UINT) * bufsize, NULL,
			sizeof(UINT) * bufsize, NULL,
			sizeof(UINT), &n
		);

		recursive_scan(d_data2, d_part2, d_flag2, k);

		/*
		clw.kernel_arg(downsweep_subarrays,
			d_data,  d_part,  d_flag,
			d_data2, d_part2, d_flag2,
			bufsize, bufsize, bufsize,
			n);
		k2 += clw.run_kernel_with_timing(downsweep_subarrays, 1, &gx, &wx);
		*/
		icl_run_kernel(downsweep_subarrays, 1, &gx, &wx, NULL, NULL, 10,
			(size_t)0, (void *)d_data,
			(size_t)0, (void *)d_part,
			(size_t)0, (void *)d_flag,
			(size_t)0, (void *)d_data2,
			(size_t)0, (void *)d_part2,
			(size_t)0, (void *)d_flag2,
			sizeof(UINT) * bufsize, NULL,
			sizeof(UINT) * bufsize, NULL,
			sizeof(UINT) * bufsize, NULL,
			sizeof(UINT), &n
		);

		/*
		clw.dev_free(d_data2);
		clw.dev_free(d_part2);
		clw.dev_free(d_flag2);
		*/
//		icl_release_buffers(3, d_data2, d_part2, d_flag2);
	}
}

void segmented_scan_init(UINT _wx, UINT n, icl_device *_dev, const char* build_options, icl_create_kernel_flag flag) {
/*
	: clw(clw), wx(wx),
	m0(0), m1(0), m2(0), m3(0),
	k0(0), k1(0), k2(0) {
*/		
	dev = _dev;
	wx = _wx;
	m = wx * 2;

	UINT k = (UINT) ceil((float)n/(float)m);
	const UINT buf_1 = sizeof(UINT)*k*m;
	d_data = icl_create_buffer(dev, CL_MEM_READ_WRITE, buf_1);
	d_part = icl_create_buffer(dev, CL_MEM_READ_WRITE, buf_1);
	d_flag = icl_create_buffer(dev, CL_MEM_READ_WRITE, buf_1);

	UINT buf_size = sizeof(UINT)*k;
	d_data2 = icl_create_buffer(dev, CL_MEM_READ_WRITE, buf_size);
	d_part2 = icl_create_buffer(dev, CL_MEM_READ_WRITE, buf_size);
	d_flag2 = icl_create_buffer(dev, CL_MEM_READ_WRITE, buf_size);

	scan_pow2 = icl_create_kernel(dev, "kernel/scan.cl", "segscan_pow2_wrapper", build_options, flag);
	scan_pad_to_pow2 = icl_create_kernel(dev, "kernel/scan.cl", "segscan_pad_to_pow2", build_options, flag);
	upsweep_subarrays = icl_create_kernel(dev, "kernel/scan.cl", "upsweep_subarrays", build_options, flag);
	downsweep_subarrays = icl_create_kernel(dev, "kernel/scan.cl", "downsweep_subarrays", build_options, flag);
}

void segmented_scan_release() {
	icl_release_kernel(scan_pow2);
	icl_release_kernel(scan_pad_to_pow2);
	icl_release_kernel(upsweep_subarrays);
	icl_release_kernel(downsweep_subarrays);

	icl_release_buffers(3, d_data, d_part, d_flag);
	icl_release_buffers(3, d_data2, d_part2, d_flag2);
}
