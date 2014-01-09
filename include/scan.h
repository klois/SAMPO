
#pragma once
#include "abms.h"

void segmented_scan_init(size_t _wx, UINT n, icl_device *dev, const char* build_options, icl_create_kernel_flag flag);
void segmented_scan_release();

void scan_host(UINT *data, UINT *flag, UINT n);
void scan(icl_buffer *data, icl_buffer *flag, UINT n);

void recursive_scan(icl_buffer *d_data, icl_buffer *d_part, icl_buffer *d_flag, UINT n);


