// System includes
#include <iostream>
#include <string>

// CUDA runtime
#include <cuda_runtime.h>

#include "vvms_cuda.h"

__global__ void himan::plugin::vvms_cuda::Calculate(const double* __restrict__ d_t,
														const double* __restrict__ d_vv,
														const double* __restrict__ d_p,
														double* __restrict__ d_vv_ms,
														options opts)
{

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < opts.N)
	{
		double P = (opts.is_constant_pressure) ? opts.p_const : d_p[idx];

		if (d_t[idx] == kFloatMissing || d_vv[idx] == kFloatMissing || P == kFloatMissing)
		{
			d_vv_ms[idx] = kFloatMissing;
		}
		else
		{
			d_vv_ms[idx] = opts.vv_ms_scale * (287 * -d_vv[idx] * (opts.t_base + d_t[idx]) / (9.80665 * P * opts.p_scale));
		}
	}
}

void himan::plugin::vvms_cuda::Process(options& opts)
{

	cudaStream_t stream;

	CUDA_CHECK(cudaStreamCreate(&stream));

	size_t memsize = opts.N * sizeof(double);

	// Allocate device arrays

	double* d_t = 0;
	double* d_p = 0;
	double* d_vv = 0;
	double* d_vv_ms = 0;

	CUDA_CHECK(cudaMalloc((void **) &d_vv_ms, memsize));
	CUDA_CHECK(cudaMalloc((void **) &d_t, sizeof(double) * memsize));
	CUDA_CHECK(cudaMalloc((void **) &d_vv, sizeof(double) * memsize));


	if (opts.t->packed_values)
	{
		opts.t->packed_values->Unpack(d_t, &stream);
		CUDA_CHECK(cudaMemcpyAsync(opts.t->values, d_t, memsize, cudaMemcpyDeviceToHost, stream));
	}
	else
	{
		CUDA_CHECK(cudaMemcpyAsync(d_t, opts.t->values, memsize, cudaMemcpyHostToDevice, stream));
	}

	if (opts.vv->packed_values)
	{
		opts.vv->packed_values->Unpack(d_vv, &stream);
		CUDA_CHECK(cudaMemcpyAsync(opts.vv->values, d_vv, memsize, cudaMemcpyDeviceToHost, stream));
	}
	else
	{
		CUDA_CHECK(cudaMemcpyAsync(d_vv, opts.vv->values, memsize, cudaMemcpyHostToDevice, stream));
	}

	if (!opts.is_constant_pressure)
	{
		CUDA_CHECK(cudaMalloc((void **) &d_p, sizeof(double) * memsize));

		if (opts.p->packed_values)
		{
			opts.p->packed_values->Unpack(d_p, &stream);
			CUDA_CHECK(cudaMemcpyAsync(opts.p->values, d_p, memsize, cudaMemcpyDeviceToHost, stream));
		}
		else
		{
			CUDA_CHECK(cudaMemcpyAsync(d_p, opts.p->values, memsize, cudaMemcpyHostToDevice, stream));
		}
	}

	// dims

	const int blockSize = 512;
	const int gridSize = opts.N/blockSize + (opts.N%blockSize == 0?0:1);

	CUDA_CHECK(cudaStreamSynchronize(stream));

	Calculate <<< gridSize, blockSize, 0, stream >>> (d_t, d_vv, d_p, d_vv_ms, opts);
	
	// block until the device has completed
	CUDA_CHECK(cudaStreamSynchronize(stream));

	CUDA_CHECK_ERROR_MSG("Kernel invocation");

	// Retrieve result from device

	CUDA_CHECK(cudaMemcpyAsync(opts.vv_ms->values, d_vv_ms, memsize, cudaMemcpyDeviceToHost, stream));

	CUDA_CHECK(cudaStreamSynchronize(stream));
	
	CUDA_CHECK(cudaFree(d_t));
	CUDA_CHECK(cudaFree(d_vv));
	CUDA_CHECK(cudaFree(d_vv_ms));

	if (d_p)
	{
		CUDA_CHECK(cudaFree(d_p));
	}

	CUDA_CHECK(cudaStreamDestroy(stream));
}
