/*
 * pcuda.h
 *
 *  Created on: Dec 19, 2012
 *	  Author: partio
 *
 * Most of the functionality is in the header -- the himan executable does not necessarily
 * find all the necessary symbols if everything is defined in source.
 * This is very irritating since this means we need to link himan executable with cuda
 * libraries.
 */

#ifndef PCUDA_H
#define PCUDA_H

#include "auxiliary_plugin.h"
#include "himan_common.h"

#ifdef HAVE_CUDA
#include <cuda_runtime_api.h>
#endif

namespace himan
{
namespace plugin
{

class pcuda : public auxiliary_plugin
{
public:
	pcuda();

	virtual ~pcuda() {};

	virtual std::string ClassName() const
	{
		return "himan::plugin::pcuda";
	}

	virtual HPPluginClass PluginClass() const
	{
		return kAuxiliary;
	}

	virtual HPVersionNumber Version() const
	{
		return HPVersionNumber(0, 1);
	}

	/**
	 * @brief Check if this server has cuda enabled devices
	 */

	bool HaveCuda() const
	{
		return !(DeviceCount() == 0);
	}
	
	int DeviceCount() const;

#ifdef HAVE_CUDA

	void Capabilities() const;
	int LibraryVersion() const;
	HPVersionNumber ComputeCapability() const;
	bool SetDevice(int deviceId) const;
	void Reset() const;

	
#endif

private:
	mutable int itsDeviceCount;
	
};

#ifndef HAVE_CUDA

inline
int pcuda::DeviceCount() const
{
	return 0;
}
#else

inline
void pcuda::Capabilities() const
{
	int devCount = DeviceCount();

	if (devCount == 0)
	{
		std::cout << "No CUDA devices found" << std::endl;
		return;
	}

	std::cout << "#---------------------------------------------------#" << std::endl;
	std::cout << "CUDA library version " << LibraryVersion() << std::endl;
	std::cout << "There are " << devCount << " CUDA device(s)" << std::endl;

	// Iterate through devices
	for (int i = 0; i < devCount; ++i)
	{
		// Get device properties
		std::cout << "CUDA Device #" << i << std::endl;

		cudaDeviceProp devProp;
		cudaGetDeviceProperties(&devProp, i);

		std::cout << "Major revision number:		 " << devProp.major << std::endl
			 << "Minor revision number:		 " << devProp.minor << std::endl
			 << "Name:						  " << devProp.name << std::endl
			 << "Total global memory:		   " << devProp.totalGlobalMem << std::endl
			 << "Total shared memory per block: " << devProp.sharedMemPerBlock << std::endl
			 << "Total registers per block:	 " << devProp.regsPerBlock << std::endl
			 << "Warp size:					 " << devProp.warpSize << std::endl
			 << "Maximum memory pitch:		  " << devProp.memPitch << std::endl
			 << "Maximum threads per block:	 " << devProp.maxThreadsPerBlock << std::endl;

		for (int i = 0; i < 3; ++i)
		{
			std::cout << "Maximum dimension " << i << " of block:  " << devProp.maxThreadsDim[i] << std::endl;
		}

		for (int i = 0; i < 3; ++i)
		{
			std::cout << "Maximum dimension " << i << " of grid:   " << devProp.maxGridSize[i] << std::endl;
		}

		std::cout << "Clock rate:					" << devProp.clockRate << std::endl
			 << "Total constant memory:		 " << devProp.totalConstMem << std::endl
			 << "Texture alignment:			 " << devProp.textureAlignment << std::endl
			 << "Concurrent copy and execution: " << (devProp.deviceOverlap ? "Yes" : "No") << std::endl
			 << "Number of multiprocessors:	 " << devProp.multiProcessorCount << std::endl
			 << "Kernel execution timeout:	  " << (devProp.kernelExecTimeoutEnabled ? "Yes" : "No") << std::endl << std::endl;

	}
	std::cout << "#---------------------------------------------------#" << std::endl;

}

inline
int pcuda::DeviceCount() const
{
	if (itsDeviceCount == kHPMissingInt)
	{
		cudaError_t err = cudaGetDeviceCount(&itsDeviceCount);

		if (err == cudaErrorNoDevice || err == cudaErrorInsufficientDriver)
		{
			// No device or no driver present

			itsDeviceCount = 0;
		}
	}

	return itsDeviceCount;
}

int pcuda::LibraryVersion() const
{
   // todo: error checking
	int ver;

	cudaDriverGetVersion(&ver);

	return ver;
}

#endif

#ifndef HIMAN_AUXILIARY_INCLUDE

// the class factory

extern "C" std::shared_ptr<himan_plugin> create()
{
	return std::shared_ptr<pcuda> (new pcuda());
}

#endif /* HIMAN_AUXILIARY_INCLUDE */

} // namespace plugin
} // namespace himan

#endif /* PCUDA_H */
