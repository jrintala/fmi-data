/**
 * @file tk2tc.cpp
 *
 * @dateNov 20, 2012
 * @author partio
 */

#include "tk2tc.h"
#include "plugin_factory.h"
#include "logger_factory.h"
#include "timer_factory.h"
#include <boost/lexical_cast.hpp>
#include <boost/thread.hpp>

#define HIMAN_AUXILIARY_INCLUDE

#include "fetcher.h"

#undef HIMAN_AUXILIARY_INCLUDE

using namespace std;
using namespace himan::plugin;

#include "tk2tc_cuda.h"
#include "cuda_helper.h"

tk2tc::tk2tc()
{
	itsClearTextFormula = "Tc = Tk - 273.15";

	itsLogger = unique_ptr<logger> (logger_factory::Instance()->GetLog("tk2tc"));
}

void tk2tc::Process(std::shared_ptr<const plugin_configuration> conf)
{

	unique_ptr<timer> aTimer;

	// Get number of threads to use

	unsigned short threadCount = ThreadCount(conf->ThreadCount());

	if (conf->StatisticsEnabled())
	{
		aTimer = unique_ptr<timer> (timer_factory::Instance()->GetTimer());
		aTimer->Start();

		conf->Statistics()->UsedThreadCount(threadCount);
		conf->Statistics()->UsedGPUCount(conf->CudaDeviceCount());
	}
	
	boost::thread_group g;

	shared_ptr<info> targetInfo = conf->Info();

	/*
	 * Set target parameter to temperature
	 * - name T-C
	 * - univ_id 4
	 * - grib2 descriptor 0'00'000
	 *
	 * We need to specify grib and querydata parameter information
	 * since we don't know which one will be the output format.
	 *
	 */

	vector<param> theParams;

	param requestedParam("T-C", 4);

	// GRIB 2
	
	requestedParam.GribDiscipline(0);
	requestedParam.GribCategory(0);
	requestedParam.GribParameter(0);

	// GRIB 1

	theParams.push_back(requestedParam);

	if (conf->OutputFileType() == kGRIB1)
	{
		StoreGrib1ParameterDefinitions(theParams, targetInfo->Producer().TableVersion());
	}

	targetInfo->Params(theParams);

	/*
	 * Create data structures.
	 */

	targetInfo->Create();

	/*
	 * Initialize parent class functions for dimension handling
	 */

	Dimension(conf->LeadingDimension());
	FeederInfo(shared_ptr<info> (new info(*targetInfo)));
	FeederInfo()->Param(requestedParam);
	
	if (conf->StatisticsEnabled())
	{
		aTimer->Stop();
		conf->Statistics()->AddToInitTime(aTimer->GetTime());

		aTimer->Start();
	}

	/*
	 * Each thread will have a copy of the target info.
	 */

	for (size_t i = 0; i < threadCount; i++)
	{

		itsLogger->Info("Thread " + boost::lexical_cast<string> (i + 1) + " starting");

		boost::thread* t = new boost::thread(&tk2tc::Run,
											 this,
											 shared_ptr<info> (new info(*targetInfo)),
											 conf,
											 i + 1);

		g.add_thread(t);

	}

	g.join_all();

	if (conf->StatisticsEnabled())
	{
		aTimer->Stop();
		conf->Statistics()->AddToProcessingTime(aTimer->GetTime());
	}

	WriteToFile(conf, targetInfo);
}

void tk2tc::Run(shared_ptr<info> myTargetInfo,
				shared_ptr<const plugin_configuration> conf,
				unsigned short threadIndex)
{

	while (AdjustLeadingDimension(myTargetInfo))
	{
		Calculate(myTargetInfo, conf, threadIndex);
	}
}

/*
 * Calculate()
 *
 * This function does the actual calculation.
 */

void tk2tc::Calculate(shared_ptr<info> myTargetInfo,
					  shared_ptr<const plugin_configuration> conf,
					  unsigned short threadIndex)
{
	shared_ptr<fetcher> aFetcher = dynamic_pointer_cast <fetcher> (plugin_factory::Instance()->Plugin("fetcher"));

	// Required source parameters

	param TParam("T-K");

	unique_ptr<logger> myThreadedLogger = std::unique_ptr<logger> (logger_factory::Instance()->GetLog("tk2tcThread #" + boost::lexical_cast<string> (threadIndex)));

	ResetNonLeadingDimension(myTargetInfo);

	myTargetInfo->FirstParam();

	bool useCudaInThisThread = conf->UseCuda() && threadIndex <= conf->CudaDeviceCount();

	while (AdjustNonLeadingDimension(myTargetInfo))
	{

		myThreadedLogger->Debug("Calculating time " + myTargetInfo->Time().ValidDateTime()->String("%Y%m%d%H%M") +
								" level " + boost::lexical_cast<string> (myTargetInfo->Level().Value()));

		// Source info for T

		shared_ptr<info> TInfo;

		try
		{
			TInfo = aFetcher->Fetch(conf,
								 myTargetInfo->Time(),
								 myTargetInfo->Level(),
								 TParam,
								 conf->UseCudaForPacking() && useCudaInThisThread);

			assert(TInfo->Param().Unit() == kK);

		}
		catch (HPExceptionType e)
		{
			switch (e)
			{
				case kFileDataNotFound:
					itsLogger->Warning("Skipping step " + boost::lexical_cast<string> (myTargetInfo->Time().Step()) + ", level " + boost::lexical_cast<string> (myTargetInfo->Level().Value()));
					myTargetInfo->Data()->Fill(kFloatMissing);

					if (conf->StatisticsEnabled())
					{
						conf->Statistics()->AddToMissingCount(myTargetInfo->Grid()->Size());
						conf->Statistics()->AddToValueCount(myTargetInfo->Grid()->Size());
					}
					
					continue;
					break;

				default:
					throw runtime_error(ClassName() + ": Unable to proceed");
					break;
			}
		}

		SetAB(myTargetInfo, TInfo);

		size_t missingCount = 0;
		size_t count = 0;

		shared_ptr<NFmiGrid> targetGrid(myTargetInfo->Grid()->ToNewbaseGrid());
		shared_ptr<NFmiGrid> TGrid(TInfo->Grid()->ToNewbaseGrid());

		bool equalGrids = (*myTargetInfo->Grid() == *TInfo->Grid());

		string deviceType;

#ifdef HAVE_CUDA
		if (useCudaInThisThread && equalGrids)
		{
	
			deviceType = "GPU";

			tk2tc_cuda::tk2tc_cuda_options opts;
			tk2tc_cuda::tk2tc_cuda_data datas;

			opts.N = TGrid->Size();

			opts.cudaDeviceIndex = static_cast <unsigned short> (threadIndex-1);

			CUDA_CHECK(cudaHostAlloc(reinterpret_cast<void**> (&datas.TC), opts.N * sizeof(double), cudaHostAllocMapped));

			if (TInfo->Grid()->DataIsPacked())
			{
				assert(TInfo->Grid()->PackedData()->ClassName() == "simple_packed");

				shared_ptr<simple_packed> s = dynamic_pointer_cast<simple_packed> (TInfo->Grid()->PackedData());

				datas.pTK = *(s);

				CUDA_CHECK(cudaHostAlloc(reinterpret_cast<void**> (&datas.TK), opts.N * sizeof(double), cudaHostAllocMapped));

				opts.pTK = true;
			}
			else
			{
				datas.TK = const_cast<double*> (TInfo->Grid()->Data()->Values());
			}
			
			tk2tc_cuda::DoCuda(opts, datas);

			myTargetInfo->Data()->Set(datas.TC, opts.N);

			SwapTo(myTargetInfo, TInfo->Grid()->ScanningMode());
			
			missingCount = opts.missingValuesCount;
			count = opts.N;

			CUDA_CHECK(cudaFreeHost(datas.TC));

			if (TInfo->Grid()->DataIsPacked())
			{
				TInfo->Data()->Set(datas.TK, opts.N);
				TInfo->Grid()->PackedData()->Clear();
				CUDA_CHECK(cudaFreeHost(datas.TK));
			}
		}
		else
#endif
		{

			deviceType = "CPU";

			assert(targetGrid->Size() == myTargetInfo->Data()->Size());

			myTargetInfo->ResetLocation();

			targetGrid->Reset();

			while (myTargetInfo->NextLocation() && targetGrid->Next())
			{

				count++;

				double T = kFloatMissing;

				InterpolateToPoint(targetGrid, TGrid, equalGrids, T);

				if (T == kFloatMissing)
				{
					missingCount++;

					myTargetInfo->Value(kFloatMissing);
					continue;
				}

				double TC = T - 273.15;

				if (!myTargetInfo->Value(TC))
				{
					throw runtime_error(ClassName() + ": Failed to set value to matrix");
				}
			}

			/*
			 * Newbase normalizes scanning mode to bottom left -- if that's not what
			 * the target scanning mode is, we have to swap the data back.
			 */

			SwapTo(myTargetInfo, kBottomLeft);

		}
	
		if (conf->StatisticsEnabled())
		{
			conf->Statistics()->AddToMissingCount(missingCount);
			conf->Statistics()->AddToValueCount(count);
		}

		/*
		 * Now we are done for this level
		 *
		 * Clone info-instance to writer since it might change our descriptor places
		 */

		myThreadedLogger->Info("Missing values: " + boost::lexical_cast<string> (missingCount) + "/" + boost::lexical_cast<string> (count));

		if (conf->FileWriteOption() != kSingleFile)
		{
			WriteToFile(conf, myTargetInfo);
		}

	}
}
