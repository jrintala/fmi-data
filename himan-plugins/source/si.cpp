/**
 * @file si.cpp
 *
 * @date Feb 13, 2014
 * @author partio
 */

#include "si.h"
#include "plugin_factory.h"
#include "logger_factory.h"
#include <boost/lexical_cast.hpp>
#include "metutil.h"
#include "NFmiGrid.h"

#define HIMAN_AUXILIARY_INCLUDE

#include "neons.h"
#include "fetcher.h"
#include "querydata.h"
#include "hitool.h"

#undef HIMAN_AUXILIARY_INCLUDE

using namespace std;
using namespace himan::plugin;

#include "NFmiQueryData.h"
#include "NFmiSoundingIndexCalculator.h"

#if 0
// nicer way to convert between std::shared_ptr and boost::shared_ptr but does not work on gcc 4.4.7 :(

boost::shared_ptr<NFmiQueryData> make_shared_ptr(std::shared_ptr<NFmiQueryData>& ptr)
{
    return boost::shared_ptr<NFmiQueryData>(ptr.get(), [ptr](NFmiQueryData*) mutable {ptr.reset();});
}

std::shared_ptr<NFmiQueryData> make_shared_ptr(boost::shared_ptr<NFmiQueryData>& ptr)
{
    return std::shared_ptr<NFmiQueryData>(ptr.get(), [ptr](NFmiQueryData*) mutable {ptr.reset();});
}
#endif

// Define a null delete that's used when std::shared_ptr is converted to boost::shared_ptr or vice-versa.
// The idea is that the new shared_ptr does not call delete when it goes out of scope since the original
// shared_ptr does that.

struct nullDeleter
{
	template <typename T>
	void operator()(T*& ptr) {}
};

si::si() : itsBottomLevel(kHPMissingInt)
{
	itsClearTextFormula = "<multiple algorithms>";

	itsLogger = unique_ptr<logger> (logger_factory::Instance()->GetLog("si"));
}

void si::Process(std::shared_ptr<const plugin_configuration> conf)
{
	compiled_plugin_base::Init(conf);

	/*
	 * Set target parameters:
	 * - name 
	 * - univ_id 
	 * - grib2 descriptor 0'00'000
	 *
	 */

	vector<param> theParams;

	theParams.push_back(param("DUMMY"));

	// GRIB 1

	SetParams(theParams);

	shared_ptr<neons> theNeons = dynamic_pointer_cast <neons> (plugin_factory::Instance()->Plugin("neons"));

	itsBottomLevel = boost::lexical_cast<int> (theNeons->ProducerMetaData(itsConfiguration->SourceProducer().Id(), "last hybrid level number"));
	itsTopLevel = boost::lexical_cast<int> (theNeons->ProducerMetaData(itsConfiguration->SourceProducer().Id(), "first hybrid level number"));

	Start();
	
}

/*
 * Calculate()
 *
 * This function does the actual calculation.
 */

void si::Calculate(shared_ptr<info> myTargetInfo, unsigned short threadIndex)
{
	auto f = dynamic_pointer_cast <fetcher> (plugin_factory::Instance()->Plugin("fetcher"));
	auto q = dynamic_pointer_cast <querydata> (plugin_factory::Instance()->Plugin("querydata"));

	// Required source parameters

	const param TParam("T-K");
	const params PParam({param("P-HPA"), param("P-PA")});
	const param RHParam("RH-PRCNT");
	const param HParam("HL-M");
	const param FFParam("FF-MS");

	unique_ptr<logger> myThreadedLogger = std::unique_ptr<logger> (logger_factory::Instance()->GetLog("siThread #" + boost::lexical_cast<string> (threadIndex)));

	ResetNonLeadingDimension(myTargetInfo);

	myTargetInfo->FirstParam();

	while (AdjustNonLeadingDimension(myTargetInfo))
	{

		myThreadedLogger->Debug("Calculating time " + static_cast<string> (myTargetInfo->Time().ValidDateTime()) +
								" level " + static_cast<string> (myTargetInfo->Level()));

		// Source infos

		shared_ptr<info> sourceInfo;

		bool haveData = true;

		for (int levelNumber = itsTopLevel; haveData && levelNumber <= itsBottomLevel; levelNumber++)
		{
			level curLevel(kHybrid, levelNumber, "HYBRID");
			
			shared_ptr<info> tempInfo;
			
			try
			{

				// Temperature

				tempInfo = f->Fetch(itsConfiguration,
									 myTargetInfo->Time(),
									 curLevel,
									 TParam);

				assert(tempInfo->Param().Unit() == kK);

				// grib-plugin does not set universal id number since it does
				// not know anything about it, but we need it in smarttools-library

				param p(tempInfo->Param());
				p.UnivId(4);
				tempInfo->SetParam(p);

				ScaleBase(tempInfo, 1, -constants::kKelvin);

				if (!sourceInfo)
				{
					sourceInfo = tempInfo;
				}
				else
				{
					sourceInfo->Merge(tempInfo);
				}

				tempInfo.reset();

				// Humidity
				
				tempInfo = f->Fetch(itsConfiguration,
									 myTargetInfo->Time(),
									 curLevel,
									 RHParam);

				p = param(tempInfo->Param());
				p.UnivId(13);
				tempInfo->SetParam(p);

				sourceInfo->Merge(tempInfo);

				tempInfo.reset();

				// Pressure

				tempInfo = f->Fetch(itsConfiguration,
									 myTargetInfo->Time(),
									 curLevel,
									 PParam);

				p = param(tempInfo->Param());
				p.UnivId(1);
				tempInfo->SetParam(p);
				
				if (tempInfo->Param().Name() == "P-PA")
				{
					ScaleBase(tempInfo, 0.01, 0);
				}
				
				sourceInfo->Merge(tempInfo);

				tempInfo.reset();

				// Height

				tempInfo = f->Fetch(itsConfiguration,
									 myTargetInfo->Time(),
									 curLevel,
									 HParam);

				p = param(tempInfo->Param());
				p.UnivId(3);
				tempInfo->SetParam(p);

				sourceInfo->Merge(tempInfo);

				tempInfo.reset();

				// Wind speed

				tempInfo = f->Fetch(itsConfiguration,
									 myTargetInfo->Time(),
									 curLevel,
									 FFParam);

				p = param(tempInfo->Param());
				p.UnivId(21);
				tempInfo->SetParam(p);

				sourceInfo->Merge(tempInfo);

				tempInfo.reset();

			}
			catch (HPExceptionType e)
			{
				switch (e)
				{
					case kFileDataNotFound:
						itsLogger->Warning("Skipping step " + boost::lexical_cast<string> (myTargetInfo->Time().Step()) + ", level " + boost::lexical_cast<string> (myTargetInfo->Level().Value()));
						myTargetInfo->Data().Fill(kFloatMissing);

						if (itsConfiguration->StatisticsEnabled())
						{
							itsConfiguration->Statistics()->AddToMissingCount(myTargetInfo->Grid()->Size());
							itsConfiguration->Statistics()->AddToValueCount(myTargetInfo->Grid()->Size());
						}
						haveData = false;
						continue;
						break;

					default:
						throw runtime_error(ClassName() + ": Unable to proceed");
						break;
				}
			}
		}

		if (!haveData)
		{
			break;
		}
		
		size_t missingCount = 0;
		size_t count = 0;

		// data read from neons does not have correct fmi producer id, copy producer
		// info from target info
		
		sourceInfo->Producer(myTargetInfo->Producer());

		// info: convert to querydata

		shared_ptr<NFmiQueryData> qdata = q->CreateQueryData(*sourceInfo, false);

#ifdef DEBUG
		ofstream in("indata.fqd");
		in << *qdata;
#endif

		// querydata: std::shared_ptr to boost::shared_ptr
		// boost::shared_ptr<NFmiQueryData> bqdata = make_shared_ptr(qdata);
		boost::shared_ptr<NFmiQueryData> bqdata (qdata.get(), nullDeleter());
		
		myThreadedLogger->Info("Calculating sounding index");
		
		// got boost::shared_ptr
		boost::shared_ptr<NFmiQueryData> bsidata = NFmiSoundingIndexCalculator::CreateNewSoundingIndexData(bqdata, "ASDF", false, 0);

#ifdef DEBUG
		ofstream out("outdata.fqd");
		out << *bsidata;
#endif

		// querydata: boost::shared_ptr to std::shared_ptr
		//shared_ptr<NFmiQueryData> sidata = make_shared_ptr(bsidata);
		shared_ptr<NFmiQueryData> sidata (bsidata.get(), nullDeleter());
		
		// querydata: convert to info
		myTargetInfo = q->CreateInfo(sidata);

		// Set correct target level to output data
		
		myTargetInfo->FirstLevel();
		level l(myTargetInfo->Level());
		l.Type(kHeight);

		myTargetInfo->SetLevel(l);

		string deviceType = "CPU";

		/*
		 * Newbase normalizes scanning mode to bottom left -- if that's not what
		 * the target scanning mode is, we have to swap the data back.
		 */

		for (myTargetInfo->ResetParam(); myTargetInfo->NextParam(); )
		{
			SwapTo(myTargetInfo, kBottomLeft);

			param p(myTargetInfo->Param());
			
			switch (p.UnivId())
			{
				case 4720:
					p.Name("LCL-HPA");
					break;
				
				case 4721:
					p.Name("LFC-HPA");
					break;

				case 4722:
					p.Name("EL-HPA");
					break;

				case 4723:
					p.Name("CAPE-JKG");
					break;

				case 4724:
					p.Name("CAPE-0-3");
					break;

				case 4725:
					p.Name("CIN-N");
					break;
					
				case 4726:
					p.Name("LCL-M");
					break;

				case 4727:
					p.Name("LFC-M");
					break;

				case 4728:
					p.Name("EL-M");
					break;

				case 4729:
					p.Name("CAPE1040");
					break;

				case 4730:
					p.Name("LCL-500-HPA");
					break;

				case 4731:
					p.Name("LFC-500-HPA");
					break;

				case 4732:
					p.Name("EL-500-HPA");
					break;

				case 4733:
					p.Name("CAPE-500");
					break;

				case 4734:
					p.Name("CAPE-0-3-500");
					break;

				case 4735:
					p.Name("CIN-500-N");
					break;

				case 4736:
					p.Name("LCL-500-M");
					break;

				case 4737:
					p.Name("LFC-500-M");
					break;

				case 4738:
					p.Name("LFC-500-M");
					break;

				case 4739:
					p.Name("CAPE1040-500");
					break;

				case 4740:
					p.Name("LCL-MU-HPA");
					break;

				case 4741:
					p.Name("LFC-MU-HPA");
					break;

				case 4742:
					p.Name("EL-MU-HPA");
					break;

				case 4743:
					p.Name("CAPE-MU-JKG");
					break;

				case 4744:
					p.Name("CAPE-0-3-MU");
					break;

				case 4745:
					p.Name("CIN-MU-N");
					break;

				case 4746:
					p.Name("LCL-MU-M");
					break;

				case 4747:
					p.Name("LFC-MU-M");
					break;

				case 4748:
					p.Name("EL-MU-M");
					break;

				case 4749:
					p.Name("CAPE1040-MU");
					break;
/*
				case 4750:
					p.Name("SI-N");
					break;

				case 4751:
					p.Name("LI-N");
					break;

				case 4752:
					p.Name("KINDEX-N");
					break;

				case 4753:
					p.Name("CTI-N");
					break;

				case 4754:
					p.Name("VTI-N");
					break;

				case 4755:
					p.Name("TTI-N");
					break;

				case 4770:
					p.Name("WSH-KT");
					break;

				case 4771:
					p.Name("WSH-1-KT");
					break;
*/
				case 4772:
					p.Name("HLCY-M2S2");
					break;

				case 4773:
					p.Name("HLCY-1-M2S2");
					break;

				case 4774:
					p.Name("FF1500-MS");
					break;

				case 4775:
					p.Name("TPE3-C");
					break;
					
				default:
					throw runtime_error("Unkown sounding parameter calculated: " + p.Name());
					break;
			}

			myTargetInfo->SetParam(p);
			
			for (myTargetInfo->ResetLocation(); myTargetInfo->NextLocation();)
			{
				count++;

				if (myTargetInfo->Value() == kFloatMissing)
				{
					missingCount++;
				}
			}
		}

		if (itsConfiguration->StatisticsEnabled())
		{
			itsConfiguration->Statistics()->AddToMissingCount(missingCount);
			itsConfiguration->Statistics()->AddToValueCount(count);
		}

		/*
		 * Now we are done for this level
		 *
		 * Clone info-instance to writer since it might change our descriptor places
		 */

		myThreadedLogger->Info("[" + deviceType + "] Missing values: " + boost::lexical_cast<string> (missingCount) + "/" + boost::lexical_cast<string> (count));

		if (itsConfiguration->FileWriteOption() != kSingleFile)
		{
			WriteToFile(*myTargetInfo);
		}

	}
}

void si::ScaleBase(shared_ptr<info> anInfo, double scale, double base)
{
	if (scale == 1 && base == 0)
	{
		return;
	}
	
	for (anInfo->ResetLocation(); anInfo->NextLocation(); )
	{
		double v = anInfo->Value();
		anInfo->Value(v * scale + base);
	}
}
#if 0
void si::LCL()
{
	// Case 1: Surface T, TD, P

	// Case 2,3,4: Average over wanted layer
}

void si::LCLAverage(shared_ptr<info> myTargetInfo, double fromZ, double toZ)
{

	// Fetch Z uncompressed since it is not transferred to cuda

	auto f = dynamic_pointer_cast <fetcher> (plugin_factory::Instance()->Plugin("fetcher"));

	auto HInfo = f->Fetch(itsConfiguration,
			myTargetInfo->Time(),
			level(kHeight, 0),
			param("Z-M2S2"),
			false);

	/*
	 * Find out the pressure of 0 meter height and 500 meter height. Note that
	 * the heights are adjusted with topography.
	 * 
	 * Limit search range to 0 .. 700 meters
	 */

	vector<double> lowerHeight(HInfo->SizeLocations(), 0);
	vector<double> upperHeight(HInfo->SizeLocations(), 700);

	vector<double> H0mVector = HInfo->Grid()->Data().Values();
	vector<double> H500mVector(HInfo->SizeLocations());

	for (size_t i = 0; i < H500mVector.size(); i++)
	{
		// H0mVector contains the height of ground (compared to MSL). Height can be negative
		// (maybe even in real life (Netherlands?)), but in our case we use 0 as smallest height.
		// TODO: check how it is in smarttools

		H0mVector[i] *= constants::kIg;

		if (H0mVector[i] < 0)
		{
			H0mVector[i] = 0;
		}

		H500mVector[i] = H0mVector[i] + 500.;
	}

	// 1. Get pressure values corresponding to fromZ and toZ

	auto h = dynamic_pointer_cast <hitool> (plugin_factory::Instance()->Plugin("hitool"));

	h->Configuration(itsConfiguration);
	h->Time(myTargetInfo->Time());

	params PParam = { param("P-PA"), param("P-HPA") };
	
	vector<double> fromP = h->VerticalHeight(PParam, lowerHeight, upperHeight, H0mVector);
	vector<double> toP = h->VerticalHeight(PParam, lowerHeight, upperHeight, H500mVector);

	// Integrate from 'fromP' to 'toP' using 100 Pa intervals

	vector<size_t> counts(fromP.size(), 0);
	vector<double> T_avg(fromP.size(), 0);

	do
	{
		auto T = h->VerticalValue(param("T-K"), fromP);

		assert(T.size() == T_avg.size());

		for (size_t i = 0; i < T.size(); i++)
		{
			if (T[i] != kFloatMissing)
			{
				T_avg[i] += T[i];
				counts[i] += 1;
			}

			if (fromP[i] == kFloatMissing)
			{
				continue;
			}
			
			fromP[i] -= 100;

			if (fromP[i] < toP[i])
			{
				fromP[i] = kFloatMissing;
			}
		}

		// for_each(fromP.begin(), fromP.end(), [](double& d) { d -= 100;});
		
	} while (true); //(!Finished(fromP));


}
#endif
