/**
 * @file tpot.cpp
 *
 * @brief Plugin to calculate potential temperature
 *
 * Created on: Nov 20, 2012
 * @author partio
 */

#include "tpot.h"
#include <iostream>
#include "plugin_factory.h"
#include "logger_factory.h"
#include <boost/lexical_cast.hpp>

#define HIMAN_AUXILIARY_INCLUDE

#include "fetcher.h"
#include "writer.h"
#include "neons.h"
#include "pcuda.h"

#undef HIMAN_AUXILIARY_INCLUDE

#ifdef DEBUG
#include "timer_factory.h"
#endif

using namespace std;
using namespace himan::plugin;

#ifdef HAVE_CUDA
namespace himan
{
namespace plugin
{
namespace tpot_cuda
{
void doCuda(const float* Tin, float TBase, const float* Pin, float TScale, float* TPout, size_t N, float PConst, unsigned short index);
}
}
}
#endif

const unsigned int MAX_THREADS = 1; // Max number of threads we allow

tpot::tpot() : itsUseCuda(false)
{
    itsClearTextFormula = "Tp = Tk * pow((1000/P), 0.286)"; // Poissons equation

    itsLogger = std::unique_ptr<logger> (logger_factory::Instance()->GetLog("tpot"));

}

void tpot::Process(shared_ptr<configuration> theConfiguration)
{

    shared_ptr<plugin::pcuda> c = dynamic_pointer_cast<plugin::pcuda> (plugin_factory::Instance()->Plugin("pcuda"));

    if (c->HaveCuda())
    {
        string msg = "I possess the powers of CUDA ";

        if (!theConfiguration->UseCuda())
        {
            msg += ", but I won't use them";
        }
        else
        {
            msg += ", and I'm not afraid to use them";
            itsUseCuda = true;
        }

        itsLogger->Info(msg);

    }

    // Get number of threads to use

    unsigned int theCoreCount = boost::thread::hardware_concurrency(); // Number of cores

    unsigned int theThreadCount = theCoreCount > MAX_THREADS ? MAX_THREADS : theCoreCount;

    boost::thread_group g;

    /*
     * The target information is parsed from the configuration file.
     */

    shared_ptr<info> theTargetInfo = theConfiguration->Info();

    /*
     * Get producer information from neons if whole_file_write is false.
     */

    if (!theConfiguration->WholeFileWrite())
    {
        shared_ptr<plugin::neons> n = dynamic_pointer_cast<plugin::neons> (plugin_factory::Instance()->Plugin("neons"));

        map<string,string> prodInfo = n->ProducerInfo(theTargetInfo->Producer().Id());

        if (prodInfo.size())
        {
            producer prod(theTargetInfo->Producer().Id());

            prod.Process(boost::lexical_cast<long> (prodInfo["process"]));
            prod.Centre(boost::lexical_cast<long> (prodInfo["centre"]));
            prod.Name(prodInfo["name"]);

            theTargetInfo->Producer(prod);
        }

    }

    /*
     * Set target parameter to potential temperature
     * - name TP-K
     * - univ_id 8
     * - grib2 descriptor 0'00'002
     *
     * We need to specify grib and querydata parameter information
     * since we don't know which one will be the output format.
     * (todo: we could check from theConfiguration but why bother?)
     *
     */

    vector<param> theParams;

    param theRequestedParam ("TP-K", 8);

    theRequestedParam.GribDiscipline(0);
    theRequestedParam.GribCategory(0);
    theRequestedParam.GribParameter(2);

    theParams.push_back(theRequestedParam);

    theTargetInfo->Params(theParams);

    /*
     * Create data structures.
     */

    theTargetInfo->Create();

    /*
     * Initialize thread manager
     */

    itsThreadManager = shared_ptr<util::thread_manager> (new util::thread_manager());

    itsThreadManager->Dimension(theConfiguration->LeadingDimension());
    itsThreadManager->FeederInfo(theTargetInfo->Clone());
    itsThreadManager->FeederInfo()->Param(theRequestedParam);

    /*
     * Each thread will have a copy of the target info.
     */

    vector<shared_ptr<info> > theTargetInfos;

    theTargetInfos.resize(theThreadCount);

    for (size_t i = 0; i < theThreadCount; i++)
    {

        itsLogger->Info("Thread " + boost::lexical_cast<string> (i + 1) + " starting");

        theTargetInfos[i] = theTargetInfo->Clone();

        boost::thread* t = new boost::thread(&tpot::Run,
                                             this,
                                             theTargetInfos[i],
                                             theConfiguration,
                                             i + 1);

        g.add_thread(t);

    }

    g.join_all();

    if (theConfiguration->WholeFileWrite())
    {

        shared_ptr<writer> theWriter = dynamic_pointer_cast <writer> (plugin_factory::Instance()->Plugin("writer"));

        theTargetInfo->FirstTime();

        string theOutputFile = "himan_" + theTargetInfo->Param().Name() + "_" + theTargetInfo->Time().OriginDateTime()->String("%Y%m%d%H");
        theWriter->ToFile(theTargetInfo, theConfiguration->OutputFileType(), false, theOutputFile);

    }
}

void tpot::Run(shared_ptr<info> myTargetInfo, shared_ptr<const configuration> theConfiguration, unsigned short theThreadIndex)
{

    while (itsThreadManager->AdjustLeadingDimension(myTargetInfo))
    {
        Calculate(myTargetInfo, theConfiguration, theThreadIndex);
    }

}

/*
 * Calculate()
 *
 * This function does the actual calculation.
 */

void tpot::Calculate(shared_ptr<info> myTargetInfo, shared_ptr<const configuration> theConfiguration, unsigned short theThreadIndex)
{

    shared_ptr<fetcher> theFetcher = dynamic_pointer_cast <fetcher> (plugin_factory::Instance()->Plugin("fetcher"));

    // Required source parameters

    param TParam ("T-K");
    param PParam ("P-HPA");

    unique_ptr<logger> myThreadedLogger = std::unique_ptr<logger> (logger_factory::Instance()->GetLog("tpotThread #" + boost::lexical_cast<string> (theThreadIndex)));

    itsThreadManager->ResetNonLeadingDimension(myTargetInfo);

    myTargetInfo->FirstParam();

    while (itsThreadManager->AdjustNonLeadingDimension(myTargetInfo))
    {

        myThreadedLogger->Debug("Calculating time " + myTargetInfo->Time().ValidDateTime()->String("%Y%m%d%H") +
                                " level " + boost::lexical_cast<string> (myTargetInfo->Level().Value()));

        myTargetInfo->Data()->Resize(theConfiguration->Ni(), theConfiguration->Nj());

        double PScale = 1;
        double TBase = 0;

#ifdef DEBUG
        unique_ptr<timer> t = unique_ptr<timer> (timer_factory::Instance()->GetTimer());
        t->Start();
#endif

        // Source info for T
        shared_ptr<info> TInfo = theFetcher->Fetch(theConfiguration,
                                 myTargetInfo->Time(),
                                 myTargetInfo->Level(),
                                 TParam);

        // Source info for P
        shared_ptr<info> PInfo;
        shared_ptr<NFmiGrid> PGrid;

        bool isPressureLevel = (myTargetInfo->Level().Type() == kPressure);

        if (!isPressureLevel)
        {
            // Source info for P
            PInfo = theFetcher->Fetch(theConfiguration,
                                      myTargetInfo->Time(),
                                      myTargetInfo->Level(),
                                      PParam);

            if (PInfo->Param().Unit() == kPa)
            {
                PScale = 0.01;
            }

            PGrid = PInfo->ToNewbaseGrid();
        }

#ifdef DEBUG
        t->Stop();
        itsLogger->Debug("Data fetching took " + boost::lexical_cast<string> (t->GetTime()) + " microseconds");
#endif

        if (TInfo->Param().Unit() == kC)
        {
            TBase = 273.15;
        }

        shared_ptr<NFmiGrid> targetGrid = myTargetInfo->ToNewbaseGrid();
        shared_ptr<NFmiGrid> TGrid = TInfo->ToNewbaseGrid();

        int missingCount = 0;
        int count = 0;

        assert(targetGrid->Size() == myTargetInfo->Data()->Size());

#ifdef HAVE_CUDA

        bool equalGrids = (myTargetInfo->GridAndAreaEquals(TInfo) && (isPressureLevel || myTargetInfo->GridAndAreaEquals(PInfo)));

        //if (itsUseCuda && equalGrids)
        if (itsUseCuda && equalGrids && theThreadIndex == 1)
        {
#ifdef DEBUG
            t->Start();
#endif
            size_t N = TGrid->Size();

            float* TPData = new float[N];

            if (!isPressureLevel)
            {
                tpot_cuda::doCuda(TGrid->DataPool()->Data(), TBase, PGrid->DataPool()->Data(), PScale, TPData, N, 0, theThreadIndex-1);
            }
            else
            {
                tpot_cuda::doCuda(TGrid->DataPool()->Data(), TBase, 0, 0, TPData, N, myTargetInfo->Level().Value(), theThreadIndex-1);
            }

            double *data = new double[N];

            for (size_t i = 0; i < N; i++)
            {
                data[i] = static_cast<float> (TPData[i]);

                if (data[i] == kFloatMissing)
                {
                    missingCount++;
                }

                count++;
            }

            myTargetInfo->Data()->Data(data, N);

            delete [] data;
            delete [] TPData;

#ifdef DEBUG
            t->Stop();
            itsLogger->Debug("Calculation took " + boost::lexical_cast<string> (t->GetTime()) + " microseconds on GPU");
#endif
        }
        else
        {

#else
        if (true)
        {
#endif

            myTargetInfo->ResetLocation();

            targetGrid->Reset();

#ifdef DEBUG
            t->Start();
#endif

            while (myTargetInfo->NextLocation() && targetGrid->Next())
            {
                count++;

                NFmiPoint thePoint = targetGrid->LatLon();

                double T = kFloatMissing;
                double P = kFloatMissing;

                TGrid->InterpolateToLatLonPoint(thePoint, T);

                if (isPressureLevel)
                {
                    P = myTargetInfo->Level().Value();
                }
                else
                {
                    PGrid->InterpolateToLatLonPoint(thePoint, P);
                }

                if (T == kFloatMissing || P == kFloatMissing)
                {
                    missingCount++;

                    myTargetInfo->Value(kFloatMissing);
                    continue;
                }

                double Tp = (T + TBase) * pow((1000 / (P * PScale)), 0.286);

                if (!myTargetInfo->Value(Tp))
                {
                    throw runtime_error(ClassName() + ": Failed to set value to matrix");
                }
            }

#ifdef DEBUG
            t->Stop();
            itsLogger->Debug("Calculation took " + boost::lexical_cast<string> (t->GetTime()) + " microseconds on CPU");
#endif

        }


        /*
         * Now we are done for this level
         *
         * Clone info-instance to writer since it might change our descriptor places
         * */

        myThreadedLogger->Info("Missing values: " + boost::lexical_cast<string> (missingCount) + "/" + boost::lexical_cast<string> (count));

        if (!theConfiguration->WholeFileWrite())
        {
            shared_ptr<writer> theWriter = dynamic_pointer_cast <writer> (plugin_factory::Instance()->Plugin("writer"));

            theWriter->ToFile(myTargetInfo->Clone(), theConfiguration->OutputFileType(), true);
        }
    }
}
