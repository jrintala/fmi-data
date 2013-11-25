/**
 * @file fetcher.cpp
 *
 * @date Nov 21, 2012
 * @author partio
 */

#include "fetcher.h"
#include "plugin_factory.h"
#include "logger_factory.h"
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/operations.hpp>
#include "timer_factory.h"
#include "util.h"

#define HIMAN_AUXILIARY_INCLUDE

#include "grib.h"
#include "neons.h"
#include "param.h"
#include "cache.h"
#include "querydata.h"

#undef HIMAN_AUXILIARY_INCLUDE

using namespace himan::plugin;
using namespace std;

const unsigned int SLEEPSECONDS = 10;

shared_ptr<cache> itsCache;

fetcher::fetcher()
{
	itsLogger = std::unique_ptr<logger> (logger_factory::Instance()->GetLog("fetcher"));
}

shared_ptr<himan::info> fetcher::Fetch(shared_ptr<const plugin_configuration> config,
										const forecast_time& requestedTime,
										const level& requestedLevel,
										const params& requestedParams,
										bool readPackedData)
{
	unsigned int waitedSeconds = 0;

	shared_ptr<info> ret;
	
	do
	{
		for (size_t i = 0; i < requestedParams.size(); i++)
		{
			try
			{
				ret = Fetch(config, requestedTime, requestedLevel, requestedParams[i], readPackedData, false);
				
				return ret;
			}
			catch (const HPExceptionType& e)
			{
				if (e != kFileDataNotFound)
				{
					throw e;
				}
			}
			catch (const exception& e)
			{
				throw e;
			}

		}
		if (config->FileWaitTimeout() > 0)
		{
			itsLogger->Debug("Sleeping for " + boost::lexical_cast<string> (SLEEPSECONDS) + " seconds (cumulative: " + boost::lexical_cast<string> (waitedSeconds) + ")");

			if (!config->ReadDataFromDatabase())
			{
				itsLogger->Warning("file_wait_timeout specified but file read from Neons is disabled");
			}

			sleep(SLEEPSECONDS);
		}

		waitedSeconds += SLEEPSECONDS;
	}
	while (waitedSeconds < config->FileWaitTimeout() * 60);

	string optsStr = "producer(s): ";

	for (size_t prodNum = 0; prodNum < config->SizeSourceProducers(); prodNum++)
	{
		optsStr += boost::lexical_cast<string> (config->SourceProducer(prodNum).Id()) + ",";
	}

	optsStr = optsStr.substr(0, optsStr.size()-1);

	optsStr += " origintime: " + requestedTime.OriginDateTime()->String() + ", step: " + boost::lexical_cast<string> (requestedTime.Step());

	optsStr += " param(s): ";
	
	for (size_t i = 0; i < requestedParams.size(); i++)
	{
		optsStr += requestedParams[i].Name() + ",";
	}

	optsStr = optsStr.substr(0, optsStr.size()-1);

	optsStr += " level: " + string(himan::HPLevelTypeToString.at(requestedLevel.Type())) + " " + boost::lexical_cast<string> (requestedLevel.Value());

	itsLogger->Warning("No valid data found with given search options " + optsStr);

	throw kFileDataNotFound;

}


shared_ptr<himan::info> fetcher::Fetch(shared_ptr<const plugin_configuration> config,
										const forecast_time& requestedTime,
										const level& requestedLevel,
										const param& requestedParam,
										bool readPackedData,
										bool controlWaitTime)
{

	unique_ptr<timer> t = unique_ptr<timer> (timer_factory::Instance()->GetTimer());
	
	if (config->StatisticsEnabled())
	{
		t->Start();
	}
	
	vector<shared_ptr<info>> theInfos;
	unsigned int waitedSeconds = 0;
		
	do
	{
		// Loop over all source producers if more than one specified

		for (size_t prodNum = 0; prodNum < config->SizeSourceProducers(); prodNum++)
		{
		
			producer sourceProd(config->SourceProducer(prodNum));
			
			const search_options opts (requestedTime, requestedParam, requestedLevel, sourceProd, config);

			// itsLogger->Trace("Current producer: " + sourceProd.Name());

			if (config->UseCache())
			{

				// 1. Fetch data from cache

				if (!itsCache)
				{
					itsCache = dynamic_pointer_cast<plugin::cache> (plugin_factory::Instance()->Plugin("cache"));
				}

				theInfos = FromCache(opts);

				if (theInfos.size())
				{
					itsLogger->Trace("Data found from cache");

					if (config->StatisticsEnabled())
					{
						config->Statistics()->AddToCacheHitCount(1);
					}

					break;
				}
			}

			/*
			 *  2. Fetch data from auxiliary files specified at command line
			 *
			 *  Even if file_wait_timeout is specified, auxiliary files is searched
			 *  only once.
			 */

			if (config->AuxiliaryFiles().size() && waitedSeconds == 0)
			{
				theInfos = FromFile(config->AuxiliaryFiles(), opts, true, readPackedData);

				if (theInfos.size())
				{
					itsLogger->Trace("Data found from auxiliary file(s)");

					if (config->StatisticsEnabled())
					{
						config->Statistics()->AddToCacheMissCount(1);
					}

					break;
				}
				else
				{
					itsLogger->Trace("Data not found from auxiliary file(s)");
				}
			}

			// 3. Fetch data from Neons

			vector<string> files;

			if (config->ReadDataFromDatabase())
			{
				shared_ptr<neons> n = dynamic_pointer_cast<neons> (plugin_factory::Instance()->Plugin("neons"));

				files = n->Files(opts);

				if (!files.empty())
				{
					theInfos = FromFile(files, opts, true, readPackedData);

					if (config->StatisticsEnabled())
					{
						config->Statistics()->AddToCacheMissCount(1);
					}

					break;
				}
			}
		}

		if (controlWaitTime && config->FileWaitTimeout() > 0)
		{
			itsLogger->Debug("Sleeping for " + boost::lexical_cast<string> (SLEEPSECONDS) + " seconds (cumulative: " + boost::lexical_cast<string> (waitedSeconds) + ")");

			if (!config->ReadDataFromDatabase())
			{
				itsLogger->Warning("file_wait_timeout specified but file read from Neons is disabled");
			}

			sleep(SLEEPSECONDS);
		}

		waitedSeconds += SLEEPSECONDS;
	}
	while (controlWaitTime && waitedSeconds < config->FileWaitTimeout() * 60);

	if (config->StatisticsEnabled())
	{
		t->Stop();

		config->Statistics()->AddToFetchingTime(t->GetTime());
	}
	/*
	 *  Safeguard; later in the code we do not check whether the data requested
	 *  was actually what was requested.
	 */

	if (theInfos.size() == 0)
	{
		// If this function is called from multi-param Fetch(), do not print
		// any messages yet since we might have another source param coming

		if (controlWaitTime)
		{
			string optsStr = "producer(s): ";

			for (size_t prodNum = 0; prodNum < config->SizeSourceProducers(); prodNum++)
			{
				optsStr += boost::lexical_cast<string> (config->SourceProducer(prodNum).Id()) + ",";
			}

			optsStr = optsStr.substr(0, optsStr.size()-1);

			optsStr += " origintime: " + requestedTime.OriginDateTime()->String() + ", step: " + boost::lexical_cast<string> (requestedTime.Step());
			optsStr += " param: " + requestedParam.Name();
			optsStr += " level: " + string(himan::HPLevelTypeToString.at(requestedLevel.Type())) + " " + boost::lexical_cast<string> (requestedLevel.Value());

			itsLogger->Warning("No valid data found with given search options " + optsStr);
		}

		throw kFileDataNotFound;
	}

	// assert(theConfiguration->SourceProducer() == theInfos[0]->Producer());

	assert((theInfos[0]->Level()) == requestedLevel);

	assert((theInfos[0]->Time()) == requestedTime);

	assert((theInfos[0]->Param()) == requestedParam);

	return theInfos[0];

}

vector<shared_ptr<himan::info>> fetcher::FromFile(const vector<string>& files, const search_options& options, bool readContents, bool readPackedData)
{

	vector<shared_ptr<himan::info>> allInfos ;

	for (size_t i = 0; i < files.size(); i++)
	{
		string inputFile = files[i];

		if (!boost::filesystem::exists(inputFile))
		{
			itsLogger->Error("Input file '" + inputFile + "' does not exist");
			continue;
		}

		vector<shared_ptr<himan::info>> curInfos;

		switch (util::FileType(inputFile))
		{

		case kGRIB:
		case kGRIB1:
		case kGRIB2:
		{

			curInfos = FromGrib(inputFile, options, readContents, readPackedData);
			break;

		}

		case kQueryData:
		{

			curInfos = FromQueryData(inputFile, options, readContents);
			break;
		}

		case kNetCDF:
			cout << "File is NetCDF" << endl;
			break;

		default:
			// Unknown file type, cannot proceed
			throw runtime_error("Input file is neither GRID, NetCDF nor QueryData");
			break;
		}

		allInfos.insert(allInfos.end(), curInfos.begin(), curInfos.end());

		if (curInfos.size())
		{
			if (options.configuration->UseCache())
			{
				itsCache->Insert(allInfos);
			}

			break; // We found what we were looking for
		}

	}

	return allInfos;
}

vector<shared_ptr<himan::info> > fetcher::FromCache(const search_options& options)
{
	vector<shared_ptr<himan::info>> infos = itsCache->GetInfo(options);

	return infos;
}

vector<shared_ptr<himan::info> > fetcher::FromGrib(const string& inputFile, const search_options& options, bool readContents, bool readPackedData)
{

	shared_ptr<grib> g = dynamic_pointer_cast<grib> (plugin_factory::Instance()->Plugin("grib"));

	vector<shared_ptr<info>> infos = g->FromFile(inputFile, options, readContents, readPackedData);

	return infos;
}

vector<shared_ptr<himan::info>> fetcher::FromQueryData(const string& inputFile, const search_options& options, bool readContents)
{

	shared_ptr<querydata> q = dynamic_pointer_cast<querydata> (plugin_factory::Instance()->Plugin("querydata"));

	shared_ptr<info> i = q->FromFile(inputFile, options, readContents);

	vector<shared_ptr<info>> theInfos;

	theInfos.push_back(i);

	return theInfos;
}
