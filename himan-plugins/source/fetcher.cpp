/**
 * @file fetcher.cpp
 *
 * @date Nov 21, 2012
 * @author partio
 */

#include "fetcher.h"
#include "interpolate.h"
#include "logger_factory.h"
#include "plugin_factory.h"
#include "timer_factory.h"
#include "util.h"
#include <boost/filesystem/operations.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>

#include "cache.h"
#include "csv.h"
#include "grib.h"
#include "neons.h"
#include "param.h"
#include "querydata.h"
#include "radon.h"

using namespace himan::plugin;
using namespace std;

fetcher::fetcher()
    : itsDoLevelTransform(true),
      itsDoInterpolation(true),
      itsUseCache(true),
      itsApplyLandSeaMask(false),
      itsLandSeaMaskThreshold(0.5)
{
	itsLogger = std::unique_ptr<logger>(logger_factory::Instance()->GetLog("fetcher"));
}

shared_ptr<himan::info> fetcher::Fetch(shared_ptr<const plugin_configuration> config, forecast_time requestedTime,
                                       level requestedLevel, const params& requestedParams, forecast_type requestedType,
                                       bool readPackedData)
{
	shared_ptr<info> ret;

	for (size_t i = 0; i < requestedParams.size(); i++)
	{
		try
		{
			ret = Fetch(config, requestedTime, requestedLevel, requestedParams[i], requestedType, readPackedData, true);

			return ret;
		}
		catch (const HPExceptionType& e)
		{
			if (e != kFileDataNotFound)
			{
				throw;
			}
		}
	}

	string optsStr = "producer(s): ";

	for (size_t prodNum = 0; prodNum < config->SizeSourceProducers(); prodNum++)
	{
		optsStr += boost::lexical_cast<string>(config->SourceProducer(prodNum).Id()) + ",";
	}

	optsStr = optsStr.substr(0, optsStr.size() - 1);

	optsStr += " origintime: " + requestedTime.OriginDateTime().String() + ", step: " +
	           boost::lexical_cast<string>(requestedTime.Step());

	optsStr += " param(s): ";

	for (size_t i = 0; i < requestedParams.size(); i++)
	{
		optsStr += requestedParams[i].Name() + ",";
	}

	optsStr = optsStr.substr(0, optsStr.size() - 1);

	optsStr += " level: " + string(himan::HPLevelTypeToString.at(requestedLevel.Type())) + " " +
	           boost::lexical_cast<string>(requestedLevel.Value());

	if (static_cast<int>(requestedType.Type()) > 2)
	{
		optsStr += " forecast type: " + string(himan::HPForecastTypeToString.at(requestedType.Type())) + "/" +
		           boost::lexical_cast<string>(requestedType.Value());
	}

	itsLogger->Warning("No valid data found with given search options " + optsStr);

	throw kFileDataNotFound;
}

shared_ptr<himan::info> fetcher::Fetch(shared_ptr<const plugin_configuration> config, forecast_time requestedTime,
                                       level requestedLevel, param requestedParam, forecast_type requestedType,
                                       bool readPackedData, bool suppressLogging)
{
	unique_ptr<timer> t = unique_ptr<timer>(timer_factory::Instance()->GetTimer());

	if (config->StatisticsEnabled())
	{
		t->Start();
	}

	vector<shared_ptr<info>> theInfos;

	level newLevel = requestedLevel;

	for (size_t prodNum = 0; prodNum < config->SizeSourceProducers() && theInfos.empty(); prodNum++)
	{
		producer sourceProd(config->SourceProducer(prodNum));

		if (itsDoLevelTransform && (requestedLevel.Type() != kHybrid && requestedLevel.Type() != kPressure))
		{
			newLevel = LevelTransform(sourceProd, requestedParam, requestedLevel);

			if (newLevel != requestedLevel)
			{
				itsLogger->Trace("Transform level " + string(HPLevelTypeToString.at(requestedLevel.Type())) + "/" +
				                 boost::lexical_cast<string>(requestedLevel.Value()) + " to " +
				                 HPLevelTypeToString.at(newLevel.Type()) + "/" +
				                 boost::lexical_cast<string>(newLevel.Value()) + " for producer " +
				                 boost::lexical_cast<string>(sourceProd.Id()) + ", parameter " + requestedParam.Name());
			}
		}

		search_options opts(requestedTime, requestedParam, newLevel, sourceProd, requestedType, config);

		theInfos = FetchFromProducer(opts, readPackedData);
	}

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
		if (!suppressLogging)
		{
			string optsStr = "producer(s): ";

			for (size_t prodNum = 0; prodNum < config->SizeSourceProducers(); prodNum++)
			{
				optsStr += boost::lexical_cast<string>(config->SourceProducer(prodNum).Id()) + ",";
			}

			optsStr = optsStr.substr(0, optsStr.size() - 1);

			optsStr += " origintime: " + requestedTime.OriginDateTime().String() + ", step: " +
			           boost::lexical_cast<string>(requestedTime.Step());
			optsStr += " param: " + requestedParam.Name();
			optsStr += " level: " + string(himan::HPLevelTypeToString.at(requestedLevel.Type())) + " " +
			           boost::lexical_cast<string>(requestedLevel.Value());

			if (static_cast<int>(requestedType.Type()) > 2)
			{
				optsStr += " forecast type: " + string(himan::HPForecastTypeToString.at(requestedType.Type())) + "/" +
				           boost::lexical_cast<string>(requestedType.Value());
			}

			itsLogger->Warning("No valid data found with given search options " + optsStr);
		}

		throw kFileDataNotFound;
	}

	// assert(theConfiguration->SourceProducer() == theInfos[0]->Producer());

	assert((theInfos[0]->Level()) == newLevel);

	assert((theInfos[0]->Time()) == requestedTime);

	assert((theInfos[0]->Param()) == requestedParam);

	auto baseInfo = make_shared<info>(*config->Info());
	assert(baseInfo->Dimensions().size());

	baseInfo->First();

	if (itsDoInterpolation)
	{
		if (!interpolate::Interpolate(*baseInfo, theInfos, config->UseCudaForInterpolation()))
		{
			// interpolation failed
			throw kFileDataNotFound;
		}
	}
	else
	{
		itsLogger->Trace("Interpolation disabled");
	}

	if (itsApplyLandSeaMask)
	{
		itsLogger->Trace("Applying land-sea mask with threshold " +
		                 boost::lexical_cast<string>(itsLandSeaMaskThreshold));

		itsApplyLandSeaMask = false;

		if (!ApplyLandSeaMask(config, *theInfos[0], requestedTime, requestedType))
		{
			itsLogger->Warning("Land sea mask apply failed");
		}

		itsApplyLandSeaMask = true;
	}

	/*
	 * Insert interpolated data to cache if
	 * 1. Cache is not disabled locally (itsUseCache) AND
	 * 2. Cache is not disabled globally (config->UseCache()) AND
	 * 3. Data is not packed
	 */

	if (itsUseCache && config->UseCache() && !theInfos[0]->Grid()->IsPackedData())
	{
		auto c = GET_PLUGIN(cache);
		c->Insert(*theInfos[0]);
	}

	baseInfo.reset();

	return theInfos[0];
}

vector<shared_ptr<himan::info>> fetcher::FromFile(const vector<string>& files, search_options& options,
                                                  bool readContents, bool readPackedData, bool readIfNotMatching)
{
	vector<shared_ptr<himan::info>> allInfos;

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
				curInfos = FromGrib(inputFile, options, readContents, readPackedData, readIfNotMatching);
				break;
			}

			case kQueryData:
			{
				curInfos = FromQueryData(inputFile, options, readContents);
				break;
			}

			case kNetCDF:
				itsLogger->Error("File is NetCDF");
				break;

			case kCSV:
			{
				curInfos = FromCSV(inputFile, options);
				break;
			}

			default:
				// Unknown file type, cannot proceed
				throw runtime_error("Input file is neither GRIB, NetCDF, QueryData nor CSV");
				break;
		}

		allInfos.insert(allInfos.end(), curInfos.begin(), curInfos.end());
	}

	return allInfos;
}

vector<shared_ptr<himan::info>> fetcher::FromCache(search_options& options)
{
	auto c = GET_PLUGIN(cache);

	vector<shared_ptr<himan::info>> infos = c->GetInfo(options);

	return infos;
}

vector<shared_ptr<himan::info>> fetcher::FromGrib(const string& inputFile, search_options& options, bool readContents,
                                                  bool readPackedData, bool forceCaching)
{
	auto g = GET_PLUGIN(grib);

	vector<shared_ptr<info>> infos = g->FromFile(inputFile, options, readContents, readPackedData, forceCaching);

	return infos;
}

vector<shared_ptr<himan::info>> fetcher::FromQueryData(const string& inputFile, search_options& options,
                                                       bool readContents)
{
	auto q = GET_PLUGIN(querydata);

	shared_ptr<info> i = q->FromFile(inputFile, options);

	vector<shared_ptr<info>> theInfos;

	theInfos.push_back(i);

	return theInfos;
}

vector<shared_ptr<himan::info>> fetcher::FromCSV(const string& inputFile, search_options& options)
{
	auto c = GET_PLUGIN(csv);

	auto info = c->FromFile(inputFile, options);

	vector<info_t> infos;
	infos.push_back(info);

	return infos;
}

himan::level fetcher::LevelTransform(const producer& sourceProducer, const param& targetParam,
                                     const level& targetLevel) const
{
	level sourceLevel = targetLevel;

	if (sourceProducer.TableVersion() != kHPMissingInt)
	{
		auto n = GET_PLUGIN(neons);

		string lvlName =
		    n->NeonsDB().GetGridLevelName(targetParam.Name(), targetLevel.Type(), 204, sourceProducer.TableVersion());

		if (lvlName.empty())
		{
			itsLogger->Trace("No level transformation found for param " + targetParam.Name() + " level " +
			                 HPLevelTypeToString.at(targetLevel.Type()));
			return targetLevel;
		}

		double lvlValue = targetLevel.Value();

		HPLevelType lvlType = HPStringToLevelType.at(boost::to_lower_copy(lvlName));

		if (lvlType == kGround)
		{
			lvlValue = 0;
		}

		sourceLevel = level(lvlType, lvlValue, lvlName);
	}
	else
	{
		sourceLevel = targetLevel;
	}

	return sourceLevel;
}

void fetcher::DoLevelTransform(bool theDoLevelTransform) { itsDoLevelTransform = theDoLevelTransform; }
bool fetcher::DoLevelTransform() const { return itsDoLevelTransform; }
void fetcher::DoInterpolation(bool theDoInterpolation) { itsDoInterpolation = theDoInterpolation; }
bool fetcher::DoInterpolation() const { return itsDoInterpolation; }
void fetcher::UseCache(bool theUseCache) { itsUseCache = theUseCache; }
bool fetcher::UseCache() const { return itsUseCache; }
vector<shared_ptr<himan::info>> fetcher::FetchFromProducer(search_options& opts, bool readPackedData)
{
	vector<shared_ptr<info>> ret;

	itsLogger->Trace("Current producer: " + boost::lexical_cast<string>(opts.prod.Id()));

	// itsLogger->Trace("Current producer: " + sourceProd.Name());

	if (itsUseCache && opts.configuration->UseCache())
	{
		// 1. Fetch data from cache
		ret = FromCache(opts);

		if (ret.size())
		{
			itsLogger->Trace("Data found from cache");

			if (dynamic_pointer_cast<const plugin_configuration>(opts.configuration)->StatisticsEnabled())
			{
				dynamic_pointer_cast<const plugin_configuration>(opts.configuration)
				    ->Statistics()
				    ->AddToCacheHitCount(1);
			}

			return ret;
		}
	}

	/*
	 *  2. Fetch data from auxiliary files specified at command line
	 *
	 *  Even if file_wait_timeout is specified, auxiliary files is searched
	 *  only once.
	 */

	if (!opts.configuration->AuxiliaryFiles().empty())
	{
		ret = FromFile(opts.configuration->AuxiliaryFiles(), opts, true, readPackedData,
		               !itsApplyLandSeaMask && !readPackedData);

		if (!ret.empty())
		{
			itsLogger->Trace("Data found from auxiliary file(s)");

			if (dynamic_pointer_cast<const plugin_configuration>(opts.configuration)->StatisticsEnabled())
			{
				dynamic_pointer_cast<const plugin_configuration>(opts.configuration)
				    ->Statistics()
				    ->AddToCacheMissCount(1);
			}

			return ret;
		}
		else
		{
			itsLogger->Trace("Data not found from auxiliary file(s)");
		}
	}

	// 3. Fetch data from Neons or Radon

	vector<string> files;

	if (opts.configuration->ReadDataFromDatabase())
	{
		HPDatabaseType dbtype = opts.configuration->DatabaseType();

		if (dbtype == kNeons || dbtype == kNeonsAndRadon)
		{
			// try neons first
			auto n = GET_PLUGIN(neons);

			itsLogger->Trace("Accessing Neons database");

			files = n->Files(opts);

			if (!files.empty())
			{
				ret = FromFile(files, opts, true, readPackedData);

				if (dynamic_pointer_cast<const plugin_configuration>(opts.configuration)->StatisticsEnabled())
				{
					dynamic_pointer_cast<const plugin_configuration>(opts.configuration)
					    ->Statistics()
					    ->AddToCacheMissCount(1);
				}

				return ret;
			}
		}

		if (dbtype == kRadon || dbtype == kNeonsAndRadon)
		{
			// try radon next

			auto r = GET_PLUGIN(radon);

			itsLogger->Trace("Accessing Radon database");

			files = r->Files(opts);

			if (!files.empty())
			{
				ret = FromFile(files, opts, true, readPackedData);

				if (dynamic_pointer_cast<const plugin_configuration>(opts.configuration)->StatisticsEnabled())
				{
					dynamic_pointer_cast<const plugin_configuration>(opts.configuration)
					    ->Statistics()
					    ->AddToCacheMissCount(1);
				}

				return ret;
			}
		}

		const string ref_prod = opts.prod.Name();
		const string analtime = opts.time.OriginDateTime().String("%Y%m%d%H%M%S");
		const vector<string> sourceGeoms = opts.configuration->SourceGeomNames();
		itsLogger->Warning("No geometries found for producer " + ref_prod + ", analysistime " + analtime +
		                   ", source geom name(s) '" + util::Join(sourceGeoms, ",") + "', param " + opts.param.Name());
	}

	return ret;
}

bool fetcher::ApplyLandSeaMask(shared_ptr<const plugin_configuration> config, info& theInfo,
                               forecast_time& requestedTime, forecast_type& requestedType)
{
	raw_time originTime = requestedTime.OriginDateTime();
	forecast_time firstTime(originTime, originTime);

	try
	{
		itsApplyLandSeaMask = false;

		auto lsmInfo = Fetch(config, firstTime, level(kHeight, 0), param("LC-0TO1"), requestedType, false);

		itsApplyLandSeaMask = true;

		lsmInfo->First();

		assert(*lsmInfo->Grid() == *theInfo.Grid());

		assert(itsLandSeaMaskThreshold >= -1 && itsLandSeaMaskThreshold <= 1);
		assert(itsLandSeaMaskThreshold != 0);

#ifdef HAVE_CUDA
		if (theInfo.Grid()->IsPackedData())
		{
			// We need to unpack
			util::Unpack({theInfo.Grid()});
		}
#endif

		assert(!theInfo.Grid()->IsPackedData());

		double multiplier = (itsLandSeaMaskThreshold > 0) ? 1. : -1.;

		for (lsmInfo->ResetLocation(), theInfo.ResetLocation(); lsmInfo->NextLocation() && theInfo.NextLocation();)
		{
			double lsm = lsmInfo->Value();

			if (theInfo.Value() == kFloatMissing || lsm == kFloatMissing)
			{
				continue;
			}

			if (multiplier * lsm <= itsLandSeaMaskThreshold)
			{
				theInfo.Value(kFloatMissing);
			}
		}
	}
	catch (HPExceptionType& e)
	{
		itsApplyLandSeaMask = true;
		return false;
	}

	return true;
}

bool fetcher::ApplyLandSeaMask() const { return itsApplyLandSeaMask; }
void fetcher::ApplyLandSeaMask(bool theApplyLandSeaMask) { itsApplyLandSeaMask = theApplyLandSeaMask; }
double fetcher::LandSeaMaskThreshold() const { return itsLandSeaMaskThreshold; }
void fetcher::LandSeaMaskThreshold(double theLandSeaMaskThreshold)
{
	if (theLandSeaMaskThreshold < -1 || theLandSeaMaskThreshold > 1)
	{
		itsLogger->Fatal("Invalid value for land sea mask threshold: " +
		                 boost::lexical_cast<string>(theLandSeaMaskThreshold));
		exit(1);
	}

	itsLandSeaMaskThreshold = theLandSeaMaskThreshold;
}
