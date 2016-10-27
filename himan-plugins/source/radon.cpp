/**
 * @file radon.cpp
 *
 */

#include "radon.h"
#include "logger_factory.h"
#include "plugin_factory.h"
#include "util.h"
#include <sstream>
#include <thread>

using namespace std;
using namespace himan::plugin;

const int MAX_WORKERS = 32;
static once_flag oflag;

radon::radon() : itsInit(false), itsRadonDB()
{
	itsLogger = unique_ptr<logger>(logger_factory::Instance()->GetLog("radon"));

	call_once(oflag, [&]() {
		PoolMaxWorkers(MAX_WORKERS);

		NFmiRadonDBPool::Instance()->Username("wetodb");
		NFmiRadonDBPool::Instance()->Password("3loHRgdio");
	});
}

void radon::PoolMaxWorkers(int maxWorkers) { NFmiRadonDBPool::Instance()->MaxWorkers(maxWorkers); }
vector<string> radon::Files(search_options& options)
{
	Init();

	vector<string> files;

	string analtime = options.time.OriginDateTime().String("%Y-%m-%d %H:%M:%S+00");
	string levelValue = boost::lexical_cast<string>(options.level.Value());
	string levelValue2 = "-1";

	if (options.level.Value2() != kHPMissingValue)
	{
		levelValue2 = boost::lexical_cast<string>(options.level.Value2());
	}

	string ref_prod = options.prod.Name();
	// long no_vers = options.prod.TableVersion();

	string level_name = HPLevelTypeToString.at(options.level.Type());

	vector<vector<string>> gridgeoms;
	vector<string> sourceGeoms = options.configuration->SourceGeomNames();

	if (sourceGeoms.empty())
	{
		// Get all geometries
		gridgeoms = itsRadonDB->GetGridGeoms(ref_prod, analtime);
	}
	else
	{
		for (size_t i = 0; i < sourceGeoms.size(); i++)
		{
			vector<vector<string>> geoms = itsRadonDB->GetGridGeoms(ref_prod, analtime, sourceGeoms[i]);
			gridgeoms.insert(gridgeoms.end(), geoms.begin(), geoms.end());
		}
	}

	if (gridgeoms.empty())
	{
		// No geometries found, fetcher checks this
		return files;
	}

	string forecastTypeValue = "-1";  // default, deterministic/analysis

	if (options.ftype.Type() > 2)
	{
		forecastTypeValue = boost::lexical_cast<string>(options.ftype.Value());
	}

	string forecastTypeId = boost::lexical_cast<string>(options.ftype.Type());

	if (options.time.Step() == 0 && options.ftype.Type() == 1)
	{
		// ECMWF (and maybe others) use forecast type id == 2 for analysis hour
		forecastTypeId += ",2";
	}

	for (size_t i = 0; i < gridgeoms.size(); i++)
	{
		string tablename = gridgeoms[i][1];
		string geomid = gridgeoms[i][0];

		string parm_name = options.param.Name();

		string query =
		    "SELECT param_id, level_id, level_value, forecast_period, file_location, file_server "
		    "FROM " +
		    tablename +
		    "_v "
		    "WHERE analysis_time = '" +
		    analtime +
		    "' "
		    "AND param_name = '" +
		    parm_name +
		    "' "
		    "AND level_name = upper('" +
		    level_name +
		    "') "
		    "AND level_value = " +
		    levelValue + " AND level_value2 = " + levelValue2 +
		    " "
		    "AND forecast_period = '" +
		    util::MakeSQLInterval(options.time) +
		    "' "
		    "AND geometry_id = " +
		    geomid +
		    " "
		    "AND forecast_type_id IN (" +
		    forecastTypeId +
		    ") "
		    "AND forecast_type_value = " +
		    forecastTypeValue +
		    " "
		    "ORDER BY forecast_period, level_id, level_value";

		itsRadonDB->Query(query);

		vector<string> values = itsRadonDB->FetchRow();

		if (values.empty())
		{
			continue;
		}

		itsLogger->Trace("Found data for parameter " + parm_name + " from radon geometry " + gridgeoms[i][3]);

		files.push_back(values[4]);

		break;  // discontinue loop on first positive match
	}

	return files;
}

bool radon::Save(const info& resultInfo, const string& theFileName)
{
	Init();

	stringstream query;

	if (resultInfo.Grid()->Class() != kRegularGrid)
	{
		itsLogger->Error("Only grid data can be stored to radon for now");
		return false;
	}

	/*
	 * 1. Get grid information
	 * 2. Get model information
	 * 3. Get data set information (ie model run)
	 * 4. Insert or update
	 */

	himan::point firstGridPoint = resultInfo.Grid()->FirstPoint();

	// get grib1 gridType

	int gribVersion = 1;
	int gridType = -1;

	switch (resultInfo.Grid()->Type())
	{
		case kLatitudeLongitude:
			gridType = 0;
			break;
		case kRotatedLatitudeLongitude:
			gridType = 10;
			break;
		case kStereographic:
			gridType = 5;
			break;
		case kReducedGaussian:
			gridType = 24;  // "stretched" gaussian
			break;
		case kAzimuthalEquidistant:
			gribVersion = 2;
			gridType = 110;
			break;
		default:
			throw runtime_error("Unsupported projection: " + boost::lexical_cast<string>(resultInfo.Grid()->Type()) +
			                    " " + HPGridTypeToString.at(resultInfo.Grid()->Type()));
	}

	map<string, string> geominfo = itsRadonDB->GetGeometryDefinition(
	    resultInfo.Grid()->Ni(), resultInfo.Grid()->Nj(), firstGridPoint.Y(), firstGridPoint.X(),
	    resultInfo.Grid()->Di(), resultInfo.Grid()->Dj(), gribVersion, gridType);

	if (geominfo.empty())
	{
		itsLogger->Warning("Grid geometry not found from radon");
		return false;
	}

	string geom_id = geominfo["id"];
	auto analysisTime = resultInfo.OriginDateTime().String("%Y-%m-%d %H:%M:%S+00");

	query.str("");

	query << "SELECT "
	      << "id, table_name "
	      << "FROM as_grid "
	      << "WHERE geometry_id = '" << geom_id << "'"
	      << " AND (min_analysis_time, max_analysis_time) OVERLAPS ('" << analysisTime << "'"
	      << ", '" << analysisTime << "')"
	      << " AND producer_id = " << resultInfo.Producer().Id();

	itsRadonDB->Query(query.str());

	auto row = itsRadonDB->FetchRow();

	if (row.empty())
	{
		itsLogger->Warning("Data set definition not found from radon");
		return false;
	}

	string table_name = row[1];

	query.str("");

	char host[255];
	gethostname(host, 255);

	auto paraminfo = itsRadonDB->GetParameterFromDatabaseName(resultInfo.Producer().Id(), resultInfo.Param().Name());

	if (paraminfo.empty())
	{
		itsLogger->Error("Parameter information not found from radon for parameter " + resultInfo.Param().Name() +
		                 ", producer " + boost::lexical_cast<string>(resultInfo.Producer().Id()));
		return false;
	}

	auto levelinfo = itsRadonDB->GetLevelFromGrib(resultInfo.Producer().Id(), resultInfo.Level().Type(), 1);

	if (levelinfo.empty())
	{
		itsLogger->Error("Level information not found from radon for level " +
		                 HPLevelTypeToString.at(resultInfo.Level().Type()) + ", producer " +
		                 boost::lexical_cast<string>(resultInfo.Producer().Id()));
		return false;
	}

	/*
	 * We have our own error logging for unique key violations
	 */

	// itsRadonDB->Verbose(false);
	int forecastTypeValue = -1;  // default, deterministic/analysis

	if (resultInfo.ForecastType().Type() > 2)
	{
		forecastTypeValue = static_cast<int>(resultInfo.ForecastType().Value());
	}

	double levelValue2 = (resultInfo.Level().Value2() == kHPMissingValue) ? -1 : resultInfo.Level().Value2();

	query
	    << "INSERT INTO data." << table_name
	    << " (producer_id, analysis_time, geometry_id, param_id, level_id, level_value, level_value2, forecast_period, "
	       "forecast_type_id, forecast_type_value, file_location, file_server) VALUES ("
	    << resultInfo.Producer().Id() << ", "
	    << "'" << analysisTime << "', " << geom_id << ", " << paraminfo["id"] << ", " << levelinfo["id"] << ", "
	    << resultInfo.Level().Value() << ", " << levelValue2 << ", "
	    << "'" << util::MakeSQLInterval(resultInfo.Time()) << "', "
	    << static_cast<int>(resultInfo.ForecastType().Type()) << ", " << forecastTypeValue << ","
	    << "'" << theFileName << "', "
	    << "'" << host << "')";

	try
	{
		itsRadonDB->Execute(query.str());
		query.str("");

		query << "UPDATE as_grid SET record_count = record_count+1 WHERE producer_id = " << resultInfo.Producer().Id()
		      << " AND geometry_id = " << geom_id << " AND analysis_time = '" << analysisTime << "'";

		itsRadonDB->Execute(query.str());
		itsRadonDB->Commit();
	}
	catch (const pqxx::unique_violation& e)
	{
		itsRadonDB->Rollback();

		query.str("");
		query << "UPDATE data." << table_name << " SET "
		      << "file_location = '" << theFileName << "', "
		      << "file_server = '" << host << "' WHERE "
		      << "producer_id = " << resultInfo.Producer().Id() << " AND "
		      << "analysis_time = '" << analysisTime << "' AND "
		      << "geometry_id = " << geom_id << " AND "
		      << "param_id = " << paraminfo["id"] << " AND "
		      << "level_id = " << levelinfo["id"] << " AND "
		      << "level_value = " << resultInfo.Level().Value() << " AND "
		      << "level_value2 = " << levelValue2 << " AND "
		      << "forecast_period = "
		      << "'" << util::MakeSQLInterval(resultInfo.Time()) << "' AND "
		      << "forecast_type_id = " << static_cast<int>(resultInfo.ForecastType().Type()) << " AND "
		      << "forecast_type_value = " << forecastTypeValue;

		itsRadonDB->Execute(query.str());
		itsRadonDB->Commit();
	}

	itsLogger->Trace("Saved information on file '" + theFileName + "' to radon");

	return true;
}

map<string, string> radon::Grib1ParameterName(long producer, long fmiParameterId, long codeTableVersion,
                                              long timeRangeIndicator, long levelId, double level_value)
{
	Init();

	map<string, string> paramName = itsRadonDB->GetParameterFromGrib1(producer, codeTableVersion, fmiParameterId,
	                                                                  timeRangeIndicator, levelId, level_value);
	return paramName;
}

map<string, string> radon::Grib2ParameterName(long fmiParameterId, long category, long discipline, long producer,
                                              long levelId, double level_value)
{
	Init();

	map<string, string> paramName =
	    itsRadonDB->GetParameterFromGrib2(producer, discipline, category, fmiParameterId, levelId, level_value);
	return paramName;
}
