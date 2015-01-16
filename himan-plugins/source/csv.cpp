/**
 * @file csv.cpp
 *
 * @date Nov 27, 2012
 * @author: partio
 */


#include "csv.h"
#include "logger_factory.h"
#include <fstream>
#include "regular_grid.h"
#include "irregular_grid.h"
#include <boost/foreach.hpp>
#include "csv_v3.h"

using namespace std;
using namespace himan::plugin;

typedef tuple<
	int,		// station id
	string,		// station name
	double,		// longitude
	double,		// latitude
	string,		// origintime as timestamp
	string,		// forecasttime as timestamp
	string,		// level name
	double,		// level value
	string,		// parameter name
	double>		// value
record;

typedef io::CSVReader<
		10, // column count
		io::trim_chars<' ', '\t'>,
		io::no_quote_escape<','>,
		io::throw_on_overflow,
		io::no_comment
> csv_reader;

bool GetLine(csv_reader& in, record& line)
{
	int station_id;
	string station_name;
	double longitude;
	double latitude;
	string origintime;
	string forecasttime;
	string level_name;
	double level_value;
	string parameter_name;
	double value;

	if (in.read_row(station_id, station_name, longitude, latitude, origintime, forecasttime, level_name, level_value, parameter_name, value))
	{
		line = make_tuple(station_id, station_name, longitude, latitude, origintime, forecasttime, level_name, level_value, parameter_name, value);
		return true;
	}

	return false;
}

csv::csv()
{
	itsLogger = std::unique_ptr<logger> (logger_factory::Instance()->GetLog("csv"));
}

bool csv::ToFile(info& theInfo, string& theOutputFile, HPFileWriteOption fileWriteOption)
{
	throw runtime_error("Not implemented yet");
	return true;
	
}

shared_ptr<himan::info> csv::FromFile(const string& inputFile, search_options& options)
{
	info_t ret = make_shared<info> ();
	
	csv_reader in(inputFile);

	in.read_header(io::ignore_no_column, "station_id", "station_name","longitude","latitude","origintime","forecasttime","level_name","level_value","parameter_name","value");

	if (!in.has_column("station_id"))
	{
		itsLogger->Error("CSV file does not have column station_id");
		throw kFileDataNotFound;
	}

	record line;

	vector<forecast_time> times;
	vector<param> params;
	vector<level> levels;

	// First create descriptors
	while (GetLine(in, line))
	{
		/* Validate time */

		forecast_time f(get<4>(line), get<5>(line));

		if (f != options.time)
		{
			itsLogger->Debug("Time does not match");
			itsLogger->Debug("Origin time " + static_cast<string> (options.time.OriginDateTime()) + " vs " + static_cast<string> (f.OriginDateTime()));
			itsLogger->Debug("Forecast time: " + static_cast<string> (options.time.ValidDateTime()) + " vs " + static_cast<string> (f.ValidDateTime()));

			continue;
		}

		string level_name = get<6>(line);
		boost::algorithm::to_lower(level_name);

		/* Validate level */
		level l(HPStringToLevelType.at(level_name), get<7>(line));

		if (l != options.level)
		{
			itsLogger->Debug("Level does not match");
			itsLogger->Debug(static_cast<string> (options.level) + " vs " + static_cast<string> (l));

			continue;
		}

		/* Validate param */
		
		param p(get<8>(line));
		
		if (p != options.param)
		{
			itsLogger->Debug("PAram does not match");
			itsLogger->Debug(options.param.Name() + " vs " + p.Name());

			continue;
		}

		bool found = false;

		/* Prevent duplicates */
		
		BOOST_FOREACH(const forecast_time& t, times)
		{
			if (f == t)
			{
				found = true;
				break;
			}
		}

		if (!found) times.push_back(f);

		found = false;

		BOOST_FOREACH(const level& t, levels)
		{
			if (l == t)
			{
				found = true;
				break;
			}
		}

		if (!found) levels.push_back(l);

		found = false;

		BOOST_FOREACH(const param& t, params)
		{
			if (p == t)
			{
				found = true;
				break;
			}
		}

		if (!found) params.push_back(p);
				
	}

	if (times.size() == 0 || params.size() == 0 || levels.size() == 0)
	{
		itsLogger->Error("Did not find valid data from file '" + inputFile + "'");
		throw kFileDataNotFound;
	}
	
	assert(times.size());
	assert(params.size());
	assert(levels.size());

	ret->Times(times);
	ret->Params(params);
	ret->Levels(levels);

	auto base = unique_ptr<grid> (new irregular_grid()); // placeholder
	base->Projection(kLatLonProjection);
	
	ret->Create(base.get());

	itsLogger->Debug("Read " + boost::lexical_cast<string> (times.size()) 
			+ " times, " + boost::lexical_cast<string> (levels.size())
			+ " levels and " + boost::lexical_cast<string> (params.size()) + " params from file '" + inputFile + "'");

	// Then set grids

	// The csv library used is sub-standard in that it doesn't allow rewinding of the
	// file. It does provide functions to set and get file line number, but that doesn't
	// affect the reading of the file!
	
	csv_reader in2(inputFile);
	in2.read_header(io::ignore_no_column, "station_id", "station_name","longitude","latitude","origintime","forecasttime","level_name","level_value","parameter_name","value");

	while (GetLine(in2, line))
	{
		forecast_time f(get<4>(line), get<5>(line));

		string level_name = get<6>(line);
		boost::algorithm::to_lower(level_name);

		level l(HPStringToLevelType.at(level_name), get<7>(line));
		param p(get<8>(line));
		
		station s (get<0>(line), get<1>(line), get<2>(line), get<3>(line));

		ret->Param(p);
		ret->Time(f);
		ret->Level(l);

		if (!ret->Grid())
		{
			throw runtime_error("ASDASD");
			auto g = make_shared<irregular_grid> ();
			vector<station> stations;
			stations.push_back(s);

			g->Stations(stations);
			ret->Grid(g);
		}
		else
		{
			// Add new station
			auto stats = dynamic_cast<irregular_grid*> (ret->Grid())->Stations();
			stats.push_back(s);
			dynamic_cast<irregular_grid*> (ret->Grid())->Stations(stats);

			// Add the data point
			ret->Grid()->Value(stats.size()-1, get<9> (line));
		}	
	}

	return ret;
	
}
