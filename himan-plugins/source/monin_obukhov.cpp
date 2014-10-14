/**
 * @file monin_obukhov.cpp
 *
 * Claculates the inverse of the Monin-Obukhov length.
 *
 * @date Aug 8, 2014
 * @author Tack
 */

#include <boost/lexical_cast.hpp>

#include "monin_obukhov.h"
#include "logger_factory.h"
#include "level.h"
#include "forecast_time.h"
#include "logger_factory.h"

using namespace std;
using namespace himan::plugin;

monin_obukhov::monin_obukhov()
{
	itsClearTextFormula = "1/L = -(k*g*Q)/(u*^3*T)";

	itsLogger = logger_factory::Instance()->GetLog("monin_obukhov");
}

void monin_obukhov::Process(std::shared_ptr<const plugin_configuration> conf)
{
	Init(conf);

	/*
	 * Set target parameter properties
	 * - name PARM_NAME, this name is found from neons. For example: T-K
	 * - univ_id UNIV_ID, newbase-id, ie code table 204
	 * - grib1 id must be in database
	 * - grib2 descriptor X'Y'Z, http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-2.shtml
	 *
	 */

	param theRequestedParam("MOL-M", 1204);

	// If this param is also used as a source param for other calculations

	SetParams({theRequestedParam});

	Start();
}

/*
 * Calculate()
 *
 * This function does the actual calculation.
 */

void monin_obukhov::Calculate(shared_ptr<info> myTargetInfo, unsigned short threadIndex)
{

	/*
	 * Required source parameters
	 *
	 * eg. param PParam("P-Pa"); for pressure in pascals
	 *
	 */

	const param TParam("T-K"); // ground Temperature
	const param QParam("FLSEN-JM2"); // surface heat flux
	const param U_SParam("FRVEL-MS"); // friction velocity

	// ----	

	auto myThreadedLogger = logger_factory::Instance()->GetLog("monin_obukhov Thread #" + boost::lexical_cast<string> (threadIndex));

	// Current time and level
	
	forecast_time forecastTime = myTargetInfo->Time();
	level forecastLevel = level(himan::kHeight, 0, "Height"); 
	
	myThreadedLogger->Debug("Calculating time " + static_cast<string> (*forecastTime.ValidDateTime()) + " level " + static_cast<string> (forecastLevel));

	info_t TInfo = Fetch(forecastTime, forecastLevel, TParam, false);
	info_t QInfo = Fetch(forecastTime, forecastLevel, QParam, false);
	info_t U_SInfo = Fetch(forecastTime, forecastLevel, U_SParam, false);

	// determine length of forecast step to calculate surface heat flux in W/m2
	double forecastStepSize;

	if ( itsConfiguration->SourceProducer().Id() != 199)
	{
		forecastStepSize = itsConfiguration->ForecastStep()*3600; //step size in seconds
	}
	else
	{
		forecastStepSize = itsConfiguration->ForecastStep()*60; //step size in seconds
	}

	if (!TInfo || !QInfo || !U_SInfo)
	{
		myThreadedLogger->Info("Skipping step " + boost::lexical_cast<string> (forecastTime.Step()) + ", level " + static_cast<string> (forecastLevel));

		if (itsConfiguration->StatisticsEnabled())
		{
			// When time or level is skipped, all values are missing
			itsConfiguration->Statistics()->AddToMissingCount(myTargetInfo->Data()->Size());
			itsConfiguration->Statistics()->AddToValueCount(myTargetInfo->Data()->Size());
		}

		return;

	}

	string deviceType = "CPU";

	LOCKSTEP(myTargetInfo, TInfo, QInfo, U_SInfo)
	{

		double T = TInfo->Value();
		double Q = QInfo->Value()/forecastStepSize; // divide by time step to obtain Watts
		double U_S = U_SInfo->Value();
		
		double mol(kFloatMissing);

		if (T == kFloatMissing || Q == kFloatMissing || U_S == kFloatMissing)
		{
			continue;
		}

		/* Calculation of the inverse of Monin-Obukhov length to avoid division by 0 */
		
		if (U_S != 0.0)
		{
			mol = -constants::kG * constants::kK * Q / (U_S * U_S * U_S * T);
		}
		myTargetInfo->Value(mol);

	}

	myThreadedLogger->Info("[" + deviceType + "] Missing values: " + boost::lexical_cast<string> (myTargetInfo->Data()->MissingCount()) + "/" + boost::lexical_cast<string> (myTargetInfo->Data()->Size()));

}