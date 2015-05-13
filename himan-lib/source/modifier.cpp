/**
 * @file modifier.cpp
 * @author partio
 */

#include "modifier.h"
#include "NFmiInterpolation.h"

using namespace himan;

#ifdef DEBUG
#include <boost/foreach.hpp>

void DumpVector(const std::vector<double>& vec)
{

	double min = 1e38, max = -1e38, sum = 0;
	size_t count = 0, missing = 0;

	BOOST_FOREACH(double val, vec)
	{
		if (val == kFloatMissing)
		{
			missing++;
			continue;
		}

		min = (val < min) ? val : min;
		max = (val > max) ? val : max;
		count++;
		sum += val;
	}

	double mean = std::numeric_limits<double>::quiet_NaN();

	if (count > 0)
	{
		mean = sum / static_cast<double> (count);
	}

	std::cout << "min " << min << " max " << max << " mean " << mean << " count " << count << " missing " << missing << std::endl;

}
void DumpVector(const std::vector<bool>& vec)
{
	size_t count = 0, trueval = 0;
	BOOST_FOREACH(bool val, vec)
	{
		count++;
		
		if (val) trueval++;
	}
	
	std::cout << "true values: " << trueval << " false values: " << (count-trueval) << " total: " << count << std::endl;

}
#endif

modifier::modifier()
	: itsMissingValuesAllowed(false)
	, itsFindNthValue(1) // first
	, itsModifierType(kUnknownModifierType)
	, itsHeightInMeters(true)
	, itsGridsProcessed(0)
{
}

modifier::modifier(HPModifierType theModifierType)
	: itsMissingValuesAllowed(false)
	, itsFindNthValue(1) // first
	, itsModifierType(theModifierType)
	, itsHeightInMeters(true)
	, itsGridsProcessed(0)
{
}

const std::vector<double>& modifier::Result() const
{
	assert(itsResult.size());
	
	return itsResult;
}

bool modifier::CalculationFinished() const
{
	if (itsResult.size() > 0 && static_cast<size_t> (count(itsOutOfBoundHeights.begin(), itsOutOfBoundHeights.end(), true)) == itsResult.size())
	{
		return true;
	}
	
	return false;
}

void modifier::Clear(double fillValue)
{
	std::fill(itsResult.begin(), itsResult.end(), fillValue);
}

void modifier::FindValue(const std::vector<double>& theFindValue)
{
	itsFindValue = theFindValue;
}

void modifier::LowerHeight(const std::vector<double>& theLowerHeight)
{
	itsLowerHeight = theLowerHeight;
	
	// If height limits have missing values we can't process those grid points
	for (size_t i = 0; i < itsLowerHeight.size(); i++)
	{
		if (itsLowerHeight[i] == kFloatMissing)
		{
			itsOutOfBoundHeights[i] = true;
		}
	}
#ifdef DEBUG
	DumpVector(itsLowerHeight);
	DumpVector(itsOutOfBoundHeights);
#endif
}

void modifier::UpperHeight(const std::vector<double>& theUpperHeight)
{
	itsUpperHeight = theUpperHeight;
	
	// If height limits have missing values we can't process those grid points
	for (size_t i = 0; i < itsUpperHeight.size(); i++)
	{
		if (itsUpperHeight[i] == kFloatMissing)
		{
			itsOutOfBoundHeights[i] = true;
		}
	}
#ifdef DEBUG
	DumpVector(itsLowerHeight);
	DumpVector(itsOutOfBoundHeights);
#endif
}

size_t modifier::FindNth() const
{
	return itsFindNthValue;
}

void modifier::FindNth(size_t theNth)
{
	itsFindNthValue = theNth;
}

double modifier::Value() const
{
	assert(itsIndex < itsResult.size());
	
	return itsResult[itsIndex];
}

void modifier::Value(double theValue)
{
	assert(itsIndex < itsResult.size());

	itsResult[itsIndex] = theValue;
}

void modifier::Init(const std::vector<double>& theData, const std::vector<double>& theHeights)
{
	if (itsResult.size() == 0)
	{
		assert(theData.size() == theHeights.size());

		itsResult.resize(theData.size(), kFloatMissing);
		itsOutOfBoundHeights.resize(theData.size(), false);

		// Absurd default limits if user has not specified any limits
		
		if (itsHeightInMeters)
		{	
			itsLowerHeight.resize(theData.size(), 1e38);
			itsUpperHeight.resize(theData.size(), -1e38);
		}
		else
		{
			itsLowerHeight.resize(theData.size(), -1e38);
			itsUpperHeight.resize(theData.size(), 1e38);
		}
	}
}

bool modifier::Evaluate(double theValue, double theHeight)
{

	assert(itsIndex < itsOutOfBoundHeights.size());

	if (kFloatMissing == theHeight || kFloatMissing == theValue || itsOutOfBoundHeights[itsIndex] == true)
	{
		return false;
	}
	
	double lowerLimit = itsLowerHeight[itsIndex];
	double upperLimit = itsUpperHeight[itsIndex];

	if (itsHeightInMeters)
	{
		if (theHeight >= upperLimit) //  || kFloatMissing == upperLimit || kFloatMissing == lowerLimit)
		{
			// height is above given height range OR either level value is missing: stop processing of this grid point
			itsOutOfBoundHeights[itsIndex] = true;
			return false;
		}
		else if (theHeight <= lowerLimit)
		{
			// height is below given height range, do not cancel calculation yet
			return false;
		}
	}
	else
	{
		if (theHeight <= upperLimit) //  || kFloatMissing == upperLimit || kFloatMissing == lowerLimit)
		{
			itsOutOfBoundHeights[itsIndex] = true;
			return false;
		}
		else if (theHeight >= lowerLimit)
		{
			// height is below given height range, do not cancel calculation yet
			return false;
		}
	}	

	assert((lowerLimit == kFloatMissing || upperLimit == kFloatMissing) || ((itsHeightInMeters && lowerLimit <= upperLimit) || (!itsHeightInMeters && lowerLimit >= upperLimit)));

	return true;
}

void modifier::Process(const std::vector<double>& theData, const std::vector<double>& theHeights)
{

	Init(theData, theHeights);
	
	//assert(itsResult.size() == theData.size() && itsResult.size() == theHeights.size());

	for (itsIndex = 0; itsIndex < itsResult.size(); itsIndex++)
	{
		double theValue = theData[itsIndex], theHeight = theHeights[itsIndex];

		/*
		 * Evaluate() function is separated from Calculate() because Evaluate() is the
		 * same for all classes and therefore needs to defined only once
		 */
		
		if (!Evaluate(theValue, theHeight))
		{
			continue;
		}

		Calculate(theValue, theHeight);
	}

	itsGridsProcessed++;
}

size_t modifier::HeightsCrossed() const
{
	return static_cast<size_t> (count(itsOutOfBoundHeights.begin(), itsOutOfBoundHeights.end(), true));
}

HPModifierType modifier::Type() const
{
	return itsModifierType;
}

bool modifier::HeightInMeters() const
{
	return itsHeightInMeters;
}

void modifier::HeightInMeters(bool theHeightInMeters)
{
	itsHeightInMeters = theHeightInMeters;
}

std::ostream& modifier::Write(std::ostream& file) const
{
	file << "<" << ClassName() << ">" << std::endl;

	file << "__itsMissingValuesAllowed__ " << itsMissingValuesAllowed << std::endl;
	file << "__itsFindNthValue__ " << itsFindNthValue << std::endl;
	file << "__itsIndex__ " << itsIndex << std::endl;
	file << "__itsResult__ size " << itsResult.size() << std::endl;
	file << "__itsFindValue__ size " << itsFindValue.size() << std::endl;
	file << "__itsLowerHeight__ size " << itsLowerHeight.size() << std::endl;
	file << "__itsUpperHeight__ size " << itsUpperHeight.size() << std::endl;
	file << "__itsOutOfBoundHeights__ size " << itsOutOfBoundHeights.size() << std::endl;
	
	return file;
}

/* ----------------- */


void modifier_max::Calculate(double theValue, double theHeight)
{
	if (kFloatMissing == Value() || theValue > Value())
	{
		Value(theValue);
	}
}

/* ----------------- */

void modifier_min::Calculate(double theValue, double theHeight)
{
	if (kFloatMissing == Value() || theValue < Value())
	{
		Value(theValue);
	}
}

/* ----------------- */

void modifier_maxmin::Init(const std::vector<double>& theData, const std::vector<double>& theHeights)
{
	if (itsResult.size() == 0)
	{
		assert(theData.size() == theHeights.size());

		itsResult.resize(theData.size(), kFloatMissing);
		itsMaximumResult.resize(theData.size(), kFloatMissing);
		itsOutOfBoundHeights.resize(theData.size(), false);
	}
}

const std::vector<double>& modifier_maxmin::Result() const
{
	itsResult.insert(itsResult.end(), itsMaximumResult.begin(), itsMaximumResult.end());
	return itsResult;
}

void modifier_maxmin::Calculate(double theValue, double theHeight)
{
	if (kFloatMissing == Value())
	{
		// Set min == max
		itsResult[itsIndex] = theValue;
		itsMaximumResult[itsIndex] = theValue;
	}
	else
	{
		if (theValue > itsMaximumResult[itsIndex])
		{
			itsMaximumResult[itsIndex] = theValue;
		}

		if (theValue < itsResult[itsIndex])
		{
			itsResult[itsIndex] = theValue;
		}
	}
}

/* ----------------- */

void modifier_sum::Calculate(double theValue, double theHeight)
{

	if (kFloatMissing == Value()) // First value
	{
		Value(theValue);
	}
	else
	{
		double val = Value();
		Value(theValue+val);
	}
}

/* ----------------- */

bool modifier_mean::Evaluate(double theValue, double theHeight)
{

	assert(itsIndex < itsOutOfBoundHeights.size());

	if (kFloatMissing == theHeight || kFloatMissing == theValue || itsOutOfBoundHeights[itsIndex])
	{
		return false;
	}

	return true;
}

void modifier_mean::Init(const std::vector<double>& theData, const std::vector<double>& theHeights)
{

	if (itsResult.size() == 0)
	{
		assert(theData.size() == theHeights.size());

		itsResult.resize(theData.size(), kFloatMissing);
		itsRange.resize(theData.size(),0);
		itsPreviousValue.resize(itsResult.size(), kFloatMissing);
		itsPreviousHeight.resize(itsResult.size(), kFloatMissing);
	
		itsOutOfBoundHeights.resize(itsResult.size(), false);
		
		// Absurd default limits if user has not specified any limits
		
		if (itsHeightInMeters)
		{	
			itsLowerHeight.resize(theData.size(), 1e38);
			itsUpperHeight.resize(theData.size(), -1e38);
		}
		else
		{
			itsLowerHeight.resize(theData.size(), -1e38);
			itsUpperHeight.resize(theData.size(), 1e38);
		}
	}
}

void modifier_mean::Calculate(double theValue, double theHeight)
{
	if (kFloatMissing == Value()) // First value
	{		
		Value(0);
	}

	double lowerHeight=itsLowerHeight[itsIndex];
	double upperHeight=itsUpperHeight[itsIndex];
	
	double previousValue = itsPreviousValue[itsIndex];
	double previousHeight = itsPreviousHeight[itsIndex];

	itsPreviousValue[itsIndex] = theValue;
	itsPreviousHeight[itsIndex] = theHeight;
	
	// check if averaging interval is larger then 0. Otherwise skip this gridpoint and return average value of 0.
	if (lowerHeight == upperHeight)
	{
		itsOutOfBoundHeights[itsIndex] = true;
		return;
	}

	double val = Value();	

	// value is below the lowest limit
	if ((itsHeightInMeters && previousHeight < lowerHeight && theHeight > lowerHeight)
		||
		(!itsHeightInMeters && previousHeight > lowerHeight && theHeight < lowerHeight))
	{
		double lowerValue = NFmiInterpolation::Linear(lowerHeight, previousHeight, theHeight, previousValue, theValue);
		Value((lowerValue + theValue) / 2 * (theHeight - lowerHeight) + val);
		itsRange[itsIndex] += theHeight - lowerHeight;
	}
	// value is above the highest limit
	else if ((itsHeightInMeters && previousHeight < upperHeight && theHeight > upperHeight)
		||
		(!itsHeightInMeters && previousHeight > upperHeight && theHeight < upperHeight))
	{
		double upperValue = NFmiInterpolation::Linear(upperHeight, previousHeight, theHeight, previousValue, theValue);
		Value((upperValue + previousValue) / 2 * (upperHeight - previousHeight) + val);
		itsRange[itsIndex] += upperHeight - previousHeight;
		// if upper height is passed for this grid point set OutOfBoundHeight = "true" to skip calculation of the integral in following iterations
		itsOutOfBoundHeights[itsIndex] = true;

	}
	else if (previousHeight != kFloatMissing &&
		((itsHeightInMeters && previousHeight >= lowerHeight && theHeight <= upperHeight)
		||
		(!itsHeightInMeters && previousHeight <= lowerHeight && theHeight >= upperHeight)))
	{
		Value((previousValue + theValue) / 2 * (theHeight - previousHeight) + val);
			itsRange[itsIndex] += theHeight - previousHeight;
	}
}

const std::vector<double>& modifier_mean::Result() const
{
	for (size_t i = 0; i < itsResult.size(); i++)
	{
	
		double val = itsResult[i];

		if (!IsMissingValue(val) && fabs(itsRange[i]) > 0.0)
		{
			itsResult[i] = val / itsRange[i];
		}
	}

	return itsResult;
}

bool modifier_mean::CalculationFinished() const
{
	
	if (itsResult.size() == 0) 
	{
		return false;
	}

	if (itsResult.size() > 0 && static_cast<size_t> (count(itsOutOfBoundHeights.begin(), itsOutOfBoundHeights.end(), true)) == itsResult.size())
	{
		return true;
	}

	if (itsUpperHeight.empty())
	{
		return false;
	}

	for (size_t i=0; i<itsPreviousHeight.size(); ++i)
	{
		if (itsUpperHeight[i] > itsPreviousHeight[i])
		{
			return false;
		}
	}

	return true;

}

/* ----------------- */

void modifier_count::Init(const std::vector<double>& theData, const std::vector<double>& theHeights)
{

	if (itsResult.size() == 0)
	{
		assert(theData.size() == theHeights.size());

		itsResult.resize(theData.size(), 0);

		itsPreviousValue.resize(itsResult.size(), kFloatMissing);

		itsOutOfBoundHeights.resize(itsResult.size(), false);

	}
}


void modifier_count::Calculate(double theValue, double theHeight)
{
	assert(itsFindValue.size());

	double findValue = itsFindValue[itsIndex];

	if (kFloatMissing == findValue)
	{
		itsOutOfBoundHeights[itsIndex] = true;
		return;
	}

	double previousValue = itsPreviousValue[itsIndex];

	itsPreviousValue[itsIndex] = theValue;

	// First level

	if (IsMissingValue(previousValue))
	{
		return;
	}

	/**
	 *
	 * If lower value is found and current value is above wanted value, wanted value
	 * is found.
	 *
	 * Made up example
	 *
	 * How many times does value 11 exist inside a value range
	 *
	 * Input data set:
	 *
	 * Value
	 *
	 * 10
	 * --- Value 11 is found between these levels" --
	 * 12
	 *  9
	 *  9
	 * --- Value 11 is found between these levels! --
	 * 16	 
	 * 17
	 *
	 * The answer is: two times (as far as we know).
	 */

	if ((previousValue <= findValue && theValue >= findValue) // updward trend
			||
		(previousValue >= findValue && theValue <= findValue)) // downward trend
	{
		double val = Value();
		Value() == kFloatMissing ? Value(1) : Value(val + 1);
	}	
}

/* ----------------- */

void modifier_findheight::Clear(double fillValue)
{
	std::fill(itsResult.begin(), itsResult.end(), fillValue);
	std::fill(itsPreviousValue.begin(), itsPreviousValue.end(), fillValue);
	std::fill(itsPreviousHeight.begin(), itsPreviousHeight.end(), fillValue);
	std::fill(itsFoundNValues.begin(), itsFoundNValues.end(), 0);
	std::fill(itsOutOfBoundHeights.begin(), itsOutOfBoundHeights.end(), false);
	itsValuesFound = 0;
}

bool modifier_findheight::CalculationFinished() const
{
	return (itsResult.size() && (itsValuesFound == itsResult.size() || static_cast<size_t> (count(itsOutOfBoundHeights.begin(), itsOutOfBoundHeights.end(), true)) == itsResult.size()));
}

void modifier_findheight::Init(const std::vector<double>& theData, const std::vector<double>& theHeights)
{
	if (itsResult.size() == 0)
	{
		assert(theData.size() == theHeights.size());
		assert(theData.size());

		itsResult.resize(theData.size(), kFloatMissing);
		itsPreviousValue.resize(itsResult.size(), kFloatMissing);
		itsPreviousHeight.resize(itsResult.size(), kFloatMissing);
		itsFoundNValues.resize(itsResult.size(), 0);
		itsOutOfBoundHeights.resize(itsResult.size(), false);
		
		itsValuesFound = 0;
	}
}

void modifier_findheight::Calculate(double theValue, double theHeight)
{

	assert(itsFindValue.size() && itsIndex < itsFindValue.size());

	double findValue = itsFindValue[itsIndex];
	
	if (IsMissingValue(findValue) || (itsFindNthValue > 0 && !IsMissingValue(Value())))
	{
		return;
	}

	double previousValue = itsPreviousValue[itsIndex];
	double previousHeight = itsPreviousHeight[itsIndex];

	itsPreviousValue[itsIndex] = theValue;
	itsPreviousHeight[itsIndex] = theHeight;

	if (IsMissingValue(previousValue))
	{
		return;
	}
	
	/**
	 *
	 * If lower value is found and current value is above wanted value, do the interpolation.
	 *
	 * Made up example
	 *
	 * Hight range: 120 - 125
	 * What is the height when parameter value is 15?
	 *
	 * Input data set:
	 *
	 * Height / Value
	 *
	 * 120 / 11
	 * 121 / 13
	 * 122 / 14
	 * --- Height of value 15 is found somewhere between these two levels! ---
	 * 123 / 16
	 * 124 / 19
	 * 125 / 19
	 *
	 * --> lowerValueThreshold = 14
	 * --> lowerHeightThreshold = 122
	 *
	 * --> theValue (== upperValueThreshold) = 16
	 * --> theHeight (== upperHeightThreshold) = 123
	 *
	 * Interpolate between (122,14),(123,16) to get the exact value !
	 * 
	 */

	if ((previousValue <= findValue && theValue >= findValue) || (previousValue > findValue && theValue <= findValue))
	{
		double actualHeight = NFmiInterpolation::Linear(findValue, previousValue, theValue, previousHeight, theHeight);

		if (actualHeight != kFloatMissing)
		{
			if (itsFindNthValue != 0)
			{

				itsFoundNValues[itsIndex] += 1;

				if (itsFindNthValue == itsFoundNValues[itsIndex])
				{
					Value(actualHeight);
					itsValuesFound++;
					itsOutOfBoundHeights[itsIndex] = true;
				}
			}
			else
			{
				// Search for the last value
				Value(actualHeight);
			}
		}
	}

}

/* ----------------- */

void modifier_findvalue::Clear(double fillValue)
{
	std::fill(itsResult.begin(), itsResult.end(), fillValue);
	std::fill(itsPreviousValue.begin(), itsPreviousValue.end(), fillValue);
	std::fill(itsPreviousHeight.begin(), itsPreviousHeight.end(), fillValue);
	std::fill(itsOutOfBoundHeights.begin(), itsOutOfBoundHeights.end(), false);
}

void modifier_findvalue::Init(const std::vector<double>& theData, const std::vector<double>& theHeights)
{

	if (itsResult.size() == 0)
	{
		assert(theData.size() == theHeights.size());

		itsResult.resize(theData.size(), kFloatMissing);

		itsPreviousValue.resize(itsResult.size(), kFloatMissing);
		itsPreviousHeight.resize(itsResult.size(), kFloatMissing);

		itsOutOfBoundHeights.resize(itsResult.size(), false);

		// Fake lower && upper heights

		double lowestHeight = 1e38;
		double highestHeight = -1e38;

		assert(itsFindValue.size());
		
		for (size_t i = 0; i < itsFindValue.size(); i++)
		{
			double h = itsFindValue[i];

			if (h == kFloatMissing)
			{
				itsOutOfBoundHeights[i] = true;
				continue;
			}

			if (h > highestHeight)
			{
				highestHeight = h;
			}
			if (h < lowestHeight)
			{
				lowestHeight = h;
			}
		}

		if (!itsHeightInMeters)
		{
			double tmp = lowestHeight;
			lowestHeight = highestHeight;
			highestHeight = tmp;
		}
		// Give some threshold to lowest and highest heights
		
		if (itsHeightInMeters)
		{
			lowestHeight = fmax(0, lowestHeight-500); // meters
			itsLowerHeight.resize(itsResult.size(), lowestHeight);
			itsUpperHeight.resize(itsResult.size(), highestHeight+500);
		}
		else
		{
			lowestHeight = lowestHeight+150; // hectopascals
			itsLowerHeight.resize(itsResult.size(), lowestHeight);
			itsUpperHeight.resize(itsResult.size(), highestHeight-150);
		}
		
		itsValuesFound = 0;
		
#ifdef DEBUG
		DumpVector(itsLowerHeight);
		DumpVector(itsUpperHeight);
		DumpVector(itsOutOfBoundHeights);
#endif
	} 
}

bool modifier_findvalue::CalculationFinished() const
{
	return (itsResult.size() && (itsValuesFound == itsResult.size() || static_cast<size_t> (count(itsOutOfBoundHeights.begin(), itsOutOfBoundHeights.end(), true)) == itsResult.size()));
}

void modifier_findvalue::Calculate(double theValue, double theHeight)
{
	assert(itsFindValue.size() && itsIndex < itsFindValue.size());

	double findHeight = itsFindValue[itsIndex];

	double previousValue = itsPreviousValue[itsIndex];
	double previousHeight = itsPreviousHeight[itsIndex];

	itsPreviousValue[itsIndex] = theValue;
	itsPreviousHeight[itsIndex] = theHeight;

	if (IsMissingValue(previousValue) && itsGridsProcessed == 0)
	{
		// It's possible that the height requested is below the lowest hybrid level, meaning
		// that we cannot interpolate the value. In this case clamp the value to the lowest
		// hybrid level.
		
		// Clamp threshold is set to 20 meters: if the difference between requested height
		// and lowest hybrid level is larger that this then clamping is not done and
		// kFloatMissing is the result

		double diff = fabs(theHeight - findHeight);
		if (findHeight < theHeight)
		{
			if (diff < 20)
			{
				Value(theValue);
				itsValuesFound++;
				itsOutOfBoundHeights[itsIndex] = true;
			}
		}

		// previous was missing but the level we want is above current height
		return;
	}

	/**
	 *
	 * If lower height is found and current height is above wanted height,
	 * do the interpolation.
	 *
	 * Made up example
	 *
	 * Height: 124
	 * What is the parameter value?
	 *
	 * Input data set:
	 *
	 * Height / Value
	 *
	 * 120 / 11
	 * 121 / 13
	 * 122 / 14
	 * 123 / 16
	 * --- Value of height 124 is found somewhere between these two levels! ---
	 * 126 / 19
	 * 128 / 19
	 *
	 * --> lowerValueThreshold = 16
	 * --> lowerHeightThreshold = 123
	 *
	 * --> theValue (== upperValueThreshold) = 19
	 * --> theHeight (== upperHeightThreshold) = 126
	 *
	 * Interpolate between (123,16),(126,19) to get the exact value !
	 *
	 */

	if ((previousHeight <= findHeight && theHeight >= findHeight) // upward trend
			||
		(previousHeight >= findHeight && theHeight <= findHeight)) // downward trend
	{
		
		double actualValue = NFmiInterpolation::Linear(findHeight, previousHeight, theHeight, previousValue, theValue);

		if (actualValue != kFloatMissing)
		{
			Value(actualValue);
			itsValuesFound++;
			itsOutOfBoundHeights[itsIndex] = true;
		}		
	}
}

/* ----------------- */

void modifier_integral::Init(const std::vector<double>& theData, const std::vector<double>& theHeights)
{

	if (itsResult.size() == 0)
	{
		assert(theData.size() == theHeights.size());

		itsResult.resize(theData.size(), kFloatMissing);
		itsPreviousValue.resize(itsResult.size(), kFloatMissing);
		itsPreviousHeight.resize(itsResult.size(), kFloatMissing);
	
		itsOutOfBoundHeights.resize(itsResult.size(), false);

	}
}

bool modifier_integral::Evaluate(double theValue, double theHeight)
{

	assert(itsIndex < itsOutOfBoundHeights.size());

	if (theHeight == kFloatMissing || theValue == kFloatMissing || itsOutOfBoundHeights[itsIndex])
	{
		return false;
	}

	/*
 	 * upper/lower limit check moved from evaluate function to calculate for the integration case
	 */

	return true;
}

void modifier_integral::Calculate(double theValue, double theHeight)
{
	if (IsMissingValue(Value())) // First value
	{		
		Value(0);
	}

	double lowerHeight=itsLowerHeight[itsIndex];
	double upperHeight=itsUpperHeight[itsIndex];

	double previousValue = itsPreviousValue[itsIndex];
	double previousHeight = itsPreviousHeight[itsIndex];

	itsPreviousValue[itsIndex] = theValue;
	itsPreviousHeight[itsIndex] = theHeight;

	if (previousHeight < lowerHeight && theHeight > lowerHeight)
	{
		double val = Value();
		double lowerValue = NFmiInterpolation::Linear(lowerHeight, previousHeight, theHeight, previousValue, theValue);
		Value((lowerValue + theValue) / 2 * (theHeight - lowerHeight) + val);
	}
	else if (previousHeight < upperHeight && theHeight > upperHeight)
	{
		double val = Value();
		double upperValue = NFmiInterpolation::Linear(upperHeight, previousHeight, theHeight, previousValue, theValue);
		Value((upperValue + previousValue) / 2 * (upperHeight - previousHeight) + val);
	}
	else if (!(previousHeight == kFloatMissing) && previousHeight >= lowerHeight && theHeight <= upperHeight)
	{
		double val = Value();
		Value((previousValue + theValue) / 2 * (theHeight - previousHeight) + val);
	}
}

bool modifier_integral::CalculationFinished() const
{
	if (itsResult.size() > 0 && static_cast<size_t> (count(itsOutOfBoundHeights.begin(), itsOutOfBoundHeights.end(), true)) == itsResult.size())
	{
		return true;
	}
	
	if (itsPreviousHeight > itsUpperHeight)
	{
		return true;
	}

	return false;

}

/* ----------------- */

bool modifier_plusminusarea::Evaluate(double theValue, double theHeight)
{
	assert(itsIndex < itsOutOfBoundHeights.size());
    if (theHeight == kFloatMissing || theValue == kFloatMissing || itsOutOfBoundHeights[itsIndex] == true)
	{
		return false;
	}

	/*
 	 * upper/lower limit check moved from evaluate function to calculate for the averaging case
	 */

	return true;
}

void modifier_plusminusarea::Init(const std::vector<double>& theData, const std::vector<double>& theHeights)
{

	if (itsResult.size() == 0)
	{
		assert(theData.size() == theHeights.size());

		itsPlusArea.resize(theData.size(), 0);
		itsMinusArea.resize(theData.size(), 0);
		itsPreviousValue.resize(theData.size(), kFloatMissing);
		itsPreviousHeight.resize(theData.size(), kFloatMissing);
	
		itsOutOfBoundHeights.resize(theData.size(), false);
	}
}

void modifier_plusminusarea::Process(const std::vector<double>& theData, const std::vector<double>& theHeights)
{

    Init(theData, theHeights);

	assert(itsPlusArea.size() == theData.size() && itsPlusArea.size() == theHeights.size());
	
	for (itsIndex = 0; itsIndex < itsPlusArea.size(); itsIndex++)
	{
		double theValue = theData[itsIndex], theHeight = theHeights[itsIndex];
		if (!Evaluate(theValue, theHeight))
		{
			continue;
		}
		
	    Calculate(theValue, theHeight);
	}

	itsGridsProcessed++;
}

void modifier_plusminusarea::Calculate(double theValue, double theHeight)
{
	theValue-=273.15;

	double lowerHeight=itsLowerHeight[itsIndex];
	double upperHeight=itsUpperHeight[itsIndex];

	double previousValue = itsPreviousValue[itsIndex];
	double previousHeight = itsPreviousHeight[itsIndex];

	itsPreviousValue[itsIndex] = theValue;
	itsPreviousHeight[itsIndex] = theHeight;

	// check if interval is larger then 0. Otherwise skip this gridpoint and return value of 0.
	if (lowerHeight == upperHeight)
	{
		itsOutOfBoundHeights[itsIndex] = true;
	}
	// integrate numerically with separating positive from negative area under the curve.
	// find lower bound
	else if (previousHeight < lowerHeight && theHeight > lowerHeight)
	{
		double lowerValue = NFmiInterpolation::Linear(lowerHeight, previousHeight, theHeight, previousValue, theValue);
		// zero is crossed from negative to positive: Interpolate height where zero is crossed and integrate positive and negative area separately
		if (lowerValue < 0 && theValue > 0)
		{
			double zeroHeight = NFmiInterpolation::Linear(0.0, lowerValue, theValue, lowerHeight, theHeight);
			itsMinusArea[itsIndex] += lowerValue / 2 * (zeroHeight - lowerHeight);
			itsPlusArea[itsIndex] += theValue / 2 * (theHeight - zeroHeight);
		}
		// zero is crossed from positive to negative
		else if (lowerValue > 0 && theValue < 0)
		{
			double zeroHeight = NFmiInterpolation::Linear(0.0, lowerValue, theValue, lowerHeight, theHeight);
			itsPlusArea[itsIndex] += lowerValue / 2 * (zeroHeight - lowerHeight);
			itsMinusArea[itsIndex] += theValue / 2 * (theHeight - zeroHeight);
		}
		// whole interval is in the negative area
		else if (lowerValue <= 0 && theValue <= 0)
		{
			itsMinusArea[itsIndex] += (lowerValue + theValue) / 2 * (theHeight - lowerHeight);
		}
		// whole interval is in the positive area
		else
		{
			itsPlusArea[itsIndex] += (lowerValue + theValue) / 2 * (theHeight - lowerHeight);
		}
	}
	// find upper bound
	else if (previousHeight < upperHeight && theHeight > upperHeight)
	{
		double upperValue = NFmiInterpolation::Linear(upperHeight, previousHeight, theHeight, previousValue, theValue);
		// zero is crossed from negative to positive
		if (previousValue < 0 && upperValue > 0)
		{
			double zeroHeight = NFmiInterpolation::Linear(0.0, previousValue, upperValue, previousHeight, upperHeight);
			itsMinusArea[itsIndex] += previousValue / 2 * (zeroHeight - previousHeight);
			itsPlusArea[itsIndex] += upperValue / 2 * (upperHeight - zeroHeight);
		}
		// zero is crossed from positive to negative
        else if (previousValue > 0 && upperValue < 0)
		{
			double zeroHeight = NFmiInterpolation::Linear(0.0, previousValue, upperValue, previousHeight, upperHeight);
			itsPlusArea[itsIndex] += previousValue / 2 * (zeroHeight - previousHeight);
			itsMinusArea[itsIndex] += upperValue / 2 * (upperHeight - zeroHeight);
		}
		// whole interval is in the negative area
		else if (previousValue <= 0 && upperValue <= 0)
		{
			itsMinusArea[itsIndex] += (previousValue + upperValue) / 2 * (upperHeight - previousHeight);
		}
		// whole interval is in the positive area
		else
		{
			itsPlusArea[itsIndex] += (previousValue + upperValue) / 2 * (upperHeight - previousHeight);
		}
		// if upper height is passed for this grid point set OutOfBoundHeight = "true" to skip calculation of the integral in following iterations
		itsOutOfBoundHeights[itsIndex] = true;
	}
	else if (!(previousHeight == kFloatMissing) && previousHeight >= lowerHeight && theHeight <= upperHeight)
	{
		// zero is crossed from negative to positive
		if (previousValue < 0 && theValue > 0)
		{
			double zeroHeight = NFmiInterpolation::Linear(0.0, previousValue, theValue, previousHeight, theHeight);
			itsMinusArea[itsIndex] += previousValue / 2 * (zeroHeight - previousHeight);
			itsPlusArea[itsIndex] += theValue / 2 * (theHeight - zeroHeight);
		}
		// zero is crossed from positive to negative
		else if (previousValue > 0 && theValue < 0)
		{
			double zeroHeight = NFmiInterpolation::Linear(0.0, previousValue, theValue, previousHeight, theHeight);
			itsPlusArea[itsIndex] += previousValue / 2 * (zeroHeight - previousHeight);
			itsMinusArea[itsIndex] += theValue / 2 * (theHeight - zeroHeight);
		}
		// whole interval is in the negative area
		else if (previousValue <= 0 && theValue <= 0)
		{
			itsMinusArea[itsIndex] += (previousValue + theValue) / 2 * (theHeight - previousHeight);
		}
		// whole interval is in the positive area
		else
		{
			itsPlusArea[itsIndex] += (previousValue + theValue) / 2 * (theHeight - previousHeight);
		}
	}
}

bool modifier_plusminusarea::CalculationFinished() const
{

    if (itsResult.size() == 0)
	{
		return false;
	}

	if (itsMinusArea.size() > 0 && static_cast<size_t> (count(itsOutOfBoundHeights.begin(), itsOutOfBoundHeights.end(), true)) == itsMinusArea.size())
	{
		return true;
	}

	if (itsUpperHeight.empty())
	{
		return false;
	}

    for (size_t i=0; i<itsPreviousHeight.size(); ++i)
    {
        if (itsUpperHeight[i] > itsPreviousHeight[i])
        {
            return false;
        }
    }

    return true;	
}

const std::vector<double>& modifier_plusminusarea::Result() const
{
	itsPlusArea.insert(itsPlusArea.end(), itsMinusArea.begin(), itsMinusArea.end()); //append MinusArea at the end of PlusArea 
	return itsPlusArea; // return PlusMinusArea
}
