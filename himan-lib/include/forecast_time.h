/**
 * @file forecast_time.h
 *
 * @date Dec 1, 2012
 * @author partio
 *
 */

#ifndef FORECAST_TIME_H
#define FORECAST_TIME_H

#include "logger.h"
#include "raw_time.h"
#include <stdexcept>

namespace himan
{

/**
 * @class forecast_time
 *
 * @brief Describe a forecast time: origin time and valid (forecasted) time
 */

class forecast_time
{

public:
    forecast_time();
    forecast_time(const raw_time& theOriginDateTime, const raw_time& theValidDateTime);
    forecast_time(std::shared_ptr<raw_time> theOriginDateTime, std::shared_ptr<raw_time> theValidDateTime);
    forecast_time(const std::string& theOriginDateTime,
                  const std::string& theValidDateTime,
                  const std::string& theDateMask = "%Y-%m-%d %H:%M:%S");

    ~forecast_time() {}
    forecast_time(const forecast_time& other);
    forecast_time& operator=(const forecast_time& other);

    std::string ClassName() const
    {
        return "himan::forecast_time";
    };

    HPVersionNumber Version() const
    {
        return HPVersionNumber(0, 1);
    }


    std::ostream& Write(std::ostream& file) const;

    bool operator==(const forecast_time& other);
    bool operator!=(const forecast_time& other);

    int Step() const;

    std::shared_ptr<raw_time> OriginDateTime() const;
    void OriginDateTime(std::shared_ptr<raw_time> theOriginDateTime);
    void OriginDateTime(std::string& theOriginDateTime, const std::string& theDateMask = "%Y-%m-%d %H:%M:%S");

    std::shared_ptr<raw_time> ValidDateTime() const;
    void ValidDateTime(std::shared_ptr<raw_time> theValidDateTime);
    void ValidDateTime(std::string& theValidDateTime, const std::string& theDateMask = "%Y-%m-%d %H:%M:%S");

    /**
	 *
     * @return Time step resolution
     */

    HPTimeResolution StepResolution() const;
    void StepResolution(HPTimeResolution theStepResolution);

private:
    std::unique_ptr<logger> itsLogger;

    std::shared_ptr<raw_time> itsOriginDateTime;
    std::shared_ptr<raw_time> itsValidDateTime;

    HPTimeResolution itsStepResolution;
};

inline
std::ostream& operator<<(std::ostream& file, const forecast_time& ob)
{
    return ob.Write(file);
}

} // namespace himan

#endif /* FORECAST_TIME_H */
