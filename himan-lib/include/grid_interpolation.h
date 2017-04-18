/*
 * matrix based interpolation utilities
 */

#ifndef GRID_INTERPOLATION_H
#define GRID_INTERPOLATION_H

#include "latitude_longitude_grid.h"
#include "plugin_configuration.h"
#include <cusp/coo_matrix.h>

namespace himan
{
namespace grid_interpolation
{
class interpolation_kernel
{
   public:
	interpolation_kernel();
	interpolation_kernel(latitude_longitude_grid* target, latitude_longitude_grid* source);
	interpolation_kernel(rotated_latitude_longitude_grid* target, latitude_longitude_grid* source);
	~interpolation_kernel(){};

	friend void interpolate(std::vector<double>& target, const std::vector<double>& source,
	                        const interpolation_kernel& kernel);

   private:
	cusp::coo_matrix<size_t, double, cusp::host_memory> kernel;
};

void interpolate(std::vector<double>& target, const std::vector<double>& source, const interpolation_kernel& kernel);

std::vector<double> get_weights(const himan::point& a, const himan::point& b, const himan::point& c,
                                const himan::point& d, const himan::point& p);

#ifdef DEBUG

void check_geometry(const himan::point& a, const himan::point& b, const himan::point& c, const himan::point& d,
                    const himan::point& p);

#endif
}
}

#endif /* GRID_INTERPOLATION_H */
