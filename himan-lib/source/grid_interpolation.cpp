/*
 * file: grid_interpolation.cpp
 */

#include "grid_interpolation.h"
#include "cuda_plugin_helper.h"

#ifdef DEBUG

#include <boost/geometry.hpp>
#include <boost/geometry/algorithms/convex_hull.hpp>
#include <boost/geometry/algorithms/equals.hpp>
#include <boost/geometry/algorithms/within.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>
#include <boost/geometry/geometries/polygon.hpp>

BOOST_GEOMETRY_REGISTER_BOOST_TUPLE_CS(cs::cartesian)

#endif

using namespace std;
using namespace himan;
using namespace grid_interpolation;

interpolation_kernel::interpolation_kernel(latitude_longitude_grid* target, latitude_longitude_grid* source)
{
	kernel = cusp::coo_matrix<size_t, double, cusp::host_memory>(target->Size(), source->Size(), 4 * target->Size());
	for (size_t i = 0; i < target->Size(); ++i)
	{
		point p = source->XY(target->LatLon(i));

		size_t a = source->Data().Index(static_cast<size_t>(floor(p.X())), static_cast<size_t>(ceil(p.Y())), 0);
		size_t b = source->Data().Index(static_cast<size_t>(ceil(p.X())), static_cast<size_t>(ceil(p.Y())), 0);
		size_t c = source->Data().Index(static_cast<size_t>(floor(p.X())), static_cast<size_t>(floor(p.Y())), 0);
		size_t d = source->Data().Index(static_cast<size_t>(ceil(p.X())), static_cast<size_t>(floor(p.Y())), 0);

		point P = target->LatLon(i);
		point A = source->LatLon(a);
		point B = source->LatLon(b);
		point C = source->LatLon(c);
		point D = source->LatLon(d);

		vector<double> w = himan::grid_interpolation::get_weights(A, B, C, D, P);

		kernel.row_indices[4 * i] = i;
		kernel.column_indices[4 * i] = a;
		kernel.values[4 * i] = w[0];
		kernel.row_indices[4 * i + 1] = i;
		kernel.column_indices[4 * i + 1] = b;
		kernel.values[4 * i + 1] = w[1];
		kernel.row_indices[4 * i + 2] = i;
		kernel.column_indices[4 * i + 2] = c;
		kernel.values[4 * i + 2] = w[2];
		kernel.row_indices[4 * i + 3] = i;
		kernel.column_indices[4 * i + 3] = d;
		kernel.values[4 * i + 3] = w[3];
	}
}

interpolation_kernel::interpolation_kernel(rotated_latitude_longitude_grid* target, latitude_longitude_grid* source)
    : interpolation_kernel(dynamic_cast<latitude_longitude_grid*>(target), source)
{
}

#ifdef DEBUG
void check_geometry(const himan::point& a, const himan::point& b, const himan::point& c, const himan::point& d,
                    const himan::point& p)
{
	typedef boost::tuple<double, double> boost_point;
	typedef boost::geometry::model::polygon<boost_point> polygon;

	polygon poly;
	boost::geometry::append(poly, boost_point(c.X(), c.Y()));
	boost::geometry::append(poly, boost_point(a.X(), a.Y()));
	boost::geometry::append(poly, boost_point(b.X(), b.Y()));
	boost::geometry::append(poly, boost_point(d.X(), d.Y()));
	boost::geometry::append(poly, boost_point(c.X(), c.Y()));

	polygon hull;

	// check if quadrilateral is convex
	boost::geometry::convex_hull(poly, hull);
	if (!boost::geometry::equals(poly, hull)) throw std::runtime_error("Grid points do not form convex hull");

	// check if point lies within quadrilateral
	if (!boost::geometry::within(boost_point(p.X(), p.Y()), poly))
		throw std::runtime_error("Grid point is outside of convex hull");
}
#endif
