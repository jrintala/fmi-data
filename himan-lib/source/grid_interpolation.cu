#include "cuda_plugin_helper.h"
#include "grid_interpolation.h"
#include <cusp/array1d.h>
#include <cusp/array2d.h>
#include <cusp/multiply.h>
#include <cusp/transpose.h>

#include <cusp/lapack/lapack.h>

void himan::grid_interpolation::interpolate(std::vector<double>& target,const std::vector<double>& source, const himan::grid_interpolation::interpolation_kernel& kernel)
{
        cusp::array1d<double,cusp::host_memory> target_vector = target;
        cusp::array1d<double,cusp::host_memory> source_vector = source;

        cusp::multiply(kernel.kernel, source_vector, target_vector);

        std::copy(target_vector.begin(),target_vector.end(),target.begin());
}

std::vector<double> himan::grid_interpolation::get_weights(const himan::point& a, const himan::point& b, const himan::point& c, const himan::point& d, const himan::point& p)
{
		if (p==a) return std::vector<double> {1.0, 0.0, 0.0, 0.0};
                if (p==b) return std::vector<double> {0.0, 1.0, 0.0, 0.0};
                if (p==c) return std::vector<double> {0.0, 0.0, 1.0, 0.0};
                if (p==d) return std::vector<double> {1.0, 0.0, 0.0, 1.0};

		std::vector<double> ret(4);

#ifdef DEBUG
		himan::grid_interpolation::check_geometry(a,b,c,d,p);
#endif

		// points form rectangle
		if (a.X() == c.X() && a.Y() == b.Y() && d.X() == b.X() && d.Y() == c.Y())
		{
			/*
			 * Compute a linear map from points c,b to unit square (0|0), (1|1)
			 * by solving the linear equations:
			 *
			 * xi = C0 + C1*x
			 * yi = C2 + C3*y
			 *
			 * where xi,yi denotes coordinates in unit square and x,y in latlon space.
			 *
			 * In matrix-vector form this can be written as
                         *
			 *                                          |    1     1     0     0|
			 * |xi_0 xi_1 yi_0 yi_1| =  |C0 C1 C2 C3| * |   x0    x1     0     0|
			 *                                          |    0     0     1     1|
			 *                                          |    0     0    y0    y1|
                         * with C0, C1, C2, C3 to solve for. In the following code segment this is done through manually inverting the 4x4 matrix first and then 
			 * compute coefficentes by matrix-vector multiplication.
			 */

			// Matrices and Vectors
                	cusp::array2d<double,cusp::host_memory> A(4,4,0.0);
                	cusp::array1d<double,cusp::host_memory> xi(4);
                        cusp::array1d<double,cusp::host_memory> C(4);
                        cusp::array2d<double,cusp::host_memory> x(2,4,0.0);
                        cusp::array1d<double,cusp::host_memory> phi(2);

			// This is the inverted 4x4 matrix
                	A(0,0) = b.X()/(b.X()-c.X());
                	A(0,1) = -c.X()/(b.X()-c.X());
                	A(1,0) = -1/(b.X()-c.X());
               		A(1,1) = 1/(b.X()-c.X());

                	A(2,2) = b.Y()/(b.Y()-c.Y());
                	A(2,3) = -c.Y()/(b.Y()-c.Y());
                	A(3,2) = -1/(b.Y()-c.Y());
                	A(3,3) =  1/(b.Y()-c.Y());

			// Corner coordinates of unit square
			xi[0] = 0;
			xi[1] = 1;
			xi[2] = 0;
			xi[3] = 1;

			// Compute coefficients
                	cusp::multiply(A,xi,C);

			/*
                         * Apply coordinate transformation to transform point p from latlon space to point phi in unit square.
                         *
                         *                               |  1   0|
			 * |xi_p yi_p| = |C0 C1 C2 C3| * |x_p   0|
			 *                               |  0   1|
                         *                               |  0 y_p|
			 */

			x(0,0) = 1; x(0,1) = p.X();
			x(1,2) = 1; x(1,3) = p.Y();

			cusp::multiply(x,C,phi);

			/*
			 * Calculate interpolation weights using bilinear interpolation formulae
			 */

			ret[0] = (1-phi[0])*phi[1];
			ret[1] = phi[0]*phi[1];
			ret[2] = (1-phi[0])*(1-phi[1]);
			ret[3] = phi[0]*(1-phi[1]);

			return ret;
		}
		// points form generic convex quadrilateral
		else
		{
                        /*
                         * Compute a linear map from points a,b,c,d to unit square (0|1), (1|1), (0|0) , (1|0)
                         * by solving the linear equations:
                         *
                         * xi = C0 + C1*x + C2*y + C3*x*y
                         * yi = C4 + C5*x + C6*y + C7*x*y
			 *
                         * where xi,yi denotes coordinates in unit square and x,y in latlon space.
                         *
                         * In matrix-vector form this can be written as
                         *
                         *                                          |    1     1     1     1|
                         * |xi_0 xi_1 xi_2 xi_3| =  |C0 C1 C2 C3| * |   x0    x1    x2    x3|
                         * |yi_0 yi_1 yi_2 yi_3|    |C4 C5 C6 C7|   |   y0    y1    y2    y3|
                         *                                          |x0*y0 x1*y1 x2*y2 x3*y3|
			 *
			 */

			// Matrices and Vectors
			cusp::array2d<double,cusp::host_memory> A(4,4);
			cusp::array2d<double,cusp::host_memory> xi(4,2);
                        cusp::array2d<double,cusp::host_memory> xi_trans(2,4);
                        cusp::array1d<int,cusp::host_memory> piv;
                        cusp::array1d<double,cusp::host_memory> x(4);
                        cusp::array1d<double,cusp::host_memory> phi(2);

			// Construct linear system of equations to be solved
			A(0,0) = 1;           A(1,0) = 1;           A(2,0) = 1;           A(3,0) = 1;
			A(0,1) = a.X();       A(1,1) = b.X();       A(2,1) = c.X();       A(3,1) = d.X();
			A(0,2) = a.Y();       A(1,2) = b.Y();       A(2,2) = c.Y();       A(3,2) = d.Y();
			A(0,3) = a.X()*a.Y(); A(1,3) = b.X()*b.Y(); A(2,3) = c.X()*c.Y(); A(3,3) = d.X()*d.Y();
			
                        xi(0,0) = 0; xi(1,0) = 1; xi(2,0) = 0; xi(3,0) = 1;
                        xi(0,1) = 1; xi(1,1) = 1; xi(2,1) = 0; xi(3,1) = 0;

			// Solve linear system
			cusp::lapack::gesv(A,xi,piv);

                        x[0] = 1; x[1] = p.X(); x[2] = p.Y(); x[3] = p.X()*p.Y();

			cusp::transpose(xi,xi_trans);
			cusp::multiply(xi_trans,x,phi);

                        ret[0] = (1-phi[0])*phi[1];
                        ret[1] = phi[0]*phi[1];
                        ret[2] = (1-phi[0])*(1-phi[1]);
                        ret[3] = phi[0]*(1-phi[1]);

			return ret;
		}
}
