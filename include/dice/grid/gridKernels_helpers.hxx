#ifndef __GRID_KERNEL_HELPERS_HXX__
#define __GRID_KERNEL_HELPERS_HXX__

/**
 * @file
 * @brief Helper function s for grid kernels
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"
/** \addtogroup GRID
 *   \{
 */

namespace gridKernel
{

	template <int NDIM, bool PERIODIC, bool CellCenteredValues, bool SIZE_IS_ODD,
			  class DT, class DT2, class DT3, class DT4, class IT1, class IT2>
	inline size_t getCorePixelCoords(const DT coords[NDIM], const DT2 origin[NDIM], const DT3 dim[NDIM],
									 const DT4 boxSizeInv[NDIM], const IT1 stride[NDIM], IT2 output_pos[NDIM])
	{
		static const int SIZE_IS_EVEN = !SIZE_IS_ODD;

		size_t adress = 0;
		for (int i = 0; i < NDIM; ++i)
		{
			double tmp = (coords[i] - origin[i]) * boxSizeInv[i] * dim[i];
			if (SIZE_IS_EVEN && CellCenteredValues)
				tmp -= 0.5;

			double k;
			if (PERIODIC)
			{
				if (tmp < 0)
					tmp += dim[i];
				else if (tmp >= dim[i])
					tmp -= dim[i];

				modf(tmp, &k);
			}
			else
			{
				modf(tmp, &k);
				if (tmp < 0)
					k = 0;
			}
			output_pos[i] = static_cast<IT2>(k);
			adress += stride[i] * output_pos[i];
		}
		return adress;
	}

	// NOTE: dx is the distance to the central value if SIZE is ODD and the distance to the closest
	// left value if it is EVEN.
	// => if SIZE is ODD  => -0.5 <= dx < 0.5
	// => if SIZE is EVEN =>  0   <= dx < 1
	template <int NDIM, bool PERIODIC, bool CellCenteredValues, bool SIZE_IS_ODD,
			  class DT, class DT2, class DT3, class DT4, class IT1, class IT2>
	inline size_t getCorePixelCoordsAndCoefs(const DT coords[NDIM], const DT2 origin[NDIM], const DT3 dim[NDIM],
											 const DT4 boxSizeInv[NDIM], const IT1 stride[NDIM], IT2 output_pos[NDIM],
											 double dx[NDIM])
	{
		static const int SIZE_IS_EVEN = !SIZE_IS_ODD;

		size_t adress = 0;
		for (int i = 0; i < NDIM; ++i)
		{
			double tmp = (coords[i] - origin[i]) * boxSizeInv[i] * dim[i];
			if (SIZE_IS_EVEN && CellCenteredValues)
				tmp -= 0.5;

			double k;
			if (PERIODIC)
			{
				if (tmp < 0)
				{
					tmp += dim[i];
				}
				else if (tmp >= dim[i])
					tmp -= dim[i];

				dx[i] = modf(tmp, &k);
				// dx[i][1] = modf(tmp,&k);
				// dx[i][0] = 1.0 - dx[i][1];
			}
			else
			{
				dx[i] = modf(tmp, &k);
				// dx[i][1] = modf(tmp,&k);
				// dx[i][0] = 1.0 - dx[i][1];
				if (tmp < 0)
				{
					k = 0;
					dx[i] = 0.0;
					// dx[i][0]=1.0;
					// dx[i][1]=0.0;
				}
				if (tmp + 1.0 >= dim[i])
				{
					dx[i] = 1.0;
					// dx[i][0]=1.0;
					// dx[i][1]=0.0;
				}
			}
			output_pos[i] = static_cast<IT2>(k);
			adress += stride[i] * output_pos[i];

			if (SIZE_IS_ODD && CellCenteredValues)
				dx[i] -= 0.5;
		}
		return adress;
	}

} // gridKernel namespace

/** \}*/
#include "../internal/namespace.footer"
#endif
