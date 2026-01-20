#ifndef __LOCAL_AMR_GRID_RAYTRACER_HXX__
#define __LOCAL_AMR_GRID_RAYTRACER_HXX__

#include <type_traits>
#include <limits>
#include <cmath>
#include "../tools/helpers/helpers.hxx"

#ifdef HAVE_BOOST
#include "../tools/wrappers/boostMultiprecisionFloat128.hxx"
// #include <boost/multiprecision/cpp_dec_float.hpp>
// #endif
#ifdef HAVE_GMP
#include <boost/cstdfloat.hpp>
#include <boost/multiprecision/gmp.hpp>
#endif
#endif
// #include "../tools/wrappers/boostMultiprecisionFloat128.hxx"

// #define DEBUG_ME 1

/**
 * @file
 * @brief  A raytracer for a local AMR grid
 * @author Thierry Sousbie
 */

#include "../internal/namespace.header"
/** \addtogroup AMR
 *   \{
 */

/**
 * \class LocalAmrGridRaytracerT
 * \brief A raytracer for a local AMR grid
 * \verbatim
 * The raytracer is initialized by calling one of the 'reset' functions for each new ray,
 * setting the entry point of the ray and its corresponding 'current' voxel.
 * After calling reset, a first call to advance() will compute the exit point of
 * the current voxel, as well as the next voxel the ray will enter upon leaving the
 * current one. Each following call to advance() will traverse one voxel, updating
 * the entry point, exit point, current voxel and next voxel accordingly.
 *
 * The reset function initialize the raytracer by setting a source point A for the ray and
 * destination B. The implementation guarantees that the ray will pass through the voxels
 * AMR::getVoxelAt(A) and AMR::getVoxelAt(B) even in degenerate cases. To achieve that,
 * we use exact computations when necessary and simulation of simplicity by perturbing
 * the voxel vertices coordinate Xi along dimension i (0<=i<NDIM) as follows:
 *
 *                 Xi -> Xi - pow(eps,NDIM-i+1)
 *
 * where this particular convention is taken for compatibility with the point in simplex
 * predicates. Note that in 3D, linking to Boost and GMP is necessary to achieve consistency.
 *
 * Note that for the raytracer to work, the AMR grid implementation of simulation of
 * simplicity should be consistent (i.e. same convention should be used). Specifically, if
 * the voxel containing a given point lying exactly (i.e. using exact computations) on
 * the boundary of several voxels is queried to the AMR implementation, it should return
 * the one voxel for which each coordinate of the barycenter is the highest.
 *
 * If a ray passes exactly through a vertex/edges/face, it will also contourn it in a
 * consistent way and it is guaranteed that it will at least go through the voxels whose
 * barycenter coordinates are the highest among all voxels sharing this vertex/edge/face.
 * Rays are also reversible in the sense that traversing the grid from point A to point B
 * or from point B to point A result in the exact same set of voxels in the opposite order.
 *
 * Finally, it is also guaranteed that a ray will never cross a face it is exactly parallel
 * to, even when the source and destination points are exactly on voxels faces, provided
 * the origin point's voxel is the one whose barycenter coordinates are the highest among
 * all voxels sharing this origin point on their boundary.
 * \endverbatim
 *
 * \tparam AMR the AMR grid to use
 * \tparam IF  The type of floating point number to use for computations and coordinates.
 *  The choice of this type does not influence the exactness of the ray pass and simulation
 *  of simplicity, but rather the precision of the estimated exit point coordinates.
 *  Defaults to 'long double'.
 *
 */

// TODO : length is not implemented yet ...
// TODO : check if input coord is outside of boundaries, check if ray is leaving boundaries

template <class AMR, typename IF = long double>
class LocalAmrGridRaytracerT
{
public:
	static const int NDIM = AMR::NDIM;

	typedef IF Float;
	typedef IF Coord;
	typedef typename AMR::ICoord ICoord;
	typedef typename AMR::Voxel Voxel;
	typedef typename AMR::GeometricProperties GeometricProperties;
	typedef AMR Amr;

#ifdef HAVE_BOOST
#ifdef HAVE_GMP
	typedef boost::multiprecision::mpf_float mpfloat;
	// typedef boost::multiprecision::mpf_float_50 mpfloat;
	// typedef Float128OrMore mpfloat;
#endif
#endif

	/*
	FILE *fl;
	void openFile()
	{
	  static int count =0;
	  char name[255];
	  sprintf(name,"ray_%6.6d.dat",count++);
	  fl=fopen(name,"w");
	}
	void closeFile()
	{
	  fclose(fl);
	}
	void printToFile()
	{
	  fprintf(fl,"%g %g\n",exitPoint[0],exitPoint[1]);
	}
	*/
	LocalAmrGridRaytracerT(AMR *amr_) : amr(amr_)
	{
		geometry = amr->getGeometry();
		curVoxel = NULL;
		// openFile();
	}

	~LocalAmrGridRaytracerT()
	{
		// closeFile();
	}

	/*
	template <class T>
	void computeRayVector(const T *entryPoint, const T *exitPoint, T* result)
	{
	  geometry->template getVector<T,NDIM>(entryPoint,exitPoint,result);

	  for (int i=0;i<NDIM;++i)
		{
	  T tmp=entryPoint[i]+result[i];
	  geometry->checkBoundary(tmp,i);
	  if (fabs(tmp-exitPoint[i])>amr->getEpsilon(i))
		{
		  double eps = (exitPoint[i]-entryPoint[i])*std::numeric_limits<T>::epsilon();
		  if (entryPoint[i]<exitPoint[i]) result[i]=-eps;
		  else
			result[i]=-geometry->correctCoordsDiff(entryPoint[i]-exitPoint[i],i);
		}
		}

	}
	*/

	/** \brief Initialize the raytracer to a given starting point and destination point.
	 * The direction of the ray is given by the difference of the two coordinates and
	 * although the ray will continue after destPoint upon calling advance(), it is
	 * guaranted to go through the voxel containing destPoint.
	 * Note that advance() must be called once before querying the first exitPoint.
	 *  \param sourcePoint the coordinates of the starting point of the ray
	 *  \param destPoint the point toward which the ray is going. The ray is guaranted to pass
	 *  through the voxel containing destPoint but wont stop there if you keep calling
	 *  advance().
	 *  \tparam InputTypePrecision The type that has a precision equivalent to the input
	 *  coordinates
	 *  \return The voxel containing sourcePoint
	 */
	template <class InputIterator, typename InputTypePrecision = double>
	Voxel *reset(InputIterator sourcePoint,
				 InputIterator destPoint)
	{
		Voxel *voxel = amr->getVoxelAt(sourcePoint);
		reset<InputIterator, InputTypePrecision>(sourcePoint, destPoint, voxel);
	}

	/** \brief Initialize the raytracer to a given starting point and destination point.
	 * The direction of the ray is given by the difference of the two coordinates and
	 * although the ray will continue after destPoint upon calling advance(), it is
	 * guaranted to go through the voxel containing destPoint.
	 * Note that advance() must be called once before querying the first exitPoint.
	 *  \param sourcePoint the coordinates of the starting point of the ray
	 *  \param destPoint the point toward which the ray is going. The ray is guaranted to pass
	 *  through the voxel containing destPoint but wont stop there if you keep calling
	 *  advance().
	 *  \param sourceVoxel The voxel that contains entryPoint (for speed, if you already queried it)
	 *  \tparam InputTypePrecision The type that has a precision equivalent to the input
	 *  coordinates
	 *  \return The voxel containing sourcePoint (i.e. == sourceVoxel parameter)
	 */
	template <class InputIterator, typename InputTypePrecision = double>
	Voxel *reset(InputIterator sourcePoint,
				 InputIterator destPoint,
				 Voxel *sourceVoxel)
	{
		/*
		typedef typename
		  std::remove_reference<
		typename std::remove_const<
		  decltype(*sourcePoint)
		  >::type
		  >::type InputType;
		*/
		rayDir = 0;
		for (int i = 0; i < NDIM; ++i)
		{
			origin[i] = exitPoint[i] = entryPoint[i] =
				geometry->checkBoundary(Coord(*sourcePoint), i);
			/*
			Coord tmpC=Coord(*destPoint)-Coord(*sourcePoint); // force evaluation
			direction[i] = geometry->correctCoordsDiff(tmpC,i);
			*/
			direction[i] = geometry->template correctCoordsDiff<Coord>(Coord(*destPoint) - Coord(*sourcePoint), i);

#ifdef HAVE_BOOST
#ifdef HAVE_GMP
			// We need the exact direction when multiprecision is required but this would slow
			// down reset(...) to compute it, so we will initialize it later only if necessary
			mp_exact_dir_init = true;
			rayDirPoint[i] = (*destPoint);
			initEntryPoint[i] = (*sourcePoint);
#endif
#endif

			int tmp = (direction[i] > 0);
			directionSign[i] = ((tmp << 1) - 1);
			rayDir |= tmp << i;

			++sourcePoint;
			++destPoint;
		}
		/*
		if (check)
		  {
		geometry->template checkBoundary<Coord>(exitPoint);
		for (int i=0;i<NDIM;++i)
		  geometry->template correctCoordsDiff<Coord>(direction[i],i);
		  }
		*/
		anyEpsilon = 0;
		invalidDir = 0;
		for (int i = 0; i < NDIM; ++i)
		{
			if (direction[i] != 0)
			{
				direction_inv[i] = 1.0L / direction[i];
			}
			else
			{
				direction_inv[i] = 1.0;
				invalidDir |= (1 << i);
			}
			// This is for the tested filter !
			directionEpsilon[i] = fabs(amr->getEpsilon(i) * direction_inv[i] * 0.5);

			// This is for the untested filter !
			/*
			directionEpsilon[i]=
			  std::max(fabs(rayDirPoint[i]),fabs(initEntryPoint[i]))*
			  std::numeric_limits<InputTypePrecision>::epsilon()*10;
			*/

			if (directionEpsilon[i] > anyEpsilon)
				anyEpsilon = directionEpsilon[i];
		}

		nextVoxel = sourceVoxel;
		curVoxel = sourceVoxel;

		// curVoxel->print(amr,"sourceVoxel:");
		// printf("entry = (%20.20lg,%20.20lg,%20.20lg) dir =(%20.20lg,%20.20lg,%20.20lg) dest = (%20.20lg,%20.20lg,%20.20lg)\n",
		// 	   origin[0],origin[1],origin[2],
		//  	   direction[0],direction[1],direction[2],
		// 	   origin[0]+direction[0],origin[1]+direction[1],origin[2]+direction[2]);

		// This is not set before calling advance at least once
		exitDim = -1;

		return sourceVoxel;
	}

	/** \brief Make the ray traverse exactly one voxel and store its entry and exit points
	 *  \return A pointer to the next voxel the ray will enter upon calling advance().
	 */
	Voxel *advance(/*bool forceExact=false*/)
	{
		// FIXME: BE CAREFULL WITH BOUNDARY CONDITIONS HERE, CHECK PERIODIC !!!
		// In particular, entryPoint / exitPoint should be on the same side, depending on rayDir

		for (int i = 0; i < NDIM; ++i)
			entryPoint[i] = exitPoint[i];
		// std::swap(entryPoint,exitPoint);

		// curVoxel->print(amr,"prevVoxel:");
		curVoxel = nextVoxel;
		// curVoxel->print(amr,"curVoxel:");

		Coord corner[NDIM];	   // This is the corner in the approx. direction of the ray
		Coord oppCorner[NDIM]; // This is its opposite corner
		// Each coord of this corner vertex is the coordinate along the direction
		// orthogonal to them of the NDIM faces of the voxel the ray may reach first.
		// We also get the coords of the opposite corner so that we have a full bounding
		// box

		getCornerCoords(curVoxel, corner, oppCorner);
		/*
		amr->index2CornerCoordsAndOpp(curVoxel->getIndex(),curVoxel->getLevel(),
					  corner,oppCorner,rayDir);

		// printf("corner=(%f,%f)/(%f,%f), direction : (%f,%f) / rayDir =%d, wrongDir=%d\n",
		//    	   corner[0],corner[1],oppCorner[0],oppCorner[1],
		//    	   direction[0],direction[1],rayDir,invalidDir);

		// This test should be removed by the compiler for non periodic boundary ...
		if (AMR::BOUNDARY_TYPE == BoundaryType::PERIODIC)
		  {
		// for PBC, the corners of the voxel are moved so that they are always in the
		// direction of the ray from the origin.
		for (int i=0;i<NDIM;++i)
		  {
			if ( (corner[i]<origin[i]) && (direction[i]>0) )
			  {
			corner[i]+=geometry->getBBoxSize(i);
			oppCorner[i]+=geometry->getBBoxSize(i);
			  }
			else if ( (corner[i]>origin[i]) && (direction[i]<0) )
			  {
			corner[i]-=geometry->getBBoxSize(i);
			oppCorner[i]-=geometry->getBBoxSize(i);
			  }
		  }
		  }
		*/

		// if (forceExact)
		//   return advanceExact(corner,oppCorner);
		// else
		return advanceAdaptive(corner, oppCorner);
	}

	/** \brief Get a constant pointer to the coordinates of the entry point.
	 *  \return A constant pointer to the entry point cooridnates.
	 */
	const Coord *getEntryPointConstPtr() const
	{
		return entryPoint;
	}

	/** \brief Get a constant pointer to the coordinates of the exit point.
	 *  \return A constant pointer to the exit point cooridnates.
	 */
	const Coord *getExitPointConstPtr() const
	{
		return exitPoint;
	}

	/** \brief Retrieve the entry point coordinates
	 *  \param[out] coordsOut an output iterator to store the coordinates
	 */
	template <class OutputIterator>
	void getEntryPointCoords(OutputIterator coordsOut) const
	{
		for (int i = 0; i < NDIM; ++i)
		{
			(*coordsOut) = entryPoint[i];
			++coordsOut;
		}
	}

	/** \brief Retrieve the exit point coordinates
	 *  \param[out] coordsOut an output iterator to store the coordinates
	 */
	template <class OutputIterator>
	void getExitPointCoords(OutputIterator coordsOut) const
	{
		for (int i = 0; i < NDIM; ++i)
		{
			(*coordsOut) = exitPoint[i];
			++coordsOut;
		}
	}

	/** \brief get the voxel for which the entry/exit points were computed
	 */
	Voxel *getCurVoxel() const
	{
		return curVoxel;
	}

	/** \brief get the voxel the ray will enter next
	 */
	Voxel *getNextVoxel() const
	{
		return nextVoxel;
	}

	/** \brief get the index of the dimension orthogonal to the ray's exit face.
	 *  \return {0,..,NDIM-1}
	 *  \see getExitFaceNormalSign
	 */
	int getExitFaceNormalDim() const
	{
		return exitDim;
	}

	/** \brief get the direction of the normal to the exit face
	 *  \return {+1,-1} if it is oriented positively/negatively along the normal to the
	 *  exit face
	 *  \see getExitFaceNormalDim
	 */
	Coord getExitFaceNormalSign() const
	{
		return directionSign[exitDim]; // exitDimSign;
	}

private:
	Voxel *advanceAdaptive(Coord corner[NDIM], Coord oppCorner[NDIM])
	{
		Float dMin;
		bool needExactComputation = false;

		exitDim = 0; // index of the dimension along which the ray is exiting

		/* THIS FILTER IS NOT TESTED VERY MUCH, CHECK IF IT FAILS */
		/*
		Coord direction_ref=direction[0];
		int iStart=0;
		dMin=(corner[0]-origin[0]);
		while ( (invalidDir&(1<<iStart)) )
		  {
		++iStart;
		direction_ref = direction[iStart];
		dMin=(corner[iStart]-origin[iStart]);
		exitDim=iStart;
		  }
		++iStart;

		for (int i=iStart;i<NDIM;++i)
		  {
		if (invalidDir&(1<<i)) continue;

		Float d=(corner[i]-origin[i]);
		Float d_dir = d*direction_ref;
		Float dMin_dir = dMin*direction[i];

		Float delta=(fabs(d_dir)+fabs(dMin_dir))*
		  anyEpsilon;

		if (fabs(d_dir-dMin_dir) <= delta)
		  {
			needExactComputation = true;

			// If the ray is crossing several boundaries at the same time, we choose to
			// cross in the positive direction first, and along the dimension with
			// highest index. If there is no positive direction crossing, then we cross
			// in the negative direction, along the dimension with lowest index.
			// This criterium complies with our simulation of simplicity convention,
			// DO NOT CHANGE IT !
			if (rayDir&(1<<i))
			  {
			Float tmp_dMin = d*direction_inv[i];
			// Check if we would cross the boundary, as in that case simulation
			// of simplicity may not work ...
			if (geometry->inBound(origin[i]+tmp_dMin*direction[i]+
						  directionSign[i]*amr->getEpsilon(i),i))
			  exitDim=i;
			  }
		  }
		else if ((direction_ref<0)==(direction[i]<0))
		  {
			if (d_dir < dMin_dir)
			  {
			dMin=d;
			exitDim=i;
			direction_ref = direction[i];
			  }
		  }
		else
		  {
			if (d_dir > dMin_dir)
			  {
			dMin=d;
			exitDim=i;
			direction_ref = direction[i];
			  }
		  }
		  }
		dMin /= direction_ref;
		*/
		/********* DOWN TO HERE *************/

		/* THIS FILTER IS VERIFIED ... Check speed compared to untested */
		if (invalidDir == 0)
		{
			dMin = (corner[0] - origin[0]) * direction_inv[0];

			for (int i = 1; i < NDIM; ++i)
			{
				Float d = (corner[i] - origin[i]) * direction_inv[i];

				// this happens when we go through a vertex/edge/face
				// the test should be if (d==dMin) but we compensate for the fact that we do
				// not use exact arithmetic ...
				if (fabs(d - dMin) < directionEpsilon[i])
				{
					needExactComputation = true;

					// If the ray is crossing several boundaries at the same time, we choose to
					// cross in the positive direction first, and along the dimension with
					// highest index. If there is no positive direction crossing, then we cross
					// in the negative direction, along the dimension with lowest index.
					// This criterium complies with our simulation of simplicity convention,
					// DO NOT CHANGE IT !
					if (rayDir & (1 << i))
					{
						// Check if we would cross the boundary, as in that case simulation
						// of simplicity may not work ...
						if (geometry->inBound(origin[i] + dMin * direction[i] + directionSign[i] * amr->getEpsilon(i), i))
							exitDim = i;
					}
				}
				else if (d < dMin)
				{
					dMin = d;
					exitDim = i;
				}
			}
		}
		else
		{
			dMin = amr->getMaxLength();
			for (int i = 0; i < NDIM; ++i)
			{
				if (!(invalidDir & (1 << i)))
				{
					Float d = (corner[i] - origin[i]) * direction_inv[i];

					// this happens when we go through a vertex/edge/face
					if (fabs(d - dMin) < directionEpsilon[i])
					{
						needExactComputation = true;
						// exitDim is chosen to be in the lowest index positive ray direction
						// if there is one, and the highest index negative direction if non
						// exist. This way, the ray trajectory is reversible even in degenerate
						// cases !
						// if (!(rayDir&(1<<exitDim))) exitDim=i;

						// This criterium complies with our simulation of simplicity convention
						// DO NOT CHANGE IT !
						if (rayDir & (1 << i))
						{
							if (geometry->inBound(origin[i] + dMin * direction[i] + directionSign[i] * amr->getEpsilon(i), i))
								exitDim = i;
						}
					}
					else if (d < dMin)
					{
						dMin = d;
						exitDim = i;
					}
				}
			}
		}
		/****DOWN TO HERE ! ****/

		// if (curVoxel->index == 3368694746067226624L)
		//   {
		// 	printf("Exact computation status : %d\n",(int)needExactComputation);
		//   }

		/*
	#ifdef HAVE_BOOST
	#ifdef HAVE_GMP
		if (needExactComputation) return advanceExact(corner,oppCorner);
	#endif
	#endif
		*/

		Coord neiCoords[NDIM];
		// Compute the exitPoint and prepare the coordinates for finding neighbor voxel
		for (int i = 0; i < NDIM; ++i)
		{
			neiCoords[i] = origin[i] + dMin * direction[i];
			exitPoint[i] = geometry->template checkBoundary<Coord>(neiCoords[i], i);

			// This is to prevent problems when a point is exactly on the edge along the non
			Float eps = amr->getEpsilon(i);
			if ((fabs(neiCoords[i] - corner[i]) < eps) && (i != exitDim))
			{
				// printf("DISPLACED REG\n");
				neiCoords[i] -= directionSign[i] * eps;
			}
			else if ((fabs(neiCoords[i] - oppCorner[i]) < eps) && (i != exitDim))
			{
				// printf("DISPLACED OPP\n");
				neiCoords[i] += directionSign[i] * eps;
			}
		}

		// Now move the exit point coordinates by epsilon along the direction orthogonal to
		// the face it lies on i order to retrieve the neighbor
		neiCoords[exitDim] += directionSign[exitDim] * amr->getEpsilon(exitDim);
		nextVoxel = amr->getVoxelAt(neiCoords);
		/*
		if (findVertex(exitPoint[0],exitPoint[1],exitPoint[2],0.248461,0.0859375,0.171875))
		  {
		printf("NON Degenerate exit point @ (%g %g %g) : exit@%d(%g)\n",
			   exitPoint[0],exitPoint[1],exitPoint[2],exitDim,directionSign[exitDim]);
		nextVoxel->print(amr,"A-DEST:");
		  }
		*/
#ifdef HAVE_BOOST
#ifdef HAVE_GMP

		// We may be exiting a coarse voxel right on the edge of a finer one
		// in which case we may need exact computations ...
		if ((nextVoxel != NULL) && (nextVoxel->getLevel() > curVoxel->getLevel()))
		{
			int reqLevel = nextVoxel->getLevel();
			for (int i = 0; i < NDIM; ++i)
			{
				if (i == exitDim)
					continue;

				int levelDiff = nextVoxel->getLevel() - curVoxel->getLevel();
				Float nSteps = Float(1L << levelDiff);
				Float c = corner[i];
				Float o = oppCorner[i];

				if (o < c)
					std::swap(o, c);

				Float delta = (o - c) / nSteps;
				Float ratio = (exitPoint[i] - c) / delta;
				Float ratioPlusHalf = ratio + 0.5;
				Float rep = c + delta * hlp::numericStaticCast<long>(ratioPlusHalf);
				// hlp::numericStaticCast<int>(ratio+0.5);

				if ((ratio < 0.5) || (ratio > nSteps - 0.5))
					continue;
				if (fabs(exitPoint[i] - rep) < amr->getEpsilon(i))
				{
					// We will use exact computation, no need to do it again when we are done.
					needExactComputation = false;

					// if (findVertex(exitPoint[0],exitPoint[1],exitPoint[2],0.248461,0.0859375,0.171875))
					//  printf("NEEDED\n");
					Float neiCoords2[NDIM];
					for (int j = 0; j < NDIM; ++j)
						neiCoords2[j] = neiCoords[j] - directionSign[j] * amr->getEpsilon(j) * 2;
					neiCoords2[exitDim] -= directionSign[exitDim] * amr->getEpsilon(exitDim);
					// printf("New nei coord = (%g %g %g) != (%g %g %g)\n",
					//        neiCoords2[0],neiCoords2[1],neiCoords2[2],
					//        neiCoords[0],neiCoords[1],neiCoords[2]);

					// reqLevel = nextVoxel->getLevel();
					ICoord index = amr->coords2Index(neiCoords2, reqLevel);

					// double dummy[NDIM]
					//  printf("Old Corners(%d) : (%g %g %g) (%g %g %g)\n",curVoxel->getLevel(),
					//         corner[0],corner[1],corner[2],
					//         oppCorner[0],oppCorner[1],oppCorner[2]);

					// amr->index2CornerCoordsAndOpp(index,reqLevel,corner,oppCorner,rayDir);
					getCornerCoords(index, reqLevel, corner, oppCorner);

					// printf("New Corners(%d) : (%g %g %g) (%g %g %g)\n",reqLevel,
					//        corner[0],corner[1],corner[2],
					//        oppCorner[0],oppCorner[1],oppCorner[2]);

					Voxel *tmpVoxel = advanceExact(corner, oppCorner);

					while (tmpVoxel == curVoxel)
					{
						// printf("LOOPING\n");
						neiCoords2[exitDim] += directionSign[exitDim] * amr->getEpsilon(exitDim) * 4;

						index = amr->coords2Index(neiCoords2, reqLevel);
						// amr->index2CornerCoordsAndOpp(index,reqLevel,
						// 				  corner,oppCorner,rayDir);
						getCornerCoords(index, reqLevel, corner, oppCorner);
						tmpVoxel = advanceExact(corner, oppCorner);
					}

					if (tmpVoxel == NULL)
					{
						PRINT_SRC_INFO(LOG_ERROR);
						glb::console->print<LOG_ERROR>("This should never happen, good luck debugging me ;)\n");
						exit(-1);
					}

					// printf("New exit point @ (%g %g %g) : exit@%d(%g)\n",
					//        exitPoint[0],exitPoint[1],exitPoint[2],exitDim,directionSign[exitDim]);

					nextVoxel = tmpVoxel;
					// nextVoxel->print(amr,"NEW-DEST:");

					// If this happens, we have to start over !!!!
					if (nextVoxel->getLevel() > reqLevel)
					{
						// printf("STARTING OVER (%d -> %d)!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n",
						// 	   reqLevel,nextVoxel->getLevel());

						// FIXME -> periodic boundary conditions ???
						// amr->index2CornerCoordsAndOpp(curVoxel->getIndex(),curVoxel->getLevel(),
						// 				  corner,oppCorner,rayDir);
						getCornerCoords(curVoxel, corner, oppCorner);
						reqLevel = nextVoxel->getLevel();
						i = -1;
						continue;
					}
				}
			}
		}
		/**
		if (findVertex(0.84780201735197679991,0.99999999999999988898,0,origin[0],origin[1],0,1.E-4))
		  {
		printf("EXACT COMP: %d (%e %e)(%e %e)(%d)\n",needExactComputation,
			   (double)direction_inv[0],(double)direction_inv[1],
			   (double)directionEpsilon[0],(double)directionEpsilon[1],invalidDir);
		  }
		  /**/
		if (needExactComputation)
			return advanceExact(corner, oppCorner);

#endif
#endif

		return nextVoxel;
	}

#ifdef HAVE_BOOST
#ifdef HAVE_GMP

	Voxel *advanceExact(Coord corner[NDIM], Coord oppCorner[NDIM])
	{
		// typedef boost::multiprecision::mpf_float mpfloat;

		if (mp_exact_dir_init)
		{
			for (int i = 0; i < NDIM; ++i)
			{
				mp_exact_direction[i] = hlp::numericStaticCast<mpfloat>(rayDirPoint[i]);
				mp_exact_direction[i] -= hlp::numericStaticCast<mpfloat>(initEntryPoint[i]);
				mp_exact_direction[i] = geometry->correctCoordsDiff(mp_exact_direction[i], i);
			}
			mp_exact_dir_init = false;
		}
		// bool debug = (curVoxel->index == 3368694746067226624L);
		// if (debug) printf("Using exact computations\n");

		// bool degenerate = false;
		// static long noNeed = 0;
		// static long count=0;
		// printf("NEEDED %ld/%ld exact computations\n",++count,noNeed);

		// typedef boost::multiprecision::mpf_float mpfloat;
		// mpfloat mp_dMin = amr->getMaxLength();

		// mpfloat mp_direction = direction[i];
		mpfloat mp_corner = hlp::numericStaticCast<mpfloat>(corner[0]);
		mpfloat mp_origin = hlp::numericStaticCast<mpfloat>(origin[0]);
		mpfloat mp_direction_ref = mp_exact_direction[0]; // direction[0];//mp_exact_direction[0];
		mpfloat mp_dMin = (mp_corner - mp_origin);
		exitDim = 0;
		/**
		if (debug)
		  {
		//std::streamsize precision = std::numeric_limits<cpp_dec_float_50>::digits10;
		std::cout.precision(30);

		std::cout<<"Origin : ("<<origin[0]<<","<<origin[1]<<","<<origin[2]<<")\n";
		std::cout<<"mp_Direction : ("<<mp_exact_direction[0]<<
		  ","<<mp_exact_direction[1]<<","<<mp_exact_direction[2]<<")\n";
		  }
		  /**/
		int iStart = 0;
		while ((invalidDir & (1 << iStart)))
		{
			++iStart;
			mp_corner = hlp::numericStaticCast<mpfloat>(corner[iStart]);
			mp_origin = hlp::numericStaticCast<mpfloat>(origin[iStart]);
			mp_direction_ref = mp_exact_direction[iStart]; // direction[iStart];//mp_exact_direction[iStart];
			mp_dMin = (mp_corner - mp_origin);
			exitDim = iStart;
		}
		++iStart;
		/**
		if (debug)
		  {
		std::cout<<"exitDim@start = "<<exitDim<<"\nmp_dMin["<<exitDim<<"]="<<mp_dMin<<"\n";
		//printf("exitDim@start = %d\nmp_dMin[%d]=%20.20e\n",
		//exitDim,exitDim,mp_dMin.convert_to<Coord>());
		  }
		  /**/
		for (int i = iStart; i < NDIM; ++i)
		{
			if (invalidDir & (1 << i))
				continue;

			mp_corner = hlp::numericStaticCast<mpfloat>(corner[i]);
			mp_origin = hlp::numericStaticCast<mpfloat>(origin[i]);
			mpfloat mp_direction = mp_exact_direction[i]; // direction[i];//mp_exact_direction[i];
			mpfloat mp_d = (mp_corner - mp_origin);
			// We use that so that we avoid dividing, which may not be exact even with
			// multiprecison
			mpfloat mp_d_dir = mp_d * mp_direction_ref;
			mpfloat mp_dMin_dir = mp_dMin * mp_direction;
			/**
			if (debug) {
			  std::cout<<"mp_d="<<mp_d<<"\n";
			  std::cout<<"mp_d_dir["<<i<<"]="<<mp_d_dir<<"\n";
			  std::cout<<"mp_dMin_dir["<<i<<"]="<<mp_dMin_dir<<"\n";
			  // printf("mp_d_dir[%d]=%20.20e\nmp_dMin_dir[%d]=%20.20e\n",
			  // 		  i,mp_d_dir.convert_to<Coord>(),
			  // 		  i,mp_dMin_dir.convert_to<Coord>());
			}
			/**/
			// this happens when we go through a vertex/edge/face
			if (mp_d_dir == mp_dMin_dir)
			{
				// if (debug) printf("Going through and edge/vertex/face @%d\n",i);
				//  These criteria comply with our simulation of simplicity convention
				//  DO NOT CHANGE THEM !
				if (rayDir & (1 << i))
				{
					// printf("EDGE!\n");
					// degenerate=true;
					//  That's good enough that way, no need for multiprecision ...
					Float dMin = mp_d.convert_to<Float>() * direction_inv[i];
					// printf("DEGENERATE(%d) -> %g\n",i,origin[i]+dMin*direction[i]+directionSign[i]*amr->getEpsilon(i));
					if (geometry->inBound(origin[i] + dMin * direction[i] + directionSign[i] * amr->getEpsilon(i), i))
					{
						// printf("ACCEPTED %d\n",i);
						exitDim = i;
						mp_direction_ref = mp_direction;
						mp_dMin = mp_d;
					}
				}
				// printf("NO EDGE!\n");
			}
			else if ((mp_direction_ref < 0) == (mp_direction < 0))
			{
				// printf("NEG!\n");
				// if (debug) printf("checking (mp_d_dir < mp_dMin_dir) (%20.20e/%20.20e)\n",
				// mp_direction_ref.convert_to<Coord>(),mp_direction.convert_to<Coord>());
				if (mp_d_dir < mp_dMin_dir)
				{
					// printf("YES!\n");
					// degenerate=false;
					mp_dMin = mp_d;
					exitDim = i;
					mp_direction_ref = mp_direction;
					// if (debug) printf("YES => newexit dim =%d\n",i);
				}
			}
			else
			{
				// printf("POS!\n");
				// if (debug) printf("checking (mp_d_dir > mp_dMin_dir) (%20.20e/%20.20e)\n",
				// mp_direction_ref.convert_to<Coord>(),mp_direction.convert_to<Coord>());
				if (mp_d_dir > mp_dMin_dir)
				{
					// printf("YES!\n");
					// degenerate=false;
					mp_dMin = mp_d;
					exitDim = i;
					mp_direction_ref = mp_direction;
					// if (debug) printf("YES => newexit dim =%d\n",i);
				}
			}
		}

		mp_dMin /= mp_direction_ref;

		Coord neiCoords[NDIM];
		// mpfloat mp_neiCoords[NDIM];
		//  Compute the exitPoint and prepare the coordinates for finding neighbor voxel
		for (int i = 0; i < NDIM; ++i)
		{
			mp_corner = hlp::numericStaticCast<mpfloat>(corner[i]);
			mp_origin = hlp::numericStaticCast<mpfloat>(origin[i]);
			mp_direction_ref = mp_exact_direction[i]; // direction[i];//mp_exact_direction[i];
			mpfloat mp_neiCoord = mp_origin + mp_dMin * mp_direction_ref;

			neiCoords[i] = mp_neiCoord.convert_to<Coord>();
			exitPoint[i] = geometry->template checkBoundary<Coord>(neiCoords[i], i);

			Float eps = amr->getEpsilon(i);
			if ((i != exitDim) && (fabs(neiCoords[i] - corner[i]) < eps))
			{
				neiCoords[i] -= directionSign[i] * eps;
			}
			else if ((i != exitDim) && (fabs(neiCoords[i] - oppCorner[i]) < eps))
			{
				neiCoords[i] += directionSign[i] * eps;
			}
		}

		// Now move the exit point coordinates by epsilon along the direction orthogonal to
		// the face it intersect to retrieve the neighbor
		neiCoords[exitDim] += directionSign[exitDim] * amr->getEpsilon(exitDim);
		// printf("Looking for exact next voxel @ : (%g %g %g)\n",neiCoords[0],neiCoords[1],neiCoords[2]);
		nextVoxel = amr->getVoxelAt(neiCoords);
		/*
		if (findVertex(exitPoint[0],exitPoint[1],exitPoint[2],0.248461,0.0859375,0.171875))
		  {
		printf("Degenerate exit point @ (%g %g %g) : exact(1) degenerate(%d), exit@%d(%g).\n",
			   exitPoint[0],exitPoint[1],exitPoint[2],
			   (int)degenerate,exitDim,directionSign[exitDim]);
		nextVoxel->print(amr,"E-DEST:");
		  }
		*/

		/*
		// We may be exiting a coarse voxel right on the edge of a finer one
		// in which case we may need exact computations ...
		if (nextVoxel->getLevel()>curVoxel->getLevel())
		  {
		int levelDiff = nextVoxel->getLevel() - curVoxel->getLevel();
		for (int i=0;i<NDIM;++i)
		  {
			if (i==exitDim) continue;

			mpfloat mp_oppCorner = oppCorner[i];
			mp_corner = corner[i];
			if (corner[i]>oppCorner[i])
			  std::swap(mp_corner,mp_oppCorner);

			mp_dMin = (1L<<levelDiff);
			mp_dMin /= (mp_oppCorner-mp_corner);
			mpfloat ratio = (mp_exitPoint[i]-mp_corner) * mp_dMin;
			mp_dMin = corner[i] + mp_dMin * ratio.convert_to<long>();
		  }
		  }
		*/
		return nextVoxel;
	}

#endif
#endif

private:
	bool findVertex(double x1, double y1, double z1, double vx, double vy, double vz, double tol = 2.E-5)
	{
		if ((fabs(vx - x1) <= tol) &&
			(fabs(vy - y1) <= tol) &&
			(fabs(vz - z1) <= tol))
		{
			return true;
		}
		return false;
	}

	// wrapper to amr->index2CornerCoordsAndOpp that deals with periodic boundaries
	void getCornerCoords(Voxel *voxel, Coord corner[NDIM], Coord oppCorner[NDIM])
	{
		getCornerCoords(voxel->getIndex(), voxel->getLevel(), corner, oppCorner);
	}

	// wrapper to amr->index2CornerCoordsAndOpp that deals with periodic boundaries
	void getCornerCoords(ICoord index, int level, Coord corner[NDIM], Coord oppCorner[NDIM])
	{
		amr->index2CornerCoordsAndOpp(index, level, corner, oppCorner, rayDir);
		if (AMR::BOUNDARY_TYPE == BoundaryType::PERIODIC)
		{
			// for PBC, the corners of the voxel are moved so that they are always in the
			// direction of the ray from the origin.
			for (int i = 0; i < NDIM; ++i)
			{
				if ((corner[i] < origin[i]) && (direction[i] > 0))
				{
					corner[i] += geometry->getBBoxSize(i);
					oppCorner[i] += geometry->getBBoxSize(i);
				}
				else if ((corner[i] > origin[i]) && (direction[i] < 0))
				{
					corner[i] -= geometry->getBBoxSize(i);
					oppCorner[i] -= geometry->getBBoxSize(i);
				}
			}
		}
	}

	AMR *amr;
	GeometricProperties *geometry;

	Float directionEpsilon[NDIM];
	Float anyEpsilon;

	Voxel *nextVoxel;	// The voxel the ray will enter next
	Voxel *curVoxel;	// The voxel for which entry/exit points were computed
	Float origin[NDIM]; // The point from which the ray started
	Coord entryPoint[NDIM];
	Coord exitPoint[NDIM];
	Coord direction[NDIM];	   // direction of the ray
	Coord direction_inv[NDIM]; // Inverse of direction
	Coord directionSign[NDIM];

	int rayDir;		// nth bit is 1/0 if the ray is going in the positive/negative direction
					// along dimension n
	int invalidDir; // if (invalidDir&(1<<n)) != 0 then the ray is parallel to that axis n
	int exitDim;

#ifdef HAVE_BOOST
#ifdef HAVE_GMP
	Coord rayDirPoint[NDIM];
	Coord initEntryPoint[NDIM];
	// boost::multiprecision::mpf_float mp_exact_direction[NDIM];
	mpfloat mp_exact_direction[NDIM];
	bool mp_exact_dir_init;
#endif
#endif

	// Coord exitDimSign;
};

/** \}*/

#ifdef DEBUG_ME
#undef DEBUG_ME
#endif

#include "../internal/namespace.footer"
#endif
