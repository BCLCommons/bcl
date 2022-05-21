// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

#ifndef BCL_COORD_GEOMETRIC_HASHING_H_
#define BCL_COORD_GEOMETRIC_HASHING_H_

// include the namespace header
#include "bcl_coord.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_coord_geometric_hash_storage_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class GeometricHashing
    //! @brief This class is a hashmap class
    //! @details it does two things: builds a hashmap from a given PointCloud. This PointCloud could be generated
    //! from a density map.
    //!
    //! there are different things you have to know about geometric hashing to understand whats going on
    //! in this class. So I recommend two papers: Haim J. Wolfson "Geometric Hashing: An Overview" 1997
    //! Computational Science and Engineering. You can find it here: /sb/meiler/papers/coord/wr-ghao-97.pdf
    //!
    //! we have different approaches and made decisions which should reduce the size of the hashmap and
    //! so the speed quality. /n
    //! A base of three points is used. You pass a threshold which defines the minimum and maximum length
    //! of each base. The points of base is ordered by the length of its opposite side. The length of the
    //! sides have to have different length. The shortest is in the range min to min + 1/3(max-min),
    //! the middle one is from 1/3 to 2/3 and the longest one is in the last 3/3 range.
    //! All pairwise distances of the points are calculated ( function CalculatePointPairs) and if they
    //! are in one of the ranges, they are stored in three lists. /n
    //! Then for each possible base a transformation is calculated. The middlepoint of the triangle is
    //! moved to the origin, than vector 1 is roted in the xz plane and than to the x axes. After that
    //! vector 2 is rotated in the xy plane.
    //! The coordinates of the remaining points, which are in the max threshold of each base are calculated
    //! according to this transformation. /n
    //! This coordinates are used to generate the hashkey, which contains the distance to the origin,
    //! the theta and phi angle. (function HashKey) /n
    //! For each HashKey the base is pushbacked to the hashmap, which has a Sharedpointer of Vector3D behind each Key. /n
    //! /n
    //! This building of the Hashmap takes the longest time.
    //! /n
    //! The second part of this class is the possibility to search the hashmap for a given pointcloud.
    //! function SearchTarget( ..) /n
    //! There are different things which are important: /n
    //! You have to pass a PointCloud of the target. You have to define how many random searches you want to try.
    //! You can define how many of the best counts you want to get back. You may also say, whether you
    //! want to search for multiple targets. If you expect to have more than one hit in the hashmap.
    //! The search will then also accept low counts when they yield to different positions in the pointcloud,
    //! otherwise it will only store the matches with the highest counts in the hashmap. The less populated area
    //! would get no match. /n
    //! The result of the search is a storagevector of transformationmatrices. /n
    //!
    //! @see @link example_coord_geometric_hashing.cpp @endlink
    //! @author woetzen
    //! @date 16.11.2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API GeometricHashing :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      util::ShPtr< GeometricHashStorageInterface> m_HashMap; //!< storage object for geometric hashmap

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //!default Constructor
      GeometricHashing();

      //! construct Hashmap from GeometricHashStorage
      GeometricHashing( const util::ShPtr< GeometricHashStorageInterface> &STORAGE_SP);

      //! virtual copy constructor
      GeometricHashing *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief build Hashmap from PointCloud
      //! @param FEATURE_CLOUD all points that describe the object the features encode
      void BuildHash( const PointCloud &FEATURE_CLOUD);

      //! @brief search hash map and return storage list of pairs of Transformation matrices for the highest counts
      //! @param COORDINATES the coordinates that will be searched in the hash map
      //! @param SAVE_BEST the maximal number of possible transformations returned
      //! @param SEARCH_TRIALS is the number of random searched bases, SAVEBEST is the number of transformationmatrices which are stored
      //! @param DIFFERENCE_ROT_TRANS has two values for the max deviation of a new found transformation in rotation and translation in comparison to already found transformations
      storage::List< storage::Pair< math::TransformationMatrix3D, size_t> >
      SearchTarget
      (
        const PointCloud &COORDINATES,
        const size_t &SAVE_BEST,
        const size_t &SEARCH_TRIALS,
        const storage::VectorND< 2, double> &DIFFERENCE_ROT_TRANS
      ) const;

      //! @brief returns all possible bases within the given pointcloud
      //! @param COORDINATES the coordinates in the hash map
      storage::List< storage::Pair< math::TransformationMatrix3D, linal::Vector3D> >
      ReturnPossibleBases
      (
        const PointCloud &COORDINATES
      ) const;

      //! @brief return the best transformations and their hashscore for one base within a set of coordinates
      //! @param TRANSFORMATION_CENTER the base and center to be considered within the COORDINATES
      //! @param COORDINATES the set of coordinates
      //! @param SAVE_BEST the maximal number of results to be returned
      //! @return a list of Transformations and their hashscore
      storage::List< storage::Pair< util::ShPtr< math::TransformationMatrix3D>, size_t> >
      ReturnBestCounts
      (
        const storage::Pair< math::TransformationMatrix3D, linal::Vector3D> &TRANSFORMATION_CENTER,
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const size_t &SAVE_BEST
      ) const;

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! write Hash to std::ostream using the given util::Format
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! read Hash from std::ifstream
      std::istream &Read( std::istream &ISTREAM);

    ////////////////
    // operations //
    ////////////////

    public:

      //! @brief calculate Transformation, so that 3 Vectors are the Basis for a coordinate System and their center
      //! @param VECTOR_1 is connected to VECTOR_2 by long side to VECTOR_3 by short side
      //! @param VECTOR_2 is connected to VECTOR_3 by middle side to VECTOR_1 by long side
      //! @param VECTOR_3 is connected to VECTOR_1 by short side to VECTOR_2 by middle side
      //! @return the Transformation matrix, that puts center of triangle in origin, Vector2 on x-axis and vector 3 in xy-plane and the center of the untransformed vectors
      static
      storage::Pair< math::TransformationMatrix3D, linal::Vector3D>
      CalculateBasisTransformationAndCenter
      (
        const linal::Vector3D &VECTOR_1,
        const linal::Vector3D &VECTOR_2,
        const linal::Vector3D &VECTOR_3
      );

      //! @brief return pairwise distances and order them in three different sizes
      //! @param COORDINATES coordinates
      //! @param THRESHOLD the threshold that defines the different lengths of the triangle used as a base
      //! last entry is a list of points, that contains a list of points that is within the radius of the point that is key in the map
      //! @return three maps, where the first contains for each point a set of all other points that are close by the shortest distance and so on
      static
      storage::VectorND< 3, storage::Map< util::SiPtr< const linal::Vector3D>, storage::Set< util::SiPtr< const linal::Vector3D> > > >
      CalculatePointPairs
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const storage::VectorND< 4, double> &THRESHOLD
      );

      //! @brief determine triangular bases given a map of all point pairs
      //! @param POINT_PAIRS point pairs for lower, middle an upper threshold
      //! @return a list of possible triangles
      static
      storage::List< storage::VectorND< 3, util::SiPtr< const linal::Vector3D> > >
      DeterminePossibleBases
      (
        const storage::VectorND< 3, storage::Map< util::SiPtr< const linal::Vector3D>, storage::Set< util::SiPtr< const linal::Vector3D> > > > &POINT_PAIRS
      );

      //! @brief pick a subset of random bases from a list of triangular bases which are equally distributed over different distances from the given center
      //! @param POINT_PAIRS point pairs for lower, middle an upper threshold
      //! @param NUMBER size of the resulting set
      //! @param CENTER the center to which the bases should be distributed equally distant
      //! @param NUMBER_BINS number of bin to distribute bases equally over
      //! @return List of transformations and their centers
      static
      storage::List< storage::Pair< math::TransformationMatrix3D, linal::Vector3D> >
      GatherEquallyDistributedBases
      (
        const storage::VectorND< 3, storage::Map< util::SiPtr< const linal::Vector3D>, storage::Set< util::SiPtr< const linal::Vector3D> > > > &POINT_PAIRS,
        const size_t NUMBER,
        const linal::Vector3D &CENTER,
        const size_t NUMBER_BINS
      );

      //! @brief pick a subset of random bases from a list of triangular bases which are equally distributed over different distances from the given center
      //! @param POINT_PAIRS point pairs for lower, middle an upper threshold
      //! @param NUMBER size of the resulting set
      //! @param CENTER the center to which the bases should be distributed equally distant
      //! @param NUMBER_BINS number of bin to distribute bases equally over
      //! @return List of transformations and their centers
      storage::List< storage::Pair< math::TransformationMatrix3D, linal::Vector3D> >
      GatherEquallyDistributedBases3D
      (
        const storage::VectorND< 3, storage::Map< util::SiPtr< const linal::Vector3D>, storage::Set< util::SiPtr< const linal::Vector3D> > > > &POINT_PAIRS,
        const size_t NUMBER,
        const linal::Vector3D &CENTER,
        const size_t NUMBER_BINS
      ) const;

      //! @brief convert a list of triangular bases into transformations and their triangular centers
      static
      storage::List< storage::Pair< math::TransformationMatrix3D, linal::Vector3D> >
      ConvertBasesToTransformationsAndCenter
      (
        const storage::List< storage::VectorND< 3, util::SiPtr< const linal::Vector3D> > > &BASES
      );

      //! @brief filter similar TransforamtionMatrices according to DIFFERENCE
      //! @param BESTTRANSFORMATIONS the transformations to be filtered
      //! @param KEEPBEST the number of transformations to keep
      //! @param DIFFERENCE_ROT_TRANS has two values for the max deviation of a new found transformation in rotation and translation in comparison to already found transformations
      //! @return a List that contains only transformations that are different in DIFFERENCE_ROT_TRANS
      static
      storage::List< storage::Pair< math::TransformationMatrix3D, size_t> >
      FilterSimilarTransformations
      (
        const storage::List< storage::Pair< util::ShPtr< math::TransformationMatrix3D>, size_t> > &BESTTRANSFORMATIONS,
        const size_t &KEEPBEST,
        const storage::VectorND< 2, double> &DIFFERENCE_ROT_TRANS
      );

      //! @brief for a set of coordinates, return the maximal number of bases
      //! @param COORDINATES the coordinates the number of bases is determined for
      //! @param THRESHOLD the threshold used for the triangular base
      //! @return number of bases which fullfill the given THRESHOLD
      static
      size_t
      NumberPossibleBases
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const storage::VectorND< 4, double> &THRESHOLD
      );

      //! @brief for pointpairs, return the maximal number of bases
      //! @param POINT_PAIRS point pairs for lower, middle an upper threshold
      //! @return number of bases which fullfill the given THRESHOLD
      static
      size_t
      NumberPossibleBases
      (
        const storage::VectorND< 3, storage::Map< util::SiPtr< const linal::Vector3D>, storage::Set< util::SiPtr< const linal::Vector3D> > > > &POINT_PAIRS
      );

      //! @brief determine all neighbors within given a radius
      //! @param POINT_OF_INTEREST the point for which all neighbors within the FEATURE_RADIUS is returned
      //! @param COORDINATES the coordinates that are considered
      //! @param FEATURE_RADIUS
      //! @return vector of coordinates within the FEATURE_RADIUS of the POINT_OF_INTEREST
      static
      util::SiPtrVector< const linal::Vector3D>
      DetermineNeighbors
      (
        const linal::Vector3D &POINT_OF_INTEREST,
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const double FEATURE_RADIUS
      );

      //! @brief calculate the distance distribution within a given Range
      //! @param COORDINATES the coordinates that are considered
      //! @param THRESHOLD for lower and upper boundary of histogram
      //! @param NUMBER_BINS number of bins in the histogram generated
      //! @return histogram of for all the distances
      static
      math::Histogram
      CalculateDistanceDistribution
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const storage::VectorND< 2, double> &THRESHOLD,
        const size_t NUMBER_BINS = 20
      );

      //! @brief calculate the distance distribution
      //! @param COORDINATES the coordinates that are considered
      //! @param BIN_SIZE the size of the bins
      //! @param MAX_DISTANCE the maximal distance to be considered
      //! @return histogram of for all the distances
      static
      math::Histogram
      CalculateDistanceDistribution
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const double BIN_SIZE,
        const math::Range< double> &DISTANCE_RANGE
      );

      //! @brief create three equidistant ranges from two thresholds
      //! @param THRESHOLD lower and upper boundaries for ranges
      //! @return four threshold bounding three equidistant intervals
      static
      storage::VectorND< 4, double>
      CreateEquidistantIntervals( const storage::VectorND< 2, double> &THRESHOLD);

      //! @brief create three ranges containing approximately equal number of distances
      //! @param DISTANCE_DISTRIBUTION histogram of distance distributions
      //! @return four threshold bounding three intervals with equal number of distances
      static
      storage::VectorND< 4, double>
      CreateEqualOccupiedIntervals( const math::Histogram &DISTANCE_DISTRIBUTION);

    }; // class GeometricHashing

  } // namespace coord
} // namespace bcl

#endif //BCL_COORD_GEOMETRIC_HASHING_H_
