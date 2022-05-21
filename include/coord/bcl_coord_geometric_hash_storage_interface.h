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

#ifndef BCL_COORD_GEOMETRIC_HASH_STORAGE_INTERFACE_H_
#define BCL_COORD_GEOMETRIC_HASH_STORAGE_INTERFACE_H_

// include the namespace header
#include "bcl_coord.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class GeometricHashStorageInterface
    //! @brief This is a interface class for storing hashkeys with their bases
    //!
    //! @remarks example unnecessary
    //! @author haenigc, woetzen
    //! @date 28.06.2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API GeometricHashStorageInterface :
      public util::ObjectInterface
    {

    //////////
    // data //
    //////////

      //! shift bits of one key by that amount
      static const size_t s_SingleKeyShift = 8;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! virtual constructor from all information
      virtual GeometricHashStorageInterface *Construct
      (
        const std::string &MRC_NAME,
        const double MRC_RESOLUTION,
        const double MRC_VOXELSIZE,
        const storage::VectorND< 3, size_t> &MRC_EXTENSION,
        const size_t NUMBER_POINTS,
        const double FEATURE_DISTANCE,
        const double RATIO_INTENSITY_GRADIENT,
        const storage::VectorND< 4, double> &THRESHOLD,
        const double FEATURE_RADIUS,
        const PointToKeyInterface &POINT_TO_KEY,
        const size_t MIN_NUMBER_NEIGHBORS,
        const double MIN_NEIGHBOR_DISTANCE
      ) const = 0;

    /////////////////
    // data access //
    /////////////////

      //! virtual function to return the stored number of transformation matrices
      virtual size_t NumberOfTransformationMatrices() const = 0;

      //! virtual function to return the number of different Hashkeys
      virtual size_t NumberOfHashkeys() const = 0;

      //! virtual function returning the mrc name
      virtual std::string GetMRCName() const = 0;

      //! virtual function returning the mrc_resolution
      virtual double GetMRCResolution() const = 0;

      //! virtual function returning the mrc_voxelsize
      virtual double GetMRCVoxelSize() const = 0;

      //! virtual function returning mrc_extension
      virtual storage::VectorND< 3, size_t> GetMRCExtension() const = 0;

      //! virtual function returning number points used for building hash map
      virtual size_t GetNumberPoints() const = 0;

      //! virtual function returning minimal distance between two points
      virtual double GetMinimalDistance() const = 0;

      //! virtual function returning the ratio between the intensity and gradient for a voxel to be a feature within density map
      virtual double GetRatioIntensityGradient() const = 0;

      //! virtual function returning threshold for defining proper triangles
      virtual storage::VectorND< 4, double> GetThreshold() const = 0;

      //! virtual function returning the feature radius
      virtual double GetRadius() const = 0;

      //! virtual function returning the point to key function
      virtual util::ShPtr< PointToKeyInterface> GetPointToKey() const = 0;

      //! virtual function returning the min number of neighbors for a point in the point cloud
      virtual size_t GetMinNumberNeighbors() const = 0;

      //! virtual function returning the min distance for the neighbors
      virtual double GetMinNeighborDistance() const = 0;

      //! @brief virtual function finalizing the map
      virtual void Finalize() = 0;

      //! virtual function setting read only
      virtual bool IsReadOnly() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief pure virtual function to store a set of coordinates with the according Transformatiomatrix
      //! @param COORDINATES the set of coordinates that is associated with the transformation
      //! @param TRANSFORMATIONMATRIX3D the transformation for that set of coordinates
      virtual void Store
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const math::TransformationMatrix3D &TRANSFORMATIONMATRIX3D
      ) = 0;

      //! @brief pure virtual function to return the best counting bases
      //! @param COORDINATES the coordinates to be considered
      //! @param TRANSFORMATIONMATRIX3D the transformation for that set of coordinates
      //! @param NUMBERBESTCOUNTS the max number of transformations to return
      //! @return the list of transformations with their hash score
      virtual storage::List< storage::Pair< util::ShPtr< math::TransformationMatrix3D>, size_t> >
      ReturnBestCounts
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const math::TransformationMatrix3D &TRANSFORMATIONMATRIX3D,
        const size_t NUMBERBESTCOUNTS
      ) const = 0;

      //! checks if the TRANSFORMATIONMATRIX3D is similar to any stored transformationmatrix according to DIFFERENCE
      virtual bool IsSimilarTransformation
      (
        const math::TransformationMatrix3D &TRANSFORMATIONMATRIX3D,
        const storage::VectorND< 2, double> &DIFFERENCE_ROT_TRANS
      ) const = 0;

      //! This function converts triplet of ints to a hashkey by shifting each integer by certain amount, so that in the binary representation these three integers do not overlap, given that the ints are smaller than pow(2,s_SingleKeyShift)-1
      static size_t ConvertTripletToKey( const storage::Triplet< int, int, int> &HASHKEYS);

      //! checks if the TRANSFORAMTIONS are similar to any stored transformation matrix according to DIFFERENCE_ROT_TRANS
      static bool IsSimilarTransformation
      (
        const util::ShPtrVector< math::TransformationMatrix3D> &TRANSFORAMTIONS,
        const math::TransformationMatrix3D &TRANSFORMATIONMATRIX3D,
        const storage::VectorND< 2, double> &DIFFERENCE_ROT_TRANS
      );

    }; //class GeometricHashStorageInterface

  } // namespace coord
} // namespace bcl

#endif //BCL_COORD_GEOMETRIC_HASH_STORAGE_INTERFACE_H_
