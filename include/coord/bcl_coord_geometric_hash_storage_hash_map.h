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

#ifndef BCL_COORD_GEOMETRIC_HASH_STORAGE_HASH_MAP_H_
#define BCL_COORD_GEOMETRIC_HASH_STORAGE_HASH_MAP_H_

// include the namespace header
#include "bcl_coord.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_coord_geometric_hash_storage_interface.h"
#include "bcl_coord_point_to_key_interface.h"
#include "math/bcl_math_transformation_matrix_3d.h"
#include "storage/bcl_storage_hash_map.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class GeometricHashStorageHashMap
    //! @brief This is one implementation for the GeometricHashStorage using a memory HashMap
    //!
    //! @see @link example_coord_geometric_hash_storage_hash_map.cpp @endlink
    //! @author woetzen
    //! @date 28.06.2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API GeometricHashStorageHashMap :
      public GeometricHashStorageInterface
    {
    private:

    //////////
    // data //
    //////////

      //! ShPtrVector of all Transformationmatrices
      util::ShPtrVector< math::TransformationMatrix3D>                                   m_TransformationMatrices;
      //! Hash map of size_t which is build from the keys related to each TransformationMatrix3D and stores the number of occurrences for the key in each base
      storage::HashMap< size_t, storage::Map< util::SiPtr< const math::TransformationMatrix3D>, size_t> > m_Hash;

      std::string                   m_MRCName;       //!< name of density map
      double                        m_MRCResolution; //!< resolution of density map
      double                        m_MRCVoxelSize;  //!< voxel size of density map
      storage::VectorND< 3, size_t> m_MRCExtension;  //!< extension of density map
      size_t                        m_NumberPoints;  //!< number of points representing features in density map
      double                        m_MinimalDistance; //!< minimal distance between two points in feature cloud
      double                        m_RatioIntensityGradient; //!< ratio between intensity and gradient that was used to build the map
      storage::VectorND< 4, double> m_Threshold; //!< thresholds for triangular base
      double                        m_Radius;    //!< feature radiu below which features are quantized for each triangular base
      util::ShPtr< PointToKeyInterface> m_PointToKey; //!< the point to key function used
      size_t                        m_MinNumberNeighbors; //!< minimal number of features if filtering was used
      double                        m_MinNeighborDistance; //!< minimal distance for those features

      bool                          m_Final; //!< indicate if the storage was finalized

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default private constructor
      GeometricHashStorageHashMap();

      //! @brief construct from all informations
      GeometricHashStorageHashMap
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
      );

      //! virtual constructor from all information
      GeometricHashStorageInterface *Construct
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
      ) const;

      //! @brief copy constructor
      //! @return pointer to a new GeometricHashStorageHashMap copied from this one
      GeometricHashStorageHashMap *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! return the stored number of Transformationmatrices
      size_t NumberOfTransformationMatrices() const
      {
        return m_TransformationMatrices.GetSize();
      }

      //! return the number of different Hashkeys
      size_t NumberOfHashkeys() const
      {
        return m_Hash.GetSize();
      }

      //! return the mrc name
      std::string GetMRCName() const
      {
        return m_MRCName;
      }

      //! return the mrc_resolution
      double GetMRCResolution() const
      {
        return m_MRCResolution;
      }

      //! return the mrc_voxelsize
      double GetMRCVoxelSize() const
      {
        return m_MRCVoxelSize;
      }

      //! return mrc_extension
      storage::VectorND< 3, size_t> GetMRCExtension() const
      {
        return m_MRCExtension;
      }

      //! return number points used for building hashmap
      size_t GetNumberPoints() const
      {
        return m_NumberPoints;
      }

      //! return minimal distance between two points
      double GetMinimalDistance() const
      {
        return m_MinimalDistance;
      }

      //! function returning the weight for the intensity for a voxel to be a point within electron density
      double GetRatioIntensityGradient() const
      {
        return m_RatioIntensityGradient;
      }

      //! return lower and upper threshold for triangles
      storage::VectorND< 4, double> GetThreshold() const
      {
        return m_Threshold;
      }

      //! return the feature radius
      double GetRadius() const
      {
        return m_Radius;
      }

      //! virtual function returning the point to key function
      util::ShPtr< PointToKeyInterface> GetPointToKey() const
      {
        return m_PointToKey;
      }

      //! function returning the min number of neighbors for a point in the feature cloud
      size_t GetMinNumberNeighbors() const
      {
        return m_MinNumberNeighbors;
      }

      //! function return the min distance for the neighbors
      double GetMinNeighborDistance() const
      {
        return m_MinNeighborDistance;
      }

      //! @brief virtual function finalizing the map
      void Finalize()
      {
        m_Final = true;
      }

      //! virtual function setting read only
      bool IsReadOnly() const
      {
        return m_Final;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief store a set of coordinates with the according Transformatiomatrix
      //! @param COORDINATES the set of coordinates that is associated with the transformation
      //! @param TRANSFORMATIONMATRIX3D the transformation for that set of coordinates
      void Store
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const math::TransformationMatrix3D &TRANSFORMATIONMATRIX3D
      );

      //! @brief return the best counting bases
      //! @param COORDINATES the coordinates to be considered
      //! @param TRANSFORMATIONMATRIX3D the transformation for that set of coordinates
      //! @param NUMBERBESTCOUNTS the max number of transformations to return
      //! @return the list of transformations with their hash score
      storage::List< storage::Pair< util::ShPtr< math::TransformationMatrix3D>, size_t> >
      ReturnBestCounts
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const math::TransformationMatrix3D &TRANSFORMATIONMATRIX3D,
        const size_t NUMBERBESTCOUNTS
      ) const;

      //! checks if the TRANSFORMATIONMATRIX3D is similar to any stored transformation matrix according to DIFFERENCE
      bool IsSimilarTransformation
      (
        const math::TransformationMatrix3D &TRANSFORMATIONMATRIX3D,
        const storage::VectorND< 2, double> &DIFFERENCE_ROT_TRANS
      ) const
      {
        return GeometricHashStorageInterface::IsSimilarTransformation
        (
          m_TransformationMatrices, TRANSFORMATIONMATRIX3D, DIFFERENCE_ROT_TRANS
        );
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! read Hash from io::IFStream
      std::istream &Read( std::istream &ISTREAM);

      //! write Hash to std::ostream using the given util::Format
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class Sp_CountCompare
      //! @brief compares only the second element of the pair
      //! @remarks example unnecessary
      //! @author woetzen
      //! @date 28.06.2006
      //!
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      struct Sp_CountCompare :
        public std::binary_function< std::pair< util::SiPtr< const math::TransformationMatrix3D>, size_t>, std::pair< util::SiPtr< const math::TransformationMatrix3D>, size_t>, bool>
      {
        bool operator()
        (
          const std::pair< util::SiPtr< const math::TransformationMatrix3D>, size_t> &PAIR1,
          const std::pair< util::SiPtr< const math::TransformationMatrix3D>, size_t> &PAIR2
        ) const
        {
          return PAIR1.second < PAIR2.second;
        }
      };

    }; // class GeometricHashStorageHashMap

  } // namespace coord
} // namespace bcl

#endif //BCL_COORD_GEOMETRIC_HASH_STORAGE_HASH_MAP_H_
