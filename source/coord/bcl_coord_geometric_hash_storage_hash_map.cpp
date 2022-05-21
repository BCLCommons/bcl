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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "coord/bcl_coord_geometric_hash_storage_hash_map.h"

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> GeometricHashStorageHashMap::s_Instance
    (
      GetObjectInstances().AddInstance( new GeometricHashStorageHashMap())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default private constructor
    GeometricHashStorageHashMap::GeometricHashStorageHashMap() :
      GeometricHashStorageInterface(),
      m_MRCName(),
      m_MRCResolution(),
      m_MRCVoxelSize(),
      m_MRCExtension(),
      m_NumberPoints(),
      m_MinimalDistance(),
      m_Threshold(),
      m_Radius(),
      m_PointToKey(),
      m_MinNumberNeighbors(),
      m_MinNeighborDistance(),
      m_Final( false)
    {
    }

    //! @brief construct from all informations
    GeometricHashStorageHashMap::GeometricHashStorageHashMap
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
    ) :
      GeometricHashStorageInterface(),
      m_MRCName( MRC_NAME),
      m_MRCResolution( MRC_RESOLUTION),
      m_MRCVoxelSize( MRC_VOXELSIZE),
      m_MRCExtension( MRC_EXTENSION),
      m_NumberPoints( NUMBER_POINTS),
      m_MinimalDistance( FEATURE_DISTANCE),
      m_RatioIntensityGradient( RATIO_INTENSITY_GRADIENT),
      m_Threshold( THRESHOLD),
      m_Radius( FEATURE_RADIUS),
      m_PointToKey( POINT_TO_KEY.Clone()),
      m_MinNumberNeighbors( MIN_NUMBER_NEIGHBORS),
      m_MinNeighborDistance( MIN_NEIGHBOR_DISTANCE),
      m_Final( false)
    {
    }

    //! virtual constructor from all information
    GeometricHashStorageInterface *GeometricHashStorageHashMap::Construct
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
    ) const
    {
      return new GeometricHashStorageHashMap
                 (
                   MRC_NAME,
                   MRC_RESOLUTION,
                   MRC_VOXELSIZE,
                   MRC_EXTENSION,
                   NUMBER_POINTS,
                   FEATURE_DISTANCE,
                   RATIO_INTENSITY_GRADIENT,
                   THRESHOLD,
                   FEATURE_RADIUS,
                   POINT_TO_KEY,
                   MIN_NUMBER_NEIGHBORS,
                   MIN_NEIGHBOR_DISTANCE
                 );
    }

    //! @brief copy constructor
    //! @return pointer to a new GeometricHashStorageHashMap copied from this one
    GeometricHashStorageHashMap *GeometricHashStorageHashMap::Clone() const
    {
      return new GeometricHashStorageHashMap( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GeometricHashStorageHashMap::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief store a set of coordinates with the according Transformatiomatrix
    //! @param COORDINATES the set of coordinates that is associated with the transformation
    //! @param TRANSFORMATIONMATRIX3D the transformation for that set of coordinates
    void
    GeometricHashStorageHashMap::Store
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const math::TransformationMatrix3D &TRANSFORMATIONMATRIX3D
    )
    {
      BCL_Assert( !m_Final, "unable to store into a map, that is already final");
      // store the transformation matrix
      m_TransformationMatrices.PushBack( util::ShPtr< math::TransformationMatrix3D>( TRANSFORMATIONMATRIX3D.Clone()));

      // util::SiPtr to the Last Transformationmatrix3D in this Vector
      const util::SiPtr< const math::TransformationMatrix3D> sp_transformationmatrix3D( m_TransformationMatrices.LastElement());

      // iterate over all coordinates
      for
      (
        util::SiPtrVector< const linal::Vector3D>::const_iterator
          coord_itr( COORDINATES.Begin()), coord_itr_end( COORDINATES.End());
        coord_itr != coord_itr_end;
        ++coord_itr
      )
      {
        // quantize the coordinate according to the point to key function
        const storage::Triplet< int, int, int> triplet_hash_key
        (
          m_PointToKey->operator()( linal::Vector3D( **coord_itr).Transform( TRANSFORMATIONMATRIX3D))
        );

        // generate hash key
        const size_t hash_key( ConvertTripletToKey( triplet_hash_key));

        // increase the count for that key and that base by 1
        ++m_Hash[ hash_key][ sp_transformationmatrix3D];
      }
    }

    //! @brief return the best counting bases
    //! @param COORDINATES the coordinates to be considered
    //! @param NUMBERBESTCOUNTS the max number of transformations to return
    //! @return the list of transformations with their hash score
    storage::List< storage::Pair< util::ShPtr< math::TransformationMatrix3D>, size_t> >
    GeometricHashStorageHashMap::ReturnBestCounts
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const math::TransformationMatrix3D &TRANSFORMATIONMATRIX3D,
      const size_t NUMBERBESTCOUNTS
    ) const
    {
      BCL_Assert( m_Final, "unable to count within an unfinished map");

      // store the number of occurrences for each of the keys that come from the set of coordinates
      storage::Map< size_t, size_t> keys_count;

      // iterate over all coordinates
      for
      (
        util::SiPtrVector< const linal::Vector3D>::const_iterator
          coord_itr( COORDINATES.Begin()), coord_itr_end( COORDINATES.End());
        coord_itr != coord_itr_end;
        ++coord_itr
      )
      {
        // quantize the coordinate according to the point to key function
        const storage::Triplet< int, int, int> triplet_hash_key
        (
          m_PointToKey->operator()( linal::Vector3D( **coord_itr).Transform( TRANSFORMATIONMATRIX3D))
        );

        // generate hash key
        const size_t hash_key( ConvertTripletToKey( triplet_hash_key));

        ++keys_count[ hash_key];
      }

      std::map< util::SiPtr< const math::TransformationMatrix3D>, size_t> basis_count;

      // iterate over all the keys in the target
      for
      (
        storage::Map< size_t, size_t>::const_iterator
          keys_count_itr( keys_count.Begin()), keys_count_itr_end( keys_count.End());
        keys_count_itr != keys_count_itr_end;
        ++keys_count_itr
      )
      {
        // find hash key and count the matches in the histogram
        const storage::HashMap< size_t, storage::Map< util::SiPtr< const math::TransformationMatrix3D>, size_t> >::const_iterator
          itr( m_Hash.Find( keys_count_itr->first));

        // add for each key the number of bases associated with it
        if( itr != m_Hash.End())
        {
          for
          (
            storage::Map< util::SiPtr< const math::TransformationMatrix3D>, size_t>::const_iterator
              sp_tm_itr( itr->second.Begin()), sp_tm_itr_end( itr->second.End());
            sp_tm_itr != sp_tm_itr_end;
            ++sp_tm_itr
          )
          {
            basis_count[ sp_tm_itr->first] += std::min( sp_tm_itr->second, keys_count_itr->second);
          }
        }
      }

      // these will be the best transformationmatrices sorted from highest to lowest count
      storage::List< storage::Pair< util::ShPtr< math::TransformationMatrix3D>, size_t> > transformationmatrices;

      const std::multiset< std::pair< util::SiPtr< const math::TransformationMatrix3D>, size_t>, Sp_CountCompare>
        sorted_basis_count( basis_count.begin(), basis_count.end());
      size_t count( 0);

      // insert the best NUMBERBESTCOUNTS matrices to the best matrices beginning from the highest count
      for
      (
        std::multiset
        <
          std::pair< util::SiPtr< const math::TransformationMatrix3D>, size_t>, Sp_CountCompare
        >::const_reverse_iterator
          sorted_basis_revitr( sorted_basis_count.rbegin()),
          sorted_basis_revitr_end( sorted_basis_count.rend());
        sorted_basis_revitr != sorted_basis_revitr_end && count < NUMBERBESTCOUNTS; ++sorted_basis_revitr, ++count)
      {
        transformationmatrices.PushBack
        (
          storage::Pair< util::ShPtr< math::TransformationMatrix3D>, size_t>
          (
            util::ShPtr< math::TransformationMatrix3D>( sorted_basis_revitr->first->Clone()),
            sorted_basis_revitr->second
          )
        );
      }

      //return best transformationmatrices
      return transformationmatrices;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! read Hash from io::IFStream
    std::istream &GeometricHashStorageHashMap::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_TransformationMatrices, ISTREAM);
      io::Serialize::Read( m_Hash, ISTREAM);
      io::Serialize::Read( m_MRCName, ISTREAM);
      io::Serialize::Read( m_MRCResolution, ISTREAM);
      io::Serialize::Read( m_MRCVoxelSize, ISTREAM);
      io::Serialize::Read( m_MRCExtension, ISTREAM);
      io::Serialize::Read( m_NumberPoints, ISTREAM);
      io::Serialize::Read( m_MinimalDistance, ISTREAM);
      io::Serialize::Read( m_RatioIntensityGradient, ISTREAM);
      io::Serialize::Read( m_Threshold, ISTREAM);
      io::Serialize::Read( m_Radius, ISTREAM);
      io::Serialize::Read( m_PointToKey, ISTREAM);
      io::Serialize::Read( m_MinNumberNeighbors, ISTREAM);
      io::Serialize::Read( m_MinNeighborDistance, ISTREAM);

      // set final to true
      m_Final = true;

      // return
      return ISTREAM;
    }

    //! write Hash to std::ostream using the given util::Format
    std::ostream &GeometricHashStorageHashMap::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // check that this is a final geometric hash storage
      BCL_Assert( m_Final, "writing of unfinished map not supported");

      // write member
      io::Serialize::Write( m_TransformationMatrices, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Hash, OSTREAM, INDENT)                   << '\n';
      io::Serialize::Write( m_MRCName, OSTREAM, INDENT)                << '\n';
      io::Serialize::Write( m_MRCResolution, OSTREAM, INDENT)          << '\n';
      io::Serialize::Write( m_MRCVoxelSize, OSTREAM, INDENT)           << '\n';
      io::Serialize::Write( m_MRCExtension, OSTREAM, INDENT)           << '\n';
      io::Serialize::Write( m_NumberPoints, OSTREAM, INDENT)           << '\n';
      io::Serialize::Write( m_MinimalDistance, OSTREAM, INDENT)        << '\n';
      io::Serialize::Write( m_RatioIntensityGradient, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Threshold, OSTREAM, INDENT)              << '\n';
      io::Serialize::Write( m_Radius, OSTREAM, INDENT)                 << '\n';
      io::Serialize::Write( m_PointToKey, OSTREAM, INDENT)             << '\n';
      io::Serialize::Write( m_MinNumberNeighbors, OSTREAM, INDENT)     << '\n';
      io::Serialize::Write( m_MinNeighborDistance, OSTREAM, INDENT)    << '\n';

      //return
      return OSTREAM;
    }

  } // namespace coord
} // namespace bcl

