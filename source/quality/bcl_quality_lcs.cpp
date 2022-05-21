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
#include "quality/bcl_quality_lcs.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_comparisons.h"
#include "math/bcl_math_range.h"
#include "math/bcl_math_transformation_matrix_3d.h"
#include "quality/bcl_quality_rmsd.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_logger_interface.h"
#include "util/bcl_util_si_ptr_vector.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace quality
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> LCS::s_Instance
    (
      GetObjectInstances().AddInstance( new LCS())
    );

    //! @brief returns the default Rmsd cutoff
    //! @return the default Rmsd cutoff
    double LCS::GetDefaultRmsdCutoff()
    {
      return 5.0;
    }

    //! @brief returns the default seed length
    //! @return the default seed length
    size_t LCS::GetDefaultSeedLength()
    {
      return 3;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from a RMSD cutoff and a seed length
    //! @param RMSD_CUTOFF distance cutoff
    //! @param SEED_LENGTH length of seeds
    LCS::LCS
    (
      const double RMSD_CUTOFF,
      const size_t SEED_LENGTH
    ) :
      m_RmsdCutoff( RMSD_CUTOFF),
      m_SeedLength( SEED_LENGTH)
    {
      BCL_Assert( SEED_LENGTH >= 3, "The seed length has to be at least 3 not :" + util::Format()( SEED_LENGTH));
    }

    //! @brief Clone function
    //! @return pointer to new LCS
    LCS *LCS::Clone() const
    {
      return new LCS( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &LCS::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the optimal value for that quality measurement
    //! @return the best value by which two sets of coordinates can agree
    double LCS::OptimalValue() const
    {
      return util::GetUndefined< double>();
    }

    //! @brief return the comparison function for better quality
    //! @return binary function to compare two quality measure values
    const util::BinaryFunctionInterface< double, double, bool> &LCS::GetComparisonFunction() const
    {
      return **math::Comparisons< double>::GetEnums().e_Greater;
    }

    //! @brief return rmsd cutoff
    //! @return rmsd cutoff
    double LCS::GetCutoff() const
    {
      return m_RmsdCutoff;
    }

    //! @brief get seed length
    //! @return seed length
    size_t LCS::GetSeedLength() const
    {
      return m_SeedLength;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief find a larger range by extending the given one that has a RMSD below the cutoff
    //! @param RANGE range of coordinates that are used to be as seed and to be extended
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return extended range
    math::Range< size_t> LCS::ExtendRange
    (
      const math::Range< size_t> &RANGE,
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      // store number of points and number of seed fragments
      const size_t number_points( COORDINATES.GetSize());

      // initialize rmsd object
      const RMSD rmsd( true);

      // calculate range length
      const size_t initial_range_length( RANGE.GetWidth() + 1);

      // get the initial coordinates and allocate memory for expansion
      util::SiPtrVector< const linal::Vector3D>
        coords( COORDINATES.SubSiPtrVector( RANGE.GetMin(), initial_range_length));
      coords.AllocateMemory( number_points - RANGE.GetMin());
      util::SiPtrVector< const linal::Vector3D>
        reference_coords( REFERENCE_COORDINATES.SubSiPtrVector( RANGE.GetMin(), initial_range_length));
      reference_coords.AllocateMemory( number_points - RANGE.GetMin());

      // make a copy of the range
      math::Range< size_t> longest_range( RANGE);
      math::Range< size_t> current_range( RANGE);

      // while the range can still be extended
      while( current_range.GetMax() < number_points - 1)
      {
        // increment the current range length
        current_range.SetMax( current_range.GetMax() + 1);

        // add the new coordinates
        coords.PushBack( COORDINATES( current_range.GetMax()));
        reference_coords.PushBack( REFERENCE_COORDINATES( current_range.GetMax()));

        // if the current RMSD is smaller
        if( rmsd.CalculateMeasure( coords, reference_coords) < m_RmsdCutoff)
        {
          // update the longest range
          longest_range = current_range;
        }
      }

      // return the longest range found so far
      return longest_range;
    }

    //! @brief returns the ranges of longest continuous segments that can be superimposed below cutoff for given coordinates
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return the range of longest continuous segment that can be superimposed below cutoff for given coordinates
    storage::List< math::Range< size_t> > LCS::CalculateRanges
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      // initialize variable to store the longest fragments
      storage::List< math::Range< size_t> > longest_fragments;

      // number points
      const size_t number_points( COORDINATES.GetSize());

      size_t max_frag_length( number_points), min_frag_length( 3);
      static size_t s_last_super_imp_size( size_t( math::Sqrt( max_frag_length * min_frag_length) + 1));
      size_t test_frag_length
      (
        s_last_super_imp_size + 1 < max_frag_length
        ? s_last_super_imp_size
        : size_t( math::Sqrt( max_frag_length * min_frag_length) + 1)
      );
      while( max_frag_length > min_frag_length)
      {
        bool found_good_range( false);
        for
        (
          size_t index( 0), mx_index( number_points - test_frag_length + 1);
          index < mx_index;
          ++index
        )
        {
          // range
          const math::Range< size_t> range( index, index + test_frag_length - 1);
          // fragment should be superimposable below cutoff
          if( IsGoodRange( range, COORDINATES, REFERENCE_COORDINATES))
          {
            found_good_range = true;
            break;
          }
        }
        if( found_good_range)
        {
          min_frag_length = test_frag_length;
        }
        else
        {
          max_frag_length = test_frag_length - 1;
        }
        test_frag_length = size_t( math::Sqrt( max_frag_length * min_frag_length) + 1);
      }

      s_last_super_imp_size = max_frag_length;
      // iterate over starting ranges
      for( size_t index( 0), mx_index( number_points - max_frag_length + 1); index < mx_index; ++index)
      {
        // range
        const math::Range< size_t> range( index, index + max_frag_length - 1);
        // fragment should be superimposable below cutoff
        if( IsGoodRange( range, COORDINATES, REFERENCE_COORDINATES))
        {
          longest_fragments.PushBack( range);
        }
      }

      // if message level is verbose or less
      if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
      {
        size_t ctr( 1);
        // iterate over each range
        for
        (
          storage::List< math::Range< size_t> >::const_iterator
            range_itr( longest_fragments.Begin()), range_itr_end( longest_fragments.End());
          range_itr != range_itr_end;
          ++range_itr, ++ctr
        )
        {
          // get the coordinates
          const util::SiPtrVector< const linal::Vector3D>
            coords( COORDINATES.SubSiPtrVector( range_itr->GetMin(), range_itr->GetWidth() + 1));
          const util::SiPtrVector< const linal::Vector3D>
            reference_coords( REFERENCE_COORDINATES.SubSiPtrVector( range_itr->GetMin(), range_itr->GetWidth() + 1));
          // make a copy of COORDINATES
          const math::TransformationMatrix3D transformation
          (
            RMSD::SuperimposeCoordinates( reference_coords, coords)
          );
          // make a copy of the actual coordinates so that they can be transformed and create SiPtrVector
          storage::Vector< linal::Vector3D>
            copy_coords( util::ConvertToStorageVector< linal::Vector3D, const linal::Vector3D>( COORDINATES));
          util::SiPtrVector< linal::Vector3D>
            copy_sp_coords( util::ConvertToSiPtrVector< linal::Vector3D>( copy_coords));
          // transform all the coordinates with the calculated translation
          coord::TransformCoordinates( copy_sp_coords, transformation);
          const double local_rmsd( RMSD( true).CalculateMeasure( coords, reference_coords));
          const double global_rmsd( RMSD( false).CalculateMeasure( copy_sp_coords, reference_coords));

          // print out the values
          BCL_MessageDbg
          (
            "the longest fragment # " + util::Format()( ctr) + " => " + range_itr->GetString() +
            " local_rmsd:\t" + util::Format()( local_rmsd) + " global_rmsd:\t" + util::Format()( global_rmsd)
          );
        }
      }

      // end
      return longest_fragments;
    }

    //! @brief returns the indices to the coordinates of longest continuous segments that can be superimposed below cutoff for given coordinates
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return the indices of longest continuous segments that can be superimposed below cutoff for given coordinates
    storage::List< storage::List< size_t> > LCS::CalculateIndices
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      return ConvertRangesToLists( CalculateRanges( COORDINATES, REFERENCE_COORDINATES));
    }

    //! @brief calculates LCS between given coordinates
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return LCS between COORDINATES and REFERENCE_COORDINATES
    double LCS::CalculateMeasure
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      // calculate the LCS and store the range for the first one
      const storage::List< math::Range< size_t> > longest_ranges( CalculateRanges( COORDINATES, REFERENCE_COORDINATES));
      if( longest_ranges.IsEmpty())
      {
        return 0.0;
      }

      // take first of the longest ranges, since all are equally long, they might just be at different ranges
      const math::Range< size_t> &lcs( longest_ranges.FirstElement());

      // return the length of the segment
      return lcs.GetWidth() + 1;
    }

    //! @brief calculates the transformation matrix that superimposes COORDINATES onto REFERENCE_COORDINATES
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return the transformation matrix that superimposes COORDINATES onto REFERENCE_COORDINATES
    math::TransformationMatrix3D LCS::CalculateSuperimposition
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      // calculate the LCS and store the range for the first one
      const storage::List< math::Range< size_t> > longest_ranges( CalculateRanges( COORDINATES, REFERENCE_COORDINATES));
      if( longest_ranges.IsEmpty())
      {
        return math::TransformationMatrix3D( util::UndefinedObject());
      }

      // take first of the longest ranges, since all are equally long, they might just be at different ranges
      const math::Range< size_t> &lcs( longest_ranges.FirstElement());

      // get the coordinates for the subset
      const util::SiPtrVector< const linal::Vector3D>
        coords( COORDINATES.SubSiPtrVector( lcs.GetMin(), lcs.GetWidth() + 1));
      const util::SiPtrVector< const linal::Vector3D>
        reference_coords( REFERENCE_COORDINATES.SubSiPtrVector( lcs.GetMin(), lcs.GetWidth() + 1));

      // calculate and return the transformation
      return RMSD::SuperimposeCoordinates( reference_coords, coords);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &LCS::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_RmsdCutoff, ISTREAM);
      io::Serialize::Read( m_SeedLength, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &LCS::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_RmsdCutoff, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_SeedLength, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief check if a given range superimposes below the cutoff
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return true, if coordinates within the given range are superimposable below the cutoff
    bool LCS::IsGoodRange
    (
      const math::Range< size_t> &RANGE,
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      // calculate range length
      const size_t initial_range_length( RANGE.GetWidth() + 1);

      // coordinates for that range
      const util::SiPtrVector< const linal::Vector3D>
        coords( COORDINATES.SubSiPtrVector( RANGE.GetMin(), initial_range_length));
      const util::SiPtrVector< const linal::Vector3D>
        reference_coords( REFERENCE_COORDINATES.SubSiPtrVector( RANGE.GetMin(), initial_range_length));

      // if the current RMSD is smaller
      return RMSD( true).CalculateMeasure( coords, reference_coords) < m_RmsdCutoff;
    }

    //! @brief converts a range to list of indices
    //! @param RANGE list of indices
    //! @return list that contains the indices in the range
    storage::List< size_t> LCS::ConvertRangeToIndices
    (
      const math::Range< size_t> &RANGE
    )
    {
      // initialize list to return
      storage::List< size_t> list;

      // iterate over the number
      for( size_t index( RANGE.GetMin()), index_end( RANGE.GetMax() + 1); index < index_end; ++index)
      {
        list.PushBack( index);
      }

      // end
      return list;
    }

    //! @brief converts a vector of ranges to list of list of indices
    //! @param RANGES vector of range of indices
    //! @return vector of lists that contains the indices in the range vector
    storage::List< storage::List< size_t> > LCS::ConvertRangesToLists
    (
      const storage::List< math::Range< size_t> > &RANGES
    )
    {
      // initialize list to return
      storage::List< storage::List< size_t> > list;

      // iterate over vector
      std::transform( RANGES.Begin(), RANGES.End(), std::inserter( list.InternalData(), list.Begin()), std::ptr_fun( &ConvertRangeToIndices));

      // end
      return list;
    }

  } // namespace quality
} // namespace bcl
