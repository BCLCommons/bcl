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
#include "contact/bcl_contact_types.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {

    //! @brief default constructor that constructs all Types
                     //                                                                                                                                                                                                                                                                                                                                          is_valid  dist   min  distance                          preferred distance                   tilt angle                                                                      min frag.                  sstypes
    Types::Types() : //                                                      WINDOW_RADII                                                                                                                               WINDOW_LENGTH_PAIR                                                                                                                       cutoff ssedist        range                             range                                range                                                                           intlength
      HELIX_HELIX(             AddEnum( "HELIX_HELIX",             TypeData( storage::Pair< size_t, size_t>( biol::GetSSTypes().HELIX->GetContactWindowRadius(),  biol::GetSSTypes().HELIX->GetContactWindowRadius()),  storage::Pair< size_t, size_t>( biol::GetSSTypes().HELIX->GetContactWindowLength(),  biol::GetSSTypes().HELIX->GetContactWindowLength()),  true,   8.0,   4.0, math::Range< double>( 5.0, 16.0), math::Range< double>( 8.00,  12.00), math::Range< double>( math::Angle::Radian( -60.0), math::Angle::Radian( 45.0)), 4.0, storage::Set< biol::SSType>( biol::GetSSTypes().HELIX)                           ))),
      HELIX_SHEET(             AddEnum( "HELIX_SHEET",             TypeData( storage::Pair< size_t, size_t>( biol::GetSSTypes().HELIX->GetContactWindowRadius(),  biol::GetSSTypes().STRAND->GetContactWindowRadius()), storage::Pair< size_t, size_t>( biol::GetSSTypes().HELIX->GetContactWindowLength(),  biol::GetSSTypes().STRAND->GetContactWindowLength()), true,   8.0,   4.0, math::Range< double>( 5.0, 16.0), math::Range< double>( 8.00,  13.00), math::Range< double>( math::Angle::Radian( -60.0), math::Angle::Radian(  0.0)), 4.0, storage::Set< biol::SSType>( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND)))),
      SHEET_HELIX(             AddEnum( "SHEET_HELIX",             TypeData( storage::Pair< size_t, size_t>( biol::GetSSTypes().STRAND->GetContactWindowRadius(), biol::GetSSTypes().HELIX->GetContactWindowRadius()),  storage::Pair< size_t, size_t>( biol::GetSSTypes().STRAND->GetContactWindowLength(), biol::GetSSTypes().HELIX->GetContactWindowLength()),  true,   8.0,   4.0, math::Range< double>( 5.0, 16.0), math::Range< double>( 8.00,  13.00), math::Range< double>( math::Angle::Radian( -60.0), math::Angle::Radian(  0.0)), 4.0, storage::Set< biol::SSType>( biol::GetSSTypes().STRAND, biol::GetSSTypes().HELIX)))),
      HELIX_STRAND(            AddEnum( "HELIX_STRAND",            TypeData( storage::Pair< size_t, size_t>( biol::GetSSTypes().HELIX->GetContactWindowRadius(),  biol::GetSSTypes().STRAND->GetContactWindowRadius()), storage::Pair< size_t, size_t>( biol::GetSSTypes().HELIX->GetContactWindowLength(),  biol::GetSSTypes().STRAND->GetContactWindowLength()), true,   8.0,   4.0, math::Range< double>( 7.0, 18.0), math::Range< double>( 9.00,  16.00), math::Range< double>( math::Angle::Radian( -60.0), math::Angle::Radian(  0.0)), 4.0, storage::Set< biol::SSType>( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND)))),
      STRAND_HELIX(            AddEnum( "STRAND_HELIX",            TypeData( storage::Pair< size_t, size_t>( biol::GetSSTypes().STRAND->GetContactWindowRadius(), biol::GetSSTypes().HELIX->GetContactWindowRadius()),  storage::Pair< size_t, size_t>( biol::GetSSTypes().STRAND->GetContactWindowLength(), biol::GetSSTypes().HELIX->GetContactWindowLength()),  true,   8.0,   4.0, math::Range< double>( 7.0, 18.0), math::Range< double>( 9.00,  16.00), math::Range< double>( math::Angle::Radian( -60.0), math::Angle::Radian(  0.0)), 4.0, storage::Set< biol::SSType>( biol::GetSSTypes().STRAND, biol::GetSSTypes().HELIX)))),
      STRAND_STRAND(           AddEnum( "STRAND_STRAND",           TypeData( storage::Pair< size_t, size_t>( biol::GetSSTypes().STRAND->GetContactWindowRadius(), biol::GetSSTypes().STRAND->GetContactWindowRadius()), storage::Pair< size_t, size_t>( biol::GetSSTypes().STRAND->GetContactWindowLength(), biol::GetSSTypes().STRAND->GetContactWindowLength()), true,   8.0,   3.0, math::Range< double>( 3.5,  5.5), math::Range< double>( 4.00,   5.25), math::Range< double>( math::Angle::Radian( -30.0), math::Angle::Radian(  0.0)), 4.0, storage::Set< biol::SSType>( biol::GetSSTypes().STRAND)                          ))),
      SHEET_SHEET(             AddEnum( "SHEET_SHEET",             TypeData( storage::Pair< size_t, size_t>( biol::GetSSTypes().STRAND->GetContactWindowRadius(), biol::GetSSTypes().STRAND->GetContactWindowRadius()), storage::Pair< size_t, size_t>( biol::GetSSTypes().STRAND->GetContactWindowLength(), biol::GetSSTypes().STRAND->GetContactWindowLength()), true,   8.0,   3.0, math::Range< double>( 6.0, 14.0), math::Range< double>( 8.00,  12.00), math::Range< double>( math::Angle::Radian( -60.0), math::Angle::Radian(  0.0)), 4.0, storage::Set< biol::SSType>( biol::GetSSTypes().STRAND)                          ))),
      UNDEFINED_HELIX_STRAND(  AddEnum( "UNDEFINED_HELIX_STRAND",  TypeData( GetUndefinedLengthPair(),                                                                                                                  GetUndefinedLengthPair(),                                                                                       false, GetUnknownResidueDistanceCutoff(), 4.0, GetUnknownDistanceRange()       , GetUnknownPreferredDistanceRange() , GetUnknownTiltAngleRange()                                                    , 4.0, storage::Set< biol::SSType>( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND)))),
      UNDEFINED_STRAND_HELIX(  AddEnum( "UNDEFINED_STRAND_HELIX",  TypeData( GetUndefinedLengthPair(),                                                                                                                  GetUndefinedLengthPair(),                                                                                       false, GetUnknownResidueDistanceCutoff(), 4.0, GetUnknownDistanceRange()       , GetUnknownPreferredDistanceRange() , GetUnknownTiltAngleRange()                                                    , 4.0, storage::Set< biol::SSType>( biol::GetSSTypes().STRAND, biol::GetSSTypes().HELIX)))),
      UNDEFINED_STRAND_STRAND( AddEnum( "UNDEFINED_STRAND_STRAND", TypeData( GetUndefinedLengthPair(),                                                                                                                  GetUndefinedLengthPair(),                                                                                       false, GetUnknownResidueDistanceCutoff(), 3.0, GetUnknownDistanceRange()       , GetUnknownPreferredDistanceRange() , GetUnknownTiltAngleRange()                                                    , 4.0, storage::Set< biol::SSType>( biol::GetSSTypes().STRAND)                          )))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Types::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief cut off distance for unknown types;
    //! @return cut off distance for unknown types;
    const double &Types::GetUnknownResidueDistanceCutoff()
    {
      static const double s_distance_cutoff( 8.0);

      // end
      return s_distance_cutoff;
    }

    //! @brief return distance range for unknown types
    //! @return distance range for unknown types
    const math::Range< double> &Types::GetUnknownDistanceRange()
    {
      static const math::Range< double> s_distance_range( 5.0, 14.0);

      // end
      return s_distance_range;
    }

    //! @brief return preferred distance range for unknown types
    //! @return preferred distance range for unknown types
    const math::Range< double> &Types::GetUnknownPreferredDistanceRange()
    {
      static const math::Range< double> s_distance_range( 5.0, 14.0);

      // end
      return s_distance_range;
    }

    //! @brief return tilt angle range for unknown types
    //! @return tilt angle range for unknown types
    const math::Range< double> &Types::GetUnknownTiltAngleRange()
    {
      static const math::Range< double> s_tilt_angle_range( math::Angle::Radian( -30.0), math::Angle::Radian( 30.0));

      // end
      return s_tilt_angle_range;
    }

    //! @brief return undefined pair length
    //! @return undefined pair length
    const storage::Pair< size_t, size_t> &Types::GetUndefinedLengthPair()
    {
      static const storage::Pair< size_t, size_t> s_undefined_pair
      (
        util::GetUndefined< size_t>(), util::GetUndefined< size_t>()
      );

      // end
      return s_undefined_pair;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief returns the reverse pair for the contact type
    //! @param TYPE Type of interest
    //! @return the reverse pair for the contact type
    const Type &Types::Reverse( const Type &TYPE) const
    {
      if( TYPE == HELIX_SHEET)
      {
        return SHEET_HELIX;
      }
      else if( TYPE == SHEET_HELIX)
      {
        return HELIX_SHEET;
      }
      if( TYPE == HELIX_STRAND)
      {
        return STRAND_HELIX;
      }
      else if( TYPE == STRAND_HELIX)
      {
        return HELIX_STRAND;
      }
      else if( TYPE == UNDEFINED_HELIX_STRAND)
      {
        return UNDEFINED_STRAND_HELIX;
      }
      else if( TYPE == UNDEFINED_STRAND_HELIX)
      {
        return UNDEFINED_HELIX_STRAND;
      }

      return TYPE;
    }

    //! @brief Given two SSEGeometryInterfaces, merges their SSTypes to form a ContactType
    //! @param SSE_A first SSEGeometryInterface
    //! @param SSE_B second SSEGeometryInterface
    //! @return Contact type formed by SSE_A and SSE_B
    const Type &Types::TypeFromSSTypes
    (
      const assemble::SSEGeometryInterface &SSE_A,
      const assemble::SSEGeometryInterface &SSE_B
    ) const
    {
      // HELIX_HELIX
      if( SSE_A.GetType() == biol::GetSSTypes().HELIX && SSE_B.GetType() == biol::GetSSTypes().HELIX)
      {
        return HELIX_HELIX;
      }
      // HELIX_SHEET
      if( SSE_A.GetType() == biol::GetSSTypes().HELIX && SSE_B.GetType() == biol::GetSSTypes().STRAND)
      {
        return HELIX_SHEET;
      }
      // HELIX_SHEET
      if( SSE_A.GetType() == biol::GetSSTypes().STRAND && SSE_B.GetType() == biol::GetSSTypes().HELIX)
      {
        return SHEET_HELIX;
      }
      // STRAND_STRAND ( or SHEET_SHEET)
      if( SSE_A.GetType() == biol::GetSSTypes().STRAND && SSE_B.GetType() == biol::GetSSTypes().STRAND)
      {
        return STRAND_STRAND;
      }

      return e_Undefined;
    }

    //! @brief returns a map with distance ranges for valid contact types
    //! @return a map with distance ranges for valid contact types
    const storage::Map< Type, math::Range< double> > &Types::GetValidDistanceRanges() const
    {
      // initialize static storage for the distance cutoffs
      static const storage::Map< Type, math::Range< double> > s_distance_ranges( CollectValidDistanceRanges());

      // end
      return s_distance_ranges;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief generates a map with distance ranges for valid contact types
    //! @return a map with distance ranges for valid contact types
    storage::Map< Type, math::Range< double> > Types::CollectValidDistanceRanges() const
    {
      // initialize storage
      storage::Map< Type, math::Range< double> > distance_range_map;

      // iterator over standard amino acid types
      for
      (
        const_iterator type_itr( Begin()), type_itr_end( GetEnumIteratorFromIndex( s_NumberValidTypes));
        type_itr != type_itr_end;
        ++type_itr
      )
      {
        // store the cutoff distance thresholds
        distance_range_map[ *type_itr] = ( *type_itr)->GetDistanceRange();
      }
      // end
      return distance_range_map;
    }

    //! @brief retrieves Types enumerator
    //! @return Types enumerator
    const Types &GetTypes()
    {
      return Types::GetEnums();
    }

  } // namespace contact

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< contact::TypeData, contact::Types>;

  } // namespace util
} // namespace bcl
