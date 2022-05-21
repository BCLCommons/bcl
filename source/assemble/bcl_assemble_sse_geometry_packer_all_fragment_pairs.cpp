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
#include "assemble/bcl_assemble_sse_geometry_packer_all_fragment_pairs.h"

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> SSEGeometryPackerAllFragmentPairs::s_Instance
    (
      GetObjectInstances().AddInstance( new SSEGeometryPackerAllFragmentPairs())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SSEGeometryPackerAllFragmentPairs::SSEGeometryPackerAllFragmentPairs() :
      m_MinimalInterfaceLength( SSEGeometryPacking::GetDefaultFragmentMinimalInterfaceLength()),
      m_UseDistanceCutoffs( false),
      m_DistanceCutoffs()
    {
    }

    //! @brief constructor from a minimal interface length
    //! @param MINIMAL_INTERFACE_LENGTH minimal interface length
    SSEGeometryPackerAllFragmentPairs::SSEGeometryPackerAllFragmentPairs
    (
      const double MINIMAL_INTERFACE_LENGTH
    ) :
      m_MinimalInterfaceLength( MINIMAL_INTERFACE_LENGTH),
      m_UseDistanceCutoffs( false),
      m_DistanceCutoffs()
    {
    }

    //! @brief constructor from a minimal interface length
    //! @param MINIMAL_INTERFACE_LENGTH minimal interface length
    //! @param DISTANCE_CUTOFFS map of distance cutoffs for each contact type
    SSEGeometryPackerAllFragmentPairs::SSEGeometryPackerAllFragmentPairs
    (
      const double MINIMAL_INTERFACE_LENGTH,
      const storage::Map< biol::SSType, storage::Map< biol::SSType, double> > &DISTANCE_CUTOFFS
    ) :
      m_MinimalInterfaceLength( MINIMAL_INTERFACE_LENGTH),
      m_UseDistanceCutoffs( true),
      m_DistanceCutoffs( DISTANCE_CUTOFFS)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SSEGeometryPackerAllFragmentPairs
    SSEGeometryPackerAllFragmentPairs *SSEGeometryPackerAllFragmentPairs::Clone() const
    {
      return new SSEGeometryPackerAllFragmentPairs( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSEGeometryPackerAllFragmentPairs::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate the packing between all fragments of the given geometries and return it in a list of lists
    //! @param SSE_GEOMETRY_A first SSEGeometryInterface derived class of interest
    //! @param SSE_GEOMETRY_B second SSEGeometryInterface derived class of interest
    //! @return list of lists containing packing between all fragment pairings between given geometry pair
    storage::Vector< storage::List< SSEGeometryPacking> > SSEGeometryPackerAllFragmentPairs::operator()
    (
      const SSEGeometryInterface &SSE_GEOMETRY_A,
      const SSEGeometryInterface &SSE_GEOMETRY_B
    ) const
    {
      // initialize a static undefined sse geometry packing
      static const SSEGeometryPacking s_undefined_packing;

      // if both SSEs are not helix or strand or they are the same geometry
      if
      (
        !SSE_GEOMETRY_A.GetType()->IsStructured() ||
        !SSE_GEOMETRY_B.GetType()->IsStructured() ||
        &SSE_GEOMETRY_A == &SSE_GEOMETRY_B
      )
      {
        // return empty list
        return storage::Vector< storage::List< SSEGeometryPacking> >();
      }

      // get the sub-geometries for both SSEGeometries and store them
      const util::SiPtrVector< const SSEGeometryInterface> geometries_a( SSE_GEOMETRY_A.GetSSEGeometries());
      const util::SiPtrVector< const SSEGeometryInterface> geometries_b( SSE_GEOMETRY_B.GetSSEGeometries());

      // declare the list of lists SSEPackings to be returned and initialize them
      storage::Vector< storage::List< SSEGeometryPacking> > packing_list( geometries_a.GetSize());

      // initialize iterators for the packing list
      storage::Vector< storage::List< SSEGeometryPacking> >::iterator pack_itr( packing_list.Begin());
      const storage::Vector< storage::List< SSEGeometryPacking> >::const_iterator pack_itr_end( packing_list.End());

      // determine the max cutoff for this type
      const double cutoff_distance
      (
        m_UseDistanceCutoffs ?
          m_DistanceCutoffs.Find( SSE_GEOMETRY_A.GetType())->second.Find( SSE_GEOMETRY_B.GetType())->second :
          0.0
      );

      // iterate over the first fragments list while also iterating over the packing list
      for
      (
        util::SiPtrVector< const SSEGeometryInterface>::const_iterator
          fragment_itr_a( geometries_a.Begin()), fragment_itr_a_end( geometries_a.End());
        fragment_itr_a != fragment_itr_a_end && pack_itr != pack_itr_end;
        ++fragment_itr_a, ++pack_itr
      )
      {
        // iterate over second list
        for
        (
          util::SiPtrVector< const SSEGeometryInterface>::const_iterator
            fragment_itr_b( geometries_b.Begin()), fragment_itr_b_end( geometries_b.End());
          fragment_itr_b != fragment_itr_b_end;
          ++fragment_itr_b
        )
        {
          // list the fragments
          BCL_MessageDbg
          (
            "\t\tLooking at " + ( *fragment_itr_a)->GetIdentification() + " vs " +
            ( *fragment_itr_b)->GetIdentification()
          )

          // if using distance cutoffs
          if( m_UseDistanceCutoffs)
          {
            // calculate the distance between centers of geometries
            const double center_distance
            (
              linal::Distance( ( *fragment_itr_a)->GetCenter(), ( *fragment_itr_b)->GetCenter())
            );

            // if the distance is greater than the cutoff for that contact type
            if( center_distance > cutoff_distance)
            {
              // insert undefined SSEPacking
              pack_itr->PushBack( s_undefined_packing);

              // move to next fragment pair
              continue;
            }
          }

          SSEGeometryPacking pack( **fragment_itr_a, **fragment_itr_b, m_MinimalInterfaceLength);
//          if( pack.GetDistance() > pack.GetContactType()->GetDistanceRange().GetMax())
//          {
//            BCL_MessageStd
//            (
//              "Pack Dist: " + util::Format()( pack.GetDistance())
//              + " contact type: " + pack.GetContactType().GetName()
//              + " max dist: " + util::Format()( pack.GetContactType()->GetDistanceRange().GetMax())
//            );
//          }
          // insert a packing for these two fragments into the list
          pack_itr->PushBack( pack);
//          BCL_MessageDbg( "\t\t\twith packing:" + pack_itr->LastElement().GetIdentification());

          // if the packing is defined and the distance is less than the critical distance give a warning
          if
          (
            pack_itr->LastElement().GetContactType().IsDefined() &&
            pack_itr->LastElement().GetDistance() <= SSEGeometryPacking::s_CriticalDistance
          )
          {
            BCL_MessageVrb
            (
              "distance_check_frag " + pack_itr->LastElement().GetIdentification()
            );
          }
        }
      }

      // return list of packings
      return packing_list;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEGeometryPackerAllFragmentPairs::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_MinimalInterfaceLength, ISTREAM);
      io::Serialize::Read( m_UseDistanceCutoffs, ISTREAM);
      io::Serialize::Read( m_DistanceCutoffs, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SSEGeometryPackerAllFragmentPairs::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_MinimalInterfaceLength, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_UseDistanceCutoffs, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DistanceCutoffs, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return map of maps for upper distance ranges for sstype pairs
    //! @return map of maps for upper distance ranges for sstype pairs
    const storage::Map< biol::SSType, storage::Map< biol::SSType, double> > &
    SSEGeometryPackerAllFragmentPairs::GetDistanceCutoffMap()
    {
      // initialize map
      static storage::Map< biol::SSType, storage::Map< biol::SSType, double> > s_cutoff_map;

      // if map is not initialized yet then initialize it
      if( s_cutoff_map.IsEmpty())
      {
        // initialize static overheads for helix and strand fragments
        const double helix_overhead
        (
          biol::GetSSTypes().HELIX->GetFragmentLength() * biol::GetSSTypes().HELIX->GetRiseInZPerResidue() / 2.0
        );
        const double strand_overhead
        (
          biol::GetSSTypes().STRAND->GetFragmentLength() * biol::GetSSTypes().STRAND->GetRiseInZPerResidue() / 2.0
        );

        // store the minimal fragment interface length
        const double interface_length( SSEGeometryPacking::GetDefaultFragmentMinimalInterfaceLength());

        // insert values
        s_cutoff_map[ biol::GetSSTypes().HELIX][ biol::GetSSTypes().HELIX] =
          contact::GetTypes().HELIX_HELIX->GetDistanceRange().GetMax() + 2 * helix_overhead - interface_length;
        s_cutoff_map[ biol::GetSSTypes().HELIX][ biol::GetSSTypes().STRAND] =
          contact::GetTypes().HELIX_STRAND->GetDistanceRange().GetMax() + helix_overhead + strand_overhead - interface_length;
        s_cutoff_map[ biol::GetSSTypes().STRAND][ biol::GetSSTypes().HELIX] =
          contact::GetTypes().HELIX_STRAND->GetDistanceRange().GetMax() + helix_overhead + strand_overhead - interface_length;
        s_cutoff_map[ biol::GetSSTypes().STRAND][ biol::GetSSTypes().STRAND] =
          contact::GetTypes().SHEET_SHEET->GetDistanceRange().GetMax() + 2 * strand_overhead - interface_length;
      }

      // end
      return s_cutoff_map;
    }

    //! @brief return map of maps for upper clash distance ranges for sstype pairs
    //! @return map of maps for upper clash distance ranges for sstype pairs
    const storage::Map< biol::SSType, storage::Map< biol::SSType, double> > &
    SSEGeometryPackerAllFragmentPairs::GetClashDistanceCutoffMap()
    {
      // initialize map
      static storage::Map< biol::SSType, storage::Map< biol::SSType, double> > s_cutoff_map;

      // if map is not initialized yet then initialize it
      if( s_cutoff_map.IsEmpty())
      {
        // initialize overheads for helix and strand fragments
        const double helix_overhead
        (
          biol::GetSSTypes().HELIX->GetFragmentLength() * biol::GetSSTypes().HELIX->GetRiseInZPerResidue() / 2.0
        );
        const double strand_overhead
        (
          biol::GetSSTypes().STRAND->GetFragmentLength() * biol::GetSSTypes().STRAND->GetRiseInZPerResidue() / 2.0
        );

        // insert values
        s_cutoff_map[ biol::GetSSTypes().HELIX][ biol::GetSSTypes().HELIX] =
          contact::GetTypes().HELIX_HELIX->GetMinimalSSEDistance() + 2 * helix_overhead;
        s_cutoff_map[ biol::GetSSTypes().HELIX][ biol::GetSSTypes().STRAND] =
          contact::GetTypes().HELIX_SHEET->GetMinimalSSEDistance() + helix_overhead + strand_overhead;
        s_cutoff_map[ biol::GetSSTypes().STRAND][ biol::GetSSTypes().HELIX] =
          contact::GetTypes().HELIX_SHEET->GetMinimalSSEDistance() + helix_overhead + strand_overhead;
        s_cutoff_map[ biol::GetSSTypes().STRAND][ biol::GetSSTypes().STRAND] =
          contact::GetTypes().SHEET_SHEET->GetMinimalSSEDistance() + 2 * strand_overhead;
      }

      // end
      return s_cutoff_map;
    }

  } // namespace assemble
} // namespace bcl
