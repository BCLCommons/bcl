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
#include "assemble/bcl_assemble_domain_interface.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_compare.h"
#include "assemble/bcl_assemble_sse_geometry.h"
#include "biol/bcl_biol_atom.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  /////////////////
  // data access //
  /////////////////

    //! @brief returns all SSEs in domain of given SSTYPE in a util::SiPtrVector
    //! @param SS_TYPE specific SSTYPE
    //! @return all SSEs in domain of given SSTYPE in a util::SiPtrVector
    util::SiPtrVector< const SSE> DomainInterface::GetSSEs( const biol::SSType &SS_TYPE) const
    {
      // initialize vector of sses
      util::SiPtrVector< const SSE> sses;

      // get and store all the SSEs
      const util::SiPtrVector< const SSE> all_sses_vector( GetSSEs());

      //loop over all SSElements
      for
      (
        util::SiPtrVector< const SSE>::const_iterator sse_itr( all_sses_vector.Begin()),
          sse_itr_end( all_sses_vector.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // if SSE behind sse_itr is of given SSTYPE
        if( ( *sse_itr)->GetType() == SS_TYPE)
        {
          // insert it into the list of SSEs to be returned
          sses.PushBack( util::SiPtr< const SSE>( *sse_itr));
        }
      }

      // end
      return sses;
    }

    //! @brief returns all SSEs in domain of given SSTYPEs in a util::SiPtrVector
    //! @param SS_TYPES set of sstypes
    //! @return all SSEs in domain of given SSTYPEs in a util::SiPtrVector
    util::SiPtrVector< const SSE> DomainInterface::GetSSEs( const storage::Set< biol::SSType> &SS_TYPES) const
    {
      // initialize vector of sses
      util::SiPtrVector< const SSE> sses;

      // get and store all the SSEs
      const util::SiPtrVector< const SSE> all_sses_vector( GetSSEs());

      //loop over all SSElements
      for
      (
        util::SiPtrVector< const SSE>::const_iterator sse_itr( all_sses_vector.Begin()),
          sse_itr_end( all_sses_vector.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // if SSE behind sse_itr is of given SSTYPE
        if( SS_TYPES.Contains( ( *sse_itr)->GetType()))
        {
          // insert it into the list of SSEs to be returned
          sses.PushBack( util::ToSiPtr( **sse_itr));
        }
      }

      // end
      return sses;
    }

    //! @brief return total number of sses
    //! @return total number of sses
    size_t DomainInterface::GetNumberSSEs() const
    {
      return GetSSEs().GetSize();
    }

    //! @brief return number of SSE of specified SSTYPE
    //! @param SS_TYPE specific SSTYPE
    //! @return number of SSE of specified SSTYPE
    size_t DomainInterface::GetNumberSSEs( const biol::SSType &SS_TYPE) const
    {
      // initialize count
      size_t sum( 0);

      // get and store the SSEs
      const util::SiPtrVector< const SSE> sse_vector( GetSSEs());

      //loop over all SSElements
      for
      (
        util::SiPtrVector< const SSE>::const_iterator sse_itr( sse_vector.Begin()),
          sse_itr_end( sse_vector.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        sum += ( ( *sse_itr)->GetType() == SS_TYPE);
      }

      // end
      return sum;
    }

    //! @brief returns all atoms in domain as SiPtrVector
    //! @return all atoms in domain as SiPtrVector
    util::SiPtrVector< const biol::Atom> DomainInterface::GetAtoms() const
    {
      // initialize storage for atoms
      util::SiPtrVector< const biol::Atom> atoms;

      // get all SSEs first
      const util::SiPtrVector< const SSE> sse_vector( GetSSEs());

      //loop over all SSElements
      for
      (
        util::SiPtrVector< const SSE>::const_iterator sse_itr( sse_vector.Begin()),
          sse_itr_end( sse_vector.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        atoms.Append( ( *sse_itr)->GetAtoms());
      }

      return atoms;
    }

    //! @brief returns all atoms of specified ATOM_TYPES in domain as SiPtrVector
    //! @param ATOM_TYPES Set of AtomTypes of interest
    //! @return all atoms of specified types in domain as SiPtrVector
    util::SiPtrVector< const biol::Atom>
    DomainInterface::GetAtoms
    (
      const storage::Set< biol::AtomType> &ATOM_TYPES
    ) const
    {
      //util::SiPtrVector of atoms
      util::SiPtrVector< const biol::Atom> atoms;

      // get all SSEs first
      const util::SiPtrVector< const SSE> sse_vector( GetSSEs());

      //loop over all SSElements
      for
      (
        util::SiPtrVector< const SSE>::const_iterator sse_itr( sse_vector.Begin()),
          sse_itr_end( sse_vector.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // append the atoms of current SSE's atoms
        atoms.Append( ( *sse_itr)->GetAtoms( ATOM_TYPES));
      }

      //return
      return atoms;
    }

    //! @brief returns coordinates for all atoms in domain as SiPtrVector
    //! @return coordinates for all atoms in domain as SiPtrVector
    util::SiPtrVector< const linal::Vector3D> DomainInterface::GetAtomCoordinates() const
    {
      // util::SiPtrVector of positions
      util::SiPtrVector< const linal::Vector3D> positions;

      // get all SSEs first
      const util::SiPtrVector< const SSE> sse_vector( GetSSEs());

      //loop over all SSElements
      for
      (
        util::SiPtrVector< const SSE>::const_iterator sse_itr( sse_vector.Begin()),
          sse_itr_end( sse_vector.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        positions.Append( ( *sse_itr)->GetAtomCoordinates());
      }

      //return
      return positions;
    }

    //! @brief returns coordinates for all atoms of specified ATOM_TYPES in domain as SiPtrVector
    //! @param ATOM_TYPES Set of AtomTypes of interest
    //! @return coordinates for all atoms of specified types in domain as SiPtrVector
    util::SiPtrVector< const linal::Vector3D>
    DomainInterface::GetAtomCoordinates
    (
      const storage::Set< biol::AtomType> &ATOM_TYPES
    ) const
    {
      //util::SiPtrVector of positions
      util::SiPtrVector< const linal::Vector3D> positions;

      // get all SSEs first
      const util::SiPtrVector< const SSE> sse_vector( GetSSEs());

      //loop over all SSElements
      for
      (
        util::SiPtrVector< const SSE>::const_iterator sse_itr( sse_vector.Begin()),
          sse_itr_end( sse_vector.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // append the positions of current SSE's atoms
        positions.Append( ( *sse_itr)->GetAtomCoordinates( ATOM_TYPES));
      }

      //return
      return positions;
    }

    //! @brief returns the center of the domain
    //! @return the center of the domain
    linal::Vector3D DomainInterface::GetCenter() const
    {
      return coord::CenterOfMass( GetAtomCoordinates());
    }

    //! @brief return all the amino acids in all SSEs
    //! @return all the amino acids in all SSEs
    util::SiPtrVector< const biol::AABase> DomainInterface::GetAminoAcids() const
    {
      util::SiPtrVector< const biol::AABase> aminoacids;

      // get all SSEs first
      const util::SiPtrVector< const SSE> sse_vector( GetSSEs());

      // iterate over all SSEs in this domain
      for
      (
        util::SiPtrVector< const SSE>::const_iterator sse_itr( sse_vector.Begin()),
          sse_itr_end( sse_vector.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // add the amino acids for this SSE
        aminoacids.Append( util::SiPtrVector< const biol::AABase>( ( *sse_itr)->GetData()));
      }

      return aminoacids;
    }

    //! @brief create and return SSEGeometries for all SSEs in this Domain
    //! @return ShPtrVector of SSEGeometries corresponding to SSEs in this domain
    util::ShPtrVector< SSEGeometry> DomainInterface::GetSSEGeometries() const
    {
      // initialize storage to return geometries
      util::ShPtrVector< SSEGeometry> geometries;

      // get all SSEs first
      const util::SiPtrVector< const SSE> sse_vector( GetSSEs());

      // iterate over all SSEs
      for
      (
        util::SiPtrVector< const SSE>::const_iterator sse_itr( sse_vector.Begin()),
          sse_itr_end( sse_vector.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // construct a geometry from this SSE and add it to the vector
        geometries.PushBack( util::ShPtr< SSEGeometry>( new SSEGeometry( **sse_itr)));
      }

      // end
      return geometries;
    }

    //! @brief create and return SSEGeometries for all SSEs in this Domain
    //! @param SS_TYPE specific SSTYPE
    //! @return ShPtrVector of SSEGeometries corresponding to SSEs in this domain
    util::ShPtrVector< SSEGeometry> DomainInterface::GetSSEGeometries( const biol::SSType &SS_TYPE) const
    {
      // initialize storage to return geometries
      util::ShPtrVector< SSEGeometry> geometries;

      // get all SSEs first
      const util::SiPtrVector< const SSE> sse_vector( GetSSEs());

      for
      (
        util::SiPtrVector< const SSE>::const_iterator sse_itr( sse_vector.Begin()),
          sse_itr_end( sse_vector.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // if this SSE is of the requested SS_TYPE
        if( ( *sse_itr)->GetType() == SS_TYPE)
        {
          // construct and append a new SSEGeometry from this SSE
          geometries.PushBack( util::ShPtr< SSEGeometry>( new SSEGeometry( **sse_itr)));
        }
      }

      // end
      return geometries;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief checks if domain already contains THIS_SSE
    //! @param THIS_SSE SSE of interest
    //! @return whether domain already contains THIS_SSE
    bool DomainInterface::DoesContain( const SSE &THIS_SSE) const
    {
      // get all the SSEs
      const util::SiPtrVector< const SSE> sse_vector( GetSSEs());

      // search for this sse in m_Data
      util::SiPtrVector< const SSE>::const_iterator itr
      (
        std::find_if
        (
          sse_vector.Begin(),
          sse_vector.End(),
          SSECompare( THIS_SSE)
        )
      );

      // return if it's found
      return ( itr != sse_vector.End());
    }

    //! @brief returns true if the domain has no SSEs in it
    //! @return true if the domain has no SSEs in it
    bool DomainInterface::IsEmpty() const
    {
      return GetSSEs().IsEmpty();
    }

    //! @brief checks if domain already contains this THIS_SSE or any overlapping SSE with THIS_SSE
    //! @param THIS_SSE SSE to be searched for
    //! @return if domain already contains this THIS_SSE or any overlapping SSE with THIS_SSE
    bool DomainInterface::DoesContainOverlapping( const SSE &THIS_SSE) const
    {
      // get all the SSEs
      const util::SiPtrVector< const SSE> sse_vector( GetSSEs());

      // search for this sse or any overlapping sse in m_Data
      util::SiPtrVector< const SSE>::const_iterator itr
      (
        std::find_if
        (
          sse_vector.Begin(),
          sse_vector.End(),
          SSECompareOverlap( THIS_SSE)
        )
      );

      // return if it's found
      return ( itr != sse_vector.End());
    }

    //! @brief returns SiPtrList of sses that have short loops ( at most MAX_LOOP_LENGTH) to provided TARGET_SSE
    //! @param TARGET_SSE SSE for which short loop connecting SSEs are being search
    //! @param MAX_LOOP_LENGTH maximum loop length that is defined as short loop
    //! @return SiPtrList of sses that have short loops to provided TARGET_SSE
    util::SiPtrList< const SSE> DomainInterface::GetSSEsWithShortLoops
    (
      const SSE &TARGET_SSE,
      const size_t MAX_LOOP_LENGTH
    ) const
    {
      // initialize the sse_list to be returned
      util::SiPtrList< const SSE> sse_list;

      // get the vector of all SSEs
      const util::SiPtrVector< const SSE> sse_vector( GetSSEs());

      // iterate over sses stored in the domain
      for
      (
        util::SiPtrVector< const SSE>::const_iterator sse_itr( sse_vector.Begin()),
          sse_itr_end( sse_vector.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // calculate the loop length of this sse behind sse_itr to TARGET_SSE
        const size_t this_distance( biol::CalculateSequenceDistance( **sse_itr, TARGET_SSE));

        // if the loop length is not undefined or is not larger than the MAX_LOOP_LENGTH
        if( this_distance != util::GetUndefinedSize_t() && this_distance <= MAX_LOOP_LENGTH)
        {
          // insert this sse into the list
          sse_list.PushBack( *sse_itr);
        }
      }

      // end
      return sse_list;
    }

    //! @brief selects from provided SSE_LIST, sses that have short loops ( <=MAX_LOOP_LENGTH) to SSEs in Domain
    //! @param SSE_LIST list of SSEs on which the selection will be done
    //! @param MAX_LOOP_LENGTH maximum loop length that is defined as short loop
    //! @return Subset of SSE_LIST that has short loops to SSEs in this domain
    util::SiPtrList< const SSE> DomainInterface::GetSSEsWithShortLoops
    (
      const util::SiPtrList< const SSE> &SSE_LIST,
      const size_t MAX_LOOP_LENGTH
    ) const
    {
      // initialize the list to be returned
      util::SiPtrList< const SSE> eligible_sses;

      // iterate over every SSE in this list
      for
      (
        util::SiPtrList< const SSE>::const_iterator sse_itr( SSE_LIST.Begin()), sse_itr_end( SSE_LIST.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // if this domain has at least one sse with a short loop to the SSE behind sse_itr
        if( !GetSSEsWithShortLoops( **sse_itr, MAX_LOOP_LENGTH).IsEmpty())
        {
          // insert this sse into eligible_sses
          eligible_sses.PushBack( *sse_itr);
        }
      }

      // return
      return eligible_sses;
    }

    //! @brief returns pairs of sses that have short loops (at most MAX_LOOP_LENGTH) between each other
    //! @param MAX_LOOP_LENGTH maximum loop length that is defined as short loop
    //! @return list of pairs of sses that have short loops (at most MAX_LOOP_LENGTH) between
    //! @return each other
    storage::List< storage::VectorND< 2, util::SiPtr< const SSE> > > DomainInterface::GetSSEsWithShortLoops
    (
      const size_t MAX_LOOP_LENGTH
    ) const
    {
      // initialize the sse_list to be returned
      storage::List< storage::VectorND< 2, util::SiPtr< const SSE> > > sse_list;

      // get the vector of all SSEs
      const util::SiPtrVector< const SSE> sse_vector( GetSSEs());

      // iterate over sses stored in the domain
      for
      (
        util::SiPtrVector< const SSE>::const_iterator sse_itr_a( sse_vector.Begin()),
          sse_itr_end( sse_vector.End());
        sse_itr_a != sse_itr_end;
        ++sse_itr_a
      )
      {
        // iterate over other sses in the sequence
        for
        (
          util::SiPtrVector< const SSE>::const_iterator sse_itr_b( sse_vector.Begin());
          sse_itr_b != sse_itr_end;
          ++sse_itr_b
        )
        {
          // if same SSE skip
          if( sse_itr_a == sse_itr_b)
          {
            continue;
          }

          // calculate the loop length of this sse behind sse_itr to TARGET_SSE
          const size_t this_distance( biol::CalculateSequenceDistance( **sse_itr_a, **sse_itr_b));

          // if the loop length is not undefined or is not larger than the MAX_LOOP_LENGTH
          if( this_distance != util::GetUndefinedSize_t() && this_distance <= MAX_LOOP_LENGTH)
          {
            // insert this sse pair into the list
            sse_list.PushBack( storage::VectorND< 2, util::SiPtr< const SSE> >( *sse_itr_a, *sse_itr_b));
          }
        }
      }

      // end
      return sse_list;
    }

    //! @brief returns the SSE before and SSE after given TARGET_SSE in this Domain
    //! @param TARGET_SSE SSE of interest
    //! @return the SSE before and SSE after given TARGET_SSE in this Domain
    storage::VectorND< 2, util::SiPtr< const SSE> > DomainInterface::GetNeighborSSEs
    (
      const SSE &TARGET_SSE
    ) const
    {
      // initialize the list to be returned
      storage::VectorND< 2, util::SiPtr< const SSE> > neighbor_sses;

      // get the vector of all SSEs
      const util::SiPtrVector< const SSE> sse_vector( GetSSEs());

      storage::Pair< util::SiPtr< const SSE>, int> left_dist( util::SiPtr< const SSE>(), 99999);
      storage::Pair< util::SiPtr< const SSE>, int> right_dist( util::SiPtr< const SSE>(), 99999);

      // iterate over SSEs
      for( util::SiPtrVector< const SSE>::const_iterator sse_itr( sse_vector.Begin()), sse_itr_end( sse_vector.End()); sse_itr != sse_itr_end; ++sse_itr)
      {
        // distances
        const int distance_left( TARGET_SSE.GetFirstAA()->GetSeqID() - ( *sse_itr)->GetLastAA()->GetSeqID());

        // if TARGET_SSE comes before SSE pointed by and is closer
        if( distance_left > 0 && distance_left < left_dist.Second())
        {
          left_dist.First()  = *sse_itr;
          left_dist.Second() = distance_left;
          continue;
        }

        const int distance_right( ( *sse_itr)->GetFirstAA()->GetSeqID() - TARGET_SSE.GetLastAA()->GetSeqID());

        // if TARGET_SSE comes after SSE pointed by and is closer
        if( distance_right > 0 && distance_right < right_dist.Second())
        {
          right_dist.First()  = *sse_itr;
          right_dist.Second() = distance_right;
        }
      }

      neighbor_sses.First() = left_dist.First();
      neighbor_sses.Second() = right_dist.First();

      // return
      return neighbor_sses;
    }

    //! @brief get the directly adjacent sses by seqid
    //! @param TARGET_SSE SSE of interest
    //! @return the SSE before and SSE after given TARGET_SSE in this Domain by seqid
    storage::VectorND< 2, util::SiPtr< const SSE> > DomainInterface::GetAdjacentSSEs
    (
      const SSE &TARGET_SSE
    ) const
    {
      // initialize the list to be returned
      storage::VectorND< 2, util::SiPtr< const SSE> > neighbor_sses( GetNeighborSSEs( TARGET_SSE));

      // check if neighbor sses distance is 1
      if
      (
           neighbor_sses.First().IsDefined()
        && ( TARGET_SSE.GetFirstAA()->GetSeqID() - neighbor_sses.First()->GetLastAA()->GetSeqID()) != 1
      )
      {
        neighbor_sses.First() = util::SiPtr< const SSE>();
      }

      if
      (
           neighbor_sses.Second().IsDefined()
        && ( neighbor_sses.Second()->GetFirstAA()->GetSeqID() - TARGET_SSE.GetLastAA()->GetSeqID()) != 1
      )
      {
        neighbor_sses.Second() = util::SiPtr< const SSE>();
      }

      // end
      return neighbor_sses;
    }

    //! @brief returns SiPtrList of SSEs that overlap with the given TARGET_SSE
    //! @param TARGET_SSE SSE for which overlapping SSEs are searched for
    //! @return SiPtrList of SSEs that overlap with the given TARGET_SSE
    util::SiPtrList< const SSE> DomainInterface::GetOverlappingSSEs
    (
      const SSE &TARGET_SSE
    ) const
    {
      // get the vector of all SSEs
      const util::SiPtrVector< const SSE> sse_vector( GetSSEs());

      // search for overlapping sses with SP_SSE using equal_range
      std::pair
      <
        util::SiPtrVector< const SSE>::const_iterator,
        util::SiPtrVector< const SSE>::const_iterator
      > itr_pair
      (
        std::equal_range
        (
          sse_vector.Begin(),
          sse_vector.End(),
          TARGET_SSE,
          SSELessThanNoOverlap()
        )
      );

      // if such an sse is found
      if( itr_pair.first != sse_vector.End())
      {
        // return the overlapping sses
        return util::SiPtrList< const SSE>( itr_pair.first, itr_pair.second);
      }

      // if no match is found return empty
      return util::SiPtrList< const SSE>();
    }

    //! @brief get Identification of this Domain
    //! @return string with identification
    std::string DomainInterface::GetIdentification() const
    {
      // get the vector of all SSEs
      const util::SiPtrVector< const SSE> sse_vector( GetSSEs());

      // initialize string to hold the identification
      std::string identification( "#SSEs: " + util::Format()( sse_vector.GetSize()) + "\n");

      // iterate over SSEs
      for
      (
        util::SiPtrVector< const SSE>::const_iterator sse_itr( sse_vector.Begin()),
          sse_itr_end( sse_vector.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // appned the identification of this SSE to the identification of domain
        identification += ( *sse_itr)->GetIdentification() + " \n";
      }

      // end
      return identification;
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace assemble
} // namespace bcl
