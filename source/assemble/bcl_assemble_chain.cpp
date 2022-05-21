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
#include "assemble/bcl_assemble_chain.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_factory_conformation.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "fold/bcl_fold_mutate_domain_merge_consecutive_ss_types.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Chain::s_Instance
    (
      GetObjectInstances().AddInstance( new Chain())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Chain::Chain() :
      m_Sequence(),
      m_Data()
    {
    }

    //! @brief construct from util::ShPtr to SEQUENCE and ShPtrVector of SSEs
    //! @param SP_SEQUENCE util::ShPtr to sequence
    //! @param SSE_VECTOR ShPtrVector of SSEs
    Chain::Chain( const util::ShPtr< biol::AASequence> &SP_SEQUENCE, const util::ShPtrVector< SSE> &SSE_VECTOR) :
      m_Sequence( SP_SEQUENCE),
      m_Data()
    {
      Insert( SSE_VECTOR);
    }

    //! @brief  construct from util::ShPtr to SEQUENCE and DOMAIN
    //! @param SP_SEQUENCE util::ShPtr to sequence
    //! @param SOURCE_DOMAIN Domain which contains SSEs
    Chain::Chain( const util::ShPtr< biol::AASequence> &SP_SEQUENCE, const Domain &SOURCE_DOMAIN) :
      m_Sequence( SP_SEQUENCE),
      m_Data()
    {
      // insert the SSEs
      Insert( SOURCE_DOMAIN);
    }

    //! @brief construct from util::ShPtr to SEQUENCE with no SSE information
    //! @param SP_SEQUENCE util::ShPtr to sequence
    Chain::Chain( const util::ShPtr< biol::AASequence> &SP_SEQUENCE) :
      m_Sequence( SP_SEQUENCE),
      m_Data()
    {
    }

    //! @brief copy constructor
    //! @param CHAIN_RHS Chain to be copied
    Chain::Chain( const Chain &CHAIN_RHS) :
      m_Sequence( CHAIN_RHS.m_Sequence),
      m_Data( CHAIN_RHS.m_Data)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new copy of this Chain
    Chain *Chain::Clone() const
    {
      return new Chain( *this);
    }

    //! @brief hardcopy
    Chain *Chain::HardCopy() const
    {
      // make a new chain with just the sequence
      Chain *hard_copy( new Chain( m_Sequence.HardCopy()));

      // iterate over all sses
      for
      (
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::const_iterator
          sse_itr( m_Data.Begin()), sse_itr_end( m_Data.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // insert a hardcopy of the shared pointer to sse into the new domain
        hard_copy->Insert( sse_itr->HardCopy());
      }

      // end
      return hard_copy;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @brief the class name as const ref std::string
    const std::string &Chain::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the identification of the sequence and sses
    //! @return the identification of the sequence and sses
    std::string Chain::GetIdentification() const
    {
      // initialize string to be sequence identification
      std::string identification( m_Sequence->GetSequenceIdentification() + "\n");

      // iterate through the sses
      for
      (
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::const_iterator sse_itr( m_Data.Begin()),
          sse_itr_end( m_Data.End());
        sse_itr != sse_itr_end; ++sse_itr
      )
      {
        // add the sse identification to the string
        identification += "\t" + ( *sse_itr)->GetIdentification() + "\n";
      }

      // end
      return identification;
    }

    //! @brief set chainID to CHAINID
    //! @param CHAINID chainID to be set to
    void Chain::SetChainID( const char CHAINID)
    {
      // set chain id for the sequence after hardcopying the AAData
      m_Sequence = util::ShPtr< biol::AASequence>( m_Sequence->HardCopy());
      m_Sequence->SetChainID( CHAINID);

      // initialize set of SSEs
      storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap> new_sses;

      // iterate through the SSE data
      storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::iterator sse_itr;
      while( ( sse_itr = m_Data.Begin()) != m_Data.End())
      {
        // hardcopy the SSE
        util::ShPtr< SSE> current_sse( ( *sse_itr)->HardCopy());

        // remove it from the data
        m_Data.RemoveElement( sse_itr);

        // set the chain id on the hardcopy and insert it into the set
        current_sse->SetChainID( CHAINID);
        new_sses.Insert( current_sse);
      }

      // update the SSE data with the new SSEs
      m_Data = new_sses;
    }

    //! @brief get SiPtrVector of SSEs
    //! @return SiPtrVector of SSEs
    util::SiPtrVector< const SSE> Chain::GetSSEs() const
    {
      // construct a SiPtrVector of const SSEs and return it
      return util::SiPtrVector< const SSE>( m_Data.Begin(), m_Data.End());
    }

    //! @brief returns all SSEs in domain of given SSTYPE in a util::SiPtrVector
    //! @param SS_TYPE specific SSTYPE
    //! @return all SSEs in domain of given SSTYPE in a util::SiPtrVector
    util::SiPtrVector< const SSE> Chain::GetSSEs( const biol::SSType &SS_TYPE) const
    {
      // initialize vector of sses
      util::SiPtrVector< const SSE> sses;

      //loop over all SSElements
      for
      (
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::const_iterator sse_itr( m_Data.Begin()),
          sse_itr_end( m_Data.End());
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

      return sses;
    }

    //! @brief returns SiPtrVector of all SSEs of given SS_TYPES in the set
    //! @param SS_TYPES set of SSTypes of interest
    //! @brief returns SiPtrVector of all SSEs of given SS_TYPES in the set
    util::SiPtrVector< const SSE> Chain::GetSSEs( const storage::Set< biol::SSType> &SS_TYPES) const
    {
      // initialize vector of sses
      util::SiPtrVector< const SSE> sses;

      //loop over all SSElements
      for
      (
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::const_iterator sse_itr( m_Data.Begin()),
          sse_itr_end( m_Data.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // if SSE behind sse_itr is of any given SSTYPES
        if( SS_TYPES.Find( ( *sse_itr)->GetType()) != SS_TYPES.End())
        {
          // insert it into the list of SSEs to be returned
          sses.PushBack( util::SiPtr< const SSE>( *sse_itr));
        }
      }

      return sses;
    }

    //! @brief find and to return the ShPtr for the given SSE
    //! @param SSE_TO_SEARCH SSE of interest
    //! @return ShPtr to corresponding SSE, otherwise an empty ShPtr
    const util::ShPtr< SSE> &Chain::FindSSE( const SSE &SSE_TO_SEARCH) const
    {
      // static undefined SSE ptr
      static const util::ShPtr< SSE> s_undefined_sse_ptr;

      // search for this sse in m_Data
      storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::iterator itr
      (
        std::find_if
        (
          m_Data.Begin(),
          m_Data.End(),
          SSECompare( SSE_TO_SEARCH)
        )
      );

      // if such an sse is found
      if( itr != m_Data.End())
      {
        // return the ShPtr
        return *itr;
      }

      // otherwise return empty ShPtr
      return s_undefined_sse_ptr;
    }

    //! @brief return number of SSE of specified SSTYPE
    //! @param SS_TYPE specific SSTYPE
    //! @return number of SSE of specified SSTYPE
    size_t Chain::GetNumberSSE( const biol::SSType &SS_TYPE) const
    {
      // initialize count
      size_t number_sses( 0);

      //loop over all SSElements
      for
      (
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::const_iterator sse_itr( m_Data.Begin()),
          sse_itr_end( m_Data.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // if of correct type increment
        number_sses += ( ( *sse_itr)->GetType() == SS_TYPE);
      }

      // end
      return number_sses;
    }

    //! @brief returns the geometric center of the object
    //! @return the geometric center of the object
    linal::Vector3D Chain::GetCenter() const
    {
      return coord::CenterOfMass( GetAtomCoordinates());
    }

    //! @brief returns the number of amino acids in the chain
    //! @return the number of amino acids in the chain
    size_t Chain::GetNumberAAs() const
    {
      // initialize size of chain
      size_t size( 0);

      // iterate over SSEs and get their size
      for
      (
         auto sse_itr( m_Data.Begin()), sse_itr_end( m_Data.End());
         sse_itr != sse_itr_end;
         ++sse_itr
      )
      {
        size += ( *sse_itr)->GetData().GetSize();
      }

      return ( size);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief translates the coordinates of all SSEs by the supplied translation vector
    //! @param TRANSLATION_VECTOR_3D Translation vector to be applied
    void Chain::Translate( const linal::Vector3D &TRANSLATION_VECTOR_3D)
    {
      //loop over all secondary structure elements and transform them
      for
      (
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::iterator sse_itr( m_Data.Begin()),
          sse_itr_end( m_Data.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // make a copy of the ShPtr
        util::ShPtr< SSE> this_sse( *sse_itr);

        // translate this_sse
        this_sse->Translate( TRANSLATION_VECTOR_3D);
      }
    }

    //! @brief transforms the coordinates of all SSEs according to given transformation matrix
    //! @param TRANSFORMATION_MATRIX_3D transformation matrix to be applied
    void Chain::Transform( const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D)
    {
      //loop over all secondary structure elements and transform them
      for
      (
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::iterator sse_itr( m_Data.Begin()),
          sse_itr_end( m_Data.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // make a copy of the ShPtr
        util::ShPtr< SSE> this_sse( *sse_itr);

        // transform this_sse
        this_sse->Transform( TRANSFORMATION_MATRIX_3D);
      }
    }

    //! @brief rotate the SSE by a given rotation matrix
    //! @param ROTATION_MATRIX_3D rotation matrix to be applied
    void Chain::Rotate( const math::RotationMatrix3D &ROTATION_MATRIX_3D)
    {
      //loop over all secondary structure elements and transform them
      for
      (
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::iterator sse_itr( m_Data.Begin()),
          sse_itr_end( m_Data.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // make a copy of the ShPtr
        util::ShPtr< SSE> this_sse( *sse_itr);

        // transform this_sse
        this_sse->Rotate( ROTATION_MATRIX_3D);
      }
    }

    //! @brief insert the given NEW_SSE into this chain
    //! @param NEW_SSE SSE to be inserted
    //! @return whether insertion was successful
    bool Chain::Insert( const util::ShPtr< SSE> &NEW_SSE)
    {
      //check that chain ID of SSE fits this chain
      if( NEW_SSE->GetChainID() != m_Sequence->GetChainID())
      {
        BCL_MessageCrt
        (
          std::string( "impossible to insert SSE of different chain: this: ") + m_Sequence->GetChainID() +
          " != SSE: " + NEW_SSE->GetChainID()
        );
        return false;
      }

      // call insert function of the set, it will take care if there are overlap and return false
      return m_Data.Insert( NEW_SSE).second;
    }

    //! @brief insert the SSEs in the given NEW_SSE_VECTOR into this chain
    //! @param NEW_SSE_VECTOR Vector of SSEs to be inserted
    //! @return whether insertion was successful
    bool Chain::Insert( const util::ShPtrVector< SSE> &NEW_SSE_VECTOR)
    {
      // initialize status
      bool success( true);

      // iterate over the given SSEs
      for
      (
        util::ShPtrVector< SSE>::const_iterator sse_itr( NEW_SSE_VECTOR.Begin()), sse_itr_end( NEW_SSE_VECTOR.End());
        sse_itr != sse_itr_end; ++sse_itr
      )
      {
        success &= Insert( *sse_itr);
      }

      // end
      return success;
    }

    //! @brief insert the SSEs from the given domain
    //! @param NEW_DOMAIN Domain from which SSEs should be inserted
    //! @return whether insertion was successful
    bool Chain::Insert( const Domain &NEW_DOMAIN)
    {
      // initialize status
      bool success( true);

      // iterate over the given SSEs
      for
      (
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::const_iterator
          sse_itr( NEW_DOMAIN.GetData().Begin()), sse_itr_end( NEW_DOMAIN.GetData().End());
        sse_itr != sse_itr_end; ++sse_itr
      )
      {
        success &= Insert( *sse_itr);
      }

      // end
      return success;
    }

    //! @brief replace the given SP_SSE with already existing one
    //! @param SP_SSE ShPtr pointing to the SSE to be replaced
    bool Chain::Replace( const util::ShPtr< SSE> &SP_SSE)
    {
      // search for this sse in m_Data
      storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::iterator itr
      (
        std::find_if
        (
          m_Data.Begin(),
          m_Data.End(),
          SSECompare( *SP_SSE)
        )
      );

      // if such an sse is found
      if( itr != m_Data.End())
      {
        // remove this sse
        m_Data.RemoveElement( itr);

        // insert the new sse
        m_Data.Insert( SP_SSE);

        // return
        return true;
      }

      // if not found do nothing and return false
      return false;
    }

    //! @brief replace all SSEs that overlap with SP_SSE with SP_SSE
    //! @param SP_SSE ShPtr to SSE to be inserted
    //! @return whether replacement succeeded
    bool Chain::ReplaceWithOverlapping( const util::ShPtr< SSE> &SP_SSE)
    {
      // search for overlapping sses with SP_SSE using equal_range
      std::pair
      <
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::iterator,
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::iterator
      > itr_pair
      (
        std::equal_range
        (
          m_Data.Begin(),
          m_Data.End(),
          SP_SSE,
          SSELessThanNoOverlap()
        )
      );

      // if such an sse is found
      if( itr_pair.first != m_Data.End())
      {
        // remove the sses in the range
        m_Data.Erase( itr_pair.first, itr_pair.second);

        // insert the new sse
        m_Data.Insert( SP_SSE);

        // return
        return true;
      }

      BCL_MessageCrt
      (
        "ReplaceWithOverlapping has been called with a SSE that is not overlapping!"
      );

      // if not found do nothing and return false
      return false;
    }

    //! @brief remove given SSELEMENT from the domain
    //! @param SSELEMENT SSE to be removed
    bool Chain::Remove( const SSE &SSELEMENT)
    {
      // search for this sse in m_Data
      storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::iterator itr
      (
        std::find_if
        (
          m_Data.Begin(),
          m_Data.End(),
          SSECompare( SSELEMENT)
        )
      );

      // if such an sse is found
      if( itr != m_Data.End())
      {
        // remove this sse
        m_Data.RemoveElement( itr);

        // return
        return true;
      }

      // if not found do nothing and return false
      return false;
    }

    //! @brief sets positions of all SSEs to ideal conformation w/wo superimposing with prior coordinates
    //! @param KEEP_POSITION flag to indicate whether to original body information of SSEs should be reserved
    void Chain::SetToIdealConformation( const bool KEEP_POSITION)
    {
      //set each sse to ideal conformation
      for
      (
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::iterator sse_itr( m_Data.Begin()),
          sse_itr_end( m_Data.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // make a non-const ShPtr to this SSE since Set does not allow non-const iterators
        util::ShPtr< SSE> this_sse( *sse_itr);

        if( KEEP_POSITION)
        {
          this_sse->SetToIdealConformationInPlace();
        }
        else
        {
          this_sse->SetToIdealConformationAtOrigin();
        }
      }
    }

    //! @brief chop all SSEs of that model in pieces of the sizes defined by MIN_SSE_LENGTHS
    //! @param MIN_SSE_LENGTHS VectorND of sizes that defined min size for each SSType
    void Chain::ChopSSEs( const storage::VectorND< 3, size_t> &MIN_SSE_LENGTHS)
    {
      // instatiate set of chopped_sse
      storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap> chopped_sses;

      //iterate over all sselements
      for
      (
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::const_iterator sse_itr( m_Data.Begin()),
          sse_itr_end( m_Data.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // chop this sse and insert into new set of sses
        chopped_sses.InsertElements( ( *sse_itr)->Chop( MIN_SSE_LENGTHS( ( *sse_itr)->GetType())));
      }

      //set m_Model to new chopped elements
      m_Data = chopped_sses;
    }

    //! @brief checks if domain already contains THIS_SSE
    //! @param THIS_SSE SSE of interest
    //! @return whether domain already contains THIS_SSE
    bool Chain::DoesContain( const SSE &THIS_SSE) const
    {
      // search for this sse in m_Data
      storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::const_iterator itr
      (
        std::find_if
        (
          m_Data.Begin(),
          m_Data.End(),
          SSECompare( THIS_SSE)
        )
      );

      // return if it's found
      return ( itr != m_Data.End());
    }

    //! @brief checks if domain already contains this THIS_SSE or any overlapping SSE with THIS_SSE
    //! @param THIS_SSE SSE to be searched for
    //! @return if domain already contains this THIS_SSE or any overlapping SSE with THIS_SSE
    bool Chain::DoesContainOverlapping( const SSE &THIS_SSE) const
    {
      // search for this sse or any overlapping sse in m_Data
      storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::const_iterator itr
      (
        std::find_if
        (
          m_Data.Begin(),
          m_Data.End(),
          SSECompareOverlap( THIS_SSE)
        )
      );

      // return if it's found
      return ( itr != m_Data.End());
    }

    //! @brief join following ( progressing sequence id) SSEs of given SS_TYPE into one SSE
    //! @param SS_TYPE SSType of interest
    //! @param TEST_PEPTIDE_BOND join only, if SSEs are connected by peptide bond
    void Chain::Join( const biol::SSType &SS_TYPE, const bool TEST_PEPTIDE_BOND)
    {
      if( GetNumberSSE( SS_TYPE) < 2)
      {
        return;
      }

      // new SSEs after joining
      storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap> joined_sses;

      //set each SSE to ideal conformation
      storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::iterator
      sse_itr( m_Data.Begin()), sse_itr_next( m_Data.Begin()), sse_itr_end( m_Data.End());
      ++sse_itr_next;

      // store the current SSE
      util::ShPtr< SSE> current_sse( ( *sse_itr));

      // initialize a boolean
      bool last_sse_joinable( false);

      // iterate until anyone of the iterators hit the end
      while( sse_itr_next != sse_itr_end && sse_itr != sse_itr_end)
      {
        // if this SSE is of specified type as well as the next SSE and they have no loop residues in between
        if
        (
              current_sse->GetType() == SS_TYPE
           && current_sse->GetType() == ( *sse_itr_next)->GetType()
           && ( current_sse->GetLastAA()->GetSeqID() + 1) == ( *sse_itr_next)->GetFirstAA()->GetSeqID()
           && ( !TEST_PEPTIDE_BOND || biol::AABase::AreAminoAcidsPeptideBonded( *current_sse->GetLastAA(), *( *sse_itr_next)->GetFirstAA(), true))
        )
        {
          // update the sequence of the current SSE by appending the sequence of the next SSE that is connected
          current_sse->AppendSequence( ( **sse_itr_next), false);

          // set the flag
          last_sse_joinable = true;

          // move to next SSE
          ++sse_itr_next;
        }
        // if next SSE is not connected
        else
        {
          // but the previous pair of SSEs were connected
          if( last_sse_joinable)
          {
            // update the geometry
            current_sse->SetGeometry();
            // set the flag back
            last_sse_joinable = false;
          }

          // insert the SSE joined or not into the new list of SSEs
          joined_sses.Insert( current_sse);

          // update the iterators
          current_sse = *sse_itr_next;
          sse_itr = sse_itr_next;
          ++sse_itr_next;
        }
      }

      // if the last SSE was connected
      if( last_sse_joinable)
      {
        // update geometry
        current_sse->SetGeometry();
      }
      // insert the last SSE into new list of SSEs
      joined_sses.Insert( current_sse);

      // update the data
      m_Data = joined_sses;
    }

    //! @brief filters the current chain by given minimum SSE sizes
    //! @param MIN_SSE_SIZES minimum SSE sizes to filter the chain by
    void Chain::FilterByMinSSESizes( const storage::Map< biol::SSType, size_t> &MIN_SSE_SIZES)
    {
      // create a new set to store sses
      storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap> filtered_sses;

      // iterate over the SSEs
      for
      (
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::iterator
          sse_itr( m_Data.Begin()), sse_itr_end( m_Data.End());
        sse_itr != sse_itr_end; ++sse_itr
      )
      {
        storage::Map< biol::SSType, size_t>::const_iterator
          type_itr( MIN_SSE_SIZES.Find( ( *sse_itr)->GetType()));

        // if this SSE is one of the specified types in the given map and
        // if this SSE has a size larger than/equal to the specified size
        if
        (
          type_itr != MIN_SSE_SIZES.End() && ( *sse_itr)->GetSize() >= type_itr->second
        )
        {
          // insert it into filtered SSEs
          filtered_sses.Insert( *sse_itr);
        }
      }

      // update the SSE set with the filtered one
      m_Data = filtered_sses;
    }

    //! @brief AddLoops generates loop SSEs for the Chain and assigns them undefined coordinates
    //! @param UNDEFINED_COORDINATES create loop with undefined coordinates, or with coordinates form the member sequence
    //! @param MERGE_CONSECUTIVE_SSES merge consecutive SSEs of given type
    //! @param SS_TYPE the sstype of the SSEs to be merged
    void Chain::AddLoops
    (
      const bool UNDEFINED_COORDINATES,
      const bool MERGE_CONSECUTIVE_SSES,
      const biol::SSType &SS_TYPE
    )
    {
      // create ShPtrVector iterator "seq_itr" and "seq_itr_end" for iterating over the aa sequence of this chain
      biol::AASequence::const_iterator seq_itr( m_Sequence->Begin()), seq_itr_end( m_Sequence->End());

      // check that there is a least one sse
      if( m_Data.IsEmpty())
      {
        // create ShPtr to SSE "new_sse"
        util::ShPtr< SSE> new_sse;

        if( !UNDEFINED_COORDINATES)
        {
          new_sse = util::ShPtr< SSE>( new SSE( *m_Sequence, biol::GetSSTypes().COIL));
        }
        else
        {
          // create ShPtr to SSE "new_sse" and initialize with a new SSE of type COIL
          new_sse = util::ShPtr< SSE>( new SSE( biol::GetSSTypes().COIL));

          // set the chain ID of "new_sse"
          new_sse->SetChainID( GetChainID());

          // set the fasta header of "new_sse"
          new_sse->SetFastaHeader( m_Sequence->GetFastaHeader());

          // copy all amino acids, without coordinates
          for( ; seq_itr != seq_itr_end; ++seq_itr)
          {
            new_sse->PushBack( util::ShPtr< biol::AABase>( ( *seq_itr)->Empty( ( *seq_itr)->GetData())));
          }
        }

        m_Data.Insert( new_sse);

        return;
      }

      // create Set of SSEs "new_sses" to hold the SSEs which are created
      storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap> new_sses;

      // loop over all SSEs in this chain
      for
      (
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::const_iterator
          sse_itr( m_Data.Begin()), sse_itr_end( m_Data.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // create size_t "sse_begin_seq_id" initialize with the SeqID of the first aa in the SSE denoted by "sse_itr"
        const int sse_begin_seq_id( ( *sse_itr)->GetFirstAA()->GetSeqID());

        // create ShPtr to SSE "new_sse" and initialize with a new SSE of type COIL
        util::ShPtr< SSE> new_sse( new SSE( biol::GetSSTypes().COIL));

        // set the chain ID of "new_sse"
        new_sse->SetChainID( GetChainID());

        // set the fasta header of "new_sse"
        new_sse->SetFastaHeader( m_Sequence->GetFastaHeader());

        // while seq ID of aa denoted by "seq_itr" has not reached "sse_begin_seq_id" and the end of the aa sequence
        // has not been reached
        while( ( *seq_itr)->GetSeqID() != sse_begin_seq_id && seq_itr != seq_itr_end)
        {
          // create a new amino acid of the corresponding aa class (should set coordinates to undefined)
          util::ShPtr< biol::AABase> new_aa;

          if( UNDEFINED_COORDINATES)
          {
            new_aa = util::ShPtr< biol::AABase>( ( *seq_itr)->Empty( ( *seq_itr)->GetData()));
          }
          else
          {
            new_aa = util::ShPtr< biol::AABase>( ( *seq_itr)->Clone());
          }

          // add the new aa to "new_sse"
          new_sse->PushBack( new_aa);

          // move "seq_itr" to the next aa
          ++seq_itr;
        }

        // true if "new_sse" has some aas in it
        if( new_sse->GetSize())
        {
          // insert "new_sse" into "new_sses"
          new_sses.Insert( new_sse);
        }

        // move "seq_itr" to the place of the next non-sse region
        const size_t sse_size( ( *sse_itr)->GetSize());
        for( size_t sse_index( 0); sse_index != sse_size && seq_itr != seq_itr_end; ++sse_index)
        {
          ++seq_itr;
        }
      }

      // we have to make sure to add the loop after the last SSE
      const size_t last_sse_end_id( ( *m_Data.ReverseBegin())->GetLastAA()->GetSeqID());
      const size_t last_aa_id( m_Sequence->GetLastAA()->GetSeqID());

      // if there are any residues at the end of the sequence that are not in the last SSE
      if( last_sse_end_id < last_aa_id)
      {
        // create ShPtr to SSE "new_sse" and initialize with a new SSE of type COIL
        util::ShPtr< SSE> new_sse( new SSE( biol::GetSSTypes().COIL));

        // set the chain ID of "new_sse"
        new_sse->SetChainID( GetChainID());

        // set the fasta header of "new_sse"
        new_sse->SetFastaHeader( m_Sequence->GetFastaHeader());

        // iterate until you hit the end of amino acids in the sequecne
        while( seq_itr != seq_itr_end)
        {
          // create a new amino acid of the corresponding aa class (should set coordinates to 0.000)
          util::ShPtr< biol::AABase> new_aa;

          if( UNDEFINED_COORDINATES)
          {
            new_aa = util::ShPtr< biol::AABase>( ( *seq_itr)->Empty( ( *seq_itr)->GetData()));
          }
          else
          {
            new_aa = util::ShPtr< biol::AABase>( ( *seq_itr)->Clone());
          }

          // add the new aa to "new_sse"
          new_sse->PushBack( new_aa);

          // move "seq_itr" to the next aa
          ++seq_itr;
        }

        // push the last loop
        new_sses.Insert( new_sse);
      }

      // if consecutive SSEs are to be merged, do it with the mutate
      if( MERGE_CONSECUTIVE_SSES)
      {
        const fold::MutateDomainMergeConsecutiveSSTypes mutate( SS_TYPE);
        Domain new_domain( new_sses);

        new_sses = mutate( new_sses).GetArgument()->GetData();
      }

      // add all the new sses in "new_sses" to this chain
      for
      (
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::const_iterator
          sse_itr( new_sses.Begin()), sse_itr_end( new_sses.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        Insert( *sse_itr);
      }
    }

    //! @brief let the aadata of the each sse point to the corresponding data in the chain - determined by pdbID
    void Chain::ConnectSSEToChainData()
    {
      // set which will hold new sses which have been connected to the data in the chain
      storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap> new_sses;

      // iterate through the sses of this chain
      for
      (
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::const_iterator
          sse_itr( m_Data.Begin()), sse_itr_end( m_Data.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // get a sequence copied from the sse
        biol::AASequence sequence( **sse_itr);

        // connect the sequence from the sse to the data in the chain sequence
        m_Sequence->ConnectAADataByPdbID( sequence);

        // make a new sse out of the connected sequence and the type of sse that sse_itr is
        const util::ShPtr< SSE> new_sse( new SSE( sequence, ( *sse_itr)->GetType()));

        // add the new sse to the set of new sses - don't use replace here or segmentation fault results
        new_sses.Insert( new_sse);
      }

      // set the m_Data to the set of new sses
      m_Data = new_sses;
    }

    //! @brief Get SSE hash string to aid in identifying similar chains
    std::string Chain::GetSSEHashString() const
    {
      std::string hash;
      // iterate through the chains
      for
      (
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::const_iterator
          sse_itr( m_Data.Begin()), sse_itr_end( m_Data.End());
        sse_itr != sse_itr_end;
      )
      {
        hash += ( *sse_itr)->GetHashString();
        if( ++sse_itr != sse_itr_end)
        {
          hash += ',';
        }
      }
      return hash;
    }

    //! @brief Replace SSEs with those drawn from the pool
    void Chain::AdoptSSEsMaintainCoordinates( const util::SiPtrVector< const SSE> &SSES)
    {
      m_Data.Reset();

      // create ShPtrVector iterator "seq_itr" and "seq_itr_end" for iterating over the aa sequence of this chain
      biol::AASequence::const_iterator seq_itr( m_Sequence->Begin()), seq_itr_end( m_Sequence->End());

      for( auto itr_sses( SSES.Begin()), itr_sses_end( SSES.End()); itr_sses != itr_sses_end; ++itr_sses)
      {
        const SSE &sse( **itr_sses);
        for( auto first_seq_id( sse.GetFirstMember()->GetSeqID()); ( *seq_itr)->GetSeqID() < first_seq_id; ++seq_itr)
        {
        }
        biol::AASequence sequence;
        for
        (
          auto last_seq_id( sse.GetLastMember()->GetSeqID());
          seq_itr != seq_itr_end && ( *seq_itr)->GetSeqID() <= last_seq_id;
          ++seq_itr
        )
        {
          sequence.PushBack( *seq_itr);
        }
        this->Insert( util::ShPtr< SSE>( new SSE( sequence, sse.GetType())));
      }
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief equal operator
    //! @param CHAIN_RHS Chain to be assigned
    //! @return this chain after being assigned to given CHAIN_RHS
    Chain &Chain::operator =( const Chain &CHAIN_RHS)
    {
      // set members
      m_Sequence = CHAIN_RHS.m_Sequence;
      m_Data = CHAIN_RHS.m_Data;

      // return
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read Chain from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Chain::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Sequence, ISTREAM);
      io::Serialize::Read( m_Data, ISTREAM);

      //end
      return ISTREAM;
      }

    //! @brief write Chain to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &Chain::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Sequence, OSTREAM, INDENT);
      io::Serialize::Write( m_Data, OSTREAM, INDENT);

      //end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief construct a chain from a sequence with SSEs
    //! the phi and psi angles of the backbone are calculated and depending on the position in the ramachandran plot,
    //! sses are constructed
    //! @author woetzen
    //! @param SP_SEQUENCE ShPtr to sequence- at least the backbone atoms, otherwise they will be indentified as loops
    //! @return Chain with given sequence and identified sses
    Chain ConstructChainWithSSEsFromConformation( const util::ShPtr< biol::AASequence> &SP_SEQUENCE)
    {
      // create an sse pool using SSEFactoryConformation
      const SSEPool sse_pool( SSEFactoryConformation()( *SP_SEQUENCE));

      // collection of identified sses
      const Domain identified_sses
      (
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>( sse_pool.Begin(), sse_pool.End())
      );

      // construct chain from sequence and domain
      return Chain( SP_SEQUENCE, identified_sses);
    }

  ///////////////////
  // ChainLessThan //
  ///////////////////

    //! @brief return whether one chain  is less than another
    //! @param CHAIN_LHS first chain
    //! @param CHAIN_RHS second chain
    //! @return whether one chain is less than another
    bool ChainLessThan::operator()( const Chain &CHAIN_LHS, const Chain &CHAIN_RHS) const
    {
      // compare chain ids
      return CHAIN_LHS.GetChainID() < CHAIN_RHS.GetChainID();
    }

    //! @brief return whether one chain is less than another
    //! @param PTR_CHAIN_LHS first chain
    //! @param PTR_CHAIN_RHS second chain
    //! @return whether one chain is less than another
    bool ChainLessThan::operator()
    (
      const util::PtrInterface< Chain> &PTR_CHAIN_LHS,
      const util::PtrInterface< Chain> &PTR_CHAIN_RHS
    ) const
    {
      return operator ()( *PTR_CHAIN_LHS, *PTR_CHAIN_RHS);
    }

  } // namespace assemble
} // namespace bcl
