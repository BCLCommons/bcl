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
#include "assemble/bcl_assemble_domain.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Domain::s_Instance
    (
      GetObjectInstances().AddInstance( new Domain())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Domain::Domain() :
      m_Data(),
      m_Topology()
    {
    }

    //! @brief construct from a Set of ShPtr to SSEs
    //! @param SSE_SET Set of ShPtr to SSEs
    Domain::Domain( const storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap> &SSE_SET) :
      m_Data( SSE_SET),
      m_Topology()
    {
    }

    //! @brief construct from a Set of ShPtr to SSEs and a corresponding topology
    //! @param SSE_SET Set of ShPtr to SSEs
    //! @param SP_TOPOLOGY ShPtr to corresponding Topology
    Domain::Domain
    (
      const storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap> &SSE_SET,
      const util::ShPtr< Topology> &SP_TOPOLOGY
    ) :
      m_Data( SSE_SET),
      m_Topology( SP_TOPOLOGY)
    {
    }

    //! @brief construct from util::ShPtrVector of SSE
    //! @param SSE_VECTOR ShPtrVector of SSEs
    Domain::Domain( const util::ShPtrVector< SSE> &SSE_VECTOR) :
      m_Data(),
      m_Topology()
    {
      // iterate over given vector of SSEs and insert them into the set m_Data
      for
      (
        util::ShPtrVector< SSE>::const_iterator sse_itr( SSE_VECTOR.Begin()),
          sse_itr_end( SSE_VECTOR.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        Insert( *sse_itr);
      }
    }

    //! @brief construct from util::ShPtrVector of SSE and a corresponding topology
    //! @param SSE_VECTOR ShPtrVector of SSEs
    //! @param SP_TOPOLOGY ShPtr to corresponding Topology
    Domain::Domain
    (
      const util::ShPtrVector< SSE> &SSE_VECTOR,
      const util::ShPtr< Topology> &SP_TOPOLOGY
    ) :
      m_Data(),
      m_Topology( SP_TOPOLOGY)
    {
      // iterate over given vector of SSEs and insert them into the set m_Data
      for
      (
        util::ShPtrVector< SSE>::const_iterator sse_itr( SSE_VECTOR.Begin()),
          sse_itr_end( SSE_VECTOR.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        Insert( *sse_itr);
      }
    }

    //! @brief copy constructor for domain
    //! @param DOMAIN_RHS Domain to be copied
    Domain::Domain( const Domain &DOMAIN_RHS) :
      m_Data( DOMAIN_RHS.m_Data),
      m_Topology( DOMAIN_RHS.m_Topology)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new Domain instance copied from this one
    Domain *Domain::Clone() const
    {
      return new Domain( *this);
    }

    //! @brief virtual hard copy constructor
    //! @return pointer to a new Domain instance hard-copied from this one
    Domain *Domain::HardCopy() const
    {
      // make hard_ccopies of the SSEs
      storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap> new_sses;

      // iterate over all sses
      for
      (
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::const_iterator sse_itr( m_Data.Begin()),
          sse_itr_end( m_Data.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // insert a hardcopy of the shared pointer to sse into the new domain
        new_sses.Insert( sse_itr->HardCopy());
      }

      // if a topology was defined
      if( m_Topology.IsDefined())
      {
        // create a new topology
        // this topology will only have the geometry since
        // the graph would need to point to new SSEs
        // and there is no guarantee they won't be changed at this point
        util::ShPtr< Topology> sp_new_topology( new Topology( m_Topology->GetOrientation()));

        // now create a domain topology with the SSEs and the topology
        Domain *hard_copy( new Domain( new_sses, sp_new_topology));
        // and return
        return hard_copy;
      }
      // if no topology
      else
      {
        // create a new domain only from the SSEs
        Domain *hard_copy( new Domain( new_sses));
        // and return
        return hard_copy;
      }
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Domain::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return all chain ids within that domain
    //! @return set of all chain ids
    storage::Set< char> Domain::GetChainIds() const
    {
      storage::Set< char> chain_ids;
      for
      (
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::const_iterator sse_itr( m_Data.Begin()),
          sse_itr_end( m_Data.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        chain_ids.Insert( ( *sse_itr)->GetChainID());
      }

      return chain_ids;
    }

    //! @brief return number of SSE of specified SSTYPE
    //! @param SS_TYPE specific SSTYPE
    //! @return number of SSE of specified SSTYPE
    size_t Domain::GetNumberSSE( const biol::SSType &SS_TYPE) const
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

    //! @brief returns all SSEs in domain in a util::SiPtrVector
    //! @return all SSEs in domain in a util::SiPtrVector
    util::SiPtrVector< const SSE> Domain::GetSSEs() const
    {
      // construct a SiPtrVector of const SSEs and return it
      return util::SiPtrVector< const SSE>( m_Data.Begin(), m_Data.End());
    }

    //! @brief returns all SSEs in domain of given SSTYPE in a util::SiPtrVector
    //! @param SS_TYPE specific SSTYPE
    //! @return all SSEs in domain of given SSTYPE in a util::SiPtrVector
    util::SiPtrVector< const SSE> Domain::GetSSEs( const biol::SSType &SS_TYPE) const
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

    //! @brief returns the number of amino acids in the chain
    //! @return the number of amino acids in the chain
    size_t Domain::GetNumberAAs() const
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

      return size;
    }

    //! @brief concatenates sequences of all sses fills in the gaps with unknown AAs and returns this AASequence
    //! @return AASequence created by concatenating sequences of all sses and filling in the gaps with unknown AAs
    util::ShPtr< biol::AASequence> Domain::CreateSequenceFromSSEs() const
    {
      BCL_Assert( !m_Data.IsEmpty(), "It is not possible to create Sequence from 0 SSEs!!!");

      // initialize this sequence
      util::ShPtr< biol::AASequence> this_sequence( new biol::AASequence());
      const char chain_id( ( *m_Data.Begin())->GetChainID());

      // set the chainid to the chainid of the first sse
      this_sequence->SetChainID( chain_id);
      int last_seqid( 1);

      // iterate over sses
      for
      (
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::const_iterator sse_itr( m_Data.Begin()),
          sse_itr_end( m_Data.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // make sure the chainid of this SSE is same with the previous ones
        BCL_Assert
        (
          ( *sse_itr)->GetChainID() == chain_id,
          "This sse's chainid(" + util::Format()( ( *sse_itr)->GetChainID()) +
          ") does not match the sequence's chainid(" + util::Format()( chain_id) + ")!!"
        );

        // iterate over skipped residues from the end of last sse to the beginning of this sse
        for( int seq_id( last_seqid), begin_sse( ( *sse_itr)->GetFirstAA()->GetSeqID()); seq_id < begin_sse; ++seq_id)
        {
          util::ShPtr< biol::AAData> sp_aa_data
          (
            new biol::AAData
            (
              biol::GetAATypes().e_Undefined,
              seq_id,
              util::GetUndefined< int>(), // pdb id is unknown, since the amino acids added to the loop are not actually read in from a pdb
              biol::AAData::s_DefaultPdbICode,
              chain_id
            )
          );
          util::ShPtr< biol::AABase> new_aa( ( *( *sse_itr)->GetFirstAA()->GetAAClass())->Empty( sp_aa_data));
          // insert unknown residues with the corresponding seqid
          this_sequence->PushBack( new_aa);
        }

        // now append the sequence of this sse
        this_sequence->AppendSequence( **sse_itr);

        // update the last_seqid
        last_seqid = ( *sse_itr)->GetLastAA()->GetSeqID() + 1;
      }

      // return the created sequence
      return this_sequence;
    }

    //! @brief returns the geometric center of the object
    //! @return the geometric center of the object
    linal::Vector3D Domain::GetCenter() const
    {
      // if the topology is defined
      if( m_Topology.IsDefined())
      {
        return m_Topology->GetCenter();
      }

      // otherwise return center
      return coord::CenterOfMass( GetAtomCoordinates());
    }

    //! @brief return the orientation of the object
    //! @return orientation
    linal::Vector3D Domain::GetAxis( const coord::Axis &AXIS) const
    {
      // if the topology is defined
      if( m_Topology.IsDefined())
      {
        return m_Topology->GetAxis( AXIS);
      }

      // otherwise return center
      linal::Vector3D axis;
      axis( AXIS) = 1.0;

      return axis;
    }

    //! @brief return the orientation and Position as TransformationMatrix3D
    //! @return TransformationMatrix3D that defines orientation and position
    const math::TransformationMatrix3D Domain::GetOrientation() const
    {
      // if the topology is defined
      if( m_Topology.IsDefined())
      {
        return m_Topology->GetOrientation();
      }

      // otherwise return center
      return math::TransformationMatrix3D( GetCenter());
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief translates the coordinates of all SSEs by the supplied translation vector
    //! @param TRANSLATION_VECTOR_3D Translation vector to be applied
    void Domain::Translate( const linal::Vector3D &TRANSLATION_VECTOR_3D)
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
    void Domain::Transform( const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D)
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
    void Domain::Rotate( const math::RotationMatrix3D &ROTATION_MATRIX_3D)
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

    //! @brief pushback a util::ShPtr< SSE> to m_Data
    //! @param SSELEMENT ShPtr to SSE to be inserted
    //! @return whether insertion succeeded
    bool Domain::Insert( const util::ShPtr< SSE> &SSELEMENT)
    {
      // call insert function of the set, it will take care if there are overlap and return false
      return m_Data.Insert( SSELEMENT).second;
    }

    //! @brief pushback a vector SSEs
    //! @param SSE_VECTOR Vector of SSEs to be added
    //! @return whether insertion succeeded
    bool Domain::Insert( const util::ShPtrVector< SSE> &SSE_VECTOR)
    {
      // initalize boolean success
      bool success( true);

      // iterate over the  SSEs
      for
      (
        util::ShPtrVector< SSE>::const_iterator sse_itr( SSE_VECTOR.Begin()), sse_itr_end( SSE_VECTOR.End());
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
    bool Domain::Insert( const Domain &NEW_DOMAIN)
    {
      // initalize boolean success
      bool success( true);

      // iterate over the  SSEs
      for
      (
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>::const_iterator
          sse_itr( NEW_DOMAIN.GetData().Begin()), sse_itr_end( NEW_DOMAIN.GetData().End());
        sse_itr != sse_itr_end; ++sse_itr
      )
      {
        success &= m_Data.Insert( *sse_itr).second;
      }

      // end
      return success;
    }

    //! @brief replace the given SP_SSE with already existing one
    //! @param SP_SSE ShPtr pointing to the SSE to be replaced
    bool Domain::Replace( const util::ShPtr< SSE> &SP_SSE)
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
    bool Domain::ReplaceWithOverlapping( const util::ShPtr< SSE> &SP_SSE)
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
        "ReplaceWithOverlapping has been called with a sse that is not overlapping!"
      );

      // if not found do nothing and return false
      return false;
    }

    //! @brief remove given SSELEMENT from the domain
    //! @param SSELEMENT SSE to be removed
    bool Domain::Remove( const SSE &SSELEMENT)
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
    void Domain::SetToIdealConformation( const bool KEEP_POSITION)
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

    //! @brief chop all sselements of that model in pieces of the sizes defined by MIN_SSE_LENGTHS
    //! @param MIN_SSE_LENGTHS VectorND of sizes that defined min size for each SSType
    void Domain::ChopSSEs( const storage::VectorND< 3, size_t> &MIN_SSE_LENGTHS)
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

    //! @brief find and to return the ShPtr for the given SSE
    //! @param SSE_TO_SEARCH SSE of interest
    //! @return ShPtr to corresponding SSE, otherwise an empty ShPtr
    const util::ShPtr< SSE> &Domain::FindSSE( const SSE &SSE_TO_SEARCH) const
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

    //! @brief checks if domain already contains THIS_SSE
    //! @param THIS_SSE SSE of interest
    //! @return whether domain already contains THIS_SSE
    bool Domain::DoesContain( const SSE &THIS_SSE) const
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
    bool Domain::DoesContainOverlapping( const SSE &THIS_SSE) const
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
    void Domain::Join( const biol::SSType &SS_TYPE)
    {
      if( GetNumberSSE( SS_TYPE) < 2)
      {
        return;
      }

      // new sses after joining
      storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap> joined_sses;

      //set each sse to ideal conformation
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
          && ( current_sse->GetLastAA()->GetSeqID() + 1) == ( *sse_itr_next)->GetFirstAA()->GetSeqID())
        {
          // update the sequence of the current SSE by appending the sequence of the next SSE that is connected
          current_sse->AppendSequence( ( **sse_itr_next));

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
        // update geoemtry
        current_sse->SetGeometry();
      }
      // insert the last SSE into new list of SSEs
      joined_sses.Insert( current_sse);

      // update the data
      m_Data = joined_sses;
    }

    //! @brief filters the current chain by given minimum SSE sizes
    //! @param MIN_SSE_SIZES minimum SSE sizes to filter the chain by
    void Domain::FilterByMinSSESizes( const storage::Map< biol::SSType, size_t> &MIN_SSE_SIZES)
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

  ///////////////
  // operators //
  ///////////////

    //! @brief equal operator
    //! @param DOMAIN_RHS Domain to be assigned to
    //! @return this domain after being assigned to DOMAIN_RHS
    Domain &Domain::operator =( const Domain &DOMAIN_RHS)
    {
      // set members
      m_Data = DOMAIN_RHS.m_Data;
      m_Topology = DOMAIN_RHS.m_Topology;

      // return
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read Domain from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Domain::Read( std::istream &ISTREAM)
    {
      //read member
      io::Serialize::Read( m_Data, ISTREAM);

      //end
      return ISTREAM;
    }

    //! @brief write Domain to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &Domain::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      //write Data
      io::Serialize::Write( m_Data, OSTREAM, INDENT);

      //end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace assemble
} // namespace bcl
