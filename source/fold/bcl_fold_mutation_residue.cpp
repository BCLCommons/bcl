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
#include "fold/bcl_fold_mutation_residue.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutationResidue::s_Instance
    (
      GetObjectInstances().AddInstance( new MutationResidue())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutationResidue::MutationResidue() :
      m_MutationResidue(),
      m_PreviousResidue(),
      m_FollowingResidue()
    {
    }

    //! @brief constructor taking member variables
    //! @param RESIDUE_TO_MUTATE the residue of interest
    //! @param RESIDUE_A residue to either side of the residue to mutate
    //! @param RESIDUE_B residue to either side of residue to mutate
    MutationResidue::MutationResidue
    (
      const util::ShPtr< biol::AABase> &RESIDUE_TO_MUTATE,
      const util::ShPtr< biol::AABase> &RESIDUE_A,
      const util::ShPtr< biol::AABase> &RESIDUE_B
    ) :
      m_MutationResidue( RESIDUE_TO_MUTATE),
      m_PreviousResidue(),
      m_FollowingResidue()
    {
      BCL_Assert( m_MutationResidue.IsDefined(), "m_MutationResidue residue is not defined");
      if( RESIDUE_TO_MUTATE.IsDefined())
      {
        BCL_MessageDbg( "residue to mutate " + RESIDUE_TO_MUTATE->GetIdentification());
      }
      if( RESIDUE_A.IsDefined())
      {
        BCL_MessageDbg( "RESIDUE_A " + RESIDUE_A->GetIdentification());
      }
      if( RESIDUE_B.IsDefined())
      {
        BCL_MessageDbg( "RESIDUE_B " + RESIDUE_B->GetIdentification());
      }
      SortResidues( RESIDUE_A, RESIDUE_B);
    }

    //! @brief constructor from iterator to residue to mutate and the map of all the residues
    //! @param RESIDUE_TO_MUTATE residue indicated to mutate
    //! @param ALL_RESIDUES all residues that could be mutated
    MutationResidue::MutationResidue
    (
      const storage::List< storage::Pair< util::ShPtr< biol::AABase>, bool> >::const_iterator RESIDUE_TO_MUTATE,
      const storage::List< storage::Pair< util::ShPtr< biol::AABase>, bool> > &ALL_RESIDUES
    ) :
      m_MutationResidue( RESIDUE_TO_MUTATE->First()),
      m_PreviousResidue(),
      m_FollowingResidue()
    {
      BCL_Assert( m_MutationResidue.IsDefined(), "mutation residue is not defined");

      if( ALL_RESIDUES.GetSize() < 2)
      {
        BCL_MessageDbg( "single residue");
        return;
      }

      // true if first residue
      if( RESIDUE_TO_MUTATE == ALL_RESIDUES.Begin())
      {
        BCL_MessageDbg( "at first residue");
        SortResidues( util::ShPtr< biol::AABase>(), ( ++ALL_RESIDUES.Begin())->First());
      }

      // true if last residue
      else if( RESIDUE_TO_MUTATE == --ALL_RESIDUES.End())
      {
        BCL_MessageDbg( "at last residue");
        SortResidues( util::ShPtr< biol::AABase>(), ( ----ALL_RESIDUES.End())->First());
      }

      // interior residue
      else
      {
        BCL_MessageDbg( "at interior residue");
        storage::List< storage::Pair< util::ShPtr< biol::AABase>, bool> >::const_iterator itr_a( RESIDUE_TO_MUTATE);
        --itr_a;
        storage::List< storage::Pair< util::ShPtr< biol::AABase>, bool> >::const_iterator itr_b( RESIDUE_TO_MUTATE);
        ++itr_b;
        BCL_MessageDbg( "surrounding residue a " + itr_a->First()->GetIdentification());
        BCL_MessageDbg( "surrounding residue b " + itr_b->First()->GetIdentification());
        BCL_MessageDbg( "RESIDUE_TO_MUTATE still the same? " + RESIDUE_TO_MUTATE->First()->GetIdentification());
        SortResidues( itr_a->First(), itr_b->First());
      }

    }

    //! @brief Clone function
    //! @return pointer to new MutationResidue
    MutationResidue *MutationResidue::Clone() const
    {
      return new MutationResidue( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutationResidue::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief GetPreviousResidue gives the residue preceding the residue of interest in sequence
    //! @return the residue preceding the residue of interest in sequence
    const util::ShPtr< biol::AABase> &MutationResidue::GetPreviousResidue() const
    {
      return m_PreviousResidue;
    }

    //! @brief GetFollowingResidue gives the residue proceding the residue of interest in sequence
    //! @return the residue proceding the residue of interest in sequence
    const util::ShPtr< biol::AABase> &MutationResidue::GetFollowingResidue() const
    {
      return m_FollowingResidue;
    }

    //! @brief GetMutationResidue gives the residue of interest
    //! @return the residue of interest
    const util::ShPtr< biol::AABase> &MutationResidue::GetMutationResidue() const
    {
      return m_MutationResidue;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutationResidue::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_MutationResidue, ISTREAM);
      io::Serialize::Read( m_PreviousResidue, ISTREAM);
      io::Serialize::Read( m_FollowingResidue, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutationResidue::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_MutationResidue, OSTREAM, INDENT);
      io::Serialize::Write( m_PreviousResidue, OSTREAM, INDENT);
      io::Serialize::Write( m_FollowingResidue, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief assigns the member variables to the proper residue
    //! @param RESIDUE_A first residue
    //! @param RESIDUE_B second residue
    void MutationResidue::SortResidues
    (
      const util::ShPtr< biol::AABase> &RESIDUE_A, const util::ShPtr< biol::AABase> &RESIDUE_B
    )
    {
      const int mut_seq_id( m_MutationResidue->GetSeqID());

      // true if neither pointer is defined, just assign them
      if( !RESIDUE_A.IsDefined() && !RESIDUE_B.IsDefined())
      {
        BCL_MessageDbg( "neither residue is defined");
        m_PreviousResidue  = RESIDUE_A;
        m_FollowingResidue = RESIDUE_B;
      }

      // true if both are defined
      else if( RESIDUE_A.IsDefined() && RESIDUE_B.IsDefined())
      {
        BCL_MessageDbg( "sorting residue " + RESIDUE_A->GetIdentification() + " and " + RESIDUE_B->GetIdentification());
        // make sure they have the same chain ids as the mutation residue
        BCL_Assert
        (
          RESIDUE_A->GetChainID() == RESIDUE_B->GetChainID() &&
          RESIDUE_A->GetChainID() == m_MutationResidue->GetChainID(),
          "chain ids don't match"
        );

        // get the seq ids
        const int seq_id_a( RESIDUE_A->GetSeqID()), seq_id_b( RESIDUE_B->GetSeqID());

        // make sure the seq ids are not the same
        BCL_Assert( seq_id_a != seq_id_b, "seq ids are the same");

        // set the residues
        m_PreviousResidue = seq_id_a < seq_id_b ? RESIDUE_A : RESIDUE_B;
        m_FollowingResidue = m_PreviousResidue == RESIDUE_A ? RESIDUE_B : RESIDUE_A;

        BCL_Assert( m_PreviousResidue->GetSeqID() < mut_seq_id, "previous residue is not less in sequence");
        BCL_Assert( m_FollowingResidue->GetSeqID() > mut_seq_id, "following residue is not greater in sequence");

      } //< end if both are defined

      // true if RESIDUE_A is defined only
      else if( RESIDUE_A.IsDefined())
      {
        BCL_MessageDbg( "sorting residue a" + RESIDUE_A->GetIdentification());
        // make sure has same chain id as mutation residue
        BCL_Assert( RESIDUE_A->GetChainID() == m_MutationResidue->GetChainID(), "chain ids don't match");

        const int seq_id_a( RESIDUE_A->GetSeqID());

        m_PreviousResidue = seq_id_a < mut_seq_id ? RESIDUE_A : RESIDUE_B;
        m_FollowingResidue = m_PreviousResidue.IsDefined() ? RESIDUE_B : RESIDUE_A;
      }

      // true if RESIDUE_B is defined only
      else if( RESIDUE_B.IsDefined())
      {
        BCL_MessageDbg( "sorting residue b" + RESIDUE_B->GetIdentification());
        // make sure has same chain id as mutation residue
        BCL_Assert( RESIDUE_B->GetChainID() == m_MutationResidue->GetChainID(), "chain ids don't match");

        const int seq_id_b( RESIDUE_B->GetSeqID());

        m_PreviousResidue = seq_id_b < mut_seq_id ? RESIDUE_B : RESIDUE_A;
        m_FollowingResidue = m_PreviousResidue.IsDefined() ? RESIDUE_A : RESIDUE_B;
      }
    }

  } // namespace fold
} // namespace bcl
