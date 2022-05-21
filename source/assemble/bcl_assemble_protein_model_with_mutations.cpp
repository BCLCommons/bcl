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
#include "assemble/bcl_assemble_protein_model_with_mutations.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_aa.h"
#include "biol/bcl_biol_atom.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ProteinModelWithMutations::s_Instance
    (
      GetObjectInstances().AddInstance( new ProteinModelWithMutations())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProteinModelWithMutations::ProteinModelWithMutations()
    {
    }

    //! @brief construct from a protein model with cache
    //! @param MODEL protein model of interest
    //! @param REQUIRE_COORDINATES whether to exclude sspred analysis methods from AAs that lac defined coordinates
    //! @param MUTATIONS mutations already made to the protein
    ProteinModelWithMutations::ProteinModelWithMutations
    (
      const ProteinModel &MODEL,
      const bool &REQUIRE_COORDINATES,
      const storage::Vector< biol::Mutation> &MUTATIONS
    ) :
      ProteinModelWithCache( MODEL, REQUIRE_COORDINATES),
      m_CurrentMutations( MUTATIONS)
    {
    }

    //! @brief copy constructor
    //! @return pointer a new ProteinModelWithMutations copied from this model
    ProteinModelWithMutations *ProteinModelWithMutations::Clone() const
    {
      return new ProteinModelWithMutations( *this);
    }

    //! @brief hard copy constructor
    //! @return a ProteinModelWithMutations with chains hard copied from that model
    ProteinModelWithMutations *ProteinModelWithMutations::HardCopy() const
    {
      return new ProteinModelWithMutations( ProteinModel::HardCopy( *this), m_RequireCoordinates, m_CurrentMutations);
    }

    //! @brief empty copy constructor
    //! @return a ProteinModelWithMutations that is empty
    ProteinModelWithMutations *ProteinModelWithMutations::Empty() const
    {
      return new ProteinModelWithMutations;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinModelWithMutations::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief test whether we have the wild type protein
    //! @return true if we have the wild type protein
    bool ProteinModelWithMutations::IsWildType() const
    {
      return m_CurrentMutations.IsEmpty();
    }

    //! @brief test whether this protein only has the given mutation
    //! @return true if this protein only has the given mutation
    bool ProteinModelWithMutations::OnlyHasMutation( const biol::Mutation &MUTATION) const
    {
      return m_CurrentMutations.GetSize() == size_t( 1) && m_CurrentMutations( 0) == MUTATION;
    }

    //! @brief get the mutations already applied to this protein
    const storage::Vector< biol::Mutation> &ProteinModelWithMutations::GetMutations() const
    {
      return m_CurrentMutations;
    }

    //! @brief Apply a mutation to the given protein
    void ProteinModelWithMutations::Mutate( const biol::Mutation &MUTATION)
    {
      ApplyMutation( MUTATION, false);
      m_CurrentMutations.PushBack( MUTATION);
      ResetCache();
    }

    //! @brief revert all mutations to wild-type
    void ProteinModelWithMutations::RevertToWildType()
    {
      for( auto itr( m_CurrentMutations.Begin()), itr_end( m_CurrentMutations.End()); itr != itr_end; ++itr)
      {
        ApplyMutation( *itr, true);
      }
      m_CurrentMutations.Reset();
      ResetCache();
    }

    //! @brief Apply a mutation to the given protein
    //! @param REVERSE whether to reverse the given mutation
    //! @note unlike Mutate, this function doesn't change m_CurrentMutations or the cache
    bool ProteinModelWithMutations::ApplyMutation( const biol::Mutation &MUTATION, const bool &REVERSE)
    {
      size_t n_applied( 0);
      const biol::AAType from_type( REVERSE ? MUTATION.GetMutantType() : MUTATION.GetNativeType());
      const biol::AAType to_type( REVERSE ? MUTATION.GetNativeType() : MUTATION.GetMutantType());

      for( auto itr( GetChains().Begin()), itr_end( GetChains().End()); itr != itr_end; ++itr)
      {
        LocatorAA locator( ( *itr)->GetChainID(), MUTATION.GetResidueNumber(), from_type, false);
        util::SiPtr< biol::AABase> base( locator.Locate( *this));
        if( base.IsDefined())
        {
          base->SetType( to_type);
          MUTATION.SetAA( *base);
        }
        base = locator.Locate( this->m_AAs);
        if( base.IsDefined())
        {
          base->SetType( to_type);
          MUTATION.SetAA( *base);
          ++n_applied;
        }
      }
      if( !n_applied)
      {
        BCL_MessageStd( "Mutation " + MUTATION.ToString() + " Could not be applied to protein");
      }
      else
      {
        this->GetChangeSignal().Emit( *this);
      }
      return n_applied;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief read ProteinModelWithCache from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProteinModelWithMutations::Read( std::istream &ISTREAM)
    {
      ProteinModelWithCache::Read( ISTREAM);
      io::Serialize::Read( m_CurrentMutations, ISTREAM);
      return ISTREAM;
    }

    //! @brief write ProteinModelWithCache to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &ProteinModelWithMutations::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      ProteinModelWithCache::Write( OSTREAM, INDENT) << '\n';
      return io::Serialize::Write( m_CurrentMutations, OSTREAM, INDENT);
    }

  } // namespace assemble
} // namespace bcl
