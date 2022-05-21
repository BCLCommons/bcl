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
#include "assemble/bcl_assemble_protein_model_with_cache.h"

// includes from bcl - sorted alphabetically
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
    const util::SiPtr< const util::ObjectInterface> ProteinModelWithCache::s_Instance
    (
      GetObjectInstances().AddInstance( new ProteinModelWithCache())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProteinModelWithCache::ProteinModelWithCache() :
      ProteinModel(),
      descriptor::SequenceInterface< biol::AABase>(),
      m_RequireCoordinates( false)
    {
      GetChangeSignal().Connect( this, &ProteinModelWithCache::UpdateAAPtrs);
    }

    //! @brief construct from a protein model
    //! @param MODEL protein model of interest
    //! @param REQUIRE_COORDINATES whether to exclude sspred analysis methods from AAs that lac defined coordinates
    ProteinModelWithCache::ProteinModelWithCache( const ProteinModel &MODEL, const bool &REQUIRE_COORDINATES) :
      ProteinModel( MODEL),
      descriptor::SequenceInterface< biol::AABase>(),
      m_RequireCoordinates( REQUIRE_COORDINATES)
    {
      UpdateAAPtrs( *this);
      GetChangeSignal().Connect( this, &ProteinModelWithCache::UpdateAAPtrs);
    }

    //! @brief copy constructor
    //! @param ORIGINAL model with cache to copy
    ProteinModelWithCache::ProteinModelWithCache( const ProteinModelWithCache &ORIGINAL) :
      ProteinModel( ORIGINAL),
      descriptor::SequenceInterface< biol::AABase>(),
      m_RequireCoordinates( ORIGINAL.m_RequireCoordinates)
    {
      UpdateAAPtrs( *this);
      GetChangeSignal().Connect( this, &ProteinModelWithCache::UpdateAAPtrs);
    }

    //! @brief virtual copy constructor
    //! @return pointer a new ProteinModelWithCache copied from this model
    ProteinModelWithCache *ProteinModelWithCache::Clone() const
    {
      return new ProteinModelWithCache( *this);
    }

    //! @brief hard copy constructor
    //! @return a ProteinModelWithCache with chains hard copied from that model
    ProteinModelWithCache *ProteinModelWithCache::HardCopy() const
    {
      return new ProteinModelWithCache( ProteinModel::HardCopy( *this), m_RequireCoordinates);
    }

    //! @brief empty copy constructor
    //! @return a ProteinModelWithCache that is empty
    ProteinModelWithCache *ProteinModelWithCache::Empty() const
    {
      return new ProteinModelWithCache();
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinModelWithCache::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the length of the sequence in question
    //! @return the length of the sequence in question
    size_t ProteinModelWithCache::GetSize() const
    {
      return m_AAs.GetSize();
    }

    //! @brief get the iterator for the sequence
    //! @return the iterator for the sequence
    iterate::Generic< const biol::AABase> ProteinModelWithCache::GetIterator() const
    {
      return iterate::Generic< const biol::AABase>( m_AAs.Begin(), m_AAs.End());
    }

    //! @brief get a non-constant iterator for the sequence
    //! @return the non-constant iterator for the sequence
    iterate::Generic< biol::AABase> ProteinModelWithCache::GetIteratorNonConst()
    {
      return iterate::Generic< biol::AABase>( m_AAs.Begin(), m_AAs.End());
    }

    //! @brief get the protein model, represented as a molecule
    //! @note this will only be constructed if when this function is called
    const chemistry::AAFragmentComplete &ProteinModelWithCache::GetChemicalRepresentation() const
    {
      if( !m_FragmentComplete.IsDefined())
      {
        m_FragmentComplete = util::ShPtr< chemistry::AAFragmentComplete>( new chemistry::AAFragmentComplete( m_AAs, true));
      }
      return *m_FragmentComplete;
    }

    //! @brief Reset the cache
    void ProteinModelWithCache::ResetCache() const
    {
      descriptor::SequenceInterface< biol::AABase>::ResetCache();
      m_FragmentComplete = util::ShPtr< chemistry::AAFragmentComplete>();
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief assignment operator
    //! @param PROTEIN_MODEL ProteinModelWithCache to be copied
    //! @return This model after all members are assigned to values from PROTEIN_MODEL
    ProteinModelWithCache &ProteinModelWithCache::operator =( const ProteinModelWithCache &PROTEIN_MODEL)
    {
      // update members
      m_RequireCoordinates = PROTEIN_MODEL.m_RequireCoordinates;
      ProteinModel::operator =( PROTEIN_MODEL);

      // end
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read ProteinModelWithCache from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProteinModelWithCache::Read( std::istream &ISTREAM)
    {
      //read data
      ProteinModel::Read( ISTREAM);

      //return
      return ISTREAM;
    }

    //! @brief write ProteinModelWithCache to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &ProteinModelWithCache::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write protein model
      return ProteinModel::Write( OSTREAM, INDENT);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief update amino acid pointers based on the new model
    void ProteinModelWithCache::UpdateAAPtrs( const ProteinModel &MODEL)
    {
      util::SiPtrVector< biol::AABase> new_aas;
      new_aas.AllocateMemory( m_AAs.GetSize());

      m_AAs.Reset();
      //loop over all SSElements
      for
      (
        util::ShPtrVector< Chain>::iterator chain_itr( GetChains().Begin()), chain_itr_end( GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        new_aas.Append( ( *chain_itr)->GetSequence()->GetMembers());
      }

      // loop over AAs, remove secondary structure predictions that are invalid if the AA does not have defined coords
      if( m_RequireCoordinates)
      {
         // SS/TM analytic methods should require defined coordinates, so remove the method should have required them
         // This is kind of a work-around for the fact that 3D structural analysis methods can sometimes report values
         // even when the coordinates are missing; and even if they don't, the BCL may have inserted some default
         // prediction value for the method
        for
        (
          util::SiPtrVector< biol::AABase>::iterator itr( new_aas.Begin()), itr_end( new_aas.End());
          itr != itr_end;
          ++itr
        )
        {
          if( !( *itr)->HasDefinedCoordinates())
          {
            ( *itr)->RemoveStructureBasedSSTMInfo();
          }
        }
      }
      // Skip residues with undefined type at the beginning
      util::SiPtrVector< biol::AABase>::const_iterator itr( new_aas.Begin()), itr_end( new_aas.End());
      for( ; itr != itr_end && ( !( *itr)->GetType().IsDefined() || !( *itr)->GetType()->IsNaturalAminoAcid()); ++itr)
      {
      }
      for( ; itr != itr_end; ++itr)
      {
        if( ( *itr)->GetType().IsDefined())
        {
          m_AAs.PushBack( *itr);
        }
      }
      // prune residues at the end of the sequence with undefined residue types. These are typically expression tags
      // or, in the case of membrane proteins, dummy residues inserted by OPM
      while( !m_AAs.IsEmpty() && !m_AAs.LastElement()->GetType()->IsNaturalAminoAcid())
      {
        m_AAs.PopBack();
      }

      ResetCache();
    }

  } // namespace assemble
} // namespace bcl
