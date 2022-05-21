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
#include "biol/bcl_biol_protein_mutation_set.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ProteinMutationSet::s_Instance
    (
      GetObjectInstances().AddInstance( new ProteinMutationSet())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProteinMutationSet::ProteinMutationSet()
    {
    }

    //! @brief construct from a protein model
    //! @param MODEL protein model of interest
    //! @param REQUIRE_COORDINATES whether to exclude sspred analysis methods from AAs that lac defined coordinates
    ProteinMutationSet::ProteinMutationSet
    (
      const assemble::ProteinModel &MODEL,
      const bool &REQUIRE_COORDINATES,
      const storage::Vector< Mutation> &POSSIBLE_MUTATIONS
    ) :
      descriptor::SequenceInterface< Mutation>(),
      m_MutantModel( assemble::ProteinModel::HardCopy( MODEL), REQUIRE_COORDINATES),
      m_PossibleMutations( POSSIBLE_MUTATIONS)
    {
    }

    //! @brief copy constructor
    //! @param ORIGINAL model with cache to copy
    ProteinMutationSet::ProteinMutationSet( const ProteinMutationSet &ORIGINAL) :
      descriptor::SequenceInterface< Mutation>( ORIGINAL),
      m_MutantModel( assemble::ProteinModel::HardCopy( ORIGINAL.m_MutantModel), ORIGINAL.m_MutantModel.GetRequiresCoordinates()),
      m_PossibleMutations( ORIGINAL.m_PossibleMutations)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer a new ProteinMutationSet copied from this model
    ProteinMutationSet *ProteinMutationSet::Clone() const
    {
      return new ProteinMutationSet( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinMutationSet::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the length of the sequence in question
    //! @return the length of the sequence in question
    size_t ProteinMutationSet::GetSize() const
    {
      return m_PossibleMutations.GetSize();
    }

    //! @brief get the iterator for the sequence
    //! @return the iterator for the sequence
    iterate::Generic< const Mutation> ProteinMutationSet::GetIterator() const
    {
      return iterate::Generic< const Mutation>( m_PossibleMutations.Begin(), m_PossibleMutations.End());
    }

    //! @brief get a non-constant iterator for the sequence
    //! @return the non-constant iterator for the sequence
    iterate::Generic< Mutation> ProteinMutationSet::GetIteratorNonConst()
    {
      return iterate::Generic< Mutation>( m_PossibleMutations.Begin(), m_PossibleMutations.End());
    }

    //! @brief Reset the cache
    void ProteinMutationSet::ResetCache() const
    {
      descriptor::SequenceInterface< Mutation>::ResetCache();
      m_MutantModel.ResetCache();
    }

    //! @brief get a particular mutant protein model
    //! @return the mutated protein model
    const assemble::ProteinModelWithMutations &ProteinMutationSet::GetMutant( const Mutation &MUTATION) const
    {
      if( m_MutantModel.OnlyHasMutation( MUTATION))
      {
        return m_MutantModel;
      }
      m_MutantModel.RevertToWildType();
      m_MutantModel.Mutate( MUTATION);
      return m_MutantModel;
    }

    //! @brief get a particular mutant protein model
    //! @return the mutated protein model
    const assemble::ProteinModelWithMutations &ProteinMutationSet::GetNativeType() const
    {
      if( !m_MutantModel.IsWildType())
      {
        m_MutantModel.RevertToWildType();
      }
      return m_MutantModel;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief assignment operator
    //! @param PROTEIN_MODEL ProteinMutationSet to be copied
    //! @return This model after all members are assigned to values from PROTEIN_MODEL
    ProteinMutationSet &ProteinMutationSet::operator =( const ProteinMutationSet &PROTEIN_MODEL)
    {
      // update members
      m_MutantModel = PROTEIN_MODEL.m_MutantModel;
      m_PossibleMutations = PROTEIN_MODEL.m_PossibleMutations;

      // end
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read ProteinMutationSet from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProteinMutationSet::Read( std::istream &ISTREAM)
    {
      //read data
      io::Serialize::Read( m_MutantModel, ISTREAM);
      io::Serialize::Read( m_PossibleMutations, ISTREAM);

      //return
      return ISTREAM;
    }

    //! @brief write ProteinMutationSet to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &ProteinMutationSet::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write protein model
      io::Serialize::Write( m_MutantModel, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PossibleMutations, OSTREAM, INDENT);
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace biol
} // namespace bcl
