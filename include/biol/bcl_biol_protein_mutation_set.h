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

#ifndef BCL_BIOL_PROTEIN_MUTATION_SET_H_
#define BCL_BIOL_PROTEIN_MUTATION_SET_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_biol_mutation.h"
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "assemble/bcl_assemble_protein_model_with_mutations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinMutationSet
    //! @brief class that represents a biological protein structure as a model within the bcl.
    //! @details An adaptor for a protein model to a descriptor::SequenceInterface
    //!
    //! @see @link example_biol_protein_mutation_set.cpp @endlink
    //! @author mendenjl
    //! @date Jan 17, 2019
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinMutationSet :
      public descriptor::SequenceInterface< Mutation>
    {

    protected:

    //////////
    // data //
    //////////

      //! the protein model; changeable to allow mutations to be created on-the-fly
      mutable assemble::ProteinModelWithMutations m_MutantModel;

      //! mutations that can be applied to the protein model
      storage::Vector< Mutation> m_PossibleMutations;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ProteinMutationSet();

      //! @brief construct from a protein model
      //! @param MODEL protein model of interest
      //! @param REQUIRE_COORDINATES whether to exclude sspred analysis methods from AAs that lac defined coordinates
      ProteinMutationSet
      (
        const assemble::ProteinModel &MODEL,
        const bool &REQUIRE_COORDINATES,
        const storage::Vector< Mutation> &MUTATIONS = storage::Vector< Mutation>()
      );

      //! @brief copy constructor
      //! @param ORIGINAL model with cache to copy
      ProteinMutationSet( const ProteinMutationSet &ORIGINAL);

      //! @brief copy constructor
      //! @return pointer a new ProteinMutationSet copied from this model
      virtual ProteinMutationSet *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the length of the sequence in question
      //! @return the length of the sequence in question
      size_t GetSize() const;

      //! @brief get the iterator for the sequence
      //! @return the iterator for the sequence
      iterate::Generic< const Mutation> GetIterator() const;

      //! @brief get a non-constant iterator for the sequence
      //! @return the non-constant iterator for the sequence
      iterate::Generic< Mutation> GetIteratorNonConst();

      //! @brief get a particular mutant protein model
      //! @return the mutated protein model
      const assemble::ProteinModelWithMutations &GetMutant( const Mutation &MUTATION) const;

      //! @brief get a particular mutant protein model
      //! @return the mutated protein model
      const assemble::ProteinModelWithMutations &GetNativeType() const;

      //! @brief Reset the cache
      virtual void ResetCache() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief assignment operator
      //! @param PROTEIN_MODEL ProteinMutationSet to be copied
      //! @return This model after all members are assigned to values from PROTEIN_MODEL
      virtual ProteinMutationSet &operator =( const ProteinMutationSet &PROTEIN_MODEL);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read ProteinMutationSet from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write ProteinMutationSet to std::ostream
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class ProteinMutationSet

  } // namespace biol
} // namespace bcl

#endif //BCL_BIOL_PROTEIN_MUTATION_SET_H_
