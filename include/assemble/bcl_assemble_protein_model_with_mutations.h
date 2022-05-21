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

#ifndef BCL_ASSEMBLE_PROTEIN_MODEL_WITH_MUTATIONS_H_
#define BCL_ASSEMBLE_PROTEIN_MODEL_WITH_MUTATIONS_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_assemble_protein_model_with_cache.h"
#include "biol/bcl_biol_mutation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinModelWithMutations
    //! @brief class that represents a biological protein structure as a model within the bcl.
    //! @details An adaptor for a protein model to a descriptor::SequenceInterface
    //!
    //! @see @link example_assemble_protein_model_with_mutations.cpp @endlink
    //! @author mendenjl
    //! @date Jan 17, 2019
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinModelWithMutations :
      public ProteinModelWithCache
    {

    protected:

    //////////
    // data //
    //////////

      storage::Vector< biol::Mutation> m_CurrentMutations;

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
      ProteinModelWithMutations();

      //! @brief construct from a protein model with cache
      //! @param MODEL protein model of interest
      //! @param REQUIRE_COORDINATES whether to exclude sspred analysis methods from AAs that lac defined coordinates
      //! @param MUTATIONS mutations already made to the protein
      ProteinModelWithMutations
      (
        const ProteinModel &MODEL,
        const bool &REQUIRE_COORDINATES,
        const storage::Vector< biol::Mutation> &MUTATIONS = storage::Vector< biol::Mutation>()
      );

      //! @brief copy constructor
      //! @return pointer a new ProteinModelWithMutations copied from this model
      virtual ProteinModelWithMutations *Clone() const;

      //! @brief hard copy constructor
      //! @return a ProteinModelWithMutations with chains hard copied from that model
      virtual ProteinModelWithMutations *HardCopy() const;

      //! @brief empty copy constructor
      //! @return a ProteinModelWithMutations that is empty
      virtual ProteinModelWithMutations *Empty() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief test whether we have the wild type protein
      //! @return true if we have the wild type protein
      bool IsWildType() const;

      //! @brief test whether this protein only has the given mutation
      //! @return true if this protein only has the given mutation
      bool OnlyHasMutation( const biol::Mutation &MUTATION) const;

      //! @brief get the mutations already applied to this protein
      const storage::Vector< biol::Mutation> &GetMutations() const;

      //! @brief Apply a mutation to the given protein
      void Mutate( const biol::Mutation &MUTATION);

      //! @brief revert all mutations to wild-type
      void RevertToWildType();

    ///////////////
    // operators //
    ///////////////

    protected:

      //! @brief read ProteinModelWithCache from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write ProteinModelWithCache to std::ostream
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! @brief Apply a mutation to the given protein
      //! @param REVERSE whether to reverse the given mutation
      //! @note unlike Mutate, this function doesn't change m_CurrentMutations or the cache
      bool ApplyMutation( const biol::Mutation &MUTATION, const bool &REVERSE);

    }; // class ProteinModelWithMutations

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_PROTEIN_MODEL_WITH_MUTATIONS_H_
