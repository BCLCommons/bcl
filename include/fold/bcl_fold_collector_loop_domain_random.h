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

#ifndef BCL_FOLD_COLLECTOR_LOOP_DOMAIN_RANDOM_H_
#define BCL_FOLD_COLLECTOR_LOOP_DOMAIN_RANDOM_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_fold_mutation_residue.h"
#include "find/bcl_find_collector_interface.h"
#include "random/bcl_random_distribution_interface.h"
#include "random/bcl_random_uniform_distribution.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CollectorLoopDomainRandom
    //! @brief is for randomly collecting a desired number of mutation residues from a loop domain.
    //!
    //! @see @link example_fold_collector_loop_domain_random.cpp @endlink
    //! @author alexanns
    //! @date Sep 6, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CollectorLoopDomainRandom :
      public find::CollectorInterface< storage::List< MutationResidue>, LoopDomain>
    {

    private:

    //////////
    // data //
    //////////

      //! the number of residues that should be collected for mutation
      size_t m_NumberMutationResiduesToCollect;

      //! the random number generator used to select residues for mutation
      const random::DistributionInterface &m_RandomNumberGenerator;

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
      CollectorLoopDomainRandom();

      //! @brief constructor taking parameter
      //! @param NUMBER_MUTATION_RESIDUES_TO_COLLECT the number of residues that should be collected for mutation
      //! @param RANDOM_NUMBER_GENERATOR the random number generator used to select residues for mutation
      CollectorLoopDomainRandom
      (
        const size_t NUMBER_MUTATION_RESIDUES_TO_COLLECT,
        const random::DistributionInterface &RANDOM_NUMBER_GENERATOR = random::GetGlobalRandom()
      );

      //! @brief Clone function
      //! @return pointer to new CollectorLoopDomainRandom
      CollectorLoopDomainRandom *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Collect taking an ARGUMENT and returning a list of MutationResidues
      //! @param LOOP_DOMAIN LoopDomain of interest
      //! @return storage::List< MutationResidue> that results from collecting from "LOOP_DOMAIN"
      storage::List< MutationResidue> Collect( const LoopDomain &LOOP_DOMAIN) const;

    //////////////////////
    // input and output //
    //////////////////////

    }; // class CollectorLoopDomainRandom

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_COLLECTOR_LOOP_DOMAIN_RANDOM_H_ 
