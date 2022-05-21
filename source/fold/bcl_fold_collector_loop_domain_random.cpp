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
#include "fold/bcl_fold_collector_loop_domain_random.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_collector_loop_domain_all_non_rigid.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> CollectorLoopDomainRandom::s_Instance
    (
      util::Enumerated< find::CollectorInterface< storage::List< MutationResidue>, LoopDomain> >::AddInstance
      (
        new CollectorLoopDomainRandom()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    CollectorLoopDomainRandom::CollectorLoopDomainRandom() :
      m_NumberMutationResiduesToCollect(),
      m_RandomNumberGenerator( random::GetGlobalRandom())
    {
    }

    //! @brief constructor taking parameter
    //! @param NUMBER_MUTATION_RESIDUES_TO_COLLECT the number of residues that should be collected for mutation
    //! @param RANDOM_NUMBER_GENERATOR the random number generator used to select residues for mutation
    CollectorLoopDomainRandom::CollectorLoopDomainRandom
    (
      const size_t NUMBER_MUTATION_RESIDUES_TO_COLLECT,
      const random::DistributionInterface &RANDOM_NUMBER_GENERATOR
    ) :
      m_NumberMutationResiduesToCollect( NUMBER_MUTATION_RESIDUES_TO_COLLECT),
      m_RandomNumberGenerator( RANDOM_NUMBER_GENERATOR)
    {
    }

    //! @brief Clone function
    //! @return pointer to new CollectorLoopDomainRandom
    CollectorLoopDomainRandom *CollectorLoopDomainRandom::Clone() const
    {
      return new CollectorLoopDomainRandom( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &CollectorLoopDomainRandom::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &CollectorLoopDomainRandom::GetAlias() const
    {
      static const std::string s_alias( "CollectorLoopDomainRandom");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer CollectorLoopDomainRandom::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Collects random loop domains.");
      serializer.AddInitializer
      (
        "number residues",
        "number of residues to collect for mutation",
        io::Serialization::GetAgent( &m_NumberMutationResiduesToCollect)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief virtual Collect taking an ARGUMENT and returning a list of MutationResidues
    //! @param LOOP_DOMAIN LoopDomain of interest
    //! @return storage::List< MutationResidue> that results from collecting from "LOOP_DOMAIN"
    storage::List< MutationResidue> CollectorLoopDomainRandom::Collect( const LoopDomain &LOOP_DOMAIN) const
    {
      // collect all the possible mutation residues
      const storage::List< MutationResidue> all_mutation_residues
      (
        CollectorLoopDomainAllNonRigid().Collect( LOOP_DOMAIN)
      );

      // make siptr vector which will hold all possible mutation residues
      util::SiPtrList< const MutationResidue> all_mutation_residues_vector
      (
        all_mutation_residues.Begin(), all_mutation_residues.End()
      );

      // create list to hold the final list of selected mutation residues
      storage::List< MutationResidue> selected_mutation_residues;

      // fill "selected_mutation_residues" with MutationResidues randomly selected from "all_mutation_residues_vector"
      while( !all_mutation_residues.IsEmpty() && selected_mutation_residues.GetSize() < m_NumberMutationResiduesToCollect)
      {
        // get a random iterator
        util::SiPtrList< const MutationResidue>::iterator selected_itr
        (
          m_RandomNumberGenerator.Iterator
          (
            all_mutation_residues_vector.Begin(),
            all_mutation_residues_vector.End(),
            all_mutation_residues_vector.GetSize()
          )
        );

        // add the mutation residue
        selected_mutation_residues.PushBack( **selected_itr);
        all_mutation_residues_vector.Remove( selected_itr);
      }

      // return the list of selected mutation residues
      return selected_mutation_residues;
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace fold
} // namespace bcl
