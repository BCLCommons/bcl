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

#ifndef BCL_FOLD_MUTATE_PROTEIN_MODEL_LOOP_DOMAIN_H_
#define BCL_FOLD_MUTATE_PROTEIN_MODEL_LOOP_DOMAIN_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "find/bcl_find_collector_interface.h"
#include "math/bcl_math_mutate_interface.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateProteinModelLoopDomain
    //! @brief is for mutating the Loop Domain of a protein model. The method for locating
    //! the loop domain is behind an interface. Also, the method for mutating the LoopDomain is behind an interface.
    //!
    //! @see @link example_fold_mutate_protein_model_loop_domain.cpp @endlink
    //! @author alexanns, fischea
    //! @date Sep 7, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateProteinModelLoopDomain :
      public math::MutateInterface< assemble::ProteinModel>
    {

    private:

    //////////
    // data //
    //////////

      //! the method for locating a loop domain in a protein model
      util::Implementation< find::CollectorInterface< util::ShPtrList< LoopDomain>, assemble::DomainInterface> > m_DomainCollector;

      //! the method for mutating a loop domain
      util::Implementation< math::MutateInterface< LoopDomain> > m_DomainMutate;

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
      MutateProteinModelLoopDomain();

      //! @brief constructor taking member variable types
      //! @param LOOP_DOMAIN_COLLECTOR the method for collecting loop domains in a protein model
      //! @param DOMAIN_MUTATE the method for mutating a loop domain
      MutateProteinModelLoopDomain
      (
        const util::ShPtr< find::CollectorInterface< util::ShPtrList< LoopDomain>, assemble::DomainInterface> > &LOOP_DOMAIN_COLLECTOR,
        const util::ShPtr< math::MutateInterface< LoopDomain> > &DOMAIN_MUTATE
      );

      //! @brief Clone function
      //! @return pointer to new MutateProteinModelLoopDomainDihedral
      MutateProteinModelLoopDomain *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking an ARGUMENT and returning a mutated object of t_ArgumentType
      //! @param PROTEIN_MODEL Argument of interest
      //! @return MutateResult that results from mutating to the argument
      math::MutateResult< assemble::ProteinModel> operator()
      (
        const assemble::ProteinModel &PROTEIN_MODEL
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    }; // class MutateProteinModelLoopDomain

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_PROTEIN_MODEL_LOOP_DOMAIN_H_ 
