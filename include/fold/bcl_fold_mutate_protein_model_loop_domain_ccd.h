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

#ifndef BCL_FOLD_MUTATE_PROTEIN_MODEL_LOOP_DOMAIN_CCD_H_
#define BCL_FOLD_MUTATE_PROTEIN_MODEL_LOOP_DOMAIN_CCD_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_cyclic_coordinate_descent.h"
#include "find/bcl_find_collector_interface.h"
#include "find/bcl_find_locator_interface.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_range.h"
#include "random/bcl_random_distribution_interface.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateProteinModelLoopDomainCCD
    //! @brief is for mutating a protein model such that a residue in a loop domain
    //! has its dihedral angles changed in order to try to close a gap in the sequence.
    //!
    //! @see @link example_fold_mutate_protein_model_loop_domain_ccd.cpp @endlink
    //! @author alexanns, fischea
    //! @date Sep 6, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateProteinModelLoopDomainCCD :
      public math::MutateInterface< assemble::ProteinModel>
    {

    private:

    //////////
    // data //
    //////////

      //! the collector object to find possible loop domains to mutate
      util::ShPtr< find::CollectorInterface< util::ShPtrList< LoopDomain>, assemble::DomainInterface> > m_DomainCollector;

      //! the locator object to locate the loop domain to mutate
      util::ShPtr< find::LocatorInterface< util::ShPtr< LoopDomain>, assemble::DomainInterface> > m_DomainLocator;

      //! the random number generator that should be used
      const random::DistributionInterface &m_RandomNumberGenerator;

      //! range for a random fraction to multiply with the suggested rotation angle
      math::Range< double> m_RandomFraction;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutateProteinModelLoopDomainCCD();

      //! @brief constructor taking member parameters
      //! @param DOMAIN_COLLECTOR the collector to collect possible loop domains for mutating
      //! @param RANDOM_NUMBER_GENERATOR the random number generator that should be used
      //! @param RANDOM_FRACTION_RANGE range a random fraction is drawn from and multiplied with suggested rotation
      MutateProteinModelLoopDomainCCD
      (
        const util::ShPtr< find::CollectorInterface< util::ShPtrList< LoopDomain>, assemble::DomainInterface> > &DOMAIN_COLLECTOR,
        const random::DistributionInterface &RANDOM_NUMBER_GENERATOR,
        const math::Range< double> &RANDOM_FRACTION_RANGE
      );

      //! @brief constructor taking member parameters
      //! @param DOMAIN_LOCATOR the locator to locate loop domain for mutating
      //! @param RANDOM_NUMBER_GENERATOR the random number generator that should be used
      //! @param RANDOM_FRACTION_RANGE range a random fraction is drawn from and multiplied with suggested rotation
      MutateProteinModelLoopDomainCCD
      (
        const util::ShPtr< find::LocatorInterface< util::ShPtr< LoopDomain>, assemble::DomainInterface> > &DOMAIN_LOCATOR,
        const random::DistributionInterface &RANDOM_NUMBER_GENERATOR,
        const math::Range< double> &RANDOM_FRACTION_RANGE
      );

      //! @brief Clone function
      //! @return pointer to new MutateProteinModelLoopDomainCCD
      MutateProteinModelLoopDomainCCD *Clone() const;

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

    ////////////////
    // operations //
    ////////////////

      //! @brief Set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read of this class
      //! @param ERROR_STREAM stream with which to write errors
      //! @return return code indicating success or failure
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    ///////////////
    // operators //
    ///////////////

      //! @brief operator taking an ARGUMENT and returning a mutated object of t_ArgumentType
      //! @param PROTEIN_MODEL Argument of interest
      //! @return MutateResult that results from mutating to the argument
      math::MutateResult< assemble::ProteinModel> operator()
      (
        const assemble::ProteinModel &PROTEIN_MODEL
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    }; // class MutateProteinModelLoopDomainCCD

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_PROTEIN_MODEL_LOOP_DOMAIN_CCD_H_
