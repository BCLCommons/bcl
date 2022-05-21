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
#include "assemble/bcl_assemble_protein_model.h"
#include "fold/bcl_fold_mutate_protein_model_fix_loop_closure_wrapper.h"
#include "fold/bcl_fold_mutate_protein_model_sse_pair_clash.h"
#include "fold/bcl_fold_mutate_protein_model_sse_pair_fix_loop_closure.h"
#include "fold/bcl_fold_mutates.h"
#include "math/bcl_math_mutate_combine.h"
#include "util/bcl_util_enumerate.hpp"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Mutates::Mutates()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Mutates::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief adds a mutate to the enumerated mutates
    //! @param SP_MUTATE ShPtr to Mutate to add
    Mutates::EnumType Mutates::AddMutate( const util::ShPtr< math::MutateInterface< assemble::ProteinModel> > &SP_MUTATE)
    {
      return AddEnum( SP_MUTATE->GetScheme(), SP_MUTATE);
    }

    //! @brief adds a mutate to the enumerated mutates
    //! @param MUTATE Mutate to add
    Mutates::EnumType Mutates::AddMutate( const math::MutateInterface< assemble::ProteinModel> &MUTATE)
    {
      // construct ShPtr
      util::ShPtr< math::MutateInterface< assemble::ProteinModel> > sp_mutate( MUTATE.Clone());

      return AddEnum( sp_mutate->GetScheme(), sp_mutate);
    }

    //! @brief function for adding a new enum
    //! @param NAME name of the current enum
    //! @param OBJECT object to be enumerated
    Mutates::EnumType &Mutates::AddEnum
    (
      const std::string &NAME,
      const util::ShPtr< math::MutateInterface< assemble::ProteinModel> > &OBJECT
    )
    {
      util::ShPtr< math::MutateInterface< assemble::ProteinModel> > combined
      (
        new MutateProteinModelFixLoopClosureWrapper( *OBJECT, OBJECT->GetScheme())
      );
      return util::Enumerate< util::ShPtr< math::MutateInterface< assemble::ProteinModel> >, Mutates>::AddEnum( NAME, combined);
    }

    //! @brief construct on access function for all Mutates
    //! @return reference to only instances of Mutates
    Mutates &GetMutates()
    {
      return Mutates::GetEnums();
    }

  } // namespace fold

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< math::MutateInterface< assemble::ProteinModel> >, fold::Mutates>;

  } // namespace util
} // namespace bcl
