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

#ifndef BCL_FOLD_MUTATES_H_
#define BCL_FOLD_MUTATES_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "math/bcl_math_mutate_interface.h"
#include "util/bcl_util_enumerate.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Mutates
    //! @brief enumerator for all mutates used in folding
    //! @details Enumerate derived class that allows enumeration of all mutates used in different fold protocols
    //!
    //! @remarks example unnecessary
    //! @author karakam
    //! @date Mar 27, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Mutates :
      public util::Enumerate< util::ShPtr< math::MutateInterface< assemble::ProteinModel> >, Mutates>
    {
      friend class util::Enumerate< util::ShPtr< math::MutateInterface< assemble::ProteinModel> >, Mutates>;

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Mutates();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief adds a mutate to the enumerated mutates
      //! @param SP_MUTATE Shptr to Mutate to add
      EnumType AddMutate( const util::ShPtr< math::MutateInterface< assemble::ProteinModel> > &SP_MUTATE);

      //! @brief adds a mutate to the enumerated mutates
      //! @param MUTATE Mutate to add
      EnumType AddMutate( const math::MutateInterface< assemble::ProteinModel> &MUTATE);

    private:

      //! @brief function for adding a new enum
      //! @param NAME name of the current enum
      //! @param OBJECT object to be enumerated
      EnumType &AddEnum
      (
        const std::string &NAME,
        const util::ShPtr< math::MutateInterface< assemble::ProteinModel> > &OBJECT
      );

    }; // class Mutates

    //! @brief construct on access function for all Mutates
    //! @return reference to only instances of Mutates
    BCL_API
    Mutates &GetMutates();

  } // namespace fold

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< math::MutateInterface< assemble::ProteinModel> >, fold::Mutates>;

  } // namespace util
} // namespace bcl

#endif // BCL_FOLD_MUTATES_H_
