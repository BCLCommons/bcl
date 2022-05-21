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

#ifndef BCL_BIOL_AA_CLASSES_H_
#define BCL_BIOL_AA_CLASSES_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_biol_aa_base.h"
#include "util/bcl_util_enumerate.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AAClasses
    //! @brief Enumerator class for different aa classes in bcl library, designed to prevent templates usage
    //! @details By enumerating different aa classes, upper classes such as AASequence, Chain, etc. can be designed without
    //! templates. All 3 AAClasses AA, AACaCb, AABackBone are derived from AABase and these upper classes only know
    //! about AABase interface. This enumerator allows these classes to decide deduce what kind of AABase derived class
    //! can be found behind AABase interfaces
    //!
    //! @see @link example_biol_aa_classes.cpp @endlink
    //! @author karakam
    //! @date 10.10.2005
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AAClasses :
      public util::Enumerate< util::ShPtr< AABase>, AAClasses>
    {
      friend class util::Enumerate< util::ShPtr< AABase>, AAClasses>;
    public:

    //////////
    // data //
    //////////

      // declare all aminoacids classes
      const AAClass e_AA;         //!< AA
      const AAClass e_AACaCb;     //!< AACaCb
      const AAClass e_AABackBone; //!< AABackBone
      const AAClass e_AAComplete; //!< AAComplete

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct all AAClasses
      AAClasses();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    }; // class AATypes

    //! @brief construct on access function for all AAClasses
    //! @return reference to only instances of AAClasses
    BCL_API
    const AAClasses &GetAAClasses();

  } // namespace biol

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< biol::AABase>, biol::AAClasses>;

  } // namespace util
} // namespace bcl

#endif //BCL_BIOL_AA_CLASSES_H_
