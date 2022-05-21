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

#ifndef BCL_BIOL_ENVIRONMENT_TYPES_H_
#define BCL_BIOL_ENVIRONMENT_TYPES_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_biol_environment_type_data.h"
#include "util/bcl_util_enumerate.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EnvironmentTypes
    //! @brief This class enumerates environment types relevant to proteins.
    //! @details Every enumerator has a distinct EnvironmentTypeData behind it. Various functions enable access to
    //! reduced types (without the gaps) and the usage of two letter codes
    //!
    //! @see @link example_biol_environment_types.cpp @endlink
    //! @author karakam, woetzen
    //! @date 08/27/09
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API EnvironmentTypes :
      public util::Enumerate< EnvironmentTypeData, EnvironmentTypes>
    {
      friend class util::Enumerate< EnvironmentTypeData, EnvironmentTypes>;
    public:

    //////////
    // data //
    //////////

      // declare all environment types
      const EnvironmentType e_MembraneCore;          //!< Core region
      const EnvironmentType e_GapCoreTransition;     //!< Gap region between core and transition
      const EnvironmentType e_Transition;            //!< Transition region
      const EnvironmentType e_GapTransitionSolution; //!< Gap region between transition and solution
      const EnvironmentType e_Solution;              //!< Solution region
      const EnvironmentType e_MembraneInside;        //!< Membrane region, moving towards cytosol
      const EnvironmentType e_MembraneOutside;       //!< Membrane region, moving out of the cytosol
      const EnvironmentType e_SolutionInside;        //!< Solution region, inside the membrane
      const EnvironmentType e_SolutionOutside;       //!< Solution region, outside the membrane

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct all EnvironmentTypes
      EnvironmentTypes();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief function to deduce EnvironmentType from two letter code
      //! @param TWO_LETTER_CODE two letter code descriptor for environment type of interest
      //! @return EnvironmentType specified by the given TWO_LETTER_CODE
      const EnvironmentType &EnvironmentTypeFromTwoLetterCode( const std::string &TWO_LETTER_CODE) const;

      //! @brief returns number of reduced types
      //! @return number of reduced types
      size_t GetNumberReducedTypes() const;

      //! @brief returns a vector that contains the 3 reduced types, e_MembraneCore, e_Transition and e_Solution
      //! @return a vector that contains the 3 reduced types, e_MembraneCore, e_Transition and e_Solution
      const storage::Vector< EnvironmentType> &GetReducedTypes() const;

    }; // class EnvironmentTypes

    //! @brief construct on access function for all EnvironmentTypes
    //! @return reference to only instance of EnvironmentTypes enum
    BCL_API const EnvironmentTypes &GetEnvironmentTypes();

  } // namespace biol

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< biol::EnvironmentTypeData, biol::EnvironmentTypes>;

  } // namespace util
} // namespace bcl

#endif //BCL_BIOL_ENVIRONMENT_TYPES_H_
