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

#ifndef BCL_CHEMISTRY_CONFORMATION_COMPARISON_INTERFACE_H_
#define BCL_CHEMISTRY_CONFORMATION_COMPARISON_INTERFACE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_interface.h"
#include "util/bcl_util_binary_function_interface.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_serializable_interface.h"
#include "util/bcl_util_sh_ptr.h"

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConformationComparisonInterface
    //! @brief This class is designed to be used for determining and comparing 3D structures for molecules.
    //!
    //! @see @link example_chemistry_conformation_comparison_interface.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date Mar 05, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConformationComparisonInterface :
      public virtual util::SerializableInterface,
      public util::BinaryFunctionInterface< ConformationInterface, ConformationInterface, double>
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! virtual copy constructor
      virtual ConformationComparisonInterface *Clone() const = 0;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief prepare the class for comparing conformations in the given ensemble
      //! @param ENSEMBLE the ensemble to prepare to compare
      virtual void PrepareEnsemble( const FragmentEnsemble &ENSEMBLE) const;

      //! @brief prepare the class for comparing conformations in the given ensemble
      //! @param ENSEMBLE the ensemble to prepare to compare
      virtual void Prepare( const ConformationInterface &ENSEMBLE) const
      {
      }

      //! @brief get a bool that indicates whether to ignore atom / bond types in ConformationsAreComparable, below
      //! @note this is necessary when comparing parts of molecules, rather than complete ligands
      static bool &GetIgnoreAtomAndBondTypesWhenDeterminingComparability();

      //! @brief Get a flag to turn on/off strict atom type checking
      static const util::ShPtr< command::FlagInterface> &GetDisableStrictAtomBondTypeCheckingFlag();

      //! @brief determine whether the conformations represent identical, aligned, constitutions
      //! @param MOLECULE_A, MOLECULE_B the conformations to check for being identical and aligned
      //! @return true iff the conformations represent identical, aligned, constitutions
      static bool ConformationsAreComparable
      (
        const ConformationInterface &MOLECULE_A,
        const ConformationInterface &MOLECULE_B
      );

    };

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_CONFORMATION_COMPARISON_INTERFACE_H_
