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

#ifndef BCL_DENSITY_FIT_PROTEIN_MINIMIZERS_H_
#define BCL_DENSITY_FIT_PROTEIN_MINIMIZERS_H_

// include the namespace header
#include "bcl_density.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_density_fit_protein_minimizer_interface.h"
#include "util/bcl_util_enumerate.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FitProteinMinimizers
    //! @brief TODO: add a brief comment
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_density_fit_protein_minimizers.cpp @endlink
    //! @author woetzen
    //! @date Mar 11, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FitProteinMinimizers :
      public util::Enumerate< util::ShPtr< FitProteinMinimizerInterface>, FitProteinMinimizers>
    {
      friend class util::Enumerate< util::ShPtr< FitProteinMinimizerInterface>, FitProteinMinimizers>;
    public:

    //////////
    // data //
    //////////

      FitProteinMinimizer e_MC;     //!< Monte carlo minimizer
      FitProteinMinimizer e_Powell; //!< powell pseudo gradient

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FitProteinMinimizers();

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

    }; // class FitProteinMinimizers

    //! @brief construct on access function for all FitProteinMinimizers
    //! @return reference to only instances of FitProteinMinimizers
    BCL_API
    FitProteinMinimizers &GetFitProteinMinimizers();

  } // namespace density

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< density::FitProteinMinimizerInterface>, density::FitProteinMinimizers>;

  } // namespace util
} // namespace bcl

#endif // BCL_DENSITY_FIT_PROTEIN_MINIMIZERS_H_ 
