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
#include "density/bcl_density_fit_protein_minimizers.h"

// includes from bcl - sorted alphabetically
#include "density/bcl_density_fit_protein_minimizer_mc.h"
#include "density/bcl_density_fit_protein_minimizer_powell.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FitProteinMinimizers::FitProteinMinimizers() :
      e_MC(     AddEnum( "MC"    , util::ShPtr< FitProteinMinimizerInterface>( new FitProteinMinimizerMC()))),
      e_Powell( AddEnum( "Powell", util::ShPtr< FitProteinMinimizerInterface>( new FitProteinMinimizerPowell())))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FitProteinMinimizers::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief construct on access function for all FitProteinMinimizers
    //! @return reference to only instances of FitProteinMinimizers
    FitProteinMinimizers &GetFitProteinMinimizers()
    {
      return FitProteinMinimizers::GetEnums();
    }

  } // namespace density

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< density::FitProteinMinimizerInterface>, density::FitProteinMinimizers>;

  } // namespace util
} // namespace bcl
