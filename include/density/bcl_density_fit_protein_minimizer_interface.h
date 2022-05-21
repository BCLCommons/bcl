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

#ifndef BCL_DENSITY_FIT_PROTEIN_MINIMIZER_INTERFACE_H_
#define BCL_DENSITY_FIT_PROTEIN_MINIMIZER_INTERFACE_H_

// include the namespace header
#include "bcl_density.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "function/bcl_function_binary_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FitProteinMinimizerInterface
    //! @brief TODO: add a brief comment
    //! @details TODO: add an detailed description to this class
    //!
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date Mar 11, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FitProteinMinimizerInterface :
      public function::BinaryInterface< const assemble::ProteinModel, const Map, assemble::ProteinModel>
    {

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief set the resolution of the density map
      //! @param RESOLUTION density map and simulation resolution
      virtual void SetResolution( const double RESOLUTION) = 0;

      //! @brief set max translation and rotation
      //! @param MAX_TRANSLATION max translation in any direction for a single iteration
      //! @param MAX_ROTATION max rotation in radians in any direction for a single iteration
      virtual void SetMaxTranslationAndRotation( const double MAX_TRANSLATION, const double MAX_ROTATION) = 0;

      //! @brief set the max number of iterations for minimization
      //! @param MAX_NUMBER_ITERATIONS maximum number of iterations for minimization
      virtual void SetMaxIterations( const size_t MAX_NUMBER_ITERATIONS) = 0;

      //! @brief set protein agreement measure to be used
      //! @param AGREEMENT protein agreement enumerator
      virtual void SetProteinAgreement( const ProteinAgreement &AGREEMENT) = 0;

      //! @brief simulator to use
      //! @param DENSITY_SIMULATOR simulator enumerator
      virtual void SetSimulator( const Simulator &DENSITY_SIMULATOR) = 0;

    }; // class FitProteinMinimizerInterface

  } // namespace density
} // namespace bcl

#endif // BCL_DENSITY_FIT_PROTEIN_MINIMIZER_INTERFACE_H_ 
