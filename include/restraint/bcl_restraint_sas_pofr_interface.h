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

#ifndef BCL_RESTRAINT_SAS_POFR_INTERFACE_H_
#define BCL_RESTRAINT_SAS_POFR_INTERFACE_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_restraint_sas_experimental_and_calculated_density.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SasPofRInterface
    //! @brief Calculates protein model pair-wise distances and bins them in the same range as the experimental data
    //!
    //! @author putnamdk
    //! @date Jun 19, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SasPofRInterface :
      public math::FunctionInterfaceSerializable< assemble::ProteinModel, SasExperimentalAndCalculatedDensity>
    {

    private:

    //////////
    // data //
    //////////

      //! Experimental data
        util::Implementation< SasDensityData> m_ExpDensity;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone the SasDebyeInterface
      //! @return pointer to new SasDebyeInterface
      virtual SasPofRInterface *Clone() const = 0;

      //! @brief get the experimental saxs profile
      //! @return pointer to the experimental saxs profile
      const util::Implementation< SasDensityData> &GetExperimentalDensity() const
      {
        return m_ExpDensity;
      }

      //! @brief set the experimental saxs profile
      //! @param EXPERIMENTAL_DATA the new experimental data
      void SetExperimentalDensity( const util::Implementation< SasDensityData> &EXPERIMENTAL_DENSITY)
      {
        m_ExpDensity = EXPERIMENTAL_DENSITY;
      }

    }; // class SasPofRInterface

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_SAS_POFR_INTERFACE_H_
