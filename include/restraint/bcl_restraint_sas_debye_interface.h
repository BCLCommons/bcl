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

#ifndef BCL_RESTRAINT_SAS_DEBYE_INTERFACE_H_
#define BCL_RESTRAINT_SAS_DEBYE_INTERFACE_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_restraint_sas_experimental_and_calculated_data.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SasDebyeInterface
    //! @brief Calculates experimental and predicted values for a protein model
    //!
    //! @author mendenjl
    //! @date Oct 31, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SasDebyeInterface :
      public math::FunctionInterfaceSerializable< assemble::ProteinModel, SasExperimentalAndCalculatedData>
    {

    private:

    //////////
    // data //
    //////////

      //! Experimental data
      util::ShPtr< SasScatteringData> m_ExpData;

      //! Shared pointer to proposed reduced experimental data set
      util::ShPtr< storage::Vector< SasScatteringPoint> > m_ReducedExpData;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone the SasDebyeInterface
      //! @return pointer to new SasDebyeInterface
      virtual SasDebyeInterface *Clone() const = 0;

      //! @brief get the experimental saxs profile
      //! @return pointer to the experimental saxs profile
      const util::ShPtr< SasScatteringData> &GetExperimentalData() const
      {
        return m_ExpData;
      }

      //! @brief set the experimental saxs profile
      //! @param EXPERIMENTAL_DATA the new experimental data
      void SetExperimentalData( const util::ShPtr< SasScatteringData> &EXPERIMENTAL_DATA)
      {
        m_ExpData = EXPERIMENTAL_DATA;
      }

      //! @brief function to set the Data for the reduced representation of SAXS experimental data
      //! @param REDUCED_EXPERIMENTAL_DATA - shared pointer to data
      void SetReducedExperimentalData
      (
        const util::ShPtr< storage::Vector< SasScatteringPoint> > &REDUCED_EXPERIMENTAL_DATA
      )
      {
        m_ReducedExpData = REDUCED_EXPERIMENTAL_DATA;
      }

      //! @brief function to get the Data for the reduced representation of SAXS experimental data
      //! @return pointer to the experimental saxs profile
      const util::ShPtr< storage::Vector< SasScatteringPoint> > &GetReducedExperimentalData() const
      {
        return m_ReducedExpData;
      }

    }; // class SasDebyeInterface

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_SAS_DEBYE_INTERFACE_H_
