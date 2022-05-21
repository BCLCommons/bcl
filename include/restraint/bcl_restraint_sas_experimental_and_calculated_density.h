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

#ifndef BCL_RESTRAINT_SAS_EXPERIMENTAL_AND_CALCULATED_DENSITY_H_
#define BCL_RESTRAINT_SAS_EXPERIMENTAL_AND_CALCULATED_DENSITY_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_restraint_sas_density_data.h"
#include "bcl_restraint_sas_distance_density_point.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SasExperimentalAndCalculatedDensity
    //! @brief Stores SAS raw data
    //! @details Container class for sas distance data (r, density, error) to be used during folding runs
    //!
    //! @see @link example_restraint_saxs_experimental_and_calculated_density.cpp @endlink
    //! @author putnamdk, mendenjl
    //! @date Sep 04, 2012
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SasExperimentalAndCalculatedDensity :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! storage for experimental sas data < r-value, density, error >
      SasDensityData m_ExperimentalDensity;

      //! max dimension of experimental sas data from gnom
      double m_ExperimentalDmax;

      //! storage for calculated sas data <r-value, density, error >
      SasDensityData m_CalculatedDensity;

      //! max dimension of experimental sas data from gnom
      double m_CalculatedDmax;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SasExperimentalAndCalculatedDensity();

      //!@brief copy constructor
      SasExperimentalAndCalculatedDensity( const SasExperimentalAndCalculatedDensity &rhs);

      //! @brief constructor from given input data
      SasExperimentalAndCalculatedDensity
      (
        const SasDensityData &EXPERIMENTAL,
        const SasDensityData &CALCULATED
      );

      //! @brief Clone function
      //! @return pointer to new SasExperimentalAndCalculatedData
      SasExperimentalAndCalculatedDensity *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the experimental data
      //! @return the experimental data
      const SasDensityData &GetExperimentalDensity() const
      {
        return m_ExperimentalDensity;
      }

      //! @brief returns the calculated data
      //! @return the calculated data
      const SasDensityData &GetCalculatedDensity() const
      {
        return m_CalculatedDensity;
      }
      //! @brief returns the calculated data
      //! @return the calculated data
      SasDensityData &GetCalculatedDensity()
      {
        return m_CalculatedDensity;
      }

    /////////////////
    // operations  //
    /////////////////

      //! @brief scale the calculated and experimental data
      //! @param SCALING_FACTOR - the value to scale all of the Density values by
      void ScaleExperimentalDensity( const double &EXP_SCALING_FACTOR);

      //! @brief scale the calculated and experimental data
      //! @param SCALING_FACTOR - the value to scale all of the Density values by
      void ScaleCalculatedDensity( const double &CAL_SCALING_FACTOR);

      //! @brief align the max peak of the experimental and calculated curves
      void ShiftDensity();

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write to filename in three columns : R-value, density, experimental error
      //! @param filename to write to
      std::ostream &WriteToOstream( std::ostream &OSTREAM) const;

      //! @brief write to filename in four columns : Q-value, experimental intensity, error and calculated intensity
      //! @param filename to write to
      void WriteToFileName( const std::string &FILENAME) const;

      //! @brief write to filename in three columns : R-value, experimental density, and calculated density
      //! @param filename to write to
      void WriteToGnuplotFileName( const std::string &FILENAME) const;

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::osSasExperimentalAndCalculatedDensitytream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class SasExperimentalAndCalculatedDensity

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_SAS_EXPERIMENTAL_AND_CALCULATED_DENSITY_H_
