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

#ifndef BCL_RESTRAINT_SAS_EXPERIMENTAL_AND_CALCULATED_DATA_H_
#define BCL_RESTRAINT_SAS_EXPERIMENTAL_AND_CALCULATED_DATA_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_restraint_sas_scattering_data.h"
#include "bcl_restraint_sas_scattering_point.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SasExperimentalAndCalculatedData
    //! @brief Stores SAXS raw data
    //! @details Container class for saxs data (q, intensity, error) to be used during folding runs
    //!
    //! @see @link example_restraint_saxs_experimental_and_calculated_data.cpp @endlink
    //! @author putnamdk, mendenjl
    //! @date Sep 04, 2012
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SasExperimentalAndCalculatedData :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! storage for experimental saxs data, < q-value, intensity, error> , < r-value, density, error >
      SasScatteringData m_ExperimentalData;

      //! storage for calculated / simulated saxs data, < q-value, intensity, error> , <r-value, density, error >
      SasScatteringData m_CalculatedData;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SasExperimentalAndCalculatedData();

      //! @brief constructor from given input data
      SasExperimentalAndCalculatedData
      (
        const SasScatteringData &EXPERIMENTAL,
        const SasScatteringData &CALCULATED
      );

      //! @brief Clone function
      //! @return pointer to new SasExperimentalAndCalculatedData
      SasExperimentalAndCalculatedData *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the experimental data
      //! @return the experimental data
      const SasScatteringData &GetExperimentalData() const
      {
        return m_ExperimentalData;
      }

      //! @brief returns the calculated data
      //! @return the calculated data
      const SasScatteringData &GetCalculatedData() const
      {
        return m_CalculatedData;
      }
      //! @brief returns the calculated data
      //! @return the calculated data
      SasScatteringData &GetCalculatedData()
      {
        return m_CalculatedData;
      }

    /////////////////
    // operations  //
    /////////////////

      //! @brief scale the calculated and experimental data
      //! @param SCALING_FACTOR - the value to scale all of the Intensity values by
      void ScaleData( const double &SCALING_FACTOR);

      //! @brief scale the calculated intensity to align the experimental and calculated curves
      //! @param SCALING_FACTOR - the value to scale all of the Intensity values by
      void ScaleCalculatedData( const double &SCALING_FACTOR);

      //! @brief Set the Experimental Data values to the values of the object passed
      //! @param DATA_OBJECT - the object with the error values to copy
      void SetExperimentalData( const SasScatteringData &DATA_OBJECT);

      //! @brief Set the scale for the experimental and calculated intensities to the specified boundary
      void SetYScale( const double &Y_MAX);

      //! @brief Normalize the experimental data set by its largest intensity and then Normalize the calculated data
      //! @brief by its largest intensity value
      void NormalizeData();

      //! @brief take the data to log base 10
      void Log10();

      //! @brief transform the log10 data to the absolute scale
      void LogtoAbsolute();

      //! @brief take the derivative of the data
      void Derivative();

      //! @brief move the data to all positive numbers while maintaining identical morphology
      void SlideData( const double &Y_MIN);

      //! @brief compute the chi score between the experimental and calculated curves
      //! @return the chi score between the experimental and calculated curves
      double ComputeScoringFunction( bool USE_ERRORS) const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief read in the member data from BCL
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &ReadFromDataFile( std::istream &ISTREAM);

      //! @brief read fit file format from CRYSOL
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &ReadFromCrysolFitFile( std::istream &ISTREAM, const double &FIRST_EXPERIMENTAL_POINT);

      //! @brief read fit file format from CRYSOL
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &ReadFromFoxsFile( std::istream &ISTREAM);

      //! @brief write to std::ostream in three columns : Q-value, experimental intensity, and calculated intensity
      //! @param OSTREAM outputstream to write to
      //! @return outputstream which was written to
      std::ostream &WriteToGnuplot( std::ostream &OSTREAM) const;

      //! @brief write to filename in four columns : Q-value, experimental intensity, and calculated intensity
      //! @param filename to write to
      void WriteToGnuplotFileName( const std::string &FILENAME) const;

      //! @brief write to filename in three columns : Q-value, calculated intensity, experimental error
      //! @param filename to write to
      std::ostream &WriteToGnomeFormat( std::ostream &OSTREAM) const;

      //! @brief create filename to write and call WriteToGnomeFile
      //! @param filename to write to
      void WriteToGnomeFileName( const std::string &FILENAME) const;

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class SasExperimentalAndCalculatedData

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_SAS_EXPERIMENTAL_AND_CALCULATED_DATA_H_
