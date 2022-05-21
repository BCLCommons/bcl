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

#ifndef BCL_RESTRAINT_SAS_ANALYSIS_H_
#define BCL_RESTRAINT_SAS_ANALYSIS_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_restraint_sas_density_data.h"
#include "bcl_restraint_sas_scattering_data.h"
#include "bcl_restraint_sas_scattering_point.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "command/bcl_command_flag_interface.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SasAnalysis
    //! @brief Service Class with static functions to operate on SaxsAnalysis
    //! @details Service Class to compute Qmax, Dmax, Scale Data, Log10, Derivative, Derivative with variable delta
    //! @details This class will also be used to add additional analysis operations
    //!
    //! @see @link example_restraint_saxs_analysis.cpp @endlink
    //! @author putnamdk
    //! @date May 26, 2014
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SasAnalysis
    {

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief Compute the max dimension of the macromolecule from the p(r) curve
      //! @returns maximum pairwise distance in the macromolecule
      static double ComputeDmax( const SasDensityData &EXPERIMENTAL_DATA);

      //! @brief find maximum density point of the p(r) curve
      //! @returns the maximum value
      static double FindDensitymax( const SasDensityData &EXPERIMENTAL_DATA);

      //! @brief find x value of the maximum density point of the p(r) curve
      //! @returns the maximum value
      static double FindxDensitymax( const SasDensityData &EXPERIMENTAL_DATA);

      //! @brief scale the p(r) function to desired range
      //! @param EXPERIMENTAL_DATA - the data to scale
      //! @param SCALING_FACTOR - the value to scale all of the Intensity values by
      //! @returns scaled SasDensityData Object
      static SasDensityData ScaleDensityData( const SasDensityData &EXPERIMENTAL_DATA, const double &SCALING_FACTOR);

      //! @brief shift the calculated p(r) function to overlap the experimental calculated p(r) function
      //! @param CALCULATED_DATA - the data to scale
      //! @param SHIFT - the value to shift all of the r values by
      //! @returns scaled SasDensityData Object
      static SasDensityData ShiftDensity( const SasDensityData &CALCULATED_DATA, const double &SHIFT);

      //! @brief compute course integral of rectangles of width binsize
      //! @param DATA - the data to operate on
      //! @returns courseintegral value
      static double ComputeCourseIntegral( const SasDensityData &DATA);

      //! @brief function to calculate cumulative integral score for SAS pofr curves
      //! @param SAS_DATA experimental and calculated density distribution data
      //! @return integral score
      static double CalculatePofRIntegralScore( const SasExperimentalAndCalculatedDensity &SAS_DATA);

      //! @brief function to calculate the excess integral of the calculated and experimental pofr curves
      //! @param SAS_DATA experimental and calculated density distribution data
      //! @return excess integral score
      static double CalculatePofRExcessIntegralScore( const SasExperimentalAndCalculatedDensity &SAS_DATA);

      //! @brief function to calculate the amount of oscillation in the pofr curves
      //! @param SAS_DATscaled_data.GetCalculatedDensity()A experimental and calculated density distribution data
      //! @return oscillation score
      static double CalculatePofROscillationScore( const SasDensityData &DATA);

      //! @brief Compute the max momentum transfer vector value from the i(q) curve
      //! @returns maximum momentum transfer vector
      static double ComputeQmax( const SasScatteringData &EXPERIMENTAL_DATA);

      //! @brief scale the calculated intensity to align the experimental and calculated curves
      //! @param EXPERIMENTAL_DATA - the data to scale
      //! @param SCALING_FACTOR - the value to scale all of the Intensity values by
      //! @returns scaled SasScatteringData Object
      static SasScatteringData ScaleData
      (
        const SasScatteringData &EXPERIMENTAL_DATA,
        const double &SCALING_FACTOR
      );

      //! @brief move the data to all positive numbers while maintaining identical morphology
      //! @param SAXS_DATA_OBJECT - the data to move
      //! @param MIN_VALUE - the value to base the move from
      //! @return data set with all positive values
      static SasScatteringData SlideData( const SasScatteringData &SAXS_DATA_OBJECT, const double &MIN_VALUE);

      //! @brief find the minimum intensity value in a SasScatteringData Object
      //! @param SAXS_DATA_OBJECT - the data to search
      //! @return the minimum intensity value
      static double MinIntensity( const SasScatteringData &SAXS_DATA_OBJECT);

      //! @brief find the maximum intensity value in a SasScatteringData Object
      //! @param SAXS_DATA_OBJECT - the data to search
      //! @return the maximum intensity value
      static double MaxIntensity( const SasScatteringData &SAXS_DATA_OBJECT);

      //! @brief take the data to log base 10
      //! @param EXPERIMENTAL_DATA - the data to compute the log10 values of
      //! @returns log10 of SasScatteringData Object
      static SasScatteringData Log10( const SasScatteringData &EXPERIMENTAL_DATA);

      //! @brief transform the data from log base 10 scale to absolute scale
      //! @param DATA_OBJECT - the data to compute the log10 values of
      //! @returns data on absolute scale
      static SasScatteringData LogtoAbsolute( const SasScatteringData &DATA_OBJECT);

      //! @brief take the derivative of the data
      //! @param EXPERIMENTAL_DATA - the data to compute the derivative of
      //! @return derivative of SasScatteringDataObject
      static SasScatteringData Derivative( const SasScatteringData &EXPERIMENTAL_DATA);

      //! @brief take the derivative of the pofr data
      //! @param EXPERIMENTAL_DATA - the data to compute the derivative of
      //! @return derivative of SasDensityData object
      static SasDensityData PofrDerivative( const SasDensityData &DATA);

      //! @brief Compute Experimental error from Simulated SAXS curves.  example: Crysol file
      //! @param SIMULATED_DATA - the data to simulate the experimental error for
      //! @return data set with simulated experimental error
      static SasScatteringData AddErrors( const SasScatteringData &SIMULATED_DATA);

      //! @brief Compute Experimental error for given Data point
      //! @param Q - momentum transfer vector
      //! @param INTENISTY - Intensity value for given momentum transfer vector
      //! @return data set with simulated experimental error
      static double ComputeError( const double &Q, const double &INTENSITY);

      //! @brief Set experimental error to defined value
      //! @param SIMULATED_DATA - the data to set the experimental error for
      //! @param ERROR_VALUE - the value to set all the experimental errors to
      //! @return data set with experimental errors set to the error value
      static SasScatteringData SetErrors( const SasScatteringData &SIMULATED_DATA, const double &ERROR_VALUE);

      //! @brief function to calculate scaling factor for different SAXS intensity curves
      //! @param SAXS_DATA experimental and calculated saxs data
      //! @return scaling factor
      static double CalculateScalingWeight( const SasExperimentalAndCalculatedData &SAXS_DATA, const bool &USE_ERRORS);

      //! @brief function to calculate scaling factor for different SAXS intensity curves
      //! @param SAXS_DATA experimental and calculated saxs data
      //! @return scaling factor
      static double CalculateStovgaardScalingWeight( const SasExperimentalAndCalculatedData &SAXS_DATA);

      //! @brief function to normalize the experimental and calculated curves individually
      //! @param SAXS_DATA scattering profile to normalize.  Range is between 0 and 1
      //! @return vector of normalized intensities
      static SasScatteringData NormalizeData( const SasScatteringData &SAXS_DATA);

      //! @brief function to convert Scattering Intensities to a list of doubles
      //! @param SAXS_DATA input scattering profile to normalize.
      //! @return list of intensities of type double
      static storage::List< double> ConvertIntensityDataToList( const SasScatteringData &SAXS_DATA);

      //! @brief function to convert Density Data to a list of doubles
      //! @param DATA input pofr profile to operate on.
      //! @return list of intensities of type double
      static storage::List< double> ConvertDensityDataToList( const SasDensityData &POFR_DATA);

      //! @brief function to create Debye Implementation (GPU or CPU) based on parameters
      //! @param APPROXIMATE_LOOPS set to true to approximate loop regions
      //! @param APPROXIMATE_SIDE_CHAINS set to true to approximate side chains
      //! @param C1 Value of Exluded Volume adjustment parameter
      //! @param C2 Value of Hydration Shell adjustment parameter
      //! @param USE_CPU set to true to only use CPU implementation
      //! @param USE_SANS set to use SANS implementation, otherwise SAXS implementation
      //! @return either CPU or GPU debye implementation based on provided parameters
      static util::Implementation< SasDebyeInterface> SetDebyeImplementation
      (
        const bool &APPROXIMATE_LOOPS,
        const bool &APPROXIMATE_SIDE_CHAINS,
        const double &C1,
        const double &C2,
        const bool &USE_CPU,
        const bool &USE_SANS,
        const double &DEUTERIUM_EXCHANGE_PARAMETER
      );

    }; // class SaxsAnalysis

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_SAS_ANALYSIS_H_
