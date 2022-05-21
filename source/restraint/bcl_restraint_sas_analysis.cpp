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
#include "restraint/bcl_restraint_sas_analysis.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_cubic_spline_damped.h"
#include "math/bcl_math_statistics.h"
#include "restraint/bcl_restraint_sas_debye_interface.h"
#include "restraint/bcl_restraint_sas_experimental_and_calculated_density.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  /////////////////
  // data access //
  /////////////////

    //! @brief Compute the max dimension of the macromolecule from the p(r) curve
    //! @returns maximum pairwise distance in the macromolecule
    double SasAnalysis::ComputeDmax( const SasDensityData &EXPERIMENTAL_DATA)
    {
      // To get the last element of the Vector go the the end using the End function and move back 1 position
      storage::Vector< SasDistanceDensityPoint>::const_iterator data_itr( EXPERIMENTAL_DATA.End());
      data_itr--;

      return data_itr->GetRvalue();
    }

    //! @brief find maximum density point of the p(r) curve
    //! @returns the maximum value
    double SasAnalysis::FindDensitymax( const SasDensityData &EXPERIMENTAL_DATA)
    {
      double max_density( EXPERIMENTAL_DATA.Begin()->GetDensity());
      for
      (
        storage::Vector< SasDistanceDensityPoint>::const_iterator
          data_itr( EXPERIMENTAL_DATA.Begin()),
          data_itr_end( EXPERIMENTAL_DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        if( data_itr->GetDensity() > max_density)
        {
          max_density = data_itr->GetDensity();
        }
      }
      return max_density;
    }

    //! @brief find maximum density point of the p(r) curve
    //! @returns the maximum value
    double SasAnalysis::FindxDensitymax( const SasDensityData &EXPERIMENTAL_DATA)
    {
      double max_density( EXPERIMENTAL_DATA.Begin()->GetDensity());
      double max_r( EXPERIMENTAL_DATA.Begin()->GetRvalue());
      for
      (
        storage::Vector< SasDistanceDensityPoint>::const_iterator
          data_itr( EXPERIMENTAL_DATA.Begin()),
          data_itr_end( EXPERIMENTAL_DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        if( data_itr->GetDensity() > max_density)
        {
          max_density = data_itr->GetDensity();
          max_r = data_itr->GetRvalue();
        }
      }
      return max_r;
    }

    //! @brief scale the p(r) function to desired range
    //! @param EXPERIMENTAL_DATA - the data to scale
    //! @param SCALING_FACTOR - the value to scale all of the Intensity values by
    //! @returns scaled SasDensityData Object
    SasDensityData SasAnalysis::ScaleDensityData
    (
      const SasDensityData &EXPERIMENTAL_DATA,
      const double &SCALING_FACTOR
    )
    {
      SasDensityData scaled_data;

      scaled_data.SetBinSize( EXPERIMENTAL_DATA.GetBinSize());
      scaled_data.SetBinNumber( EXPERIMENTAL_DATA.GetBinNumber());
      scaled_data.SetDmax   ( EXPERIMENTAL_DATA.GetDmax());
      for
      (
        storage::Vector< SasDistanceDensityPoint>::const_iterator
          data_itr( EXPERIMENTAL_DATA.Begin()),
          data_itr_end( EXPERIMENTAL_DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        double r_value( data_itr->GetRvalue());
        double density( data_itr->GetDensity() * SCALING_FACTOR);
        double error( data_itr->GetError() * SCALING_FACTOR);

        SasDistanceDensityPoint scaled_point( r_value, density, error);
        scaled_data.PushBackDensity( scaled_point);
      }

      scaled_data.SetHmax( scaled_data.ComputeHmax());
      scaled_data.SetHxmax( scaled_data.ComputeHxmax());

      return scaled_data;
    }

    //! @brief shift the experimental p(r) function to overlap calculated p(r) function by aligning max peak
    //! @param EXPERIMENTAL_DATA - the data to scale
    //! @param SHIFT - the value to shift all of the r values by
    //! @returns scaled SasDensityData Object
    SasDensityData SasAnalysis::ShiftDensity( const SasDensityData &CALCULATED_DATA, const double &SHIFT)
    {
      SasDensityData shift_data;

      shift_data.SetBinSize( CALCULATED_DATA.GetBinSize());
      shift_data.SetBinNumber( CALCULATED_DATA.GetBinNumber());
      shift_data.SetDmax   ( CALCULATED_DATA.GetDmax());

      for
      (
        storage::Vector< SasDistanceDensityPoint>::const_iterator
          data_itr( CALCULATED_DATA.Begin()),
          data_itr_end( CALCULATED_DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        double r_shift_value( data_itr->GetRvalue() + SHIFT);
        double density( data_itr->GetDensity());
        double error( data_itr->GetError());

        SasDistanceDensityPoint shift_point( r_shift_value, density, error);
        shift_data.PushBackDensity( shift_point);
      }
      shift_data.SetHmax( shift_data.ComputeHmax());
      shift_data.SetHxmax( shift_data.ComputeHxmax());

      return shift_data;
    }

    //! @brief compute course integral of rectangles of width binsize
    //! @param DATA - the data to operate on
    //! @returns courseintegral value
    double SasAnalysis::ComputeCourseIntegral( const SasDensityData &DATA)
    {

      double integral( 0.0);
      double bin_size( DATA.GetBinSize());

      for
      (
        storage::Vector< SasDistanceDensityPoint>::const_iterator
          data_itr( DATA.Begin()),
          data_itr_end( DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        double density( data_itr->GetDensity());
        integral += density * bin_size;
      }

      return integral;
    }

    //! @brief function to calculate cumulative integral score for SAS pofr curves
    //! @param SAS_DATA experimental and calculated density distribution data
    //! @return integral score
    double SasAnalysis::CalculatePofRIntegralScore( const SasExperimentalAndCalculatedDensity &SAS_DATA)
    {
      SasExperimentalAndCalculatedDensity normalized_data( SAS_DATA);

      normalized_data.ScaleCalculatedDensity( 1 / SAS_DATA.GetCalculatedDensity().GetHmax());
      normalized_data.ScaleExperimentalDensity( 1 / SAS_DATA.GetExperimentalDensity().GetHmax());

      double cal_integral( ComputeCourseIntegral( normalized_data.GetCalculatedDensity()));
      double exp_integral( ComputeCourseIntegral( normalized_data.GetExperimentalDensity()));

      return exp_integral - cal_integral;
    }

    //! @brief function to calculate the excess integral of the calculated and experimental pofr curves
    //! @param SAS_DATA experimental and calculated density distribution data
    //! @return excess integral score
    double SasAnalysis::CalculatePofRExcessIntegralScore
    (
      const SasExperimentalAndCalculatedDensity &SAS_DATA
    )
    {
      double excess_integral( 0.0);

      SasExperimentalAndCalculatedDensity normalized_data( SAS_DATA);

      normalized_data.ScaleCalculatedDensity( 1 / SAS_DATA.GetCalculatedDensity().GetHmax());
      normalized_data.ScaleExperimentalDensity( 1 / SAS_DATA.GetExperimentalDensity().GetHmax());

      // Shift the calculated data
      normalized_data.ShiftDensity();

      for
      (
        storage::Vector< SasDistanceDensityPoint>::const_iterator
         cal_data_itr( normalized_data.GetCalculatedDensity().Begin()),
         exp_data_itr( normalized_data.GetExperimentalDensity().Begin()),
         cal_data_itr_end( normalized_data.GetCalculatedDensity().End());
        cal_data_itr != cal_data_itr_end;
        ++cal_data_itr, ++exp_data_itr
      )
      {
        double cal_density( cal_data_itr->GetDensity());
        double exp_density( exp_data_itr->GetDensity());

        if( cal_density > exp_density)
        {
          excess_integral += cal_density - exp_density;
        }
      }

      return excess_integral;
    }

    //! @brief function to calculate the amount of oscillation in the pofr curves
    //! @param SAS_DATA experimental and calculated density distribution data
    //! @return oscillation score
    double SasAnalysis::CalculatePofROscillationScore( const SasDensityData &DATA)
    {
      // A measure of smoothness is provided by the ratio ( ||P'|| / ||P||) / ( pi/( dmax - dmin))

      // Get the length of the distance distribution curve
      double delta( DATA.GetDmax());

      // Scale the Density Profiles
      SasDensityData scaled_data( ScaleDensityData( DATA, 1 / DATA.GetHmax()));

      // Compute the derivatives of the scaled data
      SasDensityData derivative( PofrDerivative( scaled_data));

      // Compute the Norm of the scaled pofr data
      storage::List< double> density( ConvertDensityDataToList( scaled_data));
      double density_norm( math::Statistics::Norm( density.Begin(), density.End()));

      // Compute the Norm of the derivative pofr data
      storage::List< double> der_density( ConvertDensityDataToList( derivative));
      double density_der_norm( math::Statistics::Norm( der_density.Begin(), der_density.End()));

      // Compute the oscillation score
      double oscillation_score( ( density_der_norm / density_norm) / ( math::g_Pi / delta));

      return oscillation_score;
    }

    //! @brief Compute the max momentum transfer vector value from the i(q) curve
    //! @returns maximum momentum transfer vector
    double SasAnalysis::ComputeQmax( const SasScatteringData &EXPERIMENTAL_DATA)
    {
      // To get the last element of the Vector go the the end using the End function and move back 1 position
      storage::Vector< SasScatteringPoint>::const_iterator data_itr( EXPERIMENTAL_DATA.End());
      --data_itr;

      return data_itr->GetQvalue();
    }

    //! @brief find the minimum intensity value in a SasScatteringData Object
    //! @param SAXS_DATA_OBJECT - the data to search
    //! @return the minimum intensity value
    double SasAnalysis::MinIntensity( const SasScatteringData &SAXS_DATA_OBJECT)
    {
      double min_intensity( SAXS_DATA_OBJECT.GetScatteringData().Begin()->GetIntensity());
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
          data_itr( SAXS_DATA_OBJECT.Begin()),
          data_itr_end( SAXS_DATA_OBJECT.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        if( data_itr->GetIntensity() < min_intensity)
        {
          min_intensity = data_itr->GetIntensity();
        }
      }
      return min_intensity;
    }

    //! @brief find the maximum intensity value in a SasScatteringData Object
    //! @param SAXS_DATA_OBJECT - the data to search
    //! @return the maximum intensity value
    double SasAnalysis::MaxIntensity( const SasScatteringData &SAXS_DATA_OBJECT)
    {
      double max_intensity( SAXS_DATA_OBJECT.GetScatteringData().Begin()->GetIntensity());
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
          data_itr( SAXS_DATA_OBJECT.Begin()),
          data_itr_end( SAXS_DATA_OBJECT.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        if( data_itr->GetIntensity() > max_intensity)
        {
          max_intensity = data_itr->GetIntensity();
        }
      }
      return max_intensity;
    }

    //! @brief scale the calculated intensity to align the experimental and calculated curves
    //! @param EXPERIMENTAL_DATA - the data to scale
    //! @param SCALING_FACTOR - the value to scale all of the Intensity values by
    //! @returns scaled SasScatteringData Object
    SasScatteringData SasAnalysis::ScaleData
    (
      const SasScatteringData &EXPERIMENTAL_DATA,
      const double &SCALING_FACTOR
    )
    {
      SasScatteringData scaled_data;
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
          data_itr( EXPERIMENTAL_DATA.Begin()),
          data_itr_end( EXPERIMENTAL_DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        double q_value( data_itr->GetQvalue());
        double intensity( data_itr->GetIntensity() * SCALING_FACTOR);
        double error( data_itr->GetError() * SCALING_FACTOR);

        SasScatteringPoint scaled_point( q_value, intensity, error);
        scaled_data.PushBackScattering( scaled_point);
      }
      return scaled_data;
    }

    //! @brief move the data to all positive numbers while maintaining identical morphology
    //! @brief This move will not adjust errors
    //! @param SAXS_DATA_OBJECT - the data to move
    //! @param MIN_VALUE - the value to base the move from
    //! @return data set with all positive values
    SasScatteringData SasAnalysis::SlideData( const SasScatteringData &SAXS_DATA_OBJECT, const double &MIN_VALUE)
    {
      SasScatteringData slid_data;
      double offset( math::Absolute( MIN_VALUE) + 1);

      BCL_MessageStd( "offset: " + util::Format()( offset));

      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
          data_itr( SAXS_DATA_OBJECT.Begin()),
          data_itr_end( SAXS_DATA_OBJECT.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        double q_value( data_itr->GetQvalue());
        double intensity( data_itr->GetIntensity() + offset);
        double error( data_itr->GetError());

        SasScatteringPoint scaled_point( q_value, intensity, error);
        slid_data.PushBackScattering( scaled_point);
      }
      return slid_data;
    }

    //! @brief take the data to log base 10
    SasScatteringData SasAnalysis::Log10( const SasScatteringData &EXPERIMENTAL_DATA)
    {
      SasScatteringData log_base10_data;

      // iterate over data to get values
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
        itr( EXPERIMENTAL_DATA.Begin()),
        itr_end( EXPERIMENTAL_DATA.End());
        itr != itr_end;
        ++itr
      )
      {
        // Store intensity
        double q_value( itr->GetQvalue());
        double intensity( log10( itr->GetIntensity()));
        double error( itr->GetError());

        // transform the error to a logaithmic scale.  If the error is greater than the intensity, only compute the
        // transformation based off the upper value ( should not happen),
        // otherwise take the average of the upper and lower values.

        if( itr->GetError() >= itr->GetIntensity())
        {
          error = log10( ( itr->GetIntensity() + itr->GetError()) / itr->GetIntensity());
        }
        else
        {
          error = 0.5 * log10( ( itr->GetIntensity() + itr->GetError()) / ( itr->GetIntensity() - itr->GetError()));
        }
        SasScatteringPoint log_point( q_value, intensity, error);
        log_base10_data.PushBackScattering( log_point);
      }
      return log_base10_data;
    }

    //! @brief transform the data from log base 10 scale to absolute scale
    //! @param DATA_OBJECT - the data to compute the log10 values of
    //! @returns data on absolute scale
    SasScatteringData SasAnalysis::LogtoAbsolute( const SasScatteringData &DATA_OBJECT)
    {
      SasScatteringData absolute_data;

      // iterate over data to get values
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
        itr( DATA_OBJECT.Begin()),
        itr_end( DATA_OBJECT.End());
        itr != itr_end;
        ++itr
      )
      {
        // Store intensity
        double q_value( itr->GetQvalue());
        double intensity( pow( 10, itr->GetIntensity()));
        double error( itr->GetError());

        double transform_error_1( error / 0.5);
        double transform_error_2( pow( 10, transform_error_1));

        double numerator( intensity * ( transform_error_2 - 1));
        double denominator( transform_error_2 + 1);

        double final_error( numerator / denominator);

        SasScatteringPoint absolute_point( q_value, intensity, final_error);
        absolute_data.PushBackScattering( absolute_point);
      }
      return absolute_data;
    }

    //! @brief take the pofr derivative of the data
    SasDensityData SasAnalysis::PofrDerivative( const SasDensityData &DATA)
    {

      const size_t data_size( DATA.GetDensitySize());

      // initialize math vectors with size of data set
      linal::Vector< double> data_values( data_size);

      // calculate delta
      storage::Vector< SasDistanceDensityPoint>::const_iterator delta_itr( DATA.Begin());
      storage::Vector< SasDistanceDensityPoint>::const_iterator delta_next( delta_itr);
      ++delta_next;

      const double delta( delta_next->GetRvalue() - delta_itr->GetRvalue());

      // populate math vectors with values from SAXS_DATA
      size_t pos( 0);

      for
      (
        storage::Vector< SasDistanceDensityPoint>::const_iterator
         data_itr( DATA.Begin()),
         data_itr_end( DATA.End());
        data_itr != data_itr_end;
        ++data_itr, ++pos
      )
      {
        data_values( pos) = data_itr->GetDensity();
      }

      // Create CublicSpline objects for the dataset
      math::CubicSplineDamped data_function;

      // Train the data set
      data_function.Train
      (
        delta_itr->GetRvalue(),
        delta,
        data_values
      );

      SasDensityData derivative_data;

      // Use the Cubic Spline Functions to Calculate the Derivative at provided point
      // Push the result into the Storage Containers
      for
      (
        storage::Vector< SasDistanceDensityPoint>::const_iterator
         data_itr( DATA.Begin()),
         data_itr_end( DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        double r_value( data_itr->GetRvalue());
        double density( data_function.dF( data_itr->GetRvalue()));
        double error( data_itr->GetError());

        SasDistanceDensityPoint derivative_point( r_value, density, error);
        derivative_data.PushBackDensity( derivative_point);
      }
      return derivative_data;
    }

    //! @brief take the derivative of the data using variable delta
    SasScatteringData SasAnalysis::Derivative( const SasScatteringData &EXPERIMENTAL_DATA)
    {
      const size_t data_size( EXPERIMENTAL_DATA.GetScatteringSize());

      // initialize math vectors with size of data set
      linal::Vector< double> x_values( data_size);
      linal::Vector< double> data_values( data_size);

      // populate math vectors with values from SAXS_DATA
      size_t pos( 0);

      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
          data_itr( EXPERIMENTAL_DATA.Begin()),
          data_itr_end( EXPERIMENTAL_DATA.End());
        data_itr != data_itr_end;
        ++data_itr, ++pos
      )
      {
        x_values( pos) = data_itr->GetQvalue();
        data_values( pos) = data_itr->GetIntensity();
      }

      //BCL_MessageStd( " data_values: " + util::Format()( data_values));

      // Create CublicSpline objects for the dataset
      math::CubicSplineDamped data_function;

      // Train the data set
      data_function.Train( x_values, data_values);

      SasScatteringData derivative_data;

      // Use the Cubic Spline Functions to Calculate the Derivative at provided point
      // Push the result into the Storage Containers
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
          data_itr( EXPERIMENTAL_DATA.Begin()),
          data_itr_end( EXPERIMENTAL_DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        double q_value( data_itr->GetQvalue());
        double intensity( data_function.dF( data_itr->GetQvalue()));
        double error( data_itr->GetError());

        //BCL_MessageStd( " derivative vd q_value: " + util::Format()( q_value));

        SasScatteringPoint derivative_point( q_value, intensity, error);
        derivative_data.PushBackScattering( derivative_point);
      }

      //BCL_MessageStd( " derivative vd values are: " + util::Format()( derivative_data));
      return derivative_data;
    }

    //! @brief Compute Experimental error from Simulated SAXS curves.  example: Crysol file
    //! @param SIMULATED_DATA - the data to simulate the experimental error for
    //! @return data set with simulated experimental error
    SasScatteringData SasAnalysis::AddErrors( const SasScatteringData &SIMULATED_DATA)
    {
      SasScatteringData simulated_errors;

      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
          data_itr( SIMULATED_DATA.Begin()),
          data_itr_end( SIMULATED_DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        double q_value( data_itr->GetQvalue());
        double intensity( data_itr->GetIntensity());
        double error( ComputeError( q_value, intensity));

        SasScatteringPoint simulated_error_point( q_value, intensity, error);
        simulated_errors.PushBackScattering( simulated_error_point);
      }
      return simulated_errors;
    }

    //! @brief Compute Experimental error for given Data point
    //! @param Q - momentum transfer vector
    //! @param INTENISTY - Intensity value for given momentum transfer vector
    //! @return data set with simulated experimental error
    double SasAnalysis::ComputeError( const double &Q, const double &INTENSITY)
    {
      double poisson_noise( math::Absolute( random::GetGlobalRandom().RandomPoisson( 10) / 10.0 - 1.0) + 1);
      double error( 0.15 * INTENSITY * ( Q + 0.001) * poisson_noise);
      return error;
    }

    //! @brief Set experimental error to defined value
    //! @param SIMULATED_DATA - the data to set the experimental error for
    //! @param ERROR_VALUE - the value to set all the experimental errors to
    //! @return data set with experimental errors set to the error value
    SasScatteringData SasAnalysis::SetErrors( const SasScatteringData &SIMULATED_DATA, const double &ERROR_VALUE)
    {
      SasScatteringData set_errors;
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
          data_itr( SIMULATED_DATA.Begin()),
          data_itr_end( SIMULATED_DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        double q_value( data_itr->GetQvalue());
        double intensity( data_itr->GetIntensity());
        double error( ERROR_VALUE);

        SasScatteringPoint set_error_point( q_value, intensity, error);
        set_errors.PushBackScattering( set_error_point);
      }
      return set_errors;
    }

    //! @brief function to calculate scaling factor for different SAXS intensity curves.  This version reads the error
    //! @brief associated with each q value.  The default error value is set to 1.0.
    //! @param SAXS_DATA experimental and calculated saxs data
    //! @return scaling factor
    double SasAnalysis::CalculateScalingWeight
    (
      const SasExperimentalAndCalculatedData &SAXS_DATA,
      const bool &USE_ERRORS
    )
    {
      // initialize sum
      double numerator( 0.0);
      double denominator( 0.0);

      // iterate over experimental and calculated data to get values
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
          exp_data_itr( SAXS_DATA.GetExperimentalData().Begin()),
          cal_data_itr( SAXS_DATA.GetCalculatedData().Begin()),
          exp_data_itr_end( SAXS_DATA.GetExperimentalData().End());
        exp_data_itr != exp_data_itr_end;
        ++exp_data_itr, ++cal_data_itr
      )
      {

        if( USE_ERRORS)
        {
          // if both I values are defined
          numerator   += exp_data_itr->GetIntensity() * cal_data_itr->GetIntensity() / math::Sqr( exp_data_itr->GetError());
          denominator += math::Sqr( cal_data_itr->GetIntensity()) / math::Sqr( exp_data_itr->GetError());
        }
        else
        {
          numerator   += exp_data_itr->GetIntensity() * cal_data_itr->GetIntensity();
          denominator += math::Sqr( cal_data_itr->GetIntensity());
        }
      }

      // endscaled_data.GetCalculatedDensity()
      return denominator == 0.0 ? 0.0 : numerator / denominator;
    }

    //! @brief function to calculate scaling factor for different SAXS intensity curves
    //! @param SAXS_DATA experimental and calculated saxs data
    //! @return scaling factor
    double SasAnalysis::CalculateStovgaardScalingWeight( const SasExperimentalAndCalculatedData &SAXS_DATA)
    {
      // initialize sum
      double numerator( 0.0);
      double denominator( 0.0);
      double alpha( 0.15);
      double beta( 0.30);
      double sigma( 0.0);
      double sqrSigma( 0.0);

      // iterate over experimental and calculated data to get values
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
         exp_data_itr( SAXS_DATA.GetExperimentalData().Begin()),
         cal_data_itr( SAXS_DATA.GetCalculatedData().Begin()),
         exp_data_itr_end( SAXS_DATA.GetExperimentalData().End());
        exp_data_itr != exp_data_itr_end;
        ++exp_data_itr, ++cal_data_itr
      )
      {
        // if both I values are defined
        sigma = ( exp_data_itr->GetIntensity() * ( exp_data_itr->GetQvalue() + alpha) * beta);
        sqrSigma = math::Sqr( sigma);

        numerator   += ( exp_data_itr->GetIntensity() * cal_data_itr->GetIntensity()) / sqrSigma;
        denominator += ( math::Sqr( cal_data_itr->GetIntensity())) / sqrSigma;
      }

      return denominator == 0.0 ? 0.0 : numerator / denominator;
    }

    //! @brief function to normalize the experimental and calculated curves individually
    //! @param SAXS_DATA scattering profile to normalize.  Range is between 0 and 1
    //! @return vector of normalized intensities
    SasScatteringData SasAnalysis::NormalizeData( const SasScatteringData &SAXS_DATA)
    {
      SasScatteringData normalized_data( SAXS_DATA);

      // Create list of doubles to hold raw intensity data
      storage::List< double> intensity( ConvertIntensityDataToList( SAXS_DATA));

      // Get the Norm of the experimental data to divide the errors by
      double exp_norm( math::Statistics::Norm( intensity.Begin(), intensity.End()));

      // Normalize the data set
      math::Statistics::Normalize( intensity.Begin(), intensity.End());

      // Write the Normalized data to the ScatteringData Object and scale the error by the norm value
      for
      (
        storage::Vector< SasScatteringPoint>::iterator
         norm_data_itr( normalized_data.Begin()),
         norm_data_itr_end( normalized_data.End());
         norm_data_itr != norm_data_itr_end;
        ++norm_data_itr
      )
      {
        const double normalized_value( intensity.FirstElement());
        norm_data_itr->SetIntensity( normalized_value);
        norm_data_itr->SetError( norm_data_itr->GetError() / exp_norm);
        intensity.PopFront();
      }

      //return SasScatteringData;
      return normalized_data;
    }

    //! @brief function to convert Scattering Intensities to a list of doubles
    //! @param SAXS_DATA input scattering profile to normalize.
    //! @return list of intensities of type double
    storage::List< double> SasAnalysis::ConvertIntensityDataToList( const SasScatteringData &SAXS_DATA)
    {
      // Create list of doubles to hold raw data
      storage::List< double> intensity;
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
        data_itr( SAXS_DATA.Begin()),
        data_itr_end( SAXS_DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        intensity.PushBack( data_itr->GetIntensity());
      }

      return intensity;
    }

    //! @brief function to convert Scattering Intensities to a list of doubles
    //! @param SAXS_DATA input scattering profile to normalize.
    //! @return list of intensities of type double
    storage::List< double> SasAnalysis::ConvertDensityDataToList( const SasDensityData &POFR_DATA)
    {
      // Create list of doubles to hold raw data
      storage::List< double> density;
      for
      (
        storage::Vector< SasDistanceDensityPoint>::const_iterator
        data_itr( POFR_DATA.Begin()),
        data_itr_end( POFR_DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        density.PushBack( data_itr->GetDensity());
      }

      return density;
    }

    //! @brief function to create Debye Implementation (GPU or CPU) based on parameters
    //! @param APPROXIMATE_LOOPS set to true to approximate loop regions
    //! @param APPROXIMATE_SIDE_CHAINS set to true to approximate side chains
    //! @param C1 Value of Exluded Volume adjustment parameter
    //! @param C2 Value of Hydration Shell adjustment parameter
    //! @param USE_CPU set to true to only use CPU implementation
    //! @return either CPU or GPU debye implementation based on provided parameters
    util::Implementation< SasDebyeInterface> SasAnalysis::SetDebyeImplementation
    (
      const bool &APPROXIMATE_LOOPS,
      const bool &APPROXIMATE_SIDE_CHAINS,
      const double &C1,
      const double &C2,
      const bool &USE_CPU,
      const bool &USE_SANS,
      const double &DEUTERIUM_EXCHANGE_PARAMETER
    )
    {
      // Setup Commandline Strings for either the opencl or non-opencl version of the code

      // First set up variables
      std::string opencl_parameters;
      std::string parameters;

      // Need to update opencl for SANS

      opencl_parameters =
        "OpenCLSaxsDebye(consider loops="
        + util::Format()( APPROXIMATE_LOOPS)
        + ", analytic=0, excluded volume="
        + util::Format()( C1)
        + ", hydration shell="
        + util::Format()( C2)
        + ", approximate_sidechains="
        + util::Format()( APPROXIMATE_SIDE_CHAINS)
        + " )";

      parameters =
        "SasDebye(consider loops="
        + util::Format()( APPROXIMATE_LOOPS)
        + ", analytic=0, excluded volume="
        + util::Format()( C1)
        + ", hydration shell="
        + util::Format()( C2)
        + ", approximate_sidechains="
        + util::Format()( APPROXIMATE_SIDE_CHAINS)
        + ", use_sans="
        + util::Format()( USE_SANS)
        + ", deuterium_percentage="
        + util::Format()( DEUTERIUM_EXCHANGE_PARAMETER)
        + " )";

      // Try to use OpenCL to compute the curves with the provided parameters, if that fails us the non-openCL version
      std::stringstream err_stream;
      util::Implementation< SasDebyeInterface> sas( opencl_parameters, err_stream);

      if( !sas.IsDefined())
      {
        // use the non-opencl version
        sas = parameters;
      }

      if( USE_CPU || USE_SANS)
      {
        sas = parameters;
      }

      // return the desired implementation
      return sas;
    }
  } // namespace restraint
} // namespace bcl
