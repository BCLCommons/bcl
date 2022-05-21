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
#include "restraint/bcl_restraint_saxs_data_reduction.h"

// includes from bcl - sorted alphabetically
#include "restraint/bcl_restraint_sas_analysis.h"
#include "restraint/bcl_restraint_sas_debye_interface.h"
#include "restraint/bcl_restraint_sas_transformation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SaxsDataReduction::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  /////////////////
  // operations  //
  /////////////////

    //! @brief compute the shannon numbers for saxs profile generation from experimental data
    //! @param DMAX - Maximum dimension of protein inferred from GNOM
    //! @param QMAX - Maximum momentum transfer value
    size_t SaxsDataReduction::ComputeShannonNumber( const double &DMAX, const double &QMAX)
    {
      return std::ceil( DMAX * QMAX / math::g_Pi);
    }

    //! @brief samples shannon bins n times ( n typically is 1000 - 3000) to produce noise free scattering profile
    //! @param PROTEIN_MODEL - protein model that contains the coordinates for current model
    //! @param NUMBER_OF_ITERATIONS - the number of iterations to compute ( n typically is 1000 - 3000)
    //! @param ORIGINAL_DATA - pointer to the preprocessed experimental data
    util::ShPtr< SasScatteringData> SaxsDataReduction::SasSignalRecovery
    (
      assemble::ProteinModel &PROTEIN_MODEL,
      size_t NUMBER_OF_ITERATIONS,
      const util::ShPtr< SasScatteringData> ORIGINAL_DATA,
      const double DMAX
    )
    {
      // Tasks
      // 1) partition a SAS dataset into ns equal bins for a given dmax
      // compute the number of bins on the saxs interval
      size_t number_of_bins
      (
        ComputeShannonNumber
        (
          DMAX,
          SasAnalysis::ComputeQmax( *ORIGINAL_DATA)
        )
      );

      double first_point( ORIGINAL_DATA->GetScatteringData().FirstElement().GetQvalue());
      double last_point( ORIGINAL_DATA->GetScatteringData().LastElement().GetQvalue());

      // compute the size of the bin on the interval based on the first and last point
      double binsize( ( last_point - first_point) / number_of_bins);

      double right_bin( binsize + first_point);
      int count( 0);

      storage::Vector< int> index_selected_points;

      // identify the index values that contain the desired ranges
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
         data_itr( ORIGINAL_DATA->Begin()),
         data_itr_end( ORIGINAL_DATA->End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        if( data_itr->GetQvalue() >= right_bin)
        {
          if( data_itr->GetQvalue() == right_bin)
          {
            index_selected_points.PushBack( count);
          }
          else
          {
            index_selected_points.PushBack( count - 1);
          }
          right_bin = right_bin + binsize;
        }
        ++count;
      }

      int left_index = 1;
      int right_index = index_selected_points( 0);

      // Set up containers to hold the data
      storage::Vector< SasScatteringPoint> subset;
      storage::Vector< storage::Vector< SasScatteringPoint> > reduced_data_container;
      storage::VectorND< 2, float> chi_score;
      storage::Vector< storage::VectorND< 2, float> > chi_container;

      std::string opencl_parameters
      (
        "OpenCLSaxsDebye(consider loops=0, analytic=0, excluded volume=1.0, hydration shell=0.0)"
      );
      std::string parameters
      (
        "SasDebye(consider loops=0, analytic=0, excluded volume=1.0, hydration shell=0.0)"
      );

      // try to create an opencl instance for this class, if that fails, use the non-opencl version
      std::stringstream err_stream;
      util::Implementation< SasDebyeInterface> saxs( opencl_parameters, err_stream);

      if( !saxs.IsDefined())
      {
        // Use the non-opencl version
        saxs = parameters;
      }
      saxs->SetExperimentalData( ORIGINAL_DATA);

      for( size_t iteration( 0); iteration < NUMBER_OF_ITERATIONS; ++iteration)
      {
        // Always select the first point
        subset.PushBack( ORIGINAL_DATA->GetScatteringLocation( 0));

        // randomly select points from each bin
        for( size_t i( 0); i < number_of_bins; ++i)
        {
          int index( SelectIndex( left_index, right_index));
          subset.PushBack( ORIGINAL_DATA->GetScatteringLocation( index));

          if( i != number_of_bins - 1)
          {
            left_index = right_index + 1;
            right_index = index_selected_points( i + 1);
          }
        } // end inner for loop

        reduced_data_container.PushBack( subset);
        subset.Reset();

        // 2) compute chi^2 for K iterations of shannon points
        storage::Vector< SasScatteringPoint> &data_point( reduced_data_container( iteration));

        // use the Clone to ShPtr to create a shared pointer to experimental data
        util::ShPtr< storage::Vector< SasScatteringPoint> > sp_reduced_data( util::CloneToShPtr( data_point));
        saxs->SetReducedExperimentalData( sp_reduced_data);

        // Compute the raw saxs profile
        SasExperimentalAndCalculatedData data_sets( saxs->operator()( PROTEIN_MODEL));

        // transform the profile into the desired form ( default of derivative score)
        SasExperimentalAndCalculatedData transformed_data( SasTransformation()( data_sets));

        // calculate the score
        float chi( score::SasType()( transformed_data));

        chi_score( 0) = chi;
        chi_score( 1) = iteration;

        chi_container.PushBack( chi_score);
        BCL_Message( util::Message::e_Standard, " iteration: " + util::Format()( iteration));

        // Reset indices for next iteration
        left_index = 0;
        right_index = index_selected_points( 0);
      } // end outer for loop

      // 3) chi^2 is taken as the median over k sampling rounds
      std::sort( chi_container.Begin(), chi_container.End());

      // compute the median and subtract 1 because of zero based indexing
      int median( ( NUMBER_OF_ITERATIONS / 2) - 1);

      int median_index( ( int)chi_container( median).Second());

      storage::Vector< SasScatteringPoint> &final_data( reduced_data_container( median_index));

      // Create a pointer to a new SasScatteringData Object with the reduced values

      util::ShPtr< SasScatteringData> sp_reduced_data
      (
        new SasScatteringData( final_data)
      );

      return sp_reduced_data;
    } // End Function

    //! @brief splits scatting profile into n (shannon) bins.  Selects datapoint to represent the curve from
    //!        each bin that has the least error associated with it.
    //! @param ORIGINAL_DATA - pointer to the preprocessed experimental data
    //! @param DENSITY_DATA - pointer to the transformed experimental data to P(r) domain
    util::ShPtr< SasScatteringData> SaxsDataReduction::SasSignalRecoveryEstimate
    (
      const util::ShPtr< SasScatteringData> ORIGINAL_DATA,
      const double DMAX
    )
    {
      // Tasks
      // 1) partition a SAS dataset into ns equal bins for a given dmax

      // compute the number of bins on the saxs interval
      size_t number_of_bins
      (
        ComputeShannonNumber
        (
          DMAX,
          SasAnalysis::ComputeQmax( *ORIGINAL_DATA)
        )
      );

      // compute the size of the bin on the interval

      double first_point( ORIGINAL_DATA->GetScatteringData().FirstElement().GetQvalue());
      double last_point( ORIGINAL_DATA->GetScatteringData().LastElement().GetQvalue());

      // compute the size of the bin on the interval based on the first and last point
      double binsize( ( last_point - first_point) / number_of_bins);
      double right_bin( binsize + first_point);
      int count( 0);

      storage::Vector< int> index_selected_points;

      // identify the index values that contain the desired ranges
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
        data_itr( ORIGINAL_DATA->Begin()),
        data_itr_end( ORIGINAL_DATA->End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        if( data_itr->GetQvalue() >= right_bin)
        {
          if( data_itr->GetQvalue() == right_bin)
          {
            index_selected_points.PushBack( count);
          }
          else
          {
            index_selected_points.PushBack( count - 1);
          }
          right_bin = right_bin + binsize;
        }
        ++count;
      }

      int left_index = 1;
      int right_index = index_selected_points( 0);

      // Set up containers to hold the data
      storage::Vector< SasScatteringPoint> estimated_data;

      // Always start with first point
      estimated_data.PushBack( ORIGINAL_DATA->GetScatteringLocation( 0));

      // select points with minimum error value from each bin
      for( size_t i( 0); i < index_selected_points.GetSize(); ++i)
      {
        // select points beginning with a left index of 1
        int index( SelectIndexMinError( ORIGINAL_DATA, left_index, right_index));
        estimated_data.PushBack( ORIGINAL_DATA->GetScatteringLocation( index));

        if( i != index_selected_points.GetSize() - 1)
        {
          left_index = right_index + 1;
          right_index = index_selected_points( i + 1);
        }
      } // end inner for loop

      util::ShPtr< SasScatteringData> sp_reduced_data
      (
        new SasScatteringData( estimated_data)
      );

      return sp_reduced_data;

    } // End Function

    //! @brief randomly select a point inside the given boundary
    //! @param LEFT - the left boundary of the interval
    //! @param RIGHT - the right boundary of the interval
    //! @return the randomly selected value inside the interval [ LEFT, RIGHT]
    int SaxsDataReduction::SelectIndex( int LEFT, int RIGHT)
    {
      return random::GetGlobalRandom().Random( LEFT, RIGHT);
    }

    //! @brief randomly select a point inside the given boundary
    //! @param ORIGINAL_DATA - pointer to the experimental data
    //! @param LEFT - the left boundary of the interval
    //! @param RIGHT - the right boundary of the interval
    //! @return the randomly selected value inside the interval [ LEFT, RIGHT]
    int SaxsDataReduction::SelectIndexMinError
    (
      const util::ShPtr< SasScatteringData> ORIGINAL_DATA,
      int LEFT,
      int RIGHT
    )
    {
      // Pick a point that minimizes the error in a bin
      double error( ORIGINAL_DATA->GetScatteringLocation( LEFT).GetError());

      int index( LEFT);

      for( int i( LEFT); i <= ( RIGHT); ++i)
      {
        if( ORIGINAL_DATA->GetScatteringLocation( i).GetError() < error)
        {
          error = ORIGINAL_DATA->GetScatteringLocation( i).GetError();
          index = i;
        }
        else
        {
          // Do nothing
        }
      }

      return index;
    }

  } // namespace restraint
} // namespace bcl
