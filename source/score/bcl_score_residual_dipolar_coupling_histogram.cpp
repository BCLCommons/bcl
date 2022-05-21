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
#include "score/bcl_score_residual_dipolar_coupling_histogram.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_histogram.h"
#include "math/bcl_math_running_min_max.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ResidualDipolarCouplingHistogram::s_Instance
    (
      util::Enumerated< math::FunctionInterfaceSerializable< nmr::RDCContainer, double> >::AddInstance
      (
        new ResidualDipolarCouplingHistogram()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new ResidualDipolarCoupling
    ResidualDipolarCouplingHistogram *ResidualDipolarCouplingHistogram::Clone() const
    {
       return new ResidualDipolarCouplingHistogram( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ResidualDipolarCouplingHistogram::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &ResidualDipolarCouplingHistogram::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "rdc_histogram");

      // end
      return s_default_scheme;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator() for calculating the agreement between theoretical and experimental RDCs
    //! @param RDC_CONTAINER contains both the theoretical and experimental RDCs
    //! @return double which is the agreement between the theoretical and experimental RDCs
    double ResidualDipolarCouplingHistogram::operator()( const nmr::RDCContainer &RDC_CONTAINER) const
    {
      static const size_t s_nr_bins( 25);

      // initialize min/max dataset to get histogram boundaries
      math::RunningMinMax< double> data_set_min_max;

      // iterate through the RDCContainer
      for
      (
        storage::Vector< double>::const_iterator
          exp_itr( RDC_CONTAINER.GetExperimentalValues().Begin()),
          exp_itr_end( RDC_CONTAINER.GetExperimentalValues().End()),
          calc_itr( RDC_CONTAINER.GetCalculatedlValues().Begin()),
          calc_itr_end( RDC_CONTAINER.GetCalculatedlValues().End());
        exp_itr != exp_itr_end && calc_itr != calc_itr_end; ++exp_itr, ++calc_itr
      )
      {
        data_set_min_max += *exp_itr;
        data_set_min_max += *calc_itr;
      }

      // Initialize variables for creating histograms
      const double min_value( data_set_min_max.GetMin() - 0.1);
      const double bin_size( ( data_set_min_max.GetMax() - min_value) / s_nr_bins + 0.01);

      // Create histograms for experimental and calculated RDCs
      math::Histogram exp_histogram( min_value, bin_size, s_nr_bins);
      math::Histogram calc_histogram( min_value, bin_size, s_nr_bins);

      // Calculate the histograms
      exp_histogram.CalculateHistogram( RDC_CONTAINER.GetExperimentalValues());
      calc_histogram.CalculateHistogram( RDC_CONTAINER.GetCalculatedlValues());

      // Initialize sum of square difference
      double sum_of_squares( 0.0);

      // Initialize the number of non empty bins
      size_t non_empty_bins( 0);

      // Iterate through the bins in both histograms
      for
      (
        linal::Vector< double>::const_iterator
          exp_count_itr( exp_histogram.GetHistogram().Begin()),
          exp_count_itr_end( exp_histogram.GetHistogram().End()),
          calc_count_itr( calc_histogram.GetHistogram().Begin()),
          calc_count_itr_end( calc_histogram.GetHistogram().End());
        exp_count_itr != exp_count_itr_end && calc_count_itr != calc_count_itr_end;
        ++exp_count_itr, ++calc_count_itr
      )
      {

        // Increment the number of non empty bins
        if( *exp_count_itr > 0 || *calc_count_itr > 0)
        {
          ++non_empty_bins;
        }

        // Calculate the difference then square it then add to the sum
        sum_of_squares += math::Sqr( *exp_count_itr - *calc_count_itr);
      }

      // Calculate the RMSD
      const double histogram_rmsd( math::Sqrt( sum_of_squares / double( non_empty_bins)));

      BCL_MessageStd( util::Format()( histogram_rmsd));

      // convert the rmsd to a score
      const double score( histogram_rmsd == 0.0 ? -1.0 : -1.0 / histogram_rmsd);

      BCL_MessageStd( util::Format()( score));

      // incorporate the number of rdcs available into the score by taking the natural log of that number and
      // multiplying it by the rmsd
      return log( double( RDC_CONTAINER.GetExperimentalValues().GetSize() + 1)) * score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ResidualDipolarCouplingHistogram::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Compares experimental with calculated RDC histograms");
      return parameters;
    }

  } // namespace score
} // namespace bcl
