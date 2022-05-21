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
#include "score/bcl_score_residual_dipolar_coupling_q_value.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ResidualDipolarCouplingQValue::s_Instance
    (
      util::Enumerated< math::FunctionInterfaceSerializable< nmr::RDCContainer, double> >::AddInstance
      (
        new ResidualDipolarCouplingQValue()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new ResidualDipolarCouplingQValue
    ResidualDipolarCouplingQValue *ResidualDipolarCouplingQValue::Clone() const
    {
      return new ResidualDipolarCouplingQValue( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ResidualDipolarCouplingQValue::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &ResidualDipolarCouplingQValue::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "rdc_q_value");

      // end
      return s_default_scheme;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator() for calculating the agreement between theoretical and experimental RDCs
    //! @param RDC_RESTRAINTS contains the experimental RDCs and calculated RDCs
    //! @return double which is the agreement between the theoretical and experimental RDCs
    double ResidualDipolarCouplingQValue::operator()( const nmr::RDCContainer &RDC_RESTRAINTS) const
    {
      // initialize sums
      double square_experimental_sum( 0.0);
      double square_deviation_sum( 0.0);

      // iterate over the data
      for
      (
        storage::Vector< double>::const_iterator
          exp_itr( RDC_RESTRAINTS.GetExperimentalValues().Begin()),
          exp_itr_end( RDC_RESTRAINTS.GetExperimentalValues().End()),
          calc_itr( RDC_RESTRAINTS.GetCalculatedlValues().Begin()),
          calc_itr_end( RDC_RESTRAINTS.GetCalculatedlValues().End());
        exp_itr != exp_itr_end && calc_itr != calc_itr_end; ++exp_itr, ++calc_itr
      )
      {
        // update the square deviation
        square_deviation_sum += math::Sqr( *exp_itr - *calc_itr);

        // update the square experimental sum
        square_experimental_sum += math::Sqr( *exp_itr);
      }

      // if square_experimental_sum is zero (no restraints passed or some problem with data)
      if( square_experimental_sum == 0.0)
      {
        // return 0 (basically the maximum score)
        return 0.0;
      }

      // calculate the q value
      double q_value( math::Sqrt( square_deviation_sum / square_experimental_sum));

      // convert q-value to score by subtracting 1
      return q_value - 1.0;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write detailed scheme and values to OSTREAM
    //! @param RDC_RESTRAINTS contains the experimental RDCs and calculated RDCs
    //! @param OSTREAM std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &ResidualDipolarCouplingQValue::WriteDetailedSchemeAndValues
    (
      const nmr::RDCContainer &RDC_RESTRAINTS,
      std::ostream &OSTREAM
    ) const
    {
      // write the scheme
      OSTREAM << GetScheme() << '\n';

      // initialize sums
      double square_experimental_sum( 0.0);
      double square_deviation_sum( 0.0);

      // iterate over the data
      for
      (
        storage::Vector< double>::const_iterator
          exp_itr( RDC_RESTRAINTS.GetExperimentalValues().Begin()),
          exp_itr_end( RDC_RESTRAINTS.GetExperimentalValues().End()),
          calc_itr( RDC_RESTRAINTS.GetCalculatedlValues().Begin()),
          calc_itr_end( RDC_RESTRAINTS.GetCalculatedlValues().End());
        exp_itr != exp_itr_end && calc_itr != calc_itr_end; ++exp_itr, ++calc_itr
      )
      {
        // write the experimental value, the calculated value, and the single Q value
        OSTREAM << util::Format()( *exp_itr) << ", " << util::Format()( *calc_itr) << " : "
                << util::Format()( math::Sqrt( math::Sqr( *exp_itr - *calc_itr) / math::Sqr( *exp_itr))) << '\n';

        // update the sums
        square_deviation_sum += math::Sqr( *exp_itr - *calc_itr);
        square_experimental_sum += math::Sqr( *exp_itr);
      }

      // write the Q-value
      OSTREAM << "Q-value: " << util::Format()( math::Sqrt( square_deviation_sum / square_experimental_sum)) << '\n';

      // end
      return OSTREAM;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ResidualDipolarCouplingQValue::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Calculates normalized standard deviation between computed and experimental RDCs");
      return parameters;
    }

  } // namespace score
} // namespace bcl
