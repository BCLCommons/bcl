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

// include for this class
#include "math/bcl_math_discrete_set_selector.h"
// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_polynomial.h"
#include "random/bcl_random_uniform_distribution.h"

namespace bcl
{
  namespace math
  {

    //! @brief default constructor; yields uniform probability distribution
    DiscreteSetSelector::DiscreteSetSelector() :
      m_Weights(),
      m_ProbabilityFunction()
    {
      Polynomial constant_polynomial;
      linal::Vector< double> const_vector( 1);
      const_vector( 0) = 1.0;
      constant_polynomial.SetCoefficients( const_vector);
      m_ProbabilityFunction = util::CloneToShPtr( constant_polynomial);
    }

    //! @brief constructor with a probability function and intial values
    //! @param PROBABILITY_FUNCTION the probability function to use to select from values given to the object
    //! @param VALUES the values to add initially
    DiscreteSetSelector::DiscreteSetSelector
    (
      const FunctionInterfaceSerializable< double, double> &PROBABILITY_FUNCTION,
      const linal::VectorConstInterface< double> &VALUES
    ) :
      m_Weights(),
      m_ProbabilityFunction( util::CloneToShPtr( PROBABILITY_FUNCTION))
    {
      Prepare( VALUES);
    }

    //! @brief copy constructor
    DiscreteSetSelector::DiscreteSetSelector( const DiscreteSetSelector &OTHER) :
      m_Weights( OTHER.m_Weights),
      m_ProbabilityFunction( util::CloneToShPtr( *OTHER.m_ProbabilityFunction))
    {
    }

    //! @brief assignment operator
    DiscreteSetSelector &DiscreteSetSelector::operator =( const DiscreteSetSelector &OTHER)
    {
      m_Weights = OTHER.m_Weights;
      m_ProbabilityFunction = util::CloneToShPtr( *OTHER.m_ProbabilityFunction);
      return *this;
    }

    //! @brief clone function
    //! @return a copy of this class
    DiscreteSetSelector *DiscreteSetSelector::Clone() const
    {
      return new DiscreteSetSelector( *this);
    }

    //! @brief the name of this class
    //! @return the name of the class
    const std::string &DiscreteSetSelector::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief set up the selector using the given values (resets everything beforehand)
    //! @param @VALUES a vector of states that will be converted into weights
    void DiscreteSetSelector::Prepare( const linal::VectorConstInterface< double> &VALUES)
    {
      BCL_Assert
      (
        m_ProbabilityFunction.IsDefined(),
        "Prepare was called before DiscreteSetSelector probability function was defined"
      );

      size_t n_values( VALUES.GetSize());
      m_Weights = linal::Vector< double>( n_values);
      // VALUES is empty, do nothing
      if( !n_values)
      {
        return;
      }

      // Initially set the weights vector to the raw probabilities, then scale it accoring to the sum of the values that
      // were given
      for( size_t i( 0); i < n_values; ++i)
      {
        m_Weights( i) = ( *m_ProbabilityFunction)( VALUES( i));
      }

      // if all values are zero then this should be a normalized uniform distribution, otherwise just normalize the weights
      double weights_sum( m_Weights.Sum());
      if( weights_sum == 0.0)
      {
        m_Weights = 1.0 / n_values;
      }
      else
      {
        m_Weights /= weights_sum;
      }
    }

    //! @brief Reset the class
    void DiscreteSetSelector::Reset()
    {
      m_Weights = linal::Vector< double>();
    }

    //! @brief operator() selects a state from the distrubution
    //! @param ARGUMENT a number between 0.0 and 1.0 that is used for member selection
    //! @return the index of the state that was chosen, or undefined if there were no weights selected
    size_t DiscreteSetSelector::operator()( const double &ARGUMENT) const
    {
      size_t result( util::GetUndefined< size_t>());

      double arg( ARGUMENT);
      if( !util::IsDefined( arg))
      {
        arg = random::GetGlobalRandom().Double();
      }

      if( m_Weights.GetSize())
      {
        // The objective sum should be between 0 and 1
        double sum( std::max( std::min( arg, 1.0), 0.0));

        // Add weight values until they exceed the objective sum
        for( result = 0; result < m_Weights.GetSize(); ++result)
        {
          sum -= m_Weights( result);
          if( sum <= 0.0)
          {
            break;
          }
        }
        result = std::min( result, m_Weights.GetSize() - 1);
      }
      return result;
    }

    //! @brief get the weights vector
    //! @return linal::vector containing weights
    const linal::Vector< double> &DiscreteSetSelector::GetWeights() const
    {
      return m_Weights;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DiscreteSetSelector::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return output stream which was written to
    std::ostream &DiscreteSetSelector::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace math
} // namespace bcl
