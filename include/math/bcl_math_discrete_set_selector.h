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

#ifndef BCL_MATH_DISCRETE_SET_SELECTOR_H_
#define BCL_MATH_DISCRETE_SET_SELECTOR_H_

// include the namespace header
#include "bcl_math.h"

// includes from bcl - sorted alphabetically
#include "bcl_math_function_interface_serializable.h"
#include "linal/bcl_linal_vector.h"
#include "linal/bcl_linal_vector_const_interface.h"
#include "util/bcl_util_sh_ptr.h"

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DiscreteSetSelector
    //! @brief given a vector of values and a probability distribution, this class weights each vector component
    //!        based on the component's value, and can select components from the given vector stochastically based on
    //!        the calculated weights (used for selecting from a discrete set of objects given a probability distribution)
    //!
    //! @see @link example_math_discrete_set_selector.cpp @endlink
    //! @author geanesar
    //! @date Aug 25 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DiscreteSetSelector :
      public FunctionInterfaceSerializable< double, size_t>
    {

    protected:

      //! weight vector; holds a weighting term for each state that is stored in the class
      linal::Vector< double> m_Weights;

      //! the probability function used to calculate the weighting vector (if applicable)
      util::ShPtr< FunctionInterfaceSerializable< double, double> > m_ProbabilityFunction;

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief constructor with
      DiscreteSetSelector();

      //! @brief constructor with a probability function and intial values
      //! @param PROBABILITY_FUNCTION the probability function to use to select from values given to the object
      //! @param VALUES the values to add initially
      DiscreteSetSelector
      (
        const FunctionInterfaceSerializable< double, double> &PROBABILITY_FUNCTION,
        const linal::VectorConstInterface< double> &VALUES = linal::Vector< double>()
      );

      //! @brief copy constructor
      DiscreteSetSelector( const DiscreteSetSelector &OTHER);

      //! @brief assignment operator
      DiscreteSetSelector &operator =( const DiscreteSetSelector &OTHER);

      //! @brief clone function
      //! @return a copy of this class
      DiscreteSetSelector *Clone() const;

      //! @brief the name of this class
      //! @return the name of the class
      const std::string &GetClassIdentifier() const;

      //! @brief Set the states vector
      //! @param STATE_VALUES a vector of states that will be converted into weights
      void Prepare( const linal::VectorConstInterface< double> &VALUES);

      //! @brief Reset the class, zero the weights
      void Reset();

      //! @brief operator() selects a state from the distrubution
      //! @param ARGUMENT a number between 0.0 and 1.0 that is used for member selection
      //! @return the index of the state that was chosen
      size_t operator ()( const double &ARGUMENT = util::GetUndefined< double>()) const;

      //! @brief get the weights vector
      //! @return linal::vector containing weights
      const linal::Vector< double> &GetWeights() const;

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class DiscreteSetSelector

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_DISCRETE_SET_SELECTOR_H_
