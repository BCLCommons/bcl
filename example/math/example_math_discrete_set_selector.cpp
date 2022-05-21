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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "math/bcl_math_discrete_set_selector.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_trigonometric_transition.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_discrete_set_selector.cpp
  //! @brief this example tests the implementation of DiscreteSetSelector
  //!
  //! @author geanesar
  //! @date Nov 1, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathDiscreteSetSelector :
     public ExampleInterface
   {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleMathDiscreteSetSelector
    ExampleMathDiscreteSetSelector *Clone() const
    {
      return new ExampleMathDiscreteSetSelector( *this);
    }

  //////////
  // data //
  //////////

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief run routine
    //! this is performing the execution of the example
    int Run() const
    {

      // The states to use
      linal::Vector< double> states( 5);

      states( 0) = -0.75;
      states( 1) = -0.25;
      states( 2) = 0;
      states( 3) = 0.25;
      states( 4) = 0.75;

      // Construction
      math::DiscreteSetSelector def_selector;
      
      // Copy constructor 
      // ** assignment operator is below where it can actually be tested **
      math::DiscreteSetSelector copied_selector( def_selector);

      // clone operation
      util::ShPtr< math::DiscreteSetSelector> new_selector
      (
        def_selector.Clone()
      );
      BCL_ExampleCheck( new_selector.IsDefined(), true);

      // Destruction
      new_selector = util::ShPtr< math::DiscreteSetSelector>();

      // Make sure an operator call before setting weights gives an undefined size_t
      BCL_ExampleCheck( util::IsDefined( def_selector()), false);

      // by default prepare a uniform distribution
      def_selector.Prepare( states);

      double correct_uniform_weight( 1.0 / 5.0);
      linal::Vector< double> correct_weights( size_t( 5), correct_uniform_weight);

      BCL_ExampleCheck( def_selector.GetWeights(), correct_weights);
      
      // Assignment operator
      copied_selector = def_selector;

      // construction using a math function
      math::TrigonometricTransition fxn( -0.5, 0.5, 0.0, 1.0);
      math::DiscreteSetSelector selector( fxn);
      selector.Prepare( states);

      const linal::Vector< double> &weights( selector.GetWeights());

      // Check the weights vector sums to 1.0
      BCL_ExampleCheck( weights.Sum(), 1.0);

      // update correct weights 
      correct_weights( 0) = 0.0;
      correct_weights( 1) = 0.0585787;
      correct_weights( 2) = 0.2;
      correct_weights( 3) = 0.34142;
      correct_weights( 4) = 0.4;

      // These values are just slightly greater than the minimum that is needed to achieve the index in question
      double value_returns_index[] = { 0.0, 0.01, 0.059, 0.26, 0.61};

      // Check that each weight value corresponds to what it should be
      for( size_t i( 0); i < states.GetSize(); ++i)
      {
        BCL_ExampleIndirectCheck
        ( 
          math::EqualWithinTolerance( weights( i), correct_weights( i), 0.00001), 
          true,
          "correct weight should be " + util::Format()( correct_weights( i))
        );
        BCL_ExampleIndirectCheck
        (
          selector( value_returns_index[ i]),
          i,
          "selected state should be " + util::Format()( i)
        );
      }

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

   }; // class ExampleMathDiscreteSetSelector

   const ExampleClass::EnumType ExampleMathDiscreteSetSelector::s_Instance
   (
     GetExamples().AddEnum( ExampleMathDiscreteSetSelector())
   );

} // namespace bcl
