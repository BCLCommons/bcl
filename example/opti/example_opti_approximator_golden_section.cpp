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
#include "opti/bcl_opti_approximator_golden_section.h"

// includes from bcl - sorted alphabetically
#include "opti/bcl_opti_criterion_number_iterations.h"
#include "util/bcl_util_wrapper_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opti_approximator_golden_section.cpp
  //! @brief this example tests the implementation of the golden section search algorithm
  //!
  //! @author fischea, woetzen
  //! @date Dec 18, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOptiApproximatorGoldenSection :
     public ExampleInterface
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleOptiApproximatorGoldenSection
    ExampleOptiApproximatorGoldenSection *Clone() const
    {
      return new ExampleOptiApproximatorGoldenSection( *this);
    }

  //////////
  // data //
  //////////

    class QuadraticFunction :
      public math::FunctionInterfaceSerializable< double, double>
    {

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      QuadraticFunction()
      {
      }

      //! @brief Clone function
      //! @return pointer to new QuadraticFunction
      QuadraticFunction *Clone() const
      {
        return new QuadraticFunction( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief returns the function value for the given argument quadratic function
      //! @param INPUT argument to calculate the function value for
      //! @return function value for the given argument
      double operator()( const double &INPUT) const
      {
        return INPUT * INPUT;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    }; // class QuadraticFunction

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
      storage::Pair< double, double>::s_Instance.IsDefined();

      // define borders and tolerance for the approximation
      const double border_left( -3.5);
      const double border_right( 3.0);
      const double tolerance( 0.0001);

      // create quadratic example function
      util::ShPtr< QuadraticFunction> quadratic_example_function( new QuadraticFunction());

      // check function value of quadratic function for left border
      BCL_ExampleCheck( quadratic_example_function->operator ()( border_left), 12.25);

      // check function value of quadratic function for right border
      BCL_ExampleCheck( quadratic_example_function->operator ()( border_right), 9.0);

      // create termination criterion
      opti::CriterionNumberIterations< double, double> criterion( 15);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create the approximator
      opti::ApproximatorGoldenSection< double, double> approximator
      (
        *quadratic_example_function, criterion, border_left, border_right
      );

      // clone the approximator
      util::ShPtr< opti::ApproximatorGoldenSection< double, double> > sp_approximator_clone
      (
        approximator.Clone()
      );

      BCL_Example_Check
      (
        approximator.GetTracker().GetBest()->Second() == sp_approximator_clone->GetTracker().GetBest()->Second()
        &&
        approximator.GetTracker().GetBest()->First() == sp_approximator_clone->GetTracker().GetBest()->First()
        &&
        approximator.GetTracker().GetCurrent()->Second() == sp_approximator_clone->GetTracker().GetCurrent()->Second()
        &&
        approximator.GetTracker().GetCurrent()->First() == sp_approximator_clone->GetTracker().GetCurrent()->First(),
        "members of cloned approximator don't match the members of the original approximator"
      );

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck
      (
        approximator.GetClassIdentifier(),
        ( GetStaticClassName< opti::ApproximatorGoldenSection< double, double> >())
      );

    ////////////////
    // operations //
    ////////////////

      // check ResPhi
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinAbsoluteTolerance( approximator.GetResPhi(), 0.38196, 0.0001),
        true,
        "ResPhi returns wrong result"
      );

      // approximate the local minimum
      approximator.Approximate();

      BCL_Message
      (
        util::Message::e_Standard,
        "the following minimizer could be found x,f(x) = " +
        util::Format()( approximator.GetTracker().GetBest()->Second()) + " within " +
        util::Format()( approximator.GetTracker().GetIteration()) + " iterations"
      );

      // check result pair function value
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinAbsoluteTolerance( approximator.GetTracker().GetBest()->Second(), 0.0, tolerance),
        true,
        "iterative minimizer approximation"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // write approximator
      WriteBCLObject( approximator);

      // read approximator
      opti::ApproximatorGoldenSection< double, double> approximator_read;
      ReadBCLObject( approximator_read);

      // compare read in object to written object
      BCL_Example_Check
      (
        math::EqualWithinAbsoluteTolerance
        (
          approximator.GetTracker().GetBest()->Second() - approximator_read.GetTracker().GetBest()->Second(),
          0.0,
          0.00001
        ) &&
        math::EqualWithinAbsoluteTolerance
        (
          approximator.GetTracker().GetBest()->First() - approximator_read.GetTracker().GetBest()->First(),
          0.0,
          0.00001
        ) &&
        math::EqualWithinAbsoluteTolerance
        (
          approximator.GetTracker().GetCurrent()->Second() - approximator_read.GetTracker().GetCurrent()->Second(),
          0.0,
          0.00001
        ) &&
        math::EqualWithinAbsoluteTolerance
        (
          approximator.GetTracker().GetCurrent()->First() - approximator_read.GetTracker().GetCurrent()->First(),
          0.0,
          0.00001
        )
        ,
        "best model is " + util::Format()( approximator_read.GetTracker().GetBest()) +
        " and should be " + util::Format()( approximator.GetTracker().GetBest()) +
        " current model is " + util::Format()( approximator_read.GetTracker().GetCurrent()) +
        " and should be " + util::Format()( approximator.GetTracker().GetCurrent())
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // class ExampleOptiApproximatorGoldenSection

  const ExampleClass::EnumType ExampleOptiApproximatorGoldenSection::s_Instance
  (
     GetExamples().AddEnum( ExampleOptiApproximatorGoldenSection())
  );

  //! single instance of that class
  const util::SiPtr< const util::ObjectInterface> ExampleOptiApproximatorGoldenSection::QuadraticFunction::s_Instance
  (
    GetObjectInstances().AddInstance( new ExampleOptiApproximatorGoldenSection::QuadraticFunction())
  );

} // namespace bcl
