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
#include "opti/bcl_opti_approximator_root_secant.h"

// includes from bcl - sorted alphabetically
#include "opti/bcl_opti_criterion_convergence_argument.h"
#include "util/bcl_util_wrapper_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opti_approximator_root_secant.cpp
  //! @brief this example tests the implementation of the secant root finding method implementation using a cube
  //! function
  //!
  //! @author woetzen, fischea
  //! @date Dec 17, 2012
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOptiApproximatorRootSecant :
    public ExampleInterface
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleOptiApproximatorRootSecant
    ExampleOptiApproximatorRootSecant *Clone() const
    {
      return new ExampleOptiApproximatorRootSecant( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    // function class for example_function
    class CubeExampleFunction :
      public math::FunctionInterfaceSerializable< double, double>
    {

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief returns the function value for the given argument quadratic function
      //! @param ARGUMENT argument to calculate the function value for
      //! @return function value for the given argument
      double operator()( const double &ARGUMENT) const
      {
        return ARGUMENT * ARGUMENT * ARGUMENT;
      }

      //! @brief Clone function
      //! @return pointer to a new CubeExampleFunction
      CubeExampleFunction *Clone() const
      {
        return new CubeExampleFunction( *this);
      }

      //! @brief write to ostream
      //! @param OSTREAM is the output stream
      //! @param INDENT number of indentations
      //! @return returns the output stream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // end
        return OSTREAM;
      }

      //! @brief read from istream
      //! @param ISTREAM is the input stream
      //! @return returns the input stream
      std::istream &Read( std::istream &ISTREAM)
      {
        // end
        return ISTREAM;
      }

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      std::ostream &Scheme( std::ostream &OSTREAM, const util::Format &FORMAT = util::Format()) const
      {
        OSTREAM << "x^3";
        return OSTREAM;
      }

    }; // class CubeExampleFunction

  ////////////////
  // operations //
  ////////////////

    //! @brief run routine
    //! this is performing the execution of the example
    int Run() const
    {
      // create cube example function
      CubeExampleFunction cube_example_function;

      // define borders of the approximation interval
      const double left_border( -3.5);
      const double right_border( 3.0);

      // check function value of cube function for left border
      BCL_Example_Check
      (
        cube_example_function.operator()( left_border) == -42.875,
        "function value of cube function for left border should be -42.875, but is " +
        util::Format()( cube_example_function.operator()( left_border))
      );

      // check function value of cube function for right border
      BCL_Example_Check
      (
        cube_example_function.operator()( right_border) == 27.0,
        "function value of cube function for right border should be 27.0, but is " +
        util::Format()( cube_example_function.operator()( right_border))
      );

      // create termination criterion
      opti::CriterionConvergenceArgument< double, double> criterion
      (
        1, 0.00001
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct approximator
      opti::ApproximatorRootSecant< double> approximator
      (
        cube_example_function, criterion, left_border, right_border
      );

      // clone the approximator
      util::ShPtr< opti::ApproximatorRootSecant< double> > sp_approximator_clone
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
        ( GetStaticClassName< opti::ApproximatorRootSecant< double> >())
      );

    ////////////////
    // operations //
    ////////////////

      // approximate the root of the example function
      approximator.Approximate();

      // check result pair function value
      BCL_Example_Check
      (
        math::EqualWithinAbsoluteTolerance( approximator.GetTracker().GetBest()->Second(), 0.0, 0.0001),
        "root value should be close to zero"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // write approximator
      WriteBCLObject( approximator);

      // read approximator
      opti::ApproximatorRootSecant< double> approximator_read;
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

  }; // class ExampleOptiApproximatorRootSecant

  const ExampleClass::EnumType ExampleOptiApproximatorRootSecant::s_Instance
  (
    GetExamples().AddEnum( ExampleOptiApproximatorRootSecant())
  );

  //! single instance of that class
  const util::SiPtr< const util::ObjectInterface> ExampleOptiApproximatorRootSecant::CubeExampleFunction::s_Instance
  (
    GetObjectInstances().AddInstance( new ExampleOptiApproximatorRootSecant::CubeExampleFunction())
  );

} // namespace bcl
