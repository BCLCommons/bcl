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
#include "opti/bcl_opti_approximator_root_bisect.h"

// includes from bcl - sorted alphabetically
#include "opti/bcl_opti_criterion_convergence_argument.h"
#include "util/bcl_util_wrapper_base.h"

// external includes - sorted alphabetically

namespace bcl
{
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //!
   //! @example example_opti_approximator_root_bisect.cpp
   //! @brief this example tests the implementation of the bisection method for finding roots of a function using a
   //! cube function.
   //!
   //! @author woetzen, fischea
   //! @date Dec 18, 2012
   //! @remarks status complete
   //! @remarks reviewed by nobody on
   //!
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   class ExampleOptiApproximatorRootBisect :
     public ExampleInterface
   {

  public:

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to a new ExampleOptiApproximatorRootBisect
    ExampleOptiApproximatorRootBisect *Clone() const
    {
      return new ExampleOptiApproximatorRootBisect( *this);
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

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief returns the function value for the given argument
      //! @param ARGUMENT argument to calculate the function value for
      //! @return function value for the given argument: ARGUMENT^3
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

      // borders of the interval
      const double border_left( -3.5);
      const double border_right( 3.0);

      // tolerance for the approximation result
      const double root_tolerance( 0.0001);

      // check function value of cube function for left border
      BCL_ExampleCheck( cube_example_function.operator()( border_left), -42.875);

      // check function value of cube function for right border
      BCL_ExampleCheck( cube_example_function.operator()( border_right), 27.0);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create termination criteria
      opti::CriterionConvergenceArgument< double, double> criterion( 1, root_tolerance);

      // create the approximator
      opti::ApproximatorRootBisect< double, double> approximator
      (
        cube_example_function, criterion, border_left, border_right
      );

      // clone the approximator
      util::ShPtr< opti::ApproximatorRootBisect< double, double> > sp_approximator_clone
      (
        approximator.Clone()
      );

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck
      (
        approximator.GetClassIdentifier(),
        ( GetStaticClassName< opti::ApproximatorRootBisect< double, double> >())
      );

    ///////////////
    // operators //
    ///////////////

      // approximation
      approximator.Approximate();
      const util::ShPtr< storage::Pair< double, double> > sp_result
      (
        approximator.GetTracker().GetBest()
      );

      BCL_Message
      (
        util::Message::e_Standard,
        "the following root could be found x,f(x) = " + util::Format()( *sp_result) + " within " +
        util::Format()( approximator.GetTracker().GetIteration()) + " iterations"
      );

      // check result pair function value
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinAbsoluteTolerance( sp_result->Second(), 0.0, root_tolerance),
        true,
        "iterative root approximation"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // write approximator
      WriteBCLObject( approximator);

      // read approximator
      opti::ApproximatorRootBisect< double, double> approximator_read;
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

   }; // class ExampleOptiApproximatorRootBisect

   const ExampleClass::EnumType ExampleOptiApproximatorRootBisect::s_Instance
   (
     GetExamples().AddEnum( ExampleOptiApproximatorRootBisect())
   );

   //! single instance of that class
   const util::SiPtr< const util::ObjectInterface> ExampleOptiApproximatorRootBisect::CubeExampleFunction::s_Instance
   (
     GetObjectInstances().AddInstance( new ExampleOptiApproximatorRootBisect::CubeExampleFunction())
   );

} // namespace bcl
