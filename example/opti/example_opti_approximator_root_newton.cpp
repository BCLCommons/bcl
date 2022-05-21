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
#include "opti/bcl_opti_approximator_root_newton.h"

// includes from bcl - sorted alphabetically
#include "opti/bcl_opti_criterion_combine.h"
#include "opti/bcl_opti_criterion_convergence_argument.h"
#include "opti/bcl_opti_criterion_divergence_argument.h"
#include "util/bcl_util_wrapper_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opti_approximator_root_newton.cpp
  //! @brief this example tests the implementation of newton's root finding method implementation using a cube function
  //! and its first derivative
  //!
  //! @author woetzen, fischea
  //! @date Dec 13, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOptiApproximatorRootNewton :
    public ExampleInterface
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleOptiApproximatorRootNewton
    ExampleOptiApproximatorRootNewton *Clone() const
    {
      return new ExampleOptiApproximatorRootNewton( *this);
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
      //! @return pointer to new CubeExampleFunction
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
        return OSTREAM;
      }

      //! @brief read from istream
      //! @param ISTREAM is the input stream
      //! @return returns the input stream
      std::istream &Read( std::istream &ISTREAM)
      {
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

    class QuadraticFunction :
       public math::FunctionInterfaceSerializable< double, double>
    {

    public:

    //////////
    // data //
    //////////

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
        return 3 * INPUT * INPUT;
      }

    protected:

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // end
        return OSTREAM;
      }

    }; // class QuadraticFunction

  ////////////////
  // operations //
  ////////////////

    //! @brief run routine
    //! this is performing the execution of the example
    int Run() const
    {
      opti::CriterionCombine< double, double>::s_Instance.IsDefined();
      opti::CriterionConvergenceArgument
      <
        double, double
      >::s_Instance.IsDefined();
      opti::CriterionDivergenceArgument
      <
        double, double
      >::s_Instance.IsDefined();
      storage::Pair< double, double>::s_Instance.IsDefined();

      // create cube example function
      util::ShPtr< CubeExampleFunction> sp_cube_example_function( new CubeExampleFunction());

      // create cube derivative example function
      util::ShPtr< QuadraticFunction> sp_cube_example_function_derivative( new QuadraticFunction());

      // test value for previously defined functions
      const double test_value( 5.0);
      const double root_tolerance( 0.0001);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create termination criteria
      util::ShPtr
      <
        opti::CriterionCombine< double, double>
      > sp_criterion_combine
      (
        new opti::CriterionCombine< double, double>()
      );

      sp_criterion_combine->InsertCriteria( opti::CriterionConvergenceArgument< double, double>( 1, root_tolerance));
      sp_criterion_combine->InsertCriteria( opti::CriterionDivergenceArgument< double, double>( 1, 1000.0));

      util::ShPtr< opti::ApproximatorRootNewton< double> > sp_approximator
      (
        new opti::ApproximatorRootNewton< double>
        (
          *sp_cube_example_function,
          *sp_cube_example_function_derivative,
          *sp_criterion_combine,
          test_value
        )
      );

    /////////////////
    // data access //
    /////////////////

      // check class identifier
      BCL_ExampleCheck
      (
        sp_approximator->GetClassIdentifier(),
        (
          GetStaticClassName< opti::ApproximatorRootNewton< double> >()
        )
      );

    ////////////////
    // operations //
    ////////////////

      // approximation
      sp_approximator->Approximate();

      // approximation
      util::ShPtr< storage::Pair< double, double> > result_pair
      (
        sp_approximator->GetTracker().GetBest()
      );

      BCL_Message
      (
        util::Message::e_Standard,
        "the following root could be found x,f(x) = " + util::Format()( *result_pair) + " within " +
        util::Format()( sp_approximator->GetTracker().GetIteration()) + " iterations"
      );

      // check result pair function value
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinAbsoluteTolerance( result_pair->Second(), double( 0), root_tolerance),
        true,
        "iterative root approximation"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // write approximator
      WriteBCLObject( *sp_approximator);

      // read approximator
      opti::ApproximatorRootNewton< double> approximator_read;
      ReadBCLObject( approximator_read);

      BCL_Example_Check
      (
        math::EqualWithinAbsoluteTolerance
        (
          sp_approximator->GetTracker().GetBest()->Second() - approximator_read.GetTracker().GetBest()->Second(),
          0.0,
          0.00001
        ) &&
        math::EqualWithinAbsoluteTolerance
        (
          sp_approximator->GetTracker().GetBest()->First() - approximator_read.GetTracker().GetBest()->First(),
          0.0,
          0.00001
        ) &&
        math::EqualWithinAbsoluteTolerance
        (
          sp_approximator->GetTracker().GetCurrent()->Second() - approximator_read.GetTracker().GetCurrent()->Second(),
          0.0,
          0.00001
        ) &&
        math::EqualWithinAbsoluteTolerance
        (
          sp_approximator->GetTracker().GetCurrent()->First() - approximator_read.GetTracker().GetCurrent()->First(),
          0.0,
          0.00001
        ),
        "read in object does not match the written out object"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // class ExampleOptiApproximatorRootNewton

  const ExampleClass::EnumType ExampleOptiApproximatorRootNewton::s_Instance
  (
    GetExamples().AddEnum( ExampleOptiApproximatorRootNewton())
  );

  //! single instance of that class
  const util::SiPtr< const util::ObjectInterface> ExampleOptiApproximatorRootNewton::CubeExampleFunction::s_Instance
  (
    GetObjectInstances().AddInstance( new ExampleOptiApproximatorRootNewton::CubeExampleFunction())
  );

  //! single instance of that class
  const util::SiPtr< const util::ObjectInterface> ExampleOptiApproximatorRootNewton::QuadraticFunction::s_Instance
  (
    GetObjectInstances().AddInstance( new ExampleOptiApproximatorRootNewton::QuadraticFunction())
  );

} // namespace bcl
