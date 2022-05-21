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
#include "opti/bcl_opti_approximator_root_regula_falsi.h"

// includes from bcl - sorted alphabetically
#include "opti/bcl_opti_criterion_combine.h"
#include "opti/bcl_opti_criterion_convergence_argument.h"
#include "opti/bcl_opti_criterion_number_iterations.h"
#include "util/bcl_util_wrapper_base.h"

// external includes - sorted alphabetically

namespace bcl
{
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //!
   //! @example example_opti_approximator_root_regula_falsi.cpp
   //! @brief this example tests the implementation of regula falsi (false position) root finding method implementation
   //! using a cube function
   //!
   //! @author woetzen, fischea
   //! @date Dec 13, 2012
   //! @remarks status complete
   //! @remarks reviewed by nobody on
   //!
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   class ExampleOptiApproximatorRootRegulaFalsi :
     public ExampleInterface
   {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleOptiApproximatorRootRegulaFalsi
    ExampleOptiApproximatorRootRegulaFalsi *Clone() const
    {
      return new ExampleOptiApproximatorRootRegulaFalsi( *this);
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
      opti::CriterionNumberIterations< double, double>::s_Instance.IsDefined();
      opti::CriterionCombine< double, double>::s_Instance.IsDefined();
      opti::CriterionConvergenceArgument< double, double>::s_Instance.IsDefined();

      // create cube example function
      CubeExampleFunction cube_example_function;

      // borders of the interval
      const double border_left( -3.5);
      const double border_right( 3.0);

      // tolerance for the approximation result
      const double root_tolerance( 0.0001);

      // create termination criteria
      opti::CriterionCombine< double, double> criterion_combine;
      criterion_combine.InsertCriteria( opti::CriterionConvergenceArgument< double, double>( 1, root_tolerance));
      criterion_combine.InsertCriteria( opti::CriterionNumberIterations< double, double>( 100));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create the approximator
      opti::ApproximatorRootRegulaFalsi< double, double> approximator
      (
        cube_example_function, criterion_combine, border_left, border_right
      );

      // clone the approximator
      util::ShPtr< opti::ApproximatorRootRegulaFalsi< double, double> > sp_approximator_clone
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

      // check class identifier
      BCL_ExampleCheck
      (
        approximator.GetClassIdentifier(),
        (
          GetStaticClassName< opti::ApproximatorRootRegulaFalsi< double, double> >()
        )
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

      // test the clone
      sp_approximator_clone->Approximate();
      const util::ShPtr< storage::Pair< double, double> > sp_result_clone
      (
        sp_approximator_clone->GetTracker().GetBest()
      );

      BCL_Example_Check
      (
        sp_approximator_clone->GetTracker().GetIteration() == approximator.GetTracker().GetIteration() &&
        sp_approximator_clone->GetTracker().GetBestIteration() == approximator.GetTracker().GetBestIteration() &&
        sp_approximator_clone->GetTracker().GetBest()->First() == approximator.GetTracker().GetBest()->First() &&
        sp_approximator_clone->GetTracker().GetBest()->Second() == approximator.GetTracker().GetBest()->Second() &&
        sp_approximator_clone->GetTracker().GetCurrent()->First() == approximator.GetTracker().GetCurrent()->First() &&
        sp_approximator_clone->GetTracker().GetCurrent()->First() == approximator.GetTracker().GetCurrent()->First(),
        "clone provides different result"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // write approximator
      WriteBCLObject( approximator);

      // read criterion
      opti::ApproximatorRootRegulaFalsi< double, double> approximator_read;
      ReadBCLObject( approximator_read);

      BCL_Example_Check
      (
        math::EqualWithinAbsoluteTolerance
        (
          approximator.GetTracker().GetBest()->First(),
          approximator_read.GetTracker().GetBest()->First(),
          root_tolerance
        ) &&
        math::EqualWithinAbsoluteTolerance
        (
          approximator.GetTracker().GetBest()->Second(),
          approximator_read.GetTracker().GetBest()->Second(),
          root_tolerance
        ),
        "read in object should have as best model " + util::Format()( approximator.GetTracker().GetBest()->First()) +
        " with a value of " + util::Format()( approximator.GetTracker().GetBest()->Second()) +
        " but has " + util::Format()( approximator_read.GetTracker().GetBest()->First()) +
        " with a value of " + util::Format()( approximator_read.GetTracker().GetBest()->Second())
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

   }; // class ExampleOptiApproximatorRootRegulaFalsi

   const ExampleClass::EnumType ExampleOptiApproximatorRootRegulaFalsi::s_Instance
   (
     GetExamples().AddEnum( ExampleOptiApproximatorRootRegulaFalsi())
   );

   //! single instance of that class
   const util::SiPtr< const util::ObjectInterface>
   ExampleOptiApproximatorRootRegulaFalsi::CubeExampleFunction::s_Instance
   (
     GetObjectInstances().AddInstance( new ExampleOptiApproximatorRootRegulaFalsi::CubeExampleFunction())
   );

} // namespace bcl
