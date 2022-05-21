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
#include "opti/bcl_opti_approximator_powell.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "opti/bcl_opti_criterion_number_iterations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opti_approximator_powell.cpp
  //! @brief this example tests the implementation Powell's method for finding a local minimum
  //!
  //! @author woetzen, fischea
  //! @date Dec 18, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOptiApproximatorPowell :
     public ExampleInterface
   {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleOptiApproximatorPowell
    ExampleOptiApproximatorPowell *Clone() const
    {
      return new ExampleOptiApproximatorPowell( *this);
    }

  //////////
  // data //
  //////////

    // function class for example_function
    class SquareExampleFunction :
      public math::FunctionInterfaceSerializable< linal::Vector3D, double>
    {

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
      double operator()( const linal::Vector3D &ARGUMENT) const
      {
        return math::Sqr( ARGUMENT.Norm());
      }

      //! @brief Clone function
      //! @return pointer to new SquareExampleFunction
      SquareExampleFunction *Clone() const
      {
        return new SquareExampleFunction( *this);
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
        OSTREAM << "x^2"; return OSTREAM;
      }

    }; // class SquareExampleFunction

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
      // create starting vector and tolerance for the minimum
      const linal::Vector3D start( 2.1, 3.5, 1.9);
      const double tolerance( 0.0001);

      // create example function
      util::ShPtr< math::FunctionInterfaceSerializable< linal::Vector3D, double> > sp_objective( new SquareExampleFunction());

      // create termination criteria
      opti::CriterionNumberIterations< linal::Vector3D, double> criterion( 20);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create the approximator
      opti::ApproximatorPowell< linal::Vector3D, double> approximator
      (
        *sp_objective,
        criterion,
        storage::Vector< linal::Vector3D>::Create
        (
          *coord::GetAxes().e_X, *coord::GetAxes().e_Y, *coord::GetAxes().e_Z
        ),
        start
      );

      // clone the approximator
      util::ShPtr< opti::ApproximatorPowell< linal::Vector3D, double> > sp_approximator_clone
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
        ( GetStaticClassName< opti::ApproximatorPowell< linal::Vector3D, double> >())
      );

    ////////////////
    // operations //
    ////////////////

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
      opti::ApproximatorPowell< linal::Vector3D, double> approximator_read;
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
          ( approximator.GetTracker().GetBest()->First() - approximator_read.GetTracker().GetBest()->First()).Norm(),
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
          (
            approximator.GetTracker().GetCurrent()->First() - approximator_read.GetTracker().GetCurrent()->First()
          ).Norm(),
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

   }; // class ExampleOptiApproximatorPowell

   const ExampleClass::EnumType ExampleOptiApproximatorPowell::s_Instance
   (
     GetExamples().AddEnum( ExampleOptiApproximatorPowell())
   );

   //! single instance of that class
   const util::SiPtr< const util::ObjectInterface> ExampleOptiApproximatorPowell::SquareExampleFunction::s_Instance
   (
     GetObjectInstances().AddInstance( new ExampleOptiApproximatorPowell::SquareExampleFunction())
   );

} // namespace bcl
