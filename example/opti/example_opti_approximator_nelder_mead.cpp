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
#include "opti/bcl_opti_approximator_nelder_mead.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d_operations.h"
#include "opti/bcl_opti_criterion_number_iterations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opti_approximator_nelder_mead.cpp
  //! @brief this example tests the implementation of the Nelder Mead optimization algorithm
  //!
  //! @author mendenjl, fischea
  //! @date Dec 21, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOptiApproximatorNelderMead :
     public ExampleInterface
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleOptiApproximatorNelderMead
    ExampleOptiApproximatorNelderMead *Clone() const
    {
      return new ExampleOptiApproximatorNelderMead( *this);
    }

  //////////
  // data //
  //////////

    // function class for example_function
    class SquareExampleFunction :
      public math::FunctionInterfaceSerializable< linal::Vector3D, double>
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
      double operator()( const linal::Vector3D &ARGUMENT) const
      {
        return math::Sqr( ARGUMENT.Norm());
      }

      //! @brief Clone function
      //! @return pointer to a new SquareExampleFunction
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
      opti::CriterionNumberIterations< linal::Vector3D, double>::s_Instance.IsDefined();
      storage::Pair< linal::Vector3D, double>::s_Instance.IsDefined();

      const util::Format format_double;
      const util::Format format_vector3D;
      const util::Format format_result;

      const linal::Vector3D vertex_a( 1.1, 0, 0);
      const linal::Vector3D vertex_b( 0, 1.2, 0);
      const linal::Vector3D vertex_c( 0, 0, 1.5);
      const linal::Vector3D vertex_d( -1, -1, -1);

      // build up start simplex with 4 arguments
      storage::List< linal::Vector3D> simplex;
      simplex.PushBack( vertex_a);
      simplex.PushBack( vertex_b);
      simplex.PushBack( vertex_c);
      simplex.PushBack( vertex_d);

      // create example function
      SquareExampleFunction function;

      // minimum difference
      const double min_difference( 0.0000000001);

      // create termination criterion
      opti::CriterionNumberIterations< linal::Vector3D, double> criterion( 100);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create approximator
      opti::ApproximatorNelderMead< linal::Vector3D, double> approximator
      (
        function, criterion, simplex, min_difference
      );

      // clone the approximator
      util::ShPtr< opti::ApproximatorNelderMead< linal::Vector3D, double> > sp_approximator_clone
      (
        approximator.Clone()
      );

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck
      (
        approximator.GetClassIdentifier(),
        ( GetStaticClassName< opti::ApproximatorNelderMead< linal::Vector3D, double> >())
      );

    ////////////////
    // operations //
    ////////////////

      // approximate
      approximator.Approximate();

      // get result
      util::ShPtr< storage::Pair< linal::Vector3D, double> > result( approximator.GetTracker().GetBest());

      BCL_Message( util::Message::e_Standard, "The result is: " + format_result( result));

      // check if the result is correct
      BCL_Example_Check
      (
        math::EqualWithinAbsoluteTolerance( result->Second(), 0.0, 0.0001),
        "minimized value should be close to zero"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // write approximator
      WriteBCLObject( approximator);

      // read approximator
      opti::ApproximatorNelderMead< linal::Vector3D, double> approximator_read;
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

  }; // class ExampleOptiApproximatorNelderMead

  const ExampleClass::EnumType ExampleOptiApproximatorNelderMead::s_Instance
  (
     GetExamples().AddEnum( ExampleOptiApproximatorNelderMead())
  );

  //! single instance of that class
  const util::SiPtr< const util::ObjectInterface> ExampleOptiApproximatorNelderMead::SquareExampleFunction::s_Instance
  (
    GetObjectInstances().AddInstance( new ExampleOptiApproximatorNelderMead::SquareExampleFunction())
  );

} // namespace bcl
