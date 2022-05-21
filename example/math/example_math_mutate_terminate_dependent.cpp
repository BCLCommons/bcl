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
#include "math/bcl_math_mutate_terminate_dependent.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_mutate_vector.h"
#include "opti/bcl_opti_criterion_number_iterations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    // explicit instantiation of math::MutateTerminateDependent< linal::Vector< double>, double>
    template class MutateTerminateDependent< linal::Vector< double>, double>;
  } // namespace math

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_mutate_terminate_dependent.cpp
  //!
  //! @author alexanns
  //! @date
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathMutateTerminateDependent :
    public ExampleInterface
  {
  public:

    ExampleMathMutateTerminateDependent *Clone() const
    {
      return new ExampleMathMutateTerminateDependent( *this);
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

    int Run() const
    {
      opti::CriterionNumberIterations< linal::Vector< double>, double>::s_Instance.IsDefined();
      storage::List
      <
        storage::Pair
        <
          util::ShPtr< opti::CriterionInterface< linal::Vector< double>, double> >,
          util::ShPtr< math::MutateInterface< linal::Vector< double> > >
        >
      > terminates_mutates;

      terminates_mutates.PushBack
      (
        storage::Pair
        <
          util::ShPtr< opti::CriterionInterface< linal::Vector< double>, double> >,
          util::ShPtr< math::MutateInterface< linal::Vector< double> > >
        >
        (
          util::ShPtr< opti::CriterionInterface< linal::Vector< double>, double> >
          (
            new opti::CriterionNumberIterations< linal::Vector< double>, double>( 10)
          ),
          util::ShPtr< math::MutateInterface< linal::Vector< double> > >( new math::MutateVector( 1, 1, 1, false, false))
        )
      );

      opti::Tracker< linal::Vector< double>, double> tracker;
      linal::Vector< double> argument( 1, 5.6);
      tracker.Track( util::CloneToShPtr( storage::Pair< linal::Vector< double>, double>( argument, 0.5)));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      math::MutateTerminateDependent< linal::Vector< double>, double> def_constr;

      // constructor taking parameters
      const std::string scheme( "scheme");
      math::MutateTerminateDependent< linal::Vector< double>, double> param_constr( terminates_mutates, tracker, scheme);
      {
        BCL_ExampleCheck( param_constr.GetScheme(), scheme);
        linal::Vector< double> result( *param_constr( argument).GetArgument());
        BCL_ExampleCheck( result( 0) == argument( 0), false);
      }

      // clone constructor
      util::ShPtr< math::MutateTerminateDependent< linal::Vector< double>, double> > clone_constr( param_constr.Clone());
      {
        BCL_ExampleCheck( clone_constr->GetScheme(), scheme);
        linal::Vector< double> result( *clone_constr->operator()( argument).GetArgument());
        BCL_ExampleCheck( result( 0) == argument( 0), false);
      }

    /////////////////
    // data access //
    /////////////////

      // test GetClassIdentifier
      BCL_ExampleCheck
      (
        ( GetStaticClassName< math::MutateTerminateDependent< linal::Vector< double>, double> >()),
        clone_constr->GetClassIdentifier()
      );

      // test GetScheme
      {
        BCL_ExampleCheck( param_constr.GetScheme(), scheme);
      }

    ///////////////
    // operators //
    ///////////////

    ///////////////////
    // test operator //
    ///////////////////

      // terminate and mutate critieria not met
      {
        linal::Vector< double> result( *param_constr( argument).GetArgument());
        BCL_ExampleCheck( result( 0) == argument( 0), false);
      }

      // terminate and mutate critieria met
      {
        storage::List
        <
          storage::Pair
          <
            util::ShPtr< opti::CriterionInterface< linal::Vector< double>, double> >,
            util::ShPtr< math::MutateInterface< linal::Vector< double> > >
          >
        > terminates_mutates_test;

        terminates_mutates.PushBack
        (
          storage::Pair
          <
            util::ShPtr< opti::CriterionInterface< linal::Vector< double>, double> >,
            util::ShPtr< math::MutateInterface< linal::Vector< double> > >
          >
          (
            util::ShPtr< opti::CriterionInterface< linal::Vector< double>, double> >
            (
              new opti::CriterionNumberIterations< linal::Vector< double>, double>( 0)
            ),
            util::ShPtr< math::MutateInterface< linal::Vector< double> > >( new math::MutateVector( 1, 1, 1, false, false))
          )
        );

        linal::Vector< double> argument_test( 1, 5.6);
        math::MutateTerminateDependent< linal::Vector< double>, double> test( terminates_mutates_test, tracker, scheme);
        util::ShPtr< linal::Vector< double> > result( test( argument_test).GetArgument());
        BCL_ExampleCheck( result.IsDefined(), false);
      }

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // filename
      BCL_MessageStd( "writing and reading from file");
      WriteBCLObject( param_constr);
      math::MutateTerminateDependent< linal::Vector< double>, double> read_constr;
      ReadBCLObject( read_constr);
      {
        read_constr.SetTracker( tracker);
        BCL_ExampleCheck( read_constr.GetScheme(), scheme);
        linal::Vector< double> result( *read_constr( argument).GetArgument());
        BCL_ExampleCheck( result( 0) == argument( 0), false);
      }

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMathMutateTerminateDependent

  const ExampleClass::EnumType ExampleMathMutateTerminateDependent::s_Instance
  (
    GetExamples().AddEnum( ExampleMathMutateTerminateDependent())
  );

} // namespace bcl
