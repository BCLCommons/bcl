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
#include "math/bcl_math_function_cached.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_function_cached.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathFunctionCached :
    public ExampleInterface
  {
  public:

    class Student
    {
      mutable signal::Signal1< const Student&> m_OnChangeSignal;
      mutable signal::Signal1< const Student&> m_DestructorSignal;

    public:

      signal::Signal1< const Student &> &GetOnChangeSignal() const
      {
        return m_OnChangeSignal;
      }

      signal::Signal1< const Student &> &GetDestructorSignal() const
      {
        return m_DestructorSignal;
      }

      std::string m_Name;
      int         m_Age;

      Student( const std::string &NAME, const int AGE) :
        m_Name( NAME),
        m_Age( AGE)
      {
      }

      ~Student()
      {
        m_DestructorSignal.Emit( *this);
      }

      void SetAge( const int AGE)
      {
        m_Age = AGE;
        m_OnChangeSignal.Emit( *this);
      }

    }; // class Student

    class TimeToGraduation :
      public math::FunctionInterfaceSerializable< Student, int>
    {
    private:

      mutable size_t m_NumberTimesCalled; //!< number of times the function was called

    public:
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      TimeToGraduation() :
        m_NumberTimesCalled( 0)
      {
      }

      //! @brief
      TimeToGraduation *Clone() const
      {
        return new TimeToGraduation( *this);
      }

      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      const size_t &GetNumberTimesOperatorWasCalled() const
      {
        return m_NumberTimesCalled;
      }

      int operator()( const Student &STUDENT) const
      {
        ++m_NumberTimesCalled;
        return 25 - STUDENT.m_Age;
      }

      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    }; // class TimeToGraduation

    ExampleMathFunctionCached *Clone() const
    {
      return new ExampleMathFunctionCached( *this);
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
      // sp to evaluation function
      util::ShPtr< math::FunctionInterfaceSerializable< Student, int> > sp_evaluate( new TimeToGraduation());
      util::ShPtr< TimeToGraduation> sp_time_to_graduation( sp_evaluate);

      // use another TimeToGraduation to calculate the correct time to graduation each time
      // this way sp_evaluate->GetNumberTimesOperatorWasCalled tells us how many times the operator
      // was called via the signal
      TimeToGraduation time_to_graduation_oracle;

      // student 1
      Student joe( "Joe", 20);
      Student nadja( "Nadja", 26);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // cache function from evaluation function
      math::FunctionCached< Student, int> cache_function( sp_evaluate, &Student::GetDestructorSignal);
      math::FunctionCached< Student, int> cache_function_no_invalidate( sp_evaluate, &Student::GetDestructorSignal);
      const int time_to_graduation_joe1_wrong( cache_function_no_invalidate( joe));
      const int time_to_graduation_nadja1_wrong( cache_function_no_invalidate( nadja));

      // copy constructor
      math::FunctionCached< Student, int> copy_cache_function( cache_function_no_invalidate);
      BCL_Example_Check
      (
        cache_function_no_invalidate.GetCacheSize() == 2
        && copy_cache_function.GetCacheSize() == 0
        && cache_function_no_invalidate.GetFunction() == copy_cache_function.GetFunction(),
        "copied cache function, should have cache of zero, and point to the same function"
      );

      // check that TimeToGraduation::operator() has been called twice (once for each student)
      BCL_ExampleIndirectCheck
      (
        sp_time_to_graduation->GetNumberTimesOperatorWasCalled(),
        2,
        "Cache function calls the evaluation function when it receives arguments"
      );

      // clone
      util::ShPtr< math::FunctionCached< Student, int> > sp_clone( cache_function.Clone());
      BCL_Example_Check
      (
        sp_evaluate == sp_clone->GetFunction(),
        "cloned cache function should point to evaluate function"
      );

    /////////////////
    // data access //
    /////////////////

      // class name
      BCL_MessageStd( "class name: " + sp_clone->GetClassIdentifier());
      BCL_Example_Check
      (
        ( GetStaticClassName< math::FunctionCached< Student, int> >()) == "bcl::math::FunctionCached<bcl::ExampleMathFunctionCached::Student,int>"
        && sp_clone->GetClassIdentifier() == ( GetStaticClassName< math::FunctionCached< Student, int> >()),
        "incorrect class name, should be bcl::math::FunctionCached<bcl::ExampleMathFunctionCached::Student,int> but is " +
        ( GetStaticClassName< math::FunctionCached< Student, int> >())
      );

      // access the function object
      BCL_Example_Check
      (
        sp_evaluate == cache_function.GetFunction(),
        "cache function should point to evaluate function"
      );

      // cache size
      BCL_Example_Check
      (
        cache_function.GetCacheSize() == 0,
        "cache size should be 0 but is: " + util::Format()( cache_function.GetCacheSize())
      );

    ////////////////
    // operations //
    ////////////////

      // add the onchange signal to the cache that will invalidate cached results
      const bool add1_success( cache_function.AddSignalHandlerForArgument( &Student::GetOnChangeSignal));
      BCL_Example_Check
      (
        add1_success,
        "was not able to add signal handler member function pointer for argument type Student"
      );
      const bool add2_success( cache_function.AddSignalHandlerForArgument( &Student::GetOnChangeSignal));
      BCL_Example_Check
      (
        !add2_success,
        "should not be able to successfully register the same signal handler twice"
      );

    ///////////////
    // operators //
    ///////////////

      // calculate the results directly for reference
      const int time_to_graduation_joe( time_to_graduation_oracle( joe));
      const int time_to_graduation_nadja( time_to_graduation_oracle( nadja));

      // calculate with the cache function
      const int time_to_graduation_joe1( cache_function( joe));
      const int time_to_graduation_nadja1( cache_function( nadja));

      // check that the results are correct
      BCL_Example_Check
      (
        time_to_graduation_joe == time_to_graduation_joe1
        && time_to_graduation_nadja == time_to_graduation_nadja1,
        "cache function returns different results than evaluation function. joe: " +
        util::Format()( time_to_graduation_joe) + " != " + util::Format()( time_to_graduation_joe1) +
        " nadja: " + util::Format()( time_to_graduation_nadja) + " != " + util::Format()( time_to_graduation_nadja1)
      );
      // check that the results are the same as those calculated initially
      BCL_Example_Check
      (
        time_to_graduation_joe == time_to_graduation_joe1_wrong
        && time_to_graduation_nadja == time_to_graduation_nadja1_wrong,
        "cache function returns different results than evaluation function. joe: " +
        util::Format()( time_to_graduation_joe) + " != " + util::Format()( time_to_graduation_joe1_wrong) +
        " nadja: " + util::Format()( time_to_graduation_nadja) + " != " +
        util::Format()( time_to_graduation_nadja1_wrong)
      );

      // check cache size
      BCL_ExampleCheck( cache_function.GetCacheSize(), 2);

      // check that TimeToGraduation::operator() was called 4 times (once for each student (2) X once for each cache (2))
      BCL_ExampleIndirectCheck
      (
        sp_time_to_graduation->GetNumberTimesOperatorWasCalled(),
        4,
        "Cache function only calls the evaluation function if it receives new arguments"
      );

      // calculate the second time, which should be faster
      const int time_to_graduation_joe2( cache_function( joe));
      const int time_to_graduation_nadja2( cache_function( nadja));

      // check that the results are correct
      BCL_Example_Check
      (
        time_to_graduation_joe == time_to_graduation_joe2
        && time_to_graduation_nadja == time_to_graduation_nadja2,
        "cache function returns different results than evaluation function. joe: " +
        util::Format()( time_to_graduation_joe) + " != " + util::Format()( time_to_graduation_joe2) +
        " nadja: " + util::Format()( time_to_graduation_nadja) + " != " + util::Format()( time_to_graduation_nadja2)
      );

      // check cache size
      BCL_ExampleCheck( cache_function.GetCacheSize(), 2);

      // check that operator() was not recalled, since joe and nadja haven't changed and were already in the cache
      BCL_ExampleIndirectCheck
      (
        sp_time_to_graduation->GetNumberTimesOperatorWasCalled(),
        4,
        "Cache function only calls the evaluation function if it receives new arguments"
      );

      // set age
      joe.SetAge( 21);
      const int time_to_graduation_joe3( cache_function( joe));
      const int time_to_graduation_nadja3( cache_function( nadja));

      // check that the results are correct
      BCL_Example_Check
      (
        time_to_graduation_oracle( joe) == time_to_graduation_joe3
        && time_to_graduation_nadja == time_to_graduation_nadja3,
        "cache function returns different results than evaluation function. joe: " +
        util::Format()( time_to_graduation_oracle( joe)) + " != " + util::Format()( time_to_graduation_joe) +
        " nadja: " + util::Format()( time_to_graduation_nadja) + " != " + util::Format()( time_to_graduation_nadja3)
      );

      // check cache size
      BCL_ExampleCheck( cache_function.GetCacheSize(), 2);

      // check that TimeToGraduation::operator() was called again now that joe's SetAge function was called
      BCL_ExampleIndirectCheck
      (
        sp_time_to_graduation->GetNumberTimesOperatorWasCalled(),
        5,
        "Cache function recalculation upon signal"
      );

      // cache_function_no_invalidate did not invalidate joe's result, and will return presently the incorrect result
      const int time_to_graduation_joe3_wrong( cache_function_no_invalidate( joe));
      const int time_to_graduation_nadja3_wrong( cache_function_no_invalidate( nadja));

      // check that the results are incorrect
      BCL_Example_Check
      (
        time_to_graduation_joe1_wrong == time_to_graduation_joe3_wrong
        && time_to_graduation_joe1_wrong != time_to_graduation_joe3
        && time_to_graduation_nadja1_wrong == time_to_graduation_nadja3_wrong,
        "cache function without invalidation does not return the same result, after age was set: " +
        util::Format()( time_to_graduation_joe1_wrong) + " != " + util::Format()( time_to_graduation_joe3_wrong)
      );

      // check that TimeToGraduation::operator() was not called again, since joe and nadja are the same age as the last
      // time the function was called
      BCL_ExampleIndirectCheck
      (
        sp_time_to_graduation->GetNumberTimesOperatorWasCalled(),
        5,
        "Cache function only calls the evaluation function if it receives new arguments"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // write to file
      WriteBCLObject( cache_function);
      // read into empty cache function
      math::FunctionCached< Student, int> cache_function_read
      (
        ( util::ShPtr< math::FunctionInterfaceSerializable< Student, int> >()),
        &Student::GetDestructorSignal
      );
      ReadBCLObject( cache_function_read);

      BCL_Example_Check
      (
        cache_function_read.GetFunction().IsDefined()
        && cache_function( joe) == time_to_graduation_joe3,
        "copied cache function has undefined function or returns a different result for student joe"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMathFunctionCached

  const util::SiPtr< const util::ObjectInterface> ExampleMathFunctionCached::TimeToGraduation::s_Instance
  (
    GetObjectInstances().AddInstance( new ExampleMathFunctionCached::TimeToGraduation())
  );

  const ExampleClass::EnumType ExampleMathFunctionCached::s_Instance
  (
    GetExamples().AddEnum( ExampleMathFunctionCached())
  );

} // namespace bcl
