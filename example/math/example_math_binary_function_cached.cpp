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
#include "math/bcl_math_binary_function_cached.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_binary_function_cached.cpp
  //!
  //! @author woetzen, mendenjl
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathBinaryFunctionCached :
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
      double         m_HomeWorkGrade;

      Student( const std::string &NAME, const double HOME_WORK_GRADE) :
        m_Name( NAME),
        m_HomeWorkGrade( HOME_WORK_GRADE)
      {
      }

      ~Student()
      {
        m_DestructorSignal.Emit( *this);
      }

      void SetHomeWorkGrade( const double HOME_WORK_GRADE)
      {
        m_HomeWorkGrade = HOME_WORK_GRADE;
        m_OnChangeSignal.Emit( *this);
      }

    }; // class Student

    class Teacher
    {
      mutable signal::Signal1< const Teacher&> m_OnChangeSignal;
      mutable signal::Signal1< const Teacher&> m_DestructorSignal;

    public:

      signal::Signal1< const Teacher &> &GetOnChangeSignal() const
      {
        return m_OnChangeSignal;
      }

      signal::Signal1< const Teacher &> &GetDestructorSignal() const
      {
        return m_DestructorSignal;
      }

      std::string m_Name;

      Teacher( const std::string &NAME) :
        m_Name( NAME)
      {
      }

      ~Teacher()
      {
        m_DestructorSignal.Emit( *this);
      }

      void SetName( const std::string &NAME)
      {
        m_Name = NAME;
        m_OnChangeSignal.Emit( *this);
      }

    }; // class Teacher

    class GradeExam :
      public math::BinaryFunctionInterfaceSerializable< Teacher, Student, double>
    {
      mutable size_t m_UseCounter; //!< Tracks how many times this object was called
    public:
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      GradeExam() :
        m_UseCounter( 0)
      {
      }

      //! @brief
      GradeExam *Clone() const
      {
        return new GradeExam( *this);
      }

      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      const size_t &GetUseCount() const
      {
        return m_UseCounter;
      }

      double operator()( const Teacher &TEACHER, const Student &STUDENT) const
      {
        if( TEACHER.m_Name == "Ms. Nicegrade")
        {
          return 0.25 * ( 3.5) + 0.75 * STUDENT.m_HomeWorkGrade;
        }
        if( TEACHER.m_Name == "Mr. Strictgrade")
        {
          return 0.75 * ( 3.5) + 0.25 * STUDENT.m_HomeWorkGrade;
        }
        m_UseCounter++;
        return 0.5 * ( 3.5) + 0.5 * STUDENT.m_HomeWorkGrade;
      }

      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    }; // class GradeExam

    class GradeAverage :
      public math::BinaryFunctionInterfaceSerializable< Student, Student, double>
    {
      mutable size_t m_UseCounter; //!< Tracks how many times this object was called

    public:
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      GradeAverage() :
        m_UseCounter( 0)
      {
      }

      //! @brief
      GradeAverage *Clone() const
      {
        return new GradeAverage( *this);
      }

      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      const size_t &GetUseCount() const
      {
        return m_UseCounter;
      }

      bool IsSymmetric() const
      {
        return true;
      }

      double operator()( const Student &STUDENT1, const Student &STUDENT2) const
      {
        m_UseCounter++;
        return ( STUDENT1.m_HomeWorkGrade + STUDENT2.m_HomeWorkGrade) / 2.0;
      }

      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    }; // class GradeAverage

    class GradeDifference :
      public math::BinaryFunctionInterfaceSerializable< Student, Student, double>
    {
      mutable size_t m_UseCounter; //!< Tracks how many times this object was called

    public:
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      GradeDifference() :
        m_UseCounter( 0)
      {
      }

      //! @brief
      GradeDifference *Clone() const
      {
        return new GradeDifference( *this);
      }

      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      const size_t &GetUseCount() const
      {
        return m_UseCounter;
      }

      bool IsSymmetric() const
      {
        return false;
      }

      double operator()( const Student &STUDENT1, const Student &STUDENT2) const
      {
        m_UseCounter++;
        return STUDENT1.m_HomeWorkGrade - STUDENT2.m_HomeWorkGrade;
      }

      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    }; // class GradeAverage

    ExampleMathBinaryFunctionCached *Clone() const
    {
      return new ExampleMathBinaryFunctionCached( *this);
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

    //! @brief run tests using a particular set of objects
    //! @param GRADER the grader that should be being used internally by the cache
    //! @param CACHE_FUNC the cache function
    //! @param TEACHER the teacher (1st arg of binary function)
    //! @param STUDENT the student (2nd arg of binary function)
    //! @param SHOULD_BE_CACHED whether to expect the results to already be in the cache
    void RunGradeExamChecks
    (
      const util::ShPtr< GradeExam> &GRADER,
      const math::BinaryFunctionCached< Teacher, Student, double> &CACHE_FUNC,
      const Teacher &TEACHER,
      const Student &STUDENT,
      const bool &SHOULD_BE_CACHED
    ) const
    {
      // compute the grading using the non-cache function
      const double actual_grade_value( GRADER->operator()( TEACHER, STUDENT));

      // get the initial number of function calls to the grader
      const size_t initial_gradings( GRADER->GetUseCount());

      // check that the cached function returns the same value
      BCL_ExampleIndirectCheck( CACHE_FUNC( TEACHER, STUDENT), actual_grade_value, "Cache function of GradeExam");

      BCL_ExampleIndirectCheck( CACHE_FUNC.GetCacheSize(), 1, "Cache size of GradeExam");

      // check that the function was called the expected number of times
      BCL_ExampleIndirectCheck
      (
        GRADER->GetUseCount(),
        initial_gradings + size_t( !SHOULD_BE_CACHED),
        "Calls to GradExam::operator()"
      );
    }

    //! @brief run tests that symmetry handling is working properly using a set of objects
    //! @param GRADER the grader that should be being used internally by the cache
    //! @param CACHE_FUNC the cache function
    //! @param STUDENT_A, STUDENT_B the students
    //! @param SHOULD_BE_IN_CACHE whether results should already be in cache
    template< typename t_DataType>
    void RunSymmetryTests
    (
      const util::ShPtr< t_DataType> &GRADER,
      const math::BinaryFunctionCached< Student, Student, double> &CACHE_FUNC,
      const Student &STUDENT_A,
      const Student &STUDENT_B,
      const bool &SHOULD_BE_IN_CACHE
    ) const
    {
      // compute the grading using the non-cache function
      const double actual_grade_value_a_b( GRADER->operator()( STUDENT_A, STUDENT_B));
      const double actual_grade_value_b_a( GRADER->operator()( STUDENT_B, STUDENT_A));

      // get the initial number of function calls to the grader
      const size_t initial_gradings( GRADER->GetUseCount());

      // if the cache function has the correct symmetry, check that the results are as expected
      if( CACHE_FUNC.IsSymmetric() == GRADER->IsSymmetric())
      {
        // check that the cached function returns the same value
        BCL_ExampleIndirectCheck
        (
          CACHE_FUNC( STUDENT_A, STUDENT_B),
          actual_grade_value_a_b,
          "Symmetry handling for cache on " + GRADER->GetClassIdentifier()
        );

        BCL_ExampleIndirectCheck
        (
          CACHE_FUNC( STUDENT_B, STUDENT_A),
          actual_grade_value_b_a,
          "Symmetry handling for cache on " + GRADER->GetClassIdentifier()
        );
      }
      else
      {
        // check that the cached function returns the same value
        BCL_ExampleIndirectCheck
        (
          CACHE_FUNC( STUDENT_A, STUDENT_B) == actual_grade_value_a_b
          || CACHE_FUNC( STUDENT_A, STUDENT_B) == actual_grade_value_b_a,
          true,
          "Symmetry handling for cache on " + GRADER->GetClassIdentifier()
        );

        // if the cache function has incorrect symmetry, check that the resutls are the same as the initial call
        BCL_ExampleIndirectCheck
        (
          CACHE_FUNC( STUDENT_B, STUDENT_A),
          CACHE_FUNC( STUDENT_A, STUDENT_B),
          "Symmetry handling for cache on " + GRADER->GetClassIdentifier()
        );
      }

      // determine the expected number of calls of GRADER->operator()
      const size_t expected_cache_size( CACHE_FUNC.IsSymmetric() ? 1 : 2);
      BCL_ExampleIndirectCheck
      (
        CACHE_FUNC.GetCacheSize(),
        expected_cache_size,
        "Symmetry handling for cache on " + GRADER->GetClassIdentifier()
      );
      BCL_ExampleIndirectCheck
      (
        GRADER->GetUseCount(),
        initial_gradings + ( SHOULD_BE_IN_CACHE ? 0 : expected_cache_size),
        "Symmetry handling for cache on " + GRADER->GetClassIdentifier()
      );
    }

    int Run() const
    {
      // sp to evaluation function
      util::ShPtr< GradeExam> sp_evaluate_grade_exam( new GradeExam());
      util::ShPtr< GradeAverage> sp_evaluate_grade_avg( new GradeAverage());
      util::ShPtr< GradeDifference> sp_evaluate_grade_diff( new GradeDifference());

      // students
      Student joe( "Joe", 2.1);
      Student nadja( "Nadja", 1.7);

      // Teacher
      Teacher math_teacher( "Dr. Shmoe");

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // cache function from evaluation function
      math::BinaryFunctionCached< Teacher, Student, double> cache_function_grade_exam
      (
        sp_evaluate_grade_exam, &Teacher::GetDestructorSignal, &Student::GetDestructorSignal
      );
      math::BinaryFunctionCached< Teacher, Student, double> cache_function_grade_exam_no_invalidate
      (
        sp_evaluate_grade_exam, &Teacher::GetDestructorSignal, &Student::GetDestructorSignal
      );

      // identical arguments
      math::BinaryFunctionCached< Student, Student, double> cache_function_grade_avg_sym
      (
        sp_evaluate_grade_avg, &Student::GetDestructorSignal, true
      );
      math::BinaryFunctionCached< Student, Student, double> cache_function_grade_diff
      (
        sp_evaluate_grade_diff, &Student::GetDestructorSignal, false
      );
      math::BinaryFunctionCached< Student, Student, double> cache_function_grade_diff_wrong
      (
        sp_evaluate_grade_diff, &Student::GetDestructorSignal, true
      );

      // copy constructor
      math::BinaryFunctionCached< Teacher, Student, double> copy_cache_function_grade_exam
      (
        cache_function_grade_exam_no_invalidate
      );

      BCL_Example_Check
      (
        copy_cache_function_grade_exam.GetCacheSize() == 0
        && cache_function_grade_exam_no_invalidate.GetFunction() == copy_cache_function_grade_exam.GetFunction(),
        "copied cache function, should have cache of zero, and point to the same function"
      );

      math::BinaryFunctionCached< Student, Student, double> copy_cache_function_grade_diff( cache_function_grade_diff_wrong);
      BCL_Example_Check
      (
        copy_cache_function_grade_diff.GetCacheSize() == 0
        && copy_cache_function_grade_diff.IsSymmetric() == cache_function_grade_diff_wrong.IsSymmetric()
        && copy_cache_function_grade_diff.GetFunction() == cache_function_grade_diff_wrong.GetFunction(),
        "copied cache function, should have cache of zero, and point to the same function and have the symmetry"
      );

      // clone
      util::ShPtr< math::BinaryFunctionCached< Teacher, Student, double> >
        sp_clone_cache_function_grade_exam( cache_function_grade_exam.Clone());

      util::ShPtr< math::BinaryFunctionCached< Student, Student, double> >
        sp_clone_cache_function_grade_diff( cache_function_grade_diff.Clone());

      BCL_ExampleIndirectCheck( sp_evaluate_grade_exam, sp_clone_cache_function_grade_exam->GetFunction(), "Clone");
      BCL_ExampleIndirectCheck( sp_evaluate_grade_diff, sp_clone_cache_function_grade_diff->GetFunction(), "Clone");

    /////////////////
    // data access //
    /////////////////

      // class name
      BCL_ExampleCheck
      (
        ( GetStaticClassName< math::BinaryFunctionCached< Student, Student, double> >()),
        cache_function_grade_diff.GetClassIdentifier()
      );

      // access the function object
      BCL_ExampleCheck( sp_evaluate_grade_exam, cache_function_grade_exam.GetFunction());

      // cache size
      BCL_ExampleCheck( cache_function_grade_exam.GetCacheSize(), 0);

    ////////////////
    // operations //
    ////////////////

      // add the onchange signal to the cache that will invalidate cached results
      BCL_ExampleCheck( cache_function_grade_exam.AddSignalHandlerForArgument( &Student::GetOnChangeSignal), true);
      BCL_ExampleIndirectCheck
      (
        cache_function_grade_exam.AddSignalHandlerForArgument( &Student::GetOnChangeSignal),
        false,
        "should not be able to successfully register the same signal handler twice"
      );

      // add all necessary signal handlers that the cache function has to connect to on the argument in evaluated
      cache_function_grade_exam.AddSignalHandlerForArgument( &Teacher::GetOnChangeSignal);
      cache_function_grade_avg_sym.AddSignalHandlerForArgument( &Student::GetOnChangeSignal);
      cache_function_grade_diff.AddSignalHandlerForArgument( &Student::GetOnChangeSignal);
      cache_function_grade_diff_wrong.AddSignalHandlerForArgument( &Student::GetOnChangeSignal);

    ///////////////
    // operators //
    ///////////////

      // run initial tests
      RunGradeExamChecks( sp_evaluate_grade_exam, cache_function_grade_exam, math_teacher, joe, false);

      // ensure that the caches are not storing mutual data
      RunGradeExamChecks( sp_evaluate_grade_exam, cache_function_grade_exam_no_invalidate, math_teacher, joe, false);

      // test on student-student relationships on intrinsically symmetric function
      RunSymmetryTests( sp_evaluate_grade_avg, cache_function_grade_avg_sym, joe, nadja, false);

      // test on student-student relationships with intrinsically asymmetric function
      RunSymmetryTests( sp_evaluate_grade_diff, cache_function_grade_diff, joe, nadja, false);

      // test on student-student relationships with intrinsically asymmetric function and an incorrectly setup cache
      // to demonstrate that the symmetry is being respected, even though its wrong
      RunSymmetryTests( sp_evaluate_grade_diff, cache_function_grade_diff_wrong, joe, nadja, false);

      // test that results are cached

      // run initial tests
      RunGradeExamChecks( sp_evaluate_grade_exam, cache_function_grade_exam, math_teacher, joe, true);

      // ensure that the caches are not storing mutual data
      RunGradeExamChecks( sp_evaluate_grade_exam, cache_function_grade_exam_no_invalidate, math_teacher, joe, true);

      // test on student-student relationships on intrinsically symmetric function
      RunSymmetryTests( sp_evaluate_grade_avg, cache_function_grade_avg_sym, joe, nadja, true);

      // test on student-student relationships with intrinsically asymmetric function
      RunSymmetryTests( sp_evaluate_grade_diff, cache_function_grade_diff, joe, nadja, true);

      // test on student-student relationships with intrinsically asymmetric function and an incorrectly setup cache
      // to demonstrate that the symmetry is being respected, even though its wrong
      RunSymmetryTests( sp_evaluate_grade_diff, cache_function_grade_diff_wrong, joe, nadja, true);

      const double original_grade( cache_function_grade_exam( math_teacher, joe));

      // change grade
      joe.SetHomeWorkGrade( 1.7);

      // recalculate after change

      // run initial tests
      RunGradeExamChecks( sp_evaluate_grade_exam, cache_function_grade_exam, math_teacher, joe, false);

      // ensure that without adding the signal handler for set homework grade, the results would not have been
      // recalculated
      BCL_ExampleCheck( cache_function_grade_exam_no_invalidate( math_teacher, joe), original_grade);

      // test on student-student relationships on intrinsically symmetric function
      RunSymmetryTests( sp_evaluate_grade_avg, cache_function_grade_avg_sym, joe, nadja, false);

      // test on student-student relationships with intrinsically asymmetric function
      RunSymmetryTests( sp_evaluate_grade_diff, cache_function_grade_diff, joe, nadja, false);

      // test on student-student relationships with intrinsically asymmetric function and an incorrectly setup cache
      // to demonstrate that the symmetry is being respected, even though its wrong
      RunSymmetryTests( sp_evaluate_grade_diff, cache_function_grade_diff_wrong, joe, nadja, false);

    //////////////////////
    // input and output //
    //////////////////////

      // write to file
      WriteBCLObject( cache_function_grade_exam);
      // read into empty cache function
      math::BinaryFunctionCached< Teacher, Student, double> cache_function_grade_exam_read
      (
        ( util::ShPtr< math::BinaryFunctionInterface< Teacher, Student, double> >()),
        &Teacher::GetDestructorSignal,
        &Student::GetDestructorSignal
      );
      ReadBCLObject( cache_function_grade_exam_read);

      BCL_ExampleIndirectAssert
      (
        cache_function_grade_exam_read.GetFunction().IsDefined(),
        true,
        "I/O should read a defined binary cache function"
      );

      BCL_ExampleIndirectCheck
      (
        cache_function_grade_exam_read( math_teacher, joe),
        cache_function_grade_exam( math_teacher, joe),
        "I/O"
      );

      // write to file
      WriteBCLObject( cache_function_grade_avg_sym);
      // read into empty cache function
      math::BinaryFunctionCached< Student, Student, double> cache_function_grade_avg_sym_read
      (
        ( util::ShPtr< math::BinaryFunctionInterface< Student, Student, double> >()),
        &Student::GetDestructorSignal,
        false
      );
      ReadBCLObject( cache_function_grade_avg_sym_read);

      BCL_ExampleIndirectCheck
      (
        cache_function_grade_avg_sym_read.GetFunction().IsDefined() &&
        cache_function_grade_avg_sym_read.IsSymmetric(),
        true,
        "I/O respects symmetry"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMathBinaryFunctionCached

  const util::SiPtr< const util::ObjectInterface> ExampleMathBinaryFunctionCached::GradeExam::s_Instance
  (
    GetObjectInstances().AddInstance( new ExampleMathBinaryFunctionCached::GradeExam())
  );
  const util::SiPtr< const util::ObjectInterface> ExampleMathBinaryFunctionCached::GradeAverage::s_Instance
  (
    GetObjectInstances().AddInstance( new ExampleMathBinaryFunctionCached::GradeAverage())
  );
  const util::SiPtr< const util::ObjectInterface> ExampleMathBinaryFunctionCached::GradeDifference::s_Instance
  (
    GetObjectInstances().AddInstance( new ExampleMathBinaryFunctionCached::GradeDifference())
  );

  const ExampleClass::EnumType ExampleMathBinaryFunctionCached::s_Instance
  (
    GetExamples().AddEnum( ExampleMathBinaryFunctionCached())
  );

} // namespace bcl
