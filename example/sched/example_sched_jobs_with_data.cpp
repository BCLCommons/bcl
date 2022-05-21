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
#include "sched/bcl_sched_binary_function_job_with_data.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_tertiary_function_job_with_data.h"
#include "sched/bcl_sched_thunk_job.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_time.h"

namespace bcl
{
  //! @brief Make sure JobInterface works with globals.
  //! @param DOUBLE a dummy double
  //! @param RESULT a dummy result
  void GlobalFunction( const double &DOUBLE, const int &RESULT)
  {
    BCL_MessageStd( "GlobalFunction just ran!");
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_sched_jobs_with_data.cpp
  //!
  //! @author riddeljs
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleSchedJobsWithData :
    public ExampleInterface
  {
  public:

    ExampleSchedJobsWithData *Clone() const
    {
      return new ExampleSchedJobsWithData( *this);
    }

    class Convert
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      Convert *Clone() const
      {
        return new Convert( *this);
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
    // operations //
    ///////////////

      //! @brief Test tertiary_function_job_with_data and throw in some constant reference parameters
      //! @param DOUBLE1 first double
      //! @param DOUBLE2 second double, we'll be summing these two.
      //! @param RESULT the sum of the doubles converted to an integer
      void TertiaryPerformConversion( const double &DOUBLE1, const double &DOUBLE2, int &RESULT)
      {
        BCL_MessageStd( "TertiaryPerformConversion just ran!");
        util::Time::Delay( util::Time( 1, 100000));
        RESULT = int( DOUBLE1 + DOUBLE2);
      }

      //! @brief Test tertiary_function_job_with_data and throw in some constant reference parameters
      //! @param DOUBLE1 first double
      //! @param DOUBLE2 second double, we'll be summing these two.
      //! @return the sum of the doubles converted to an integer
      int BinaryPerformConversionAndReturn( const double &DOUBLE1, const double &DOUBLE2)
      {
        BCL_MessageStd( "BinaryPerformConversionAndReturn just ran!");
        util::Time::Delay( util::Time( 1, 100000));
        return int( DOUBLE1 + DOUBLE2);
      }

      //! @brief Test binary_function_job_with_data and make sure JobInterface works with static members.
      //! @param DOUBLE the double to multiply by 2.0
      //! @param RESULT The integer representing conversion of 2.0 times DOUBLE
      static void BinaryPerformConversion( double &DOUBLE, int &RESULT)
      {
        BCL_MessageStd( "BinaryPerformConversion, doubling input, just ran!");
        util::Time::Delay( util::Time( 1, 100000));
        RESULT = int( 2.0 * DOUBLE);
      }

      //! @brief Test thunk_job and make sure JobInterface works with queries.
      void ThunkStringOutput() const
      {
        BCL_MessageStd( "ThunkStringOutput has outputted its string, this!");
      }

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM << "double to int conversion function job" << '\n';
      }

    }; // class Convert

    //a brand new class
    class Convert2
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      Convert2 *Clone() const
      {
        return new Convert2( *this);
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
    // operations //
    ///////////////

      //! @brief Test binary_function_job_with_data, function from different class
      //! @param SCHEDULER just to show JobInterface can handle BCL classes
      //! @param RESULT an integer representing 3.2 converted into an integer
      void BinaryPerformConversion( sched::SchedulerInterface &SCHEDULER, int &RESULT)
      {
        BCL_MessageStd( "BinaryPerformConversion, converting 3.2, just ran!");
        util::Time::Delay( util::Time( 1, 100000));
        RESULT = int( 3.2);
      }

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM << "double to int conversion function job" << '\n';
      }

    }; // class Convert2

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct the Convert objects;
      Convert new_convert;
      Convert2 new_convert2;

      // now construct the jobs, so prepare parameters for constructor
      // group ids don't matter yet, but will likely be used for priorities later on
      const size_t convert_group_id1( 0);

      double start_double1( 5.6);
      double start_double2( 2.2);

      int RESULT1( 9);
      int RESULT2( 18);
      int RESULT3( 27);
      int RESULT4( 36);

      // construct four jobs: a thunk, two binary functions, and a tertiary function.
      // The last a binary function from a different class

      // thunk. nothing extra is needed even though it's a query, constructor will realize this.
      sched::ThunkJob< Convert, void> ThunkConvertJob
      (
        convert_group_id1, new_convert,
        &Convert::ThunkStringOutput,
        sched::JobInterface::e_READY,
        NULL
      );

      // first binary job.  Though it is static and the syntax doesn't really need the
      // class instance, we pass it anyway since otherwise we would need it
      sched::BinaryFunctionJobWithData< double, int, void, Convert> BinaryConvertJob
      (
        convert_group_id1, new_convert,
        &Convert::BinaryPerformConversion,
        start_double1, RESULT2, sched::JobInterface::e_READY,
        NULL
      );

      // tertiary job. Note that const appears in the template arguments
      sched::TertiaryFunctionJobWithData< const double, const double, int, void, Convert> TertiaryConvertJob
      ( convert_group_id1, new_convert,
        &Convert::TertiaryPerformConversion,
        start_double1, start_double2, RESULT1, sched::JobInterface::e_READY, NULL
      );

      // another binary job
      sched::BinaryFunctionJobWithData< sched::SchedulerInterface, int, void, Convert2> BinaryConvertJob2
      (
        convert_group_id1, new_convert2,
        &Convert2::BinaryPerformConversion,
        sched::GetScheduler(), RESULT3, sched::JobInterface::e_READY,
        NULL
      );

      // first binary job. Though it is static and the syntax doesn't really need the
      // class instance, we pass it anyway since otherwise we would need it
      sched::BinaryFunctionJobWithData< const double, const int, void, Convert> GlobalJob
      (
        convert_group_id1, new_convert,
        &GlobalFunction,
        start_double1, RESULT2, sched::JobInterface::e_READY,
        NULL
      );

      // a binary job that returns a value
      sched::BinaryFunctionJobWithData< const double, const double, int, Convert> ReturnJob
      (
        convert_group_id1, new_convert,
        &Convert::BinaryPerformConversionAndReturn,
        start_double1, start_double2, sched::JobInterface::e_READY,
        &RESULT4
      );

      // now have an array of job interfaces with these five jobs
      // having an array or list of job interfaces is useful to keep track of jobs and statuses
      sched::JobInterface *Schedule[ 6];

      Schedule[0] = &ThunkConvertJob;
      Schedule[1] = &BinaryConvertJob;
      Schedule[2] = &TertiaryConvertJob;
      Schedule[3] = &BinaryConvertJob2;
      Schedule[4] = &GlobalJob;
      Schedule[5] = &ReturnJob;

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

      // run all of the jobs
      Schedule[0]->Run();
      Schedule[1]->Run();
      Schedule[2]->Run();
      Schedule[3]->Run();
      Schedule[4]->Run();
      Schedule[5]->Run();

      // verify the results of the jobs that changed any of the RESULT variables

      // first the tertiary job
      // because we are indirectly checking that the job ran, we use BCL_ExampleIndirectCheck and then
      // add a helpful message to indicate what we are actually testing
      BCL_ExampleIndirectCheck
      (
        RESULT1,
        int( start_double1 + start_double2),
        "running a 3-input job (TertiaryPerformConversion) to convert the sum of two doubles into an integer"
      );

      // check the binary job that had a return value
      BCL_ExampleIndirectCheck
      (
        RESULT4,
        int( start_double1 + start_double2),
        "running a 2-input job (BinaryPerformConversionAndReturn) that returns the sum of two doubles as an integer"
      );

      //now the binary job from Convert
      BCL_ExampleIndirectCheck
      (
        RESULT2,
        int( 2.0 * start_double1),
        "running a 2-input job (BinaryPerformConversion) to convert the sum of two doubles into an integer"
      );

      // check the binary job from Convert2
      BCL_ExampleIndirectCheck
      (
        RESULT3,
        int( 3.2),
        "running a 2-input job (BinaryPerformConversion) to convert 3.2 into an integer"
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType ExampleSchedJobsWithData_Instance;

  }; //end ExampleSchedJobsWithData

  const ExampleClass::EnumType ExampleSchedJobsWithData::ExampleSchedJobsWithData_Instance
  (
    GetExamples().AddEnum( ExampleSchedJobsWithData())
  );

} // namespace bcl
