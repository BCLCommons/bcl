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
#include "sched/bcl_sched_sum_function.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serializer.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_sched_sum_function.cpp
  //! @brief TODO: add an detailed description for this example
  //!
  //! @author woetzen
  //! @date Aug 21, 2011
  //! @remarks status incomplete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleSchedSumFunction :
    public ExampleInterface
  {
    //! @class Sleeper
    class Sleeper :
      public math::FunctionInterfaceSerializable< util::Time, double>
    {

    //////////
    // data //
    //////////

      double m_Value;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      Sleeper( const double VALUE) :
        m_Value( VALUE)
      {
      }

      //! @brief Clone function
      //! @return pointer to new Sleeper
      Sleeper *Clone() const
      {
        return new Sleeper( *this);
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

      //! @brief returns alias for referencing this class over the command line
      const std::string &GetAlias() const
      {
        static const std::string s_alias( "Sleeper");
        return s_alias;
      }

    ///////////////
    // operators //
    ///////////////

      //! @return random number after sleeping
      //! @param DURATION time to sleep
      //! @return random number
      double operator()( const util::Time &DURATION) const
      {
        util::Time::Delay( DURATION);
        return m_Value;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      io::Serializer GetSerializer() const
      {
        io::Serializer ser;
        ser.SetClassDescription( "class that returns a certain value after sleeping a specified time");
        ser.AddInitializer( "", "value to return", io::Serialization::GetAgent( &m_Value));
        return ser;
      }
    }; // class Sleeper

  public:

    ExampleSchedSumFunction *Clone() const
    {
      return new ExampleSchedSumFunction( *this);
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
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct from scheme
      sched::SumFunction< util::Time, double> sum_function;

    /////////////////
    // data access //
    /////////////////

    ////////////////
    // operations //
    ////////////////

      util::Enumerated< math::FunctionInterfaceSerializable< util::Time, double> >::AddInstance( new Sleeper( 0.0));
      sum_function.NewOperand( Sleeper( 1.0), -1.0);
      sum_function.NewOperand( Sleeper( 2.0), 1.0);
      sum_function.NewOperand( Sleeper( 3.0), -1.0);
      sum_function.NewOperand( Sleeper( 4.0), 1.0);

    ///////////////
    // operators //
    ///////////////

      const util::Time duration( 0, 50000);
      const double expected( 2.0);
      {
        util::Stopwatch stop( "sched sum function operator", util::Message::e_Standard, true);
        const double result( sum_function( duration));
        BCL_ExampleCheck( result, expected);
      }

    //////////////////////
    // input and output //
    //////////////////////

      // read and write object
      WriteBCLObject( sum_function);
      sched::SumFunction< util::Time, double> sum_function_read;
      ReadBCLObject( sum_function_read);

      const double result_read( sum_function_read( duration));
      BCL_ExampleCheck( result_read, expected);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSchedSumFunction

  const ExampleClass::EnumType ExampleSchedSumFunction::s_Instance
  (
    GetExamples().AddEnum( ExampleSchedSumFunction())
  );

} // namespace bcl
