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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class

// includes from bcl - sorted alphabetically
#include "app/bcl_app_apps.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "io/bcl_io_file.h"
#include "util/bcl_util_cpu_benchmark_whetstone.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CPUBenchmark
    //! @brief Class converts pdb file into bcl-style pdb file
    //! @details Class reads a given pdb to and converts it to a bcl-style pdb with idealized SSEs,
    //! has options to print out the fasta and bcl style pdb output for the given pdb
    //!
    //! @author woetzen
    //! @date 03/24/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CPUBenchmark :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      //! flag for benchmark duration
      util::ShPtr< command::FlagInterface> m_DurationFlag;

      //! flag for output prefix
      util::ShPtr< command::FlagInterface> m_OutputPrefixFlag;

      // instance of benchmark classes
      mutable util::CPUBenchmarkWhetstone< double> m_WhetstoneDouble;
      mutable util::CPUBenchmarkWhetstone< float> m_WhetstoneFloat;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      CPUBenchmark();

    public:

      //! @brief Clone function
      //! @return pointer to new CPUBenchmark
      CPUBenchmark *Clone() const
      {
        return new CPUBenchmark( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

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
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      static const ApplicationType CPUBenchmark_Instance;

    }; // class CPUBenchmark

    //! @brief the Main function
    //! @return error code - 0 for success
    int CPUBenchmark::Main() const
    {
      // get the duration
      const util::Time duration( m_DurationFlag->GetFirstParameter()->GetNumericalValue< size_t>(), 0);

      // run float benchmark
      BCL_MessageStd( "running floating precision benchmark with duration: " + duration.GetTimeAsHourMinuteSecondMilliSeconds());
      m_WhetstoneFloat.RunBenchmark( duration);

      // run double benchmark
      BCL_MessageStd( "running double precision benchmark with duration: " + duration.GetTimeAsHourMinuteSecondMilliSeconds());
      m_WhetstoneDouble.RunBenchmark( duration);

      // writing results
      io::OFStream write;
      const std::string float_out_name( m_OutputPrefixFlag->GetFirstParameter()->GetValue() + "whetstone_float.table");
      io::File::MustOpenOFStream( write, float_out_name);
      m_WhetstoneFloat.GetResultTable().WriteFormatted( write);
      io::File::CloseClearFStream( write);

      const std::string double_out_name( m_OutputPrefixFlag->GetFirstParameter()->GetValue() + "whetstone_double.table");
      io::File::MustOpenOFStream( write, double_out_name);
      m_WhetstoneDouble.GetResultTable().WriteFormatted( write);
      io::File::CloseClearFStream( write);

      return 0;
    } // end Main

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &CPUBenchmark::GetReadMe() const
    {
      static const std::string s_readme
      (
        "This application performs a Whetstone cpu benchmark, Given a duration (recommended is 100 seconds) it will "
        "print a table with the MOPS (10^6 operations per second) for different tests. Please read "
        "http://www.roylongbottom.org.uk/whetstone.pdf"
      );

      return s_readme;
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> CPUBenchmark::InitializeCommand() const
    {
      // initialize ShPtr to a new command object
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // duration and output
      sp_cmd->AddFlag( m_DurationFlag);
      sp_cmd->AddFlag( m_OutputPrefixFlag);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! default constructor
    CPUBenchmark::CPUBenchmark() :
      m_DurationFlag
      (
        new command::FlagStatic
        (
          "duration",
          "pass the duration for the benchmark",
          command::Parameter
          (
            "duration_time",
            "duration time in seconds - the longer, the more accurate the result will be",
            command::ParameterCheckRanged< size_t>( 1, 200),
            "100"
          )
        )
      ),
      m_OutputPrefixFlag
      (
        new command::FlagStatic
        (
          "output_prefix",
          "prefix for output table",
          command::Parameter
          (
            "prefix",
            "output prefix",
            ""
          )
        )
      )
    {
    }

    const ApplicationType CPUBenchmark::CPUBenchmark_Instance
    (
      GetAppGroups().AddAppToGroup( new CPUBenchmark(), GetAppGroups().e_Utility)
    );

  } // namespace app
} // namespace bcl
