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
#include "command/bcl_command_parameter.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_parameter_check_default.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_fixed_line_width_writer.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Parameter::s_Instance
    (
      GetObjectInstances().AddInstance( new Parameter())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Parameter::Parameter() :
      m_Name(),
      m_Description(),
      m_ParameterCheck( new ParameterCheckDefault()),
      m_WasSetInCommandline( false),
      m_DefaultGiven( false)
    {
    }

    //! @brief construct from description string, parameter check function
    //! @param NAME name of parameter
    //! @param DESCRIPTION description string of parameter
    Parameter::Parameter
    (
      const std::string &NAME,
      const std::string &DESCRIPTION
    ) :
      m_Name( NAME),
      m_Description( util::TrimString( DESCRIPTION)),
      m_ParameterCheck( new ParameterCheckDefault()),
      m_WasSetInCommandline( false),
      m_DefaultGiven( false)
    {
    }

    //! @brief construct from description string, parameter check function
    //! @param NAME name of parameter
    //! @param DESCRIPTION description string of parameter
    //! @param PARAMETER_CHECK ParameterCheckInterface derived class to be used to check parameter
    Parameter::Parameter
    (
      const std::string &NAME,
      const std::string &DESCRIPTION,
      const ParameterCheckInterface &PARAMETER_CHECK
    ) :
      m_Name( NAME),
      m_Description( util::TrimString( DESCRIPTION)),
      m_ParameterCheck( PARAMETER_CHECK.Clone()),
      m_WasSetInCommandline( false),
      m_DefaultGiven( false)
    {
    }

    //! @brief construct from name, description string and default parameter string
    //! @param NAME name of parameter
    //! @param DESCRIPTION description string of parameter
    //! @param DEFAULT_PARAMETER Default string value of the parameter
    Parameter::Parameter
    (
      const std::string &NAME,
      const std::string &DESCRIPTION,
      const std::string &DEFAULT_PARAMETER
    ) :
      m_Name( NAME),
      m_Description( util::TrimString( DESCRIPTION)),
      m_ParameterCheck( new ParameterCheckDefault()),
      m_Parameter( DEFAULT_PARAMETER),
      m_WasSetInCommandline( false),
      m_DefaultGiven( true),
      m_DefaultParameter( DEFAULT_PARAMETER)
    {
    }

    //! @brief construct from description string, parameter check, and default parameter string
    //! @param NAME name of parameter
    //! @param DESCRIPTION description string of parameter
    //! @param PARAMETER_CHECK ParameterCheckInterface derived class to be used to check parameter
    //! @param DEFAULT_PARAMETER Default string value of the parameter
    Parameter::Parameter
    (
      const std::string &NAME,
      const std::string &DESCRIPTION,
      const ParameterCheckInterface &PARAMETER_CHECK,
      const std::string &DEFAULT_PARAMETER
    ) :
      m_Name( NAME),
      m_Description( util::TrimString( DESCRIPTION)),
      m_ParameterCheck( PARAMETER_CHECK.Clone()),
      m_Parameter( DEFAULT_PARAMETER),
      m_WasSetInCommandline( false),
      m_DefaultGiven( true),
      m_DefaultParameter( DEFAULT_PARAMETER)
    {
    }

    //! @brief copy constructor
    //! @return pointer to new Parameter object
    Parameter *Parameter::Clone() const
    {
      return new Parameter( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return std::string - the class name as const ref std::string
    const std::string &Parameter::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns Name
    //! @return std::string - name
    const std::string &Parameter::GetName() const
    {
      return m_Name;
    }

    //! @brief set the name
    //! @param NAME name for the parameter
    void Parameter::SetName( const std::string &NAME)
    {
      m_Name = NAME;
    }

    //! @brief returns description
    //! @return std::string-description
    const std::string &Parameter::GetDescription() const
    {
      return m_Description;
    }

    //! @brief checks if PARAMETER is allowed and returns
    //! @param PARAMETER the parameter to set
    //! @param ERROR_STREAM the stream to which errors are written
    //! @return if the parameter was set (parameters are not set if they are not allowed parameters)
    bool Parameter::SetParameter( const std::string &PARAMETER, std::ostream &ERROR_STREAM)
    {
      // get the result of the parameter check
      io::ValidationResult result( m_ParameterCheck->Check( PARAMETER, m_Name, ERROR_STREAM));

      if( !result.IsInvalid())
      {
        // clean the parameter, if it was
        m_Parameter = PARAMETER;
        m_WasSetInCommandline = true;
      }
      return result.IsAllowed();
    } // Parameter::SetParameter

    //! @brief returns m_Parameter
    //! @return value of the parameter
    const std::string &Parameter::GetValue() const
    {
      return m_Parameter;
    }

    //! @brief returns m_DefaultParameter
    //! @return value of the default parameter
    const std::string &Parameter::GetDefaultValue() const
    {
      return m_DefaultParameter;
    }

    //! @brief set the default parameter
    //! @param PARAMETER default value for the parameter
    void Parameter::SetDefaultParameter( const std::string &PARAMETER)
    {
      if( !m_WasSetInCommandline)
      {
        m_DefaultParameter = PARAMETER;
        m_Parameter        = PARAMETER;
        m_DefaultGiven     = true;
      }
    }

    //! @brief returns if parameter was set in the command line
    //! @return m_WasSetInCommandline
    bool Parameter::GetWasSetInCommandLine() const
    {
      return m_WasSetInCommandline;
    }

    //! @brief return whether default value was given
    //! @return bool - m_DefaultGiven
    bool Parameter::GetWasDefaultGiven() const
    {
      return m_DefaultGiven;
    }

    //! @brief acces to the parameter check
    //! @return the ParameterCheck
    const ParameterCheckInterface &Parameter::GetParameterCheck() const
    {
      return *m_ParameterCheck;
    }

    //! @brief set the parameter check
    //! @param CHECK the new parameter check
    void Parameter::SetParameterCheck( const ParameterCheckInterface &CHECK)
    {
      m_ParameterCheck = util::ShPtr< ParameterCheckInterface>( CHECK.Clone());
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns if PARAMETER is allowed, i.e. the object behind the ParameterCheckInterface allows it
    //! @param PARAMETER the parameter to check
    //! @param ERROR_STREAM the stream to which errors are written
    //! @return returns true if the parameter is allowed, false otherwise
    bool Parameter::IsAllowedParameter( const std::string &PARAMETER, std::ostream &ERROR_STREAM) const
    {
      return m_ParameterCheck->IsAllowedParameter( PARAMETER, m_Name, ERROR_STREAM);
    }

    //! @brief reset into original state
    void Parameter::Reset()
    {
      m_WasSetInCommandline = false;

      // default was given, set the parameter to it
      if( m_DefaultGiven)
      {
        m_Parameter = m_DefaultParameter;

        // clean the default parameter.  This is necessary after calling reset since the clean function may have
        // runtime dependencies
        m_ParameterCheck->Clean( m_Parameter);
      }
      // otherwise clear
      else
      {
        m_Parameter.clear();
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes the help for the commandline
    //! @param OSTREAM the stream to which the help is written to
    //! @param INDENT number of indentations
    //! @return the given stream to which the help was written to
    std::ostream &Parameter::WriteHelp( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::FixedLineWidthWriter writer;
      writer.SetBclIndent( INDENT);

      if( !m_Name.empty())
      {
        // write name
        writer << '<' << m_Name << "> ";
      }

      if( !m_Description.empty())
      {
        // write description
        if( m_Description.size() < 2 * writer.GetRemainingSpaceOnLine())
        {
          writer.SetIndent( writer.GetLinePosition());
          writer << m_Description << ", ";
        }
        else
        {
          // just indent by 4 and write the description; it is too long to fit in a fully-indented line
          writer.AddIndent( 4);
          writer << m_Description << ", ";
        }
      }
      else
      {
        writer.AddIndent( 2);
      }

      // if default parameter was given, also write this
      if( !m_WasSetInCommandline && m_DefaultGiven)
      {
        if( m_Parameter.empty())
        {
          writer << "optional, ";
        }
        else
        {
          writer.WriteOnOneLine( "default: \"" + m_Parameter + "\", ");
        }
      }
      writer.PopIndent();

      // write ParameterCheck
      std::ostringstream output;
      m_ParameterCheck->WriteHelp( output, 0);
      const std::string parameter_check( util::TrimString( output.str()));

      // parameter checks that have new lines look best with an extra autoindent 4
      if( parameter_check.find( '\n') != std::string::npos)
      {
        writer.SetAutoIndentExtra( 4);
      }

      // write parameter check position;
      if( parameter_check.size() < writer.GetRemainingSpaceOnLine())
      {
        // short parameter check, write it on the same line
        writer.SetIndent( writer.GetLinePosition());
      }
      else if( parameter_check.size() < writer.GetEffectiveLineWidth() - 2)
      {
        // just indent by 2 and write the parameter check on the next line
        // it is too long to fit in a fully-indented line
        writer.AddIndent( 2);
        writer.NewLine();
      }
      else if( parameter_check.size() < 2 * writer.GetRemainingSpaceOnLine())
      {
        // parameter check will look best with its own block
        writer.SetIndent( writer.GetLinePosition());
      }
      else
      {
        // long parameter check, just write with a standard indent of 2, starting on the same
        // line as the description
        writer.AddIndent( 2);
      }

      // write out the parameter check
      writer << parameter_check;
      const std::string full_help( writer.String());
      // write the whole description, skipping terminal commas, spaces, and new lines
      OSTREAM << full_help.substr( 0, full_help.find_last_not_of( ", \n") + 1) << '\n';

      // end
      return OSTREAM;
    } // Parameter::WriteHelp

    //! @brief writes the user provided commandline
    //! @param OSTREAM the stream to which the commandline is written
    //! @return given ostream to which the commandline was written
    std::ostream &Parameter::WriteUserCommand( std::ostream &OSTREAM) const
    {
      // write description
      OSTREAM << "<" << m_Name << "> " << m_Parameter;

      if( !m_WasSetInCommandline && m_DefaultGiven)
      {
        OSTREAM << " (default)";
      }

      OSTREAM << '\n';

      // end
      return OSTREAM;
    } // Parameter::WriteUserCommand

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Parameter::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Name               , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Description        , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ParameterCheck     , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Parameter          , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_WasSetInCommandline, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DefaultGiven       , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DefaultParameter   , OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Parameter::Read( std::istream &ISTREAM)
    {
      //read members
      io::Serialize::Read( m_Name               , ISTREAM);
      io::Serialize::Read( m_Description        , ISTREAM);
      io::Serialize::Read( m_ParameterCheck     , ISTREAM);
      io::Serialize::Read( m_Parameter          , ISTREAM);
      io::Serialize::Read( m_WasSetInCommandline, ISTREAM);
      io::Serialize::Read( m_DefaultGiven       , ISTREAM);
      io::Serialize::Read( m_DefaultParameter   , ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief StringVectorFromFilenameParameter gives all the strings in a file as a Vector
    //! @param PARAMETER the name of the file
    //! @return return a Vector which has all of the strings contained within the file denoted by "FLAG"
    storage::Vector< std::string> StringVectorFromFilenameParameter( const ParameterInterface &PARAMETER)
    {
      // initialize write and read stream object
      io::IFStream read;

      // open pdb list file
      read.open( PARAMETER.GetValue().c_str());

      // make sure the file opened
      BCL_Assert( read.is_open(), "unable to open file: " + PARAMETER.GetValue());

      // create "string_vector" initialize with vector of strings returned by function util::StringListFromIStream
      const storage::Vector< std::string> string_vector( util::StringListFromIStream( read));

      // close and clear "read"
      io::File::CloseClearFStream( read);

      // return Vector filled with the strings in the pdb list
      return string_vector;
    }

  } // namespace command
} // namespace bcl
