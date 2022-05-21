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
#include "command/bcl_command_flag_static.h"

// includes from bcl - sorted alphabetically
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
    const util::SiPtr< const util::ObjectInterface> FlagStatic::s_Instance( GetObjectInstances().AddInstance( new FlagStatic()));

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FlagStatic::FlagStatic() :
      m_WasSetInCommandline( false),
      m_Signal( NULL)
    {
    }

    //! @brief construct FlagStatic from NAME, DESCRIPTION
    //! @param NAME the name of the flag as a string
    //! @param DESCRIPTION a string description of the flag
    //! @param SIGNAL optional function to call whenever calling set or reset on the flag
    FlagStatic::FlagStatic
    (
      const std::string &NAME,
      const std::string &DESCRIPTION,
      t_Signal SIGNAL
    ) :
      m_Name( NAME),
      m_Description( util::TrimString( DESCRIPTION)),
      m_WasSetInCommandline( false),
      m_Signal( SIGNAL)
    {
    }

    //! @brief construct FlagStatic from NAME, DESCRIPTION and single PARAMETER
    //! @param NAME the name of the flag as a string
    //! @param DESCRIPTION a string description of the flag
    //! @param PARAMETER the parameter associated with this flag
    //! @param SIGNAL optional function to call whenever calling set or reset on the flag
    FlagStatic::FlagStatic
    (
      const std::string &NAME,
      const std::string &DESCRIPTION,
      const ParameterInterface &PARAMETER,
      t_Signal SIGNAL
    ) :
      m_Name( NAME),
      m_Description( util::TrimString( DESCRIPTION)),
      m_WasSetInCommandline( false),
      m_ParameterList( 1, PARAMETER),
      m_Signal( SIGNAL)
    {
    }

    //! @brief construct FlagStatic from NAME, DESCRIPTION and single PARAMETER
    //! @param NAME the name of the flag as a string
    //! @param DESCRIPTION a string description of the flag
    //! @param PARAMETERS the parameters associated with this flag
    //! @param SIGNAL optional function to call whenever calling set or reset on the flag
    FlagStatic::FlagStatic
    (
      const std::string &NAME,
      const std::string &DESCRIPTION,
      const util::ShPtrVector< ParameterInterface> &PARAMETERS,
      t_Signal SIGNAL
    ) :
      m_Name( NAME),
      m_Description( util::TrimString( DESCRIPTION)),
      m_WasSetInCommandline( false),
      m_ParameterList( PARAMETERS),
      m_Signal( SIGNAL)
    {
    }

    //! @brief copy constructor
    //! @return pointer to FlagStatic object
    FlagStatic *FlagStatic::Clone() const
    {
      return new FlagStatic( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief name of class
    //! @return name of class as string
    //! @see GetStaticClassName()
    const std::string &FlagStatic::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns name of flag
    //! @return name of flag as std::string
    const std::string &FlagStatic::GetName() const
    {
      return m_Name;
    }

    //! @brief set name - only possible if there was no name given yet
    //! @param NAME - what you want to name the flag
    void FlagStatic::SetName( const std::string &NAME)
    {
      BCL_Assert( m_Name.empty(), "name was already given");
      m_Name = NAME;
    }

    //! @brief returns description of flag
    //! @return description of flag as string
    const std::string &FlagStatic::GetDescription() const
    {
      return m_Description;
    }

    //! @brief set flag to true if it was false
    void FlagStatic::SetFlag()
    {
      m_WasSetInCommandline = true;
    }

    //! @brief set flag to false if it was true
    void FlagStatic::UnsetFlag()
    {
      m_WasSetInCommandline = false;
    }

    //! @brief reset the flag
    //! @detail resets all internal parameters, then removes all internally held parameters
    void FlagStatic::ResetFlag()
    {
      m_WasSetInCommandline = false;

      // iterate over parameters in current flag and reset them all
      for
      (
        util::ShPtrVector< ParameterInterface>::iterator
          flag_param_itr( m_ParameterList.Begin()),
          flag_param_itr_end( m_ParameterList.End());
        flag_param_itr != flag_param_itr_end;
        ++flag_param_itr
      )
      {
        ( *flag_param_itr)->Reset();
      }

      // if a signal should be emitted, emit it
      if( m_Signal && FlagInterface::ShouldSignal())
      {
        ( *m_Signal)();
      }
    }

    //! @brief returns if m_IsSet in command line
    //! @return whether this flag was set in the command line
    bool FlagStatic::GetFlag() const
    {
      return m_WasSetInCommandline;
    }

    //! @brief returns the function to be called whenever this flag is updated
    //! @return the function to be called whenever this flag is updated
    FlagInterface::t_Signal FlagStatic::GetSignal() const
    {
      return m_Signal;
    }

    //! @brief returns m_ParameterList
    //! @return the parameter list as a ShPtrVector
    const util::ShPtrVector< ParameterInterface> &FlagStatic::GetParameterList() const
    {
      return m_ParameterList;
    }

    //! @brief returns m_ParameterList
    //! @return the parameter list as a util::ShPtrVector< ParameterInterface>
    util::ShPtrVector< ParameterInterface> &FlagStatic::GetParameterList()
    {
      return m_ParameterList;
    }

    //! @brief returns the first parameter
    //! @return the first parameter for this flag as a util::ShPtr< ParameterInterface>
    const util::ShPtr< ParameterInterface> &FlagStatic::GetFirstParameter() const
    {
      BCL_Assert( !m_ParameterList.IsEmpty(), "there is no first parameter for " + m_Name);
      return m_ParameterList.FirstElement();
    }

    //! @brief returns the # of required parameters for this flag
    //! @return the effective # of required parameters for this flag
    size_t FlagStatic::GetNumberRequiredParameters() const
    {
      size_t counter( m_ParameterList.GetSize());

      // find the last parameter without a default value
      for
      (
        util::ShPtrVector< ParameterInterface>::const_reverse_iterator
          itr( m_ParameterList.ReverseBegin()),
          itr_end( m_ParameterList.ReverseEnd());
        itr != itr_end;
        ++itr, --counter
      )
      {
        // if no default was given, return that #
        if( !( *itr)->GetWasDefaultGiven())
        {
          return counter;
        }
      }

      return 0;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief fill the parameter list from storage vector of strings
    //! @param PARAMETER_LIST a storage vector of strings containing the parameters as strings
    //! @param ERROR_STREAM the stream to which errors should be written
    //! @return true if successful, false otherwise
    bool FlagStatic::ReadFromList( const storage::Vector< std::string> &PARAMETER_LIST, std::ostream &ERROR_STREAM)
    {
      // check if flag was already set
      if( GetFlag())
      {
        ERROR_STREAM << "flag: \"-" << GetName() << "\" was already set" << '\n';
        return false;
      }

      // set flag
      SetFlag();

      // check to make sure they didn't pass too many parameters for the given flag
      if( PARAMETER_LIST.GetSize() > GetSize())
      {
        ERROR_STREAM << "flag: \"-" << GetName() << "\" expects " << GetSize() << " parameters but received "
          << PARAMETER_LIST.GetSize() << '\n';
        return false;
      }

      bool success_setting_parameters( true);

      // iterator on begin of given parameters after option
      std::vector< std::string>::const_iterator
        param_itr( PARAMETER_LIST.Begin()),
        param_itr_end( PARAMETER_LIST.End());

      // iterate parallel over parameters in current flag with parameters and given parameters in commandline
      for
      (
        util::ShPtrVector< ParameterInterface>::iterator
          flag_param_itr( m_ParameterList.Begin()),
          flag_param_itr_end( m_ParameterList.End());
        flag_param_itr != flag_param_itr_end && param_itr != param_itr_end;
        ++flag_param_itr, ++param_itr
      )
      {
        // set parameters in static parameter list
        success_setting_parameters &= ( *flag_param_itr)->SetParameter( *param_itr, ERROR_STREAM);
      } // for

      // check if the list of parameters is valid
      success_setting_parameters &= IsValidList( ERROR_STREAM);

      // if a signal should be emitted, emit it
      if( success_setting_parameters && m_Signal && FlagInterface::ShouldSignal())
      {
        ( *m_Signal)();
      }

      // end
      return success_setting_parameters;
    } // ReadFromList

    //! @brief checks that for every parameter was given a default value or it was passed through commandline
    //! @brief or if it is a dynamic list, if it meets the size specifications
    //! @param ERROR_STREAM the stream to which errors should be written
    //! @return true if successful, false otherwise
    bool FlagStatic::IsValidList( std::ostream &ERROR_STREAM) const
    {
      bool valid_list( true);
      size_t counter( 0);
      for
      (
        util::ShPtrVector< ParameterInterface>::const_iterator
          param_itr( m_ParameterList.Begin()),
          param_itr_end( m_ParameterList.End());
        param_itr != param_itr_end;
        ++param_itr, ++counter
      )
      {
        // if neither the parameter was set in command line, nor the parameter has a default
        if( !( *param_itr)->GetWasSetInCommandLine() && !( *param_itr)->GetWasDefaultGiven())
        {
          valid_list = false;
          ERROR_STREAM << "flag: " << m_Name;

          // if the parameter has a name, use it
          if( !( *param_itr)->GetName().empty())
          {
            ERROR_STREAM << ", parameter " << ( *param_itr)->GetName() << ":";
          }
          else if( m_ParameterList.GetSize() > size_t( 1))
          {
            // multiple parameters; this parameter unnamed, use the counter to id the parameter
            ERROR_STREAM << ", parameter #" << counter;
          }

          ERROR_STREAM << " was not given!" << '\n';
        }
      } // for

      return valid_list;
    } // IsValidList

    //! @brief adds PARAMETER to the m_ParameterList
    //! @param PARAMETER the parameter to be added
    void FlagStatic::PushBack( const util::ShPtr< ParameterInterface> &PARAMETER)
    {
      m_ParameterList.PushBack( PARAMETER);
      if( PARAMETER->GetName().empty())
      {
        m_ParameterList.LastElement()->SetName( util::Format()( m_ParameterList.GetSize()));
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes the user provided commandline
    //! @param OSTREAM the stream to which the commandline is written
    //! @return given ostream to which the commandline was written
    std::ostream &FlagStatic::WriteUserCommand( std::ostream &OSTREAM) const
    {
      // write description
      OSTREAM << "-" << m_Name << ( m_WasSetInCommandline ? " set" : " not set") << '\n';

      // write WriteUserCommand for every parameter in the list
      for
      (
        util::ShPtrVector< ParameterInterface>::const_iterator
          param_itr( m_ParameterList.Begin()),
          param_itr_end( m_ParameterList.End());
        param_itr != param_itr_end;
        ++param_itr
      )
      {
        OSTREAM << "   ";
        ( *param_itr)->WriteUserCommand( OSTREAM);
      } // for

      // end
      return OSTREAM;
    } // WriteUserCommand

    //! @brief writes the help for the commandline
    //! @param OSTREAM the stream to which the help is written to
    //! @return the given stream to which the help was written to
    std::ostream &FlagStatic::WriteHelp( std::ostream &OSTREAM) const
    {
      // write name and description
      if( !m_Name.empty() || !m_Description.empty())
      {
        io::FixedLineWidthWriter writer;
        // write out the flag name directly, wrap the description / size part, indenting to beginning of the description
        writer << '-' << m_Name << " : ";

        // write description
        if( m_Description.size() < 2 * writer.GetRemainingSpaceOnLine())
        {
          // short description
          writer.SetIndent( writer.GetLinePosition());
        }
        else
        {
          // just indent by 4 and write the description; it is too long to fit in a fully-indented line
          writer.AddIndent( 2);
        }

        writer << m_Description;
        writer.PopIndent();

        // add a newline, if one is not already present
        if( writer.GetLinePosition())
        {
          writer.NewLine();
        }
        OSTREAM << writer.String();
      }

      // static list
      for
      (
        util::ShPtrVector< ParameterInterface>::const_iterator
          param_itr( m_ParameterList.Begin()),
          param_itr_end( m_ParameterList.End());
        param_itr != param_itr_end;
        ++param_itr
      )
      {
        ( *param_itr)->WriteHelp( OSTREAM, 1);
      } // for

      //end
      return OSTREAM;
    } // WriteHelp

    //! @brief writes the usage command, complete with required parameters and flags
    //! @param OSTREAM - the outstream to write to
    //! @return std::ostream &OSTREAM - return the stream after writing
    std::ostream &FlagStatic::WriteRequiredUsage( std::ostream &OSTREAM) const
    {
      const size_t n_required_parameters( GetNumberRequiredParameters());
      if( n_required_parameters)
      {
        if( !m_Name.empty())
        {
          OSTREAM << '-' << m_Name << ' ';
        }
        size_t counter( 0);
        for
        (
          util::ShPtrVector< ParameterInterface>::const_iterator param_itr( m_ParameterList.Begin());
          counter < n_required_parameters;
          ++param_itr, ++counter
        )
        {
          OSTREAM << '<' << ( *param_itr)->GetName() << "> ";
        }
      }
      return OSTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &FlagStatic::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Name,                OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Description,         OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_WasSetInCommandline, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ParameterList,       OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FlagStatic::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Name,                ISTREAM);
      io::Serialize::Read( m_Description,         ISTREAM);
      io::Serialize::Read( m_WasSetInCommandline, ISTREAM);
      io::Serialize::Read( m_ParameterList,       ISTREAM);

      // end
      return ISTREAM;
    }

  } // namespace command
} // namespace bcl
