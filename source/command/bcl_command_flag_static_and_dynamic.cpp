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
#include "command/bcl_command_flag_static_and_dynamic.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_fixed_line_width_writer.h"
#include "util/bcl_util_data_type.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> FlagStaticAndDynamic::s_Instance
    (
      GetObjectInstances().AddInstance( new FlagStaticAndDynamic())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FlagStaticAndDynamic::FlagStaticAndDynamic() :
      m_Name(),
      m_Description(),
      m_WasSetInCommandline( false),
      m_MinNumberParameters(),
      m_MaxNumberParameters(),
      m_ParameterList(),
      m_TemplateParameter(),
      m_NumberStaticParameters( 0)
    {
    }

    //! @brief construct FlagStaticAndDynamic from NAME, DESCRIPTION, TEMPLATE_PARAMETER,
    //! @brief and optional MIN_NUMBER_PARAMETERS, MAX_NUMBER_PARAMETERS
    //! @param NAME - the name of the flag as a std::string
    //! @param DESCRIPTION - a description of the flag as a std::string
    //! @param TEMPLATE_PARAMETER
    //! @param MIN_NUMBER_PARAMETERS - the minimum number of parameters for this flag
    //! @param MAX_NUMBER_PARAMETERS - the maximum number of parameters for this flag
    //! @param SIGNAL optional function to call whenever calling set or reset on the flag
    FlagStaticAndDynamic::FlagStaticAndDynamic
    (
      const std::string &NAME,
      const std::string &DESCRIPTION,
      const ParameterInterface &TEMPLATE_PARAMETER,
      const size_t MIN_NUMBER_PARAMETERS,
      const size_t MAX_NUMBER_PARAMETERS,
      t_Signal SIGNAL
    ) :
      m_Name( NAME),
      m_Description( util::TrimString( DESCRIPTION)),
      m_WasSetInCommandline( false),
      m_MinNumberParameters( MIN_NUMBER_PARAMETERS),
      m_MaxNumberParameters( MAX_NUMBER_PARAMETERS),
      m_TemplateParameter( TEMPLATE_PARAMETER.Clone()),
      m_NumberStaticParameters( 0),
      m_Signal( SIGNAL)
    {
    }

    //! @brief copy constructor
    FlagStaticAndDynamic *FlagStaticAndDynamic::Clone() const
    {
      return new FlagStaticAndDynamic( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief name of class
    //! @return the name of the class as a std::string
    const std::string &FlagStaticAndDynamic::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns name of flag
    //! @return name of flag as std::string
    const std::string &FlagStaticAndDynamic::GetName() const
    {
      return m_Name;
    }

    //! @brief set name - only possible if there was no name given yet
    //! @param NAME - what you want to name the flag
    void FlagStaticAndDynamic::SetName( const std::string &NAME)
    {
      BCL_Assert( m_Name.empty(), "name was already given");
      m_Name = NAME;
    }

    //! @brief returns description of flag
    //! @return description of flag as string
    const std::string &FlagStaticAndDynamic::GetDescription() const
    {
      return m_Description;
    }

    //! @brief returns the function to be called whenever this flag is updated
    //! @return the function to be called whenever this flag is updated
    FlagInterface::t_Signal FlagStaticAndDynamic::GetSignal() const
    {
      return m_Signal;
    }

    //! @brief set flag to true if it was false
    void FlagStaticAndDynamic::SetFlag()
    {
      m_WasSetInCommandline = true;
    }

    //! @brief set flag to false if it was true
    void FlagStaticAndDynamic::UnsetFlag()
    {
      m_WasSetInCommandline = false;
    }

    //! @brief reset the flag
    //! @detail resets all internal parameters, then removes all internally held parameters
    void FlagStaticAndDynamic::ResetFlag()
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

      // remove any dynamic parameters
      const size_t number_dynamic_parameters( m_ParameterList.GetSize() - m_NumberStaticParameters);
      if( number_dynamic_parameters > size_t( 0))
      {
        m_ParameterList.RemoveElements( m_NumberStaticParameters, number_dynamic_parameters);
      }
      if( m_Signal && FlagInterface::ShouldSignal())
      {
        ( *m_Signal)();
      }
    }

    //! @brief returns if m_IsSet in command line
    //! @return whether this flag was set in the command line
    bool FlagStaticAndDynamic::GetFlag() const
    {
      return m_WasSetInCommandline;
    }

    //! @brief returns m_ParameterList
    //! @return the parameter list as a ShPtrVector
    const util::ShPtrVector< ParameterInterface> &FlagStaticAndDynamic::GetParameterList() const
    {
      return m_ParameterList;
    }

    //! @brief returns m_ParameterList
    //! @return the parameter list as a util::ShPtrVector< ParameterInterface>
    util::ShPtrVector< ParameterInterface> &FlagStaticAndDynamic::GetParameterList()
    {
      return m_ParameterList;
    }

    //! @brief returns the first parameter
    //! @return the first parameter for this flag as a util::ShPtr< ParameterInterface>
    const util::ShPtr< ParameterInterface> &FlagStaticAndDynamic::GetFirstParameter() const
    {
      BCL_Assert( !m_ParameterList.IsEmpty(), "there is no first parameter for " + m_Name);
      return m_ParameterList.FirstElement();
    }

    //! @brief returns the # of required parameters for this flag
    //! @return the effective # of required parameters for this flag
    size_t FlagStaticAndDynamic::GetNumberRequiredParameters() const
    {
      // if dynamic parameters are required, all static parameters must be given
      if( m_MinNumberParameters)
      {
        return m_MinNumberParameters + m_NumberStaticParameters;
      }
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
        // if no default was given, update the # required to the counter
        if( !( *itr)->GetWasDefaultGiven())
        {
          return counter;
        }
      }

      return 0;
    }

    //! @brief returns list of dynamic parameters
    //! @return the parameter list of dynamically added parameters
    util::SiPtrVector< const ParameterInterface> FlagStaticAndDynamic::GetDynamicParameterList() const
    {
      util::SiPtrVector< const ParameterInterface> dynamic_list;
      for
      (
        util::ShPtrVector< ParameterInterface>::const_iterator
          itr( m_ParameterList.Begin() + m_NumberStaticParameters), itr_end( m_ParameterList.End());
        itr != itr_end;
        ++itr
      )
      {
        // insert current parameter
        dynamic_list.PushBack( *itr);
      }

      // end
      return dynamic_list;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief fill the parameter list from storage vector of strings
    //! @param PARAMETER_LIST - the parameter list as a storage vector of strings
    //! @param ERROR_STREAM - the stream to which errors should be written
    //! @return true on success, false on failure
    //! does not allow you to add more than m_MaxNumberParameters parameters to the list
    //! @see SetFlag
    //! @see PushBack
    bool FlagStaticAndDynamic::ReadFromList( const storage::Vector< std::string> &PARAMETER_LIST, std::ostream &ERROR_STREAM)
    {
      // check if flag was already set
      if( GetFlag())
      {
        ERROR_STREAM << "flag: \"-" << GetName() << "\" was already set" << '\n';
        return false;
      }

      SetFlag();

      if( PARAMETER_LIST.IsEmpty())
      {
        return true;
      }

      bool success_setting_parameters( true);

      // iterators on commandline arguments
      storage::Vector< std::string>::const_iterator
        param_itr( PARAMETER_LIST.Begin()),
        param_itr_end( PARAMETER_LIST.End());

      // read static params
      for
      (
        util::ShPtrVector< ParameterInterface>::iterator
          itr( m_ParameterList.Begin()), itr_end( m_ParameterList.Begin() + m_NumberStaticParameters);
        itr != itr_end && param_itr != param_itr_end;
        ++itr, ++param_itr
      )
      {
        // set parameters in static parameter list
        success_setting_parameters &= ( *itr)->SetParameter( *param_itr, ERROR_STREAM);
      }

      // iterate over dynamic parameters in given commandline
      for( ; param_itr != param_itr_end; ++param_itr)
      {
        // pushback new entries in the parameterlist
        success_setting_parameters &= PushBack( *param_itr, ERROR_STREAM);
      }

      // ensure that list is long enough - but not too long
      success_setting_parameters &= IsValidList( ERROR_STREAM);

      // call the signalling function if successful and one was given
      if( success_setting_parameters && m_Signal && FlagInterface::ShouldSignal())
      {
        ( *m_Signal)();
      }

      return success_setting_parameters;
    } // ReadFromList

    //! @brief checks if there are too few items in your parameter list
    //! @param ERROR_STREAM - the stream to which errors should be written
    //! @return true if the parameter list is of a valid size, false otherwise
    //! @see MeetSizeSpecification
    bool FlagStaticAndDynamic::IsValidList( std::ostream &ERROR_STREAM) const
    {
      bool valid_list( true);
      size_t counter( 0);
      for
      (
        util::ShPtrVector< ParameterInterface>::const_iterator
          param_itr( GetParameterList().Begin()),
          param_itr_end( GetParameterList().End());
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

      // check if parameter list has sufficient size
      if( !MeetSizeSpecification())
      {
        valid_list = false;
        ERROR_STREAM << "parameter list after \"-" + GetName() + "\" has wrong size: " << GetSize()
          << " is not between: [" << m_MinNumberParameters << ".." << m_MaxNumberParameters << "]" << '\n';
      }

      return valid_list;
    }

    //! @brief returns true if the number of parameters is between m_MinNumberParameters and m_MaxNumberParameters
    //! @return true if the parameter list meets the size specification, false otherwise
    bool FlagStaticAndDynamic::MeetSizeSpecification() const
    {
      const size_t number_dynamic_parameters( GetSize() - m_NumberStaticParameters);
      return ( number_dynamic_parameters >= m_MinNumberParameters) && ( number_dynamic_parameters <= m_MaxNumberParameters);
    }

    //! @brief adds PARAMETER to the m_ParameterList
    //! @param PARAMETER - the parameter to be added to the parameter list
    //! @param ERROR_STREAM - the stream to which errors should be written
    //! @return true on success, false otherwise
    //! @see FlagInterface::PushBack
    bool FlagStaticAndDynamic::PushBack( const std::string &PARAMETER, std::ostream &ERROR_STREAM)
    {
      util::ShPtr< ParameterInterface> new_parameter( m_TemplateParameter.HardCopy());
      // try to set parameter
      if( new_parameter->SetParameter( PARAMETER, ERROR_STREAM))
      {
        // if this was an allowed parameter => pushback
        m_ParameterList.PushBack( new_parameter);
        return true;
      }
      return false;
    } // PushBack

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes the help for the command line
    //! @param OSTREAM - the stream to write to
    //! @return the stream after you wrote to it
    std::ostream &FlagStaticAndDynamic::WriteHelp( std::ostream &OSTREAM) const
    {
      // write name and description
      io::FixedLineWidthWriter writer;

      // write name and description
      if( !m_Name.empty() || !m_Description.empty())
      {
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
      }

      // add a newline, if one is not already present
      if( writer.GetLinePosition())
      {
        writer.NewLine();
      }
      OSTREAM << writer.String();

      // get the description for the sizes
      std::string size_description;
      {
        std::ostringstream description_stream;
        util::DataType::WriteSizeRequirements( description_stream, m_MinNumberParameters, m_MaxNumberParameters);
        size_description = description_stream.str();
        // erase with, if it starts with that
        if( util::StartsWith( size_description, " with"))
        {
          size_description.erase( 0, 5);
        }
        else
        {
          // any number
          size_description = " any number of ";
        }
      }

      // static list
      for
      (
        util::ShPtrVector< ParameterInterface>::const_iterator
          param_itr( GetParameterList().Begin()),
          param_itr_end( GetParameterList().End());
        param_itr != param_itr_end;
        ++param_itr
      )
      {
        ( *param_itr)->WriteHelp( OSTREAM, 1);
      } // for

      // write the template help
      OSTREAM << "  Followed by " << size_description << " of ";
      m_TemplateParameter->WriteHelp( OSTREAM, 1);

      //end
      return OSTREAM;
    }

    //! @brief writes the user provided commandline to a stream
    //! @param OSTREAM stream to write to
    //! @return the stream written to
    std::ostream &FlagStaticAndDynamic::WriteUserCommand( std::ostream &OSTREAM) const
    {
      // write description
      OSTREAM << "-" << m_Name << ( m_WasSetInCommandline ? " set" : " not set") << '\n';

      // write WriteUserCommand for every parameter in the list
      for
      (
        util::ShPtrVector< ParameterInterface>::const_iterator
          param_itr( m_ParameterList.Begin()), param_itr_end( m_ParameterList.End());
        param_itr != param_itr_end;
        ++param_itr
      )
      {
        OSTREAM << "   ";
        ( *param_itr)->WriteUserCommand( OSTREAM);
      } // for

      // end
      return OSTREAM;
    }

    //! @brief writes the usage command, complete with required parameters and flags
    //! @param OSTREAM - the outstream to write to
    //! @return std::ostream &OSTREAM - return the stream after writing
    std::ostream &FlagStaticAndDynamic::WriteRequiredUsage( std::ostream &OSTREAM) const
    {
      const size_t n_required_parameters( GetNumberRequiredParameters());
      if( n_required_parameters)
      {
        OSTREAM << '-' << m_Name << ' ';
        size_t parameter_counter( 0);
        for
        (
          util::ShPtrVector< ParameterInterface>::const_iterator
            param_itr( m_ParameterList.Begin()),
            param_itr_end( m_ParameterList.End());
          param_itr != param_itr_end && parameter_counter < n_required_parameters;
          ++param_itr, ++parameter_counter
        )
        {
          OSTREAM << '<' << ( *param_itr)->GetName() << "> ";
        }
        while( parameter_counter < n_required_parameters)
        {
          OSTREAM << '<' << m_TemplateParameter->GetName() << "> ";
        }
      }
      return OSTREAM;
    }

    //! @brief write Flag with Params to ostream
    //! @param OSTREAM - the stream to write to
    //! @return the stream after you wrote to it
    std::ostream &FlagStaticAndDynamic::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Name,                OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Description,         OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_WasSetInCommandline, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ParameterList,       OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MinNumberParameters, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MaxNumberParameters, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_TemplateParameter,   OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief read Flag with params from istream
    //! @param ISTREAM - the stream to read from
    //! @return the stream after you read from it
    std::istream &FlagStaticAndDynamic::Read( std::istream &ISTREAM)
    {
      // write members
      io::Serialize::Read( m_Name,                ISTREAM);
      io::Serialize::Read( m_Description,         ISTREAM);
      io::Serialize::Read( m_WasSetInCommandline, ISTREAM);
      io::Serialize::Read( m_ParameterList,       ISTREAM);
      io::Serialize::Read( m_MinNumberParameters, ISTREAM);
      io::Serialize::Read( m_MaxNumberParameters, ISTREAM);
      io::Serialize::Read( m_TemplateParameter,   ISTREAM);

      // end
      return ISTREAM;
    }

  } // namespace command
} // namespace bcl
