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
#include "command/bcl_command_flag_dynamic.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_parameter.h"
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

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> FlagDynamic::s_Instance
    (
      GetObjectInstances().AddInstance( new FlagDynamic())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FlagDynamic::FlagDynamic() :
      m_WasSetInCommandline( false)
    {
    }

    //! @brief construct FlagDynamic from NAME, DESCRIPTION, TEMPLATE_PARAMETER,
    //! @brief and optional MIN_NUMBER_PARAMETERS, MAX_NUMBER_PARAMETERS
    //! @param NAME - the name of the flag as a std::string
    //! @param DESCRIPTION - a description of the flag as a std::string
    //! @param TEMPLATE_PARAMETER
    //! @param MIN_NUMBER_PARAMETERS - the minimum number of parameters for this flag
    //! @param MAX_NUMBER_PARAMETERS - the maximum number of parameters for this flag
    //! @param SIGNAL optional function to call whenever calling set or reset on the flag
    FlagDynamic::FlagDynamic
    (
      const std::string &NAME,
      const std::string &DESCRIPTION,
      const ParameterInterface &TEMPLATE_PARAMETER,
      const size_t MIN_NUMBER_PARAMETERS,
      const size_t MAX_NUMBER_PARAMETERS,
      FlagInterface::t_Signal SIGNAL
    ) :
      m_Name( NAME),
      m_Description( util::TrimString( DESCRIPTION)),
      m_WasSetInCommandline( false),
      m_MinNumberParameters( MIN_NUMBER_PARAMETERS),
      m_MaxNumberParameters( MAX_NUMBER_PARAMETERS),
      m_TemplateParameters( 1, util::CloneToShPtr( TEMPLATE_PARAMETER)),
      m_Signal( SIGNAL)
    {
      BCL_Assert
      (
        MIN_NUMBER_PARAMETERS <= MAX_NUMBER_PARAMETERS,
        "min number parameters > max number parameters: "
          + util::Format()( MIN_NUMBER_PARAMETERS) + " > "
          + util::Format()( MAX_NUMBER_PARAMETERS)
      );
      if( TEMPLATE_PARAMETER.GetName().empty())
      {
        m_TemplateParameters.LastElement()->SetName( "1");
      }
    } // FlagDynamic

    //! @brief construct FlagDynamic from NAME, DESCRIPTION, TEMPLATE_PARAMETERS,
    //! @brief and optional MIN_NUMBER_PARAMETERS, MAX_NUMBER_PARAMETERS
    //! @param NAME - the name of the flag as a std::string
    //! @param DESCRIPTION - a description of the flag as a std::string
    //! @param TEMPLATE_PARAMETER
    //! @param MIN_NUMBER_PARAMETERS - the minimum number of parameters for this flag
    //! @param MAX_NUMBER_PARAMETERS - the maximum number of parameters for this flag
    //! @param SIGNAL optional function to call whenever calling set or reset on the flag
    FlagDynamic::FlagDynamic
    (
      const std::string &NAME,
      const std::string &DESCRIPTION,
      const storage::Vector< Parameter> &TEMPLATE_PARAMETERS,
      const size_t MIN_NUMBER_PARAMETERS,
      const size_t MAX_NUMBER_PARAMETERS,
      FlagInterface::t_Signal SIGNAL
    ) :
      m_Name( NAME),
      m_Description( DESCRIPTION),
      m_WasSetInCommandline( false),
      m_MinNumberParameters( MIN_NUMBER_PARAMETERS),
      m_MaxNumberParameters( MAX_NUMBER_PARAMETERS),
      m_TemplateParameters(),
      m_Signal( SIGNAL)
    {
      BCL_Assert
      (
        MIN_NUMBER_PARAMETERS <= MAX_NUMBER_PARAMETERS,
        "min number parameters > max number parameters: "
          + util::Format()( MIN_NUMBER_PARAMETERS) + " > "
          + util::Format()( MAX_NUMBER_PARAMETERS)
      );
      BCL_Assert
      (
        ( MIN_NUMBER_PARAMETERS % TEMPLATE_PARAMETERS.GetSize()) == 0,
        "min number parameters is not a multiple of the number of internal parameters"
      );
      BCL_Assert
      (
        !util::IsDefined( MAX_NUMBER_PARAMETERS) || ( MAX_NUMBER_PARAMETERS % TEMPLATE_PARAMETERS.GetSize()) == 0,
        "max number parameters is not a multiple of the number of internal parameters"
      );

      size_t counter( 1);
      for
      (
        storage::Vector< Parameter>::const_iterator
          itr( TEMPLATE_PARAMETERS.Begin()), itr_end( TEMPLATE_PARAMETERS.End());
        itr != itr_end;
        ++itr, ++counter
      )
      {
        m_TemplateParameters.PushBack( util::CloneToShPtr( *itr));
        if( itr->GetName().empty())
        {
          m_TemplateParameters.LastElement()->SetName( util::Format()( counter));
        }
      }
    } // FlagDynamic

    //! @brief virtual copy constructor
    FlagDynamic *FlagDynamic::Clone() const
    {
      return new FlagDynamic( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief name of class
    //! @return the name of the class as a std::string
    const std::string &FlagDynamic::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns name of flag
    //! @return name of flag as std::string
    const std::string &FlagDynamic::GetName() const
    {
      return m_Name;
    }

    //! @brief set name - only possible if there was no name given yet
    //! @param NAME - what you want to name the flag
    void FlagDynamic::SetName( const std::string &NAME)
    {
      BCL_Assert( m_Name.empty(), "name was already given");
      m_Name = NAME;
    }

    //! @brief returns description of flag
    //! @return description of flag as string
    const std::string &FlagDynamic::GetDescription() const
    {
      return m_Description;
    }

    //! @brief set flag to true if it was false
    void FlagDynamic::SetFlag()
    {
      m_WasSetInCommandline = true;
    }

    //! @brief set flag to false if it was true
    void FlagDynamic::UnsetFlag()
    {
      m_WasSetInCommandline = false;
    }

    //! @brief reset the flag
    //! @detail resets all internal parameters, then removes all internally held parameters
    void FlagDynamic::ResetFlag()
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

      m_ParameterList.Reset();

      // call the handling function, if it was defined
      if( m_Signal && FlagInterface::ShouldSignal())
      {
        ( *m_Signal)();
      }
    }

    //! @brief returns the function to be called whenever this flag is updated
    //! @return the function to be called whenever this flag is updated
    FlagInterface::t_Signal FlagDynamic::GetSignal() const
    {
      return m_Signal;
    }

    //! @brief returns if m_IsSet in command line
    //! @return whether this flag was set in the command line
    bool FlagDynamic::GetFlag() const
    {
      return m_WasSetInCommandline;
    }

    //! @brief returns m_ParameterList
    //! @return the parameter list as a ShPtrVector
    const util::ShPtrVector< ParameterInterface> &FlagDynamic::GetParameterList() const
    {
      return m_ParameterList;
    }

    //! @brief returns m_ParameterList
    //! @return the parameter list as a util::ShPtrVector< ParameterInterface>
    util::ShPtrVector< ParameterInterface> &FlagDynamic::GetParameterList()
    {
      return m_ParameterList;
    }

    //! @brief returns the first parameter
    //! @return the first parameter for this flag as a util::ShPtr< ParameterInterface>
    const util::ShPtr< ParameterInterface> &FlagDynamic::GetFirstParameter() const
    {
      BCL_Assert( !m_ParameterList.IsEmpty(), "there is no first parameter for " + m_Name);
      return m_ParameterList.FirstElement();
    }

    //! @brief returns the # of required parameters for this flag
    //! @return the effective # of required parameters for this flag
    size_t FlagDynamic::GetNumberRequiredParameters() const
    {
      return m_TemplateParameters.GetSize() * m_MinNumberParameters;
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
    bool FlagDynamic::ReadFromList( const storage::Vector< std::string> &PARAMETER_LIST, std::ostream &ERROR_STREAM)
    {
      // check if flag was already set
      if( GetFlag())
      {
        ERROR_STREAM << "flag: \"-" << GetName() << "\" was already set" << '\n';
        return false;
      }

      // set flag
      SetFlag();

      bool success_setting_parameters( true);

      // iterate over parameter in given commandline
      for
      (
        std::vector< std::string>::const_iterator
          param_itr( PARAMETER_LIST.Begin()),
          param_itr_end( PARAMETER_LIST.End());
        param_itr != param_itr_end;
        ++param_itr
      )
      {
        // pushback new entries in the parameterlist
        success_setting_parameters &= PushBack( *param_itr, ERROR_STREAM);
      }

      // ensure that list is long enough - but not too long
      success_setting_parameters &= IsValidList( ERROR_STREAM);

      // call the handling function, if it was defined
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
    bool FlagDynamic::IsValidList( std::ostream &ERROR_STREAM) const
    {
      bool is_valid = MeetSizeSpecification();

      // check if parameter list has sufficient size
      if( !is_valid)
      {
        if( m_ParameterList.GetSize() % m_TemplateParameters.GetSize())
        {
          ERROR_STREAM << "parameter list after \"-" + GetName() + "\" has improper size: given " << GetSize()
            << " parameters, but need a multiple of " << m_TemplateParameters.GetSize() << '\n';
        }
        else
        {
          const size_t multiplicity( m_ParameterList.GetSize() / m_TemplateParameters.GetSize());
          ERROR_STREAM << "parameter list after \"-" + GetName() + "\" has wrong size: " << multiplicity
            << " is not between: [" << m_MinNumberParameters << ".." << m_MaxNumberParameters << "]" << '\n';
        }
      }

      return is_valid;
    } // IsValidList

    //! @brief returns true if the number of parameters is between m_MinNumberParameters and m_MaxNumberParameters
    //! @return true if the parameter list meets the size specification, false otherwise
    bool FlagDynamic::MeetSizeSpecification() const
    {
      return
          ( m_ParameterList.GetSize() % m_TemplateParameters.GetSize()) == 0
          && ( GetSize() / m_TemplateParameters.GetSize() >= m_MinNumberParameters)
          && ( GetSize() / m_TemplateParameters.GetSize() <= m_MaxNumberParameters);
    }

    //! @brief adds PARAMETER to the m_ParameterList
    //! @param PARAMETER - the parameter to be added to the parameter list
    //! @param ERROR_STREAM - the stream to which errors should be written
    //! @return true on success, false otherwise
    //! @see FlagDynamic::PushBack
    bool FlagDynamic::PushBack( const std::string &PARAMETER, std::ostream &ERROR_STREAM)
    {
      // determine which template parameter this member of the flag belongs to
      const size_t parameter_id( m_ParameterList.GetSize() % m_TemplateParameters.GetSize());

      util::ShPtr< ParameterInterface> new_parameter
      (
        m_TemplateParameters( parameter_id).HardCopy()
      );
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
    std::ostream &FlagDynamic::WriteHelp( std::ostream &OSTREAM) const
    {
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
        writer << m_Description << ", ";

        // get the description for the sizes
        std::string size_description;
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

        // add the size description to the writer
        writer << "This flag can be followed by " << size_description;
        writer.PopIndent();
      }

      // add a newline, if one was not already present
      if( writer.GetLinePosition())
      {
        writer.NewLine();
      }
      OSTREAM << writer.String();

      // write the writer to the output stream

      for
      (
        util::ShPtrVector< ParameterInterface>::const_iterator
          param_itr( m_TemplateParameters.Begin()), param_itr_end( m_TemplateParameters.End());
        param_itr != param_itr_end;
        ++param_itr
      )
      {
        ( *param_itr)->WriteHelp( OSTREAM, 1);
      }

      // end
      return OSTREAM;
    } // WriteHelp

    //! @brief writes the user provided commandline
    //! @return the stream after you write to it
    std::ostream &FlagDynamic::WriteUserCommand( std::ostream &OSTREAM) const
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
    } // WriteUserCommand

    //! @brief writes the usage command, complete with required parameters and flags
    //! @param OSTREAM - the outstream to write to
    //! @return std::ostream &OSTREAM - return the stream after writing
    std::ostream &FlagDynamic::WriteRequiredUsage( std::ostream &OSTREAM) const
    {
      const size_t n_required_parameters( GetNumberRequiredParameters());
      if( n_required_parameters)
      {
        OSTREAM << '-' << m_Name << ' ';
        for( size_t counter( 0); counter < m_MinNumberParameters; ++counter)
        {
          for
          (
            util::ShPtrVector< ParameterInterface>::const_iterator
              param_itr( m_TemplateParameters.Begin()),
              param_itr_end( m_TemplateParameters.End());
            param_itr != param_itr_end;
            ++param_itr
          )
          {
            OSTREAM << '<' << ( *param_itr)->GetName() << "> ";
          }
        }
      }
      return OSTREAM;
    }

    //! @brief write Flag with Params to ostream
    //! @param OSTREAM - the stream to write to
    //! @param INDENT indentation
    //! @return the stream after you wrote to it
    std::ostream &FlagDynamic::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Name,                OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Description,         OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_WasSetInCommandline, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ParameterList,       OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MinNumberParameters, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MaxNumberParameters, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_TemplateParameters,  OSTREAM, INDENT);

      // end
      return OSTREAM;
    } // Write

    //! @brief read Flag with params from istream
    //! @param ISTREAM - the stream to read from
    //! @return the stream after you read from it
    std::istream &FlagDynamic::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Name,                ISTREAM);
      io::Serialize::Read( m_Description,         ISTREAM);
      io::Serialize::Read( m_WasSetInCommandline, ISTREAM);
      io::Serialize::Read( m_ParameterList,       ISTREAM);
      io::Serialize::Read( m_MinNumberParameters, ISTREAM);
      io::Serialize::Read( m_MaxNumberParameters, ISTREAM);
      io::Serialize::Read( m_TemplateParameters,  ISTREAM);

      // end
      return ISTREAM;
    } // Read

  } // namespace command
} // namespace bcl
