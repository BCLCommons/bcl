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
#include "app/bcl_app_group_handler.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_command.h"
#include "io/bcl_io_fixed_line_width_writer.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {
    // Length of the longest app name in any group
    size_t GroupHandler::s_LongestGroupSize( std::string( "<group>:Help").size());

    // Length of the longest app name in any group
    size_t GroupHandler::s_LongestNameSize( 8);

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default group handler
    GroupHandler::GroupHandler() :
      m_GroupName( ""),
      m_GroupDescription( "Internal applications that have not yet been categorized"),
      m_LongestNameSize( 0)
    {
    }

    //! @brief constructor from string of group name
    GroupHandler::GroupHandler( const std::string &NAME, const std::string &DESCRIPTION) :
      m_GroupName( NAME),
      m_GroupDescription( DESCRIPTION),
      m_NameToDescription(),
      m_NameToAliases(),
      m_Flags(),
      m_LongestNameSize( m_GroupName.size() + size_t( 1))
    {
      s_LongestGroupSize = std::max( s_LongestGroupSize, m_GroupName.size() + size_t( 1));
      s_LongestNameSize = std::max( s_LongestNameSize, m_GroupName.size() + size_t( 1));
    }

    //! @brief Clone function
    //! @return pointer to new GroupHandler
    GroupHandler *GroupHandler::Clone() const
    {
      return new GroupHandler( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &GroupHandler::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the group
    const std::string &GroupHandler::GetName() const
    {
      return m_GroupName;
    }

    //! @brief get the description of the group
    const std::string &GroupHandler::GetDescription() const
    {
      return m_GroupDescription;
    }

    //! @brief test whether a particular application is in this group
    //! @param APP the application of interst
    bool GroupHandler::Contains( const Interface &APP) const
    {
      if( m_GroupName.empty())
      {
        return m_NameToDescription.Has( APP.GetNameForGroup());
      }
      return m_NameToDescription.Has( APP.GetNameForGroup( m_GroupName).substr( m_GroupName.size() + 1));
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief add group flags to a given command
    //! @param COMMAND the command to add group flags to
    void GroupHandler::AddGroupFlags( command::Command &COMMAND)
    {
      COMMAND.PushBack( m_Flags);
    }

    //! @brief add a flag to the application group
    //! @param FLAG the new flag to add
    void GroupHandler::AddFlagToGroup( const util::ShPtr< command::FlagInterface> &FLAG)
    {
      m_Flags.PushBack( FLAG);
    }

    //! @brief add an application to this group
    ApplicationType &GroupHandler::AddInstance( const util::ShPtr< Interface> &APP)
    {
      std::string app_name_with_group( APP->GetNameForGroup( m_GroupName));

      // strip the appropriate group identifying characters
      std::string app_name_wo_group
      (
        m_GroupName.empty()
        ? app_name_with_group
        : app_name_with_group.substr( m_GroupName.size() + 1)
      );
      m_LongestNameSize = std::max( app_name_wo_group.size(), m_LongestNameSize);
      s_LongestNameSize = std::max( app_name_wo_group.size(), s_LongestNameSize);

      BCL_Assert
      (
        m_NameToDescription.Insert( std::make_pair( app_name_wo_group, APP->GetDescription())).second,
        "App " + app_name_wo_group + " was added twice to " + m_GroupName
      );

      // get all the application names for this app
      storage::Vector< std::string> all_app_names( APP->GetDeprecatedAppNames());

      // copy the application names, sort them, to display for formerly-known
      if( !all_app_names.IsEmpty())
      {
        storage::Vector< std::string> app_names_sorted( all_app_names);
        if
        (
          APP->GetLicensedName() != APP->GetNameForGroup()
          && APP->GetLicensedName() != app_name_with_group
          && APP->GetLicensedName() != app_name_wo_group
        )
        {
          app_names_sorted.PushBack( APP->GetLicensedName());
        }
        app_names_sorted.Sort( std::less< std::string>());
        if( app_names_sorted.GetSize())
        {
          m_NameToAliases[ app_name_wo_group] = "Formerly known as " + util::Join( ", ", app_names_sorted);
        }
        else
        {
          // insert an empty statement for the aliases string
          m_NameToAliases[ app_name_wo_group] = "";
        }
      }

      // add the actual name to the group
      if( APP->GetLicensedName() != APP->GetNameForGroup() && APP->GetLicensedName() != app_name_with_group)
      {
        all_app_names.PushBack( APP->GetLicensedName());
      }

      if( APP->GetNameForGroup() != app_name_with_group)
      {
        all_app_names.PushBack( APP->GetNameForGroup());
      }

      // add the app to the enum
      for
      (
        storage::Vector< std::string>::const_iterator itr( all_app_names.Begin()), itr_end( all_app_names.End());
        itr != itr_end;
        ++itr
      )
      {
        if( !GetApps().HaveEnumWithName( *itr))
        {
          GetApps().AddEnum( *itr, APP);
        }
      }
      return GetApps().AddEnum( app_name_with_group, APP);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes the help for the commandline
    //! @param OSTREAM the stream to which the help is written to
    //! @param INDENT the amount of indent to write before the help
    //! @param INCLUDE_MEMBERS whether to write the help for every group member
    //! @param INCLUDE_MEMBER_DESCRIPTIONS whether to write a brief description for each group member
    //! @param INCLUDE_MEMBER_ALIASES whether to write aliases (even if they are deprecated) for each group member
    //! @return the given stream to which the help was written to
    std::ostream &GroupHandler::WriteHelp
    (
      std::ostream &OSTREAM,
      const size_t INDENT,
      const size_t &APP_INDENT,
      const bool &INCLUDE_DESCRIPTION,
      const bool &INCLUDE_MEMBERS,
      const bool &INCLUDE_MEMBER_DESCRIPTIONS,
      const bool &INCLUDE_MEMBER_ALIASES
    ) const
    {
      // check whether any applications are present, if not, just return
      if( m_NameToDescription.IsEmpty() || ( m_NameToDescription.GetSize() == size_t( 1) && m_GroupName != "all"))
      {
        return OSTREAM;
      }

      // write the name, left flush
      const std::string full_group_name( m_GroupName + s_GroupDelimiter);

      io::FixedLineWidthWriter writer( 0, util::GetLogger().GetMaxLineWidth() - 4);
      const size_t name_indent( APP_INDENT + 2);
      if( INCLUDE_MEMBERS)
      {
        OSTREAM << '\n';
      }

      // handle the case where there is only one application (typically help)
      if( m_NameToDescription.GetSize() == size_t( 1))
      {
        writer << full_group_name;

        // get an iterator to the only app in this group
        storage::Map< std::string, std::string>::const_iterator itr_app( m_NameToDescription.Begin());

        // write the name, followed by the spaces
        writer << itr_app->first;

        // determine spaces until description, which should be aligned between the groups
        const size_t spaces_till_description
        (
          std::max( s_LongestGroupSize + size_t( 2), writer.GetLinePosition()) - writer.GetLinePosition()
        );
        writer << std::string( spaces_till_description, ' ');

        // set indent for the description
        writer.SetIndent( writer.GetLinePosition());
        writer << itr_app->second;
        if( writer.GetLinePosition())
        {
          writer.NewLine();
        }
        OSTREAM << writer.String();
        return OSTREAM;
      }

      writer << full_group_name;
      if( INCLUDE_DESCRIPTION)
      {
        writer << std::string( s_LongestGroupSize + size_t( 2) - full_group_name.size(), ' ');
        writer.SetIndent( writer.GetLinePosition());
        writer << "Group of " << m_GroupDescription;
        writer.PopIndent();
      }

      if( INCLUDE_MEMBERS)
      {
        // add an indent for each of the sub-apps
        if( INCLUDE_MEMBER_DESCRIPTIONS && writer.GetLinePosition())
        {
          writer.AddIndent( 2);
          writer << ". Member applications:";
          writer.NewLine();
        }
        else
        {
          if( writer.GetLinePosition())
          {
            writer.NewLine();
          }
          writer.SetIndent( s_LongestGroupSize + size_t( 2));
          writer << "Applications: ";
        }
        // output each group member
        for
        (
          storage::Map< std::string, std::string>::const_iterator
            itr( m_NameToDescription.Begin()), itr_end( m_NameToDescription.End());
          itr != itr_end;
          ++itr
        )
        {
          if( itr->first == "Help")
          {
            // skip help
            continue;
          }
          // write the name, followed by the spaces
          writer << itr->first;

          if( INCLUDE_MEMBER_DESCRIPTIONS)
          {
            writer << std::string( name_indent - itr->first.size(), ' ');

            // set indent for the description
            writer.SetIndent( writer.GetLinePosition());
            writer << itr->second;
            if( INCLUDE_MEMBER_ALIASES)
            {
              // write out the aliases for each member
              const storage::Map< std::string, std::string>::const_iterator itr_alias( m_NameToAliases.Find( itr->first));
              if( itr_alias != m_NameToAliases.End() && !itr_alias->second.empty())
              {
                writer << ". ";
                writer.WriteOnOneLine( itr_alias->second);
                if( writer.GetLinePosition())
                {
                  writer.NewLine();
                }
              }
            }
            if( writer.GetLinePosition())
            {
              writer.NewLine();
            }
            writer.PopIndent();
          }
          else
          {
            writer << ", ";
          }
        }
      }
      if( INCLUDE_MEMBERS && !INCLUDE_MEMBER_DESCRIPTIONS)
      {
        // remove any trailing , and newlines
        OSTREAM << util::RStrip( writer.String(), "\n, ");
      }
      else
      {
        // write the string directly
        OSTREAM << writer.String();
      }

      return OSTREAM;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &GroupHandler::Read( std::istream &ISTREAM)
    {
      BCL_Exit( GetClassIdentifier() + " cannot be read", -1);
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    //! @see ParameterCheckInterface::Write
    std::ostream &GroupHandler::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      BCL_Exit( GetClassIdentifier() + " cannot be written", -1);
      // return the stream
      return OSTREAM;
    }

  } // namespace app
} // namespace bcl
