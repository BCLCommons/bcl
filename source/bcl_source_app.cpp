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
#include "app/bcl_app_apps.h"

// includes from bcl - sorted alphabetically
#include "bcl_version.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "storage/bcl_storage_template_instantiations.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {

  //////////
  // data //
  //////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Apps::Apps() :
      util::Enumerate< util::ShPtr< Interface>, Apps>( false),
      m_ReleaseApplicationsNames()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Apps::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return vector of application names
    //! @return StorageVector<std::string> with application name
    const storage::Vector< std::string> &Apps::GetApplicationNames() const
    {
      static storage::Vector< std::string> s_app_names( GetEnumStrings());
      static bool s_sorted( false);

      // if this function was called before static initialization occurred, it is possible that
      // the number of applications has increased
      if( s_app_names.GetSize() != GetEnumCount())
      {
        BCL_MessageStd
        (
          "app::Apps::GetApplicationNames() should not normally be called before static initialization is performed"
        );
        s_app_names = GetEnumStrings();
        s_sorted = false;
      }

      // sort once
      if( !s_sorted)
      {
        std::sort( s_app_names.Begin(), s_app_names.End());
        s_sorted = true;
      }

      // return
      return s_app_names;
    }

    //! @brief return vector of release application names
    //! @return StorageVector<std::string> with release application name
    const storage::Vector< std::string> &Apps::GetReleaseApplicationNames() const
    {
      // return
      return m_ReleaseApplicationsNames;
    }

    //! @brief return a parameter with all application strings
    //! @return command::Parameter with check that contains all possible application name strings
    util::ShPtr< command::ParameterInterface> &Apps::GetApplicationsParameter()
    {
      // create Parameter with all application names
      static util::ShPtr< command::ParameterInterface> s_apps_parameter
      (
        new command::Parameter
        (
          "application_group>:<application",
          "name of the bcl application to be executed",
          command::ParameterCheckEnumerate< Apps>()
        )
      );

      // return
      return s_apps_parameter;
    }

    //! @brief return the executable path
    //! @return the actual executable path
    std::string &Apps::GetExecutablePath()
    {
      static std::string s_executable_path;
      return s_executable_path;
    }

    //! @brief add new Interface derived Application to Apps
    ApplicationType &Apps::AddEnum
    (
      const Interface &APPLICATION
    )
    {
      return this->AddEnum( APPLICATION.GetNameForGroup(), util::CloneToShPtr( APPLICATION));
    }

    //! @brief add the app directly
    //! @param NAME           name of the current application
    //! @param SP_APPLICATION object to be enumerated
    ApplicationType &Apps::AddEnum( const std::string &NAME, const util::ShPtr< Interface> &SP_APPLICATION)
    {
      // add release applications to the release app vector
      if( SP_APPLICATION->IsReleaseApplication())
      {
        if( m_ReleaseApplicationsNames.Find( NAME) >= m_ReleaseApplicationsNames.GetSize())
        {
          m_ReleaseApplicationsNames.PushBack( NAME);
          m_ReleaseApplicationsNames.Sort( std::less< std::string>());
        }
      }

      return util::Enumerate< util::ShPtr< Interface>, Apps>::AddEnum( NAME, SP_APPLICATION);
    }

    //! @brief writes the list of enums
    //! @param OSTREAM the stream to which the help is written to
    //! @return the given stream to which the list was written to
    //! Virtual to allow derived classes alter how the help is displayed without overriding Enum
    std::ostream &Apps::WriteList( std::ostream &OSTREAM) const
    {
      // defer the call to app groups
      return GetAppGroups().WriteList( OSTREAM);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief writes the help for the commandline
    //! @param OSTREAM - the outstream to write to
    //! @return std::ostream &OSTREAM - return the stream after writing
    std::ostream &Apps::WriteGenericHelp( std::ostream &OSTREAM)
    {
      if( !GetApps().GetApplicationsParameter()->GetWasSetInCommandLine())
      {
        OSTREAM << "BCL Help\n";
        OSTREAM << GetVersion().GetDetailedInfoString() << '\n';
        OSTREAM << "Usage: " << GetExecutablePath() << " <application_group>:<application> [Other parameters] [Flags] [@filenames]\n";
        GetApps().GetApplicationsParameter()->WriteHelp( OSTREAM, 1);

        command::Command cmd;
        command::GetAppDefaultFlags().AddRequiredCommandlineFlags( cmd);
        cmd.WriteHelp( OSTREAM);
      }
      else
      {
        const std::string app_name( GetApps().GetApplicationsParameter()->GetValue());
        ApplicationType app( GetApps().GetApplicationsParameter()->GetValue());
        util::ShPtr< Interface> app_copy( app->HardCopy());
        util::ShPtr< command::Command> sp_cmd( app_copy->InitializeCommand());

        OSTREAM << '\n' << command::Command::DefaultSectionSeparator() << '\n';
        OSTREAM << app.GetName() << " Help\n";
        OSTREAM << GetVersion().GetDetailedInfoString() << '\n';
        OSTREAM << "Usage: " << GetExecutablePath() << ' ' << app.GetName() << ' ';
        sp_cmd->WriteUsage( util::GetLogger());
        sp_cmd->WriteHelp( util::GetLogger());
      }
      OSTREAM
       << "additional arguments (including the application) can be loaded from files by passing e.g. @my_file.txt"
       << "\n  If no application is specified, attempts to read commands from " << GetDefaultCommandLineArgumentsFile()
       << '\n';
      return OSTREAM;
    }

    //! @brief get the default command line
    //! If no arguments are given to an application, this file will be opened and its contents will be used as a command
    const std::string &Apps::GetDefaultCommandLineArgumentsFile()
    {
      static const std::string s_default_command_line_arguments_file( "bcl_commands.txt");
      return s_default_command_line_arguments_file;
    }

    Apps &GetApps()
    {
      return Apps::GetEnums();
    }

  } // namespace app

  namespace util
  {
  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< app::Interface>, app::Apps>;

  } // namespace util
} // namespace bcl
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
#include "app/bcl_app.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {

    //! @brief identifier for the name space
    //! @return the name of the namespace
    const std::string &GetNamespaceIdentifier()
    {
      static const std::string *s_namespace_name( new std::string( util::ExtractNamespaceIdentifier( __PRETTY_FUNCTION__)));
      return *s_namespace_name;
    }

  } // namespace app
} // namespace bcl
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
#include "app/bcl_app_groups.h"

// includes from bcl - sorted alphabetically
#include "bcl_version.h"
#include "io/bcl_io_fixed_line_width_writer.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {

  //////////
  // data //
  //////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Groups::Groups() :
      util::Enumerate< GroupHandler, Groups>( false),
      e_Bcl(          AddGroup( GroupHandler( "bcl", "applications that use only BCL-generated input"))),
      e_Protein   (   AddGroup( GroupHandler( "protein", "applications that primarily use or create proteins (usually PDB files)"))),
      e_Sequence  (   AddGroup( GroupHandler( "sequence", "applications that primarily use AA sequences (usually FASTA files)"))),
      e_Molecule  (   AddGroup( GroupHandler( "molecule", "applications that primarily use or create molecules (usually SDF files)"))),
      e_Descriptor(   AddGroup( GroupHandler( "descriptor", "applications that primarily use or create descriptor objects"))),
      e_Model     (   AddGroup( GroupHandler( "model", "applications that primarily use or create machine learning models"))),
      e_ChemInfo  (   AddGroup( GroupHandler( "cheminfo", "applications that primarily use or create molecules using machine learning models"))),
      e_BioInfo(      AddGroup( GroupHandler( "bioinfo", "applications that primarily use AA sequences and machine learning models"))),
      e_Restraint (   AddGroup( GroupHandler( "restraint", "applications that primarily use or create restraints that guide protein folding"))),
      e_Density   (   AddGroup( GroupHandler( "density", "applications that primarily use or create density maps from sequence or protein models"))),
      e_InternalBiol( AddGroup( GroupHandler( "bioutil", "bcl development and internal analysis tools for proteins or sequences"))),
      e_InternalChem( AddGroup( GroupHandler( "molutil", "bcl development and internal analysis tools for molecules"))),
      e_Utility(      AddGroup( GroupHandler( "util", "utilities used primarily for bcl development and internal analysis"))),
      e_All(          AddGroup( GroupHandler( "all", "")))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Groups::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief add a particular group to the app
    //! @param GROUP application group of interest
    //! @return the created enum
    ApplicationGroup &Groups::AddGroup( const GroupHandler &GROUP)
    {
      return AddEnum( GROUP.GetName(), GROUP);
    }

    //! @brief add new Interface derived application to a single application group
    //! @param APPLICATION the newly-created application
    //! @param GROUP the group to add it to
    //! @return the newly created application enum
    ApplicationType &Groups::AddAppToGroup( Interface *const &APPLICATION, ApplicationGroup &GROUP)
    {
      if( GetVersion().IsLicense() && !APPLICATION->IsReleaseApplication())
      {
        static ApplicationType s_undefined_app;
        // do not add non-release apps to the enum; this keeps internal apps invisible to external users
        return s_undefined_app;
      }
      return GROUP->AddInstance( util::ShPtr< Interface>( APPLICATION));
    }

    //! @brief add new Interface derived application to several groups
    //! @param APPLICATION the newly-created application
    //! @param GROUPS the set of groups to add it to
    //! @return the newly created application enum
    ApplicationType &Groups::AddAppToGroups
    (
      Interface *const &APPLICATION,
      const storage::Vector< ApplicationGroup> &GROUPS
    )
    {
      // handle if no groups were provided
      BCL_Assert( !GROUPS.IsEmpty(), "Tried to add an application without specifying any groups");

      // create a shared pointer with the new app
      util::ShPtr< Interface> sp_app( APPLICATION);

      // handle all but the last group
      for
      (
        storage::Vector< ApplicationGroup>::const_iterator itr( GROUPS.Begin()), itr_end( GROUPS.End() - 1);
        itr != itr_end;
        ++itr
      )
      {
        ApplicationGroup( *itr)->AddInstance( sp_app);
      }

      // add it to the last-specified group
      return ApplicationGroup( GROUPS.LastElement())->AddInstance( sp_app);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief writes the list of enums
    //! @param OSTREAM the stream to which the help is written to
    //! @return the given stream to which the list was written to
    //! Virtual to allow derived classes alter how the help is displayed without overriding Enum
    std::ostream &Groups::WriteList( std::ostream &OSTREAM) const
    {
      OSTREAM << '\n';
      io::FixedLineWidthWriter writer( 0, util::GetLogger().GetMaxLineWidth() - 4);
      const size_t description_indent( GroupHandler::GetLengthLongestAppGroupName() + 2);
      writer << "<group>:";
      writer << std::string( description_indent - std::min( writer.GetLinePosition(), description_indent), ' ') << "Description" << '\n';
      writer << "<group>:Help";
      writer << std::string( description_indent - std::min( writer.GetLinePosition(), description_indent), ' ');
      writer.SetIndent( description_indent);
      writer << "prints all <application>'s and descriptions for the group";
      OSTREAM << writer.String() << std::flush;

      for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        ( *itr)->WriteHelp( OSTREAM, 0, 2, true, true, false);
        OSTREAM << std::flush;
      }
      e_Undefined->WriteHelp( OSTREAM, 0, 2, true, true, false);
      return OSTREAM;
    }

    //! @brief find all the groups that a particular app belongs to
    //! @param APP the application of interest
    //! @return all application group names that the given app belongs to
    storage::Vector< std::string> Groups::GetApplicationGroupsForApp( const Interface &APP)
    {
      storage::Vector< std::string> app_groups;
      for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        if( ( *itr)->Contains( APP) && !( *itr)->GetName().empty())
        {
          app_groups.PushBack( ( *itr)->GetName());
        }
      }
      return app_groups;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    Groups &GetAppGroups()
    {
      return Groups::GetEnums();
    }

  } // namespace app

  namespace util
  {
  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< app::GroupHandler, app::Groups>;

  } // namespace util
} // namespace bcl
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

// include the namespace header
#include "app/bcl_app.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "app/bcl_app_apps.h"
#include "command/bcl_command_command.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Help
    //! @brief prints help for a given application group if any was given
    //!
    //! @author mendenjl
    //! @date Mar 5, 2013
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Help :
      public Interface
    {

    public:

      // instantiate enumerator for Help class
      static const ApplicationType Help_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new Interface
      Help *Clone() const
      {
        return new Help( *this);
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

      //! @brief return the Command object
      util::ShPtr< command::Command> InitializeCommand() const
      {
        util::ShPtr< command::Command> sp_cmd( new command::Command);

        return sp_cmd;
      }

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const
      {
        return "Lists descriptions of all applications in group";
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const
      {
        // get the group that was chosen over the command line
        const std::string full_app_name
        (
          GetApps().GetApplicationsParameter()->GetValue()
        );

        // split it, to get the application group
        const size_t group_end_pos( full_app_name.find( GroupHandler::s_GroupDelimiter));
        const std::string group_name
        (
          group_end_pos == std::string::npos
          ? std::string()
          : full_app_name.substr( 0, group_end_pos)
        );

        // check for valid group name
        ApplicationGroup real_group( group_name);

        // check whether all group help was requested
        if( real_group == GetAppGroups().e_All)
        {
          util::GetLogger() << "All applications with descriptions\n";
          util::GetLogger() << "Format:\n" << "<group>:    Group description\n";
          const std::string app_name_str( "  <application_name>");
          util::GetLogger() << app_name_str
                            << std::string( GroupHandler::GetLengthLongestAppName() - app_name_str.size(), ' ')
                            << "Description\n"
                            << DefaultSectionSeparator();
          for
          (
            Groups::const_iterator itr( GetAppGroups().Begin()), itr_end( GetAppGroups().End());
            itr != itr_end;
            ++itr
          )
          {
            ( *itr)->WriteHelp( util::GetLogger());
          }
          GetAppGroups().e_Undefined->WriteHelp( util::GetLogger());
        }
        else
        {
          util::GetLogger() << "Members of Group " << group_name << ": " << std::endl;
          // just write the help for the given application group
          real_group->WriteHelp( util::GetLogger(), 0, real_group->GetLengthLongestAppNameThisGroup());
        }

        return 0;
      }

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

    }; // Help

    // instantiate enumerator for Help class
    const ApplicationType Help::Help_Instance
    (
      GetAppGroups().AddAppToGroups
      (
        new Help(),
        storage::Vector< ApplicationGroup>( GetAppGroups().Begin(), GetAppGroups().End())
      )
    );

  } // namespace app
} // namespace bcl
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
#include "app/bcl_app_interface.h"

// includes from bcl - sorted alphabetically
#include "app/bcl_app_group_handler.h"
#include "command/bcl_command_command.h"
#include "util/bcl_util_loggers.h"
#include "util/bcl_util_time.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {

  //////////
  // data //
  //////////

    //! e-mail address for commercial license support
    const std::string &Interface::GetSupportEmailCommercial()
    {
      static const std::string s_support_email_commercial( "bcl-support-commercial@meilerlab.org");
      return s_support_email_commercial;
    }

    //! e-mail address for academic license support
    const std::string &Interface::GetSupportEmailAcademic()
    {
      static const std::string s_support_email_academic(   "bcl-support-academic@meilerlab.org");
      return s_support_email_academic;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief the expiration interval from last compilation
    //! @return the time interval from the last compilation, the license is valid - by default, it is around 20 years
    const util::Time &Interface::GetLicenseExpirationTime() const
    {
      // every 14 months which turns out to be about 425 days
      static const util::Time s_default_license_expiration( 425 * 24 * 60 * 60, 0);

      // by default, it returns 425 days
      return s_default_license_expiration;
    }

    //! @brief return the bcl::commons name
    //! @return string for the bcl::commons name of that application
    std::string Interface::GetBCLScopedName() const
    {
      return "BCL" + this->GetClassIdentifier().substr( 8);
    }

    //! @brief return the original name of the application, as it was used in license files
    //! @return string for the bcl::commons name of that application
    //! This is necessary so that, if release application names change, licenses will continue to work
    std::string Interface::GetLicensedName() const
    {
      // by default, return the app name
      return GetNameForGroup();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns ReadMe information
    //! @return string containing information about application
    const std::string &Interface::GetReadMe() const
    {
      // read me message
      return DefaultReadMe();
    }

    //! @brief returns web text information
    //! @return text (html allowed but not required) that will be displayed on the website
    //! @note use confluence syntax, e.g. !my_image.png! to refer to an image under documentation/image
    const std::string &Interface::GetWebText() const
    {
      return GetReadMe();
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief get the name of the application, excluding the group name (if present
    //! @param GROUP_NAME name of the group to exclude from the beginning of teh applications name
    //! @return the application name if it were part of the group in question
    std::string Interface::GetNameForGroup( const std::string &GROUP_NAME) const
    {
      // generate the string that needs to be removed e.g. "const-bcl::app::"
      static const std::string s_obsolete_string( GetNamespaceIdentifier() + "::");

      // get the full app name
      std::string full_app_name_no_group( GetClassIdentifier().substr( s_obsolete_string.length()));

      // return it if group name was empty
      if( GROUP_NAME.empty())
      {
        return full_app_name_no_group;
      }

      // otherwise, check whether the application name and group name start off the same, ignoring case
      if( util::StartsWith( util::ToLower( full_app_name_no_group), util::ToLower( GROUP_NAME)))
      {
        // erase the group name from the full app name
        full_app_name_no_group.erase( 0, GROUP_NAME.size());
      }

      // return the group:full_app_name_no_group
      return GROUP_NAME + GroupHandler::s_GroupDelimiter + full_app_name_no_group;
    }

    //! @brief returns the string that can be used as default for the technical support section in the readme
    std::string Interface::DefaultTechnicalSupportString() const
    {
      return
        "Technical support is available throughout the license period.  For scientific or technical questions related "
        "to " + GetBCLScopedName() + ", commercial entities should contact us at " + GetSupportEmailCommercial() + " while "
        "academic institutions should contact us at " + GetSupportEmailAcademic() + ".  If there is a major bug fix or "
        "significant improvements to the application, entities with current licenses will be notified via email, and "
        "will be given access to the updated software.\n";
    }

    //! @brief returns the string that can be used as default for the terms of use section in the readme
    //! @return line break terminated terms of use string
    const std::string &Interface::DefaultTermsOfUseString()
    {
      static const std::string s_terms_of_use
      (
        "Licensees are bound by the license agreement in place at the time this software was obtained, unless legally "
        "modified by the parties. The following brief summary does not replace the original agreement. It is expected "
        "that the software will be used for its intended purpose and that appropriate citations and acknowledgments "
        "will be given in connection to this software.  The license provides access to the software, updates, and "
        "reasonable technical support for one year, after which time renewal will be required for continued use and "
        "support.  The license allows the licensee to have the software installed only on site(s) designated at the "
        "time the license was obtained, usually a single personal computer (PC). The license does not allow the "
        "software to be passed to third parties or to be modified in any way.\n"
      );

      return s_terms_of_use;
    }

    //! @brief returns the string that can be used as default for the installation procedure section in the readme
    //! @return line break terminate installation procedure
    const std::string &Interface::DefaultInstallationProcedure()
    {
      static const std::string s_installation_procedure
      (
        "Run the installer and follow the directions. On Linux, if a \"lib\" folder is present in the installation "
        "directory, it will be necessary to add that with its full path to the LD_LIBRARY_PATH environment variable, "
        "in order for the dynamic linker (ldd) is able to find all libraries to run the application. Please be also "
        "aware of the licenses for the external libraries, which are installed as well.\n"
        "Uninstallers are provided for some operating systems.\n"
        "If there are major bug fixes or revisions, an email will be sent out to entities with current licenses so "
        "that the newest application can be installed.\n"
      );

      return s_installation_procedure;
    }

    //! @brief returns the string seperating sections within a readme
    //! @return line break terminated separator string
    const std::string &Interface::DefaultSectionSeparator()
    {
      return command::Command::DefaultSectionSeparator();
    }

    //! @brief returns the default read me
    //! @return the default read me
    const std::string &Interface::DefaultReadMe()
    {
      static const std::string s_read_me( "no ReadMe information is provided for this application");
      return s_read_me;
    }

  } // namespace app
} // namespace bcl
