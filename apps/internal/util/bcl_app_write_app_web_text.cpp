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
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_parameter_check_allowed_non_const.h"
#include "storage/bcl_storage_template_instantiations.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class WriteAppWebText
    //! App used to print off the web-text for all release applications; used to update the bcl apps on the MeilerLab
    //! website whenever releases are made
    //!
    //! @author mendenjl
    //! @date Feb 15, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API WriteAppWebText :
      public Interface
    {

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      WriteAppWebText *Clone() const
      {
        return new WriteAppWebText( *this);
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

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const
      {
        return "Writes out all information used by the UpdateBclWebText.py script to update meilerlab.org/bclcommons";
      }

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! Main
      int Main() const;

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

    private:

      static const ApplicationType WriteAppWebText_Instance;

    }; // class WriteAppWebText

    //! Main
    int WriteAppWebText::Main() const
    {
      // Iterate over all application groups

      util::GetLogger() << "<?xml version=\"1.0\"?>\n"
                        << "<groups>\n";
      for
      (
        Groups::const_iterator itr( GetAppGroups().Begin()), itr_end( GetAppGroups().End());
        itr != itr_end;
        ++itr
      )
      {
        util::GetLogger() << "  <group>" << ( *itr)->GetName() << "</group>\n"
                          << "  <description>" << ( *itr)->GetDescription() << "</description>\n";
      }
      util::GetLogger() << "</groups>\n";

      util::GetLogger() << "<apps>\n";
      storage::Set< util::ShPtr< Interface> > all_apps;
      for( Apps::const_iterator itr( GetApps().Begin()), itr_end( GetApps().End()); itr != itr_end; ++itr)
      {
        if( !( **itr)->IsReleaseApplication())
        {
          // skip unreleased apps
          continue;
        }

        // skip apps that were already written out
        if( !all_apps.Insert( **itr).second)
        {
          continue;
        }

        const Interface &app( ***itr);

        // get all groups that this app is part of
        const storage::Vector< std::string> app_groups( GetAppGroups().GetApplicationGroupsForApp( app));

        // skip apps that are not part of any group
        if( app_groups.IsEmpty())
        {
          continue;
        }

        util::GetLogger() << "  <app>\n"
                          << "    <name>" << app.GetNameForGroup() << "</name>\n"
                          << "    <classname>" << app.GetClassIdentifier() << "<classname>\n"
                          << "    <licensed_name>" << app.GetLicensedName() << "</licensed_name>\n"
                          << "    <deprecated_names>\n";
        if( !app.GetDeprecatedAppNames().IsEmpty())
        {
          util::GetLogger() << "      <deprecated>"
                            << util::Join( "</deprecated>\n      <deprecated>", app.GetDeprecatedAppNames())
                            << "</deprecated>\n";
        }
        util::GetLogger() << "    </deprecated_names>\n"
                          << "    <description>" << app.GetDescription() << "</description>\n"
                          << "    <readme>" << app.GetReadMe() << "</readme>\n"
                          << "    <webtext>" << app.GetWebText() << "</webtext>\n";
        util::GetLogger() << "    <app_groups>\n";
        for
        (
          storage::Vector< std::string>::const_iterator
            itr_group( app_groups.Begin()), itr_group_end( app_groups.End());
          itr_group != itr_group_end;
          ++itr_group
        )
        {
          util::GetLogger() << "      <app_group>\n"
                            << "        <app_group_name>" << *itr_group << "</app_group_name>\n"
                            << "        <app_name_in_group>" << app.GetNameForGroup( *itr_group)
                            << "</app_name_in_group>\n"
                            << "      </app_group>\n";
        }
        util::GetLogger() << "    </app_groups>\n";
        util::GetLogger() << "  </app>\n";
      }
      util::GetLogger() << "</apps>\n";

      // end
      return 0;
    }

    util::ShPtr< command::Command> WriteAppWebText::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // add default bcl parameters
      command::GetAppDefaultFlags().AddRequiredCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    const ApplicationType WriteAppWebText::WriteAppWebText_Instance
    (
      GetAppGroups().AddAppToGroup( new WriteAppWebText(), GetAppGroups().e_Utility)
    );

  } // namespace app
} // namespace bcl
