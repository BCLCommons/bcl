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
