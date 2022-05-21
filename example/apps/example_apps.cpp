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
#include "app/bcl_app_apps.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_command.h"
#include "pdb/bcl_pdb_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_apps.cpp
  //! @brief tests that all apps can be run with -help and -readme
  //!
  //! @author mendenjl
  //! @date Nov 26, 2012
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleApps :
    public ExampleInterface
  {
  public:

    ExampleApps *Clone() const
    {
      return new ExampleApps( *this);
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
      // get a reference to all application names
      const storage::Vector< std::string> &app_names( app::GetApps().GetApplicationNames());

      // iterate over all apps
      for
      (
        storage::Vector< std::string>::const_iterator itr( app_names.Begin()), itr_end( app_names.End());
        itr != itr_end;
        ++itr
      )
      {
        // get the associated application
        const app::Interface &application( **app::ApplicationType( *itr));

        // do not run this application for the example application, otherwise any example flags that were set externally
        // will be ignored
        if( application.GetLicensedName() == "bcl:Examples")
        {
          continue;
        }

        BCL_MessageStd( "Trying to call initialize command for " + *itr);

        // get the initial command from the application
        util::ShPtr< command::Command> command( application.InitializeCommand());
        BCL_ExampleIndirectCheck
        (
          command.IsDefined(),
          true,
          "Check that " + *itr + ".InitializeCommand() returned a valid ShPtr"
        );

        BCL_MessageStd( "Trying to call help for " + *itr);

        std::stringstream help_output;
        command->WriteHelp( help_output);

        // for release apps, check that GetReadMe was overridden
        if( application.IsReleaseApplication())
        {
          BCL_ExampleIndirectCheck
          (
            application.GetReadMe() != app::Interface::DefaultReadMe(),
            true,
            "release app overrides GetReadMe"
          );
        }
      }

      // reset all pdb factory flags
      // this is required because several applications change them via SetDefaultParameter
      pdb::Factory::ResetFlagDefaults();

      // this
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleApps

  const ExampleClass::EnumType ExampleApps::s_Instance( GetExamples().AddEnum( ExampleApps()));

} // namespace bcl
