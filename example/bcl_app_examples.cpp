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
// include the interface for all apps
#include "app/bcl_app_apps.h"

// includes from bcl - sorted alphabetically
#include "example_application_example_helper.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_allowed_non_const.h"
#include "command/bcl_command_parameter_check_enumerate.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Examples
    //! @brief This class is the application for running all examples/test cases.
    //!
    //! @author woetzen
    //! @date 08/07/2008
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class Examples :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      //! flag to enable execution of all examples
      util::ShPtr< command::FlagInterface> m_AllFlag;

      //! flag to run Examples within a certain namespace(s)
      util::ShPtr< command::FlagInterface> m_NamespaceList;

      //! flag to run Examples while excluding certain namespace(s)
      util::ShPtr< command::FlagInterface> m_NamespaceExcludeList;

      //! flag with dynamic parameterlist for example to be executed if "-all" was not given,
      //! or not to be executed if "-all" was given
      util::ShPtr< command::FlagInterface> m_ExampleList;

      //! set of namespaces to ignore by default
      storage::Set< std::string>           m_NamespacesExcludedDefault;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Examples();

    public:

      //! @brief Clone function
      //! @return pointer to new Examples
      Examples *Clone() const
      {
        return new Examples( *this);
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
        return "Unit and integration tests";
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

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

      static const ApplicationType Examples_Instance;

    }; // Examples

    // initialize the command object
    util::ShPtr< command::Command> Examples::InitializeCommand() const
    {
      // initialize ShPtr to a new command object
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      sp_cmd->AddFlag( m_NamespaceList);
      sp_cmd->AddFlag( m_NamespaceExcludeList);

      //pushback all parameters
      sp_cmd->AddFlag( m_ExampleList);
      sp_cmd->AddFlag( m_AllFlag);
      sp_cmd->AddFlag( ExampleClass::GetExamplePathFlag());
      sp_cmd->AddFlag( ApplicationExampleHelper::GetApplicationExamplePathFlag());

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      //end
      return sp_cmd;
    }

    //! procedure to execute all examples
    int Examples::Main() const
    {
      // set the path for application examples, if it was given in the command line
      if( ApplicationExampleHelper::GetApplicationExamplePathFlag()->GetFirstParameter()->GetWasSetInCommandLine())
      {
        ApplicationExampleHelper::SetApplicationExampleFilePath
        (
          ApplicationExampleHelper::GetApplicationExamplePathFlag()->GetFirstParameter()->GetValue()
        );
      }

      // list of examples that will be started
      storage::Set< ExampleClass::EnumType> examples_to_run;

      std::string example_string;
      if( m_AllFlag->GetFlag())
      {
        if( m_ExampleList->GetFlag() || m_NamespaceList->GetFlag() || m_NamespaceExcludeList->GetFlag())
        {
          BCL_MessageCrt
          (
            "all flag selected; ignoring namespace flags and example list flag"
          );
        }
        examples_to_run.InsertElements( GetExamples().Begin(), GetExamples().End());
      }
      else if( m_ExampleList->GetFlag() || m_NamespaceList->GetFlag())
      {
        if( m_NamespaceExcludeList->GetFlag())
        {
          BCL_MessageCrt
          (
            "example or namespace list flag selected; ignoring exclude namespace flag"
          );
        }

        // add all examples given in commandline
        examples_to_run = m_ExampleList->GetObjectSet< ExampleClass::EnumType>();

        // add all examples that belong to namespaces given in command line
        if( m_NamespaceList->GetFlag())
        {
          // list of namespaces from command line
          storage::Set< std::string> namespaces( m_NamespaceList->GetObjectSet< std::string>());

          // iterate over all available examples
          for
          (
            ExampleClass::const_iterator example_itr( GetExamples().Begin()), example_itr_end( GetExamples().End());
            example_itr != example_itr_end;
            ++example_itr
          )
          {
            // if this example belongs to one of the namespace that are supposed to be executed, add that to the list
            if( namespaces.Contains( ExampleClass::FirstWord( example_itr->GetName())))
            {
              examples_to_run.Insert( *example_itr);
            }
          }
        }
      }
      // no flag was used - run all examples that are not in an excluded namespace
      else
      {
        // get the excluded namespaces from the flag; then add the examples excluded by default
        storage::Set< std::string> namespaces_excluded( m_NamespaceExcludeList->GetObjectSet< std::string>());
        for
        (
          storage::Set< std::string>::const_iterator
            itr_namespace( m_NamespacesExcludedDefault.Begin()),
            itr_namespace_end( m_NamespacesExcludedDefault.End());
          itr_namespace != itr_namespace_end;
          ++itr_namespace
        )
        {
          namespaces_excluded.Insert( *itr_namespace);
        }

        // iterate over all available examples, add those to the examples to run that are not in an excluded namespace
        for
        (
          ExampleClass::const_iterator example_itr( GetExamples().Begin()), example_itr_end( GetExamples().End());
          example_itr != example_itr_end;
          ++example_itr
        )
        {
          ExampleClass::EnumType example( *example_itr);

          // get the namespace of the example
          const std::string example_namespace( ExampleClass::FirstWord( example.GetName()));

          // if the namespace is not among the excluded namespaces, add the example to the set of examples to run
          if( !namespaces_excluded.Contains( example_namespace))
          {
            examples_to_run.Insert( example);
          }
        }
      }

      // begin executing the examples
      BCL_MessageCrt( "BCL Example | BEGIN: All Examples ====================");

      int run_return_sum( 0);
      //iterate over all instances of examples
      for
      (
        storage::Set< ExampleClass::EnumType>::const_iterator
          itr( examples_to_run.Begin()), itr_end( examples_to_run.End());
        itr != itr_end;
        ++itr
      )
      {
        //start the current example
        BCL_MessageCrt( "BCL Example | BEGIN: " + itr->GetName() + " ====================");

        // set random seed to what was given in command line
        random::GetGlobalRandom().SetSeedFromCommandlineFlag();

        // no need to check that GetClassIdentifier for example was overwritten because it is pure virtual

        // insert new row for that example into the example results
        ExampleClass::GetResults().InsertDefaultTotalRow( itr->GetName());

        // run the example
        run_return_sum += ( **itr)->Run();

        // print end
        BCL_MessageStd( "BCL Example | END  : " + itr->GetName() + " ====================");
      }

      BCL_MessageCrt( "BCL Example | END  : All Examples ====================");

      // output a description for the table
      BCL_MessageStd( "The results table");

      ExampleClass::GetResults().WriteResults( util::GetLogger());

      //end
      return ExampleClass::GetResults().GetNumberErrors() + run_return_sum;
    }

    //! @brief default constructor
    Examples::Examples() :
      m_AllFlag( new command::FlagStatic( "all", "execute all examples")),
      m_NamespaceList
      (
        new command::FlagDynamic
        (
          "namespace",
          "this is a list of namespaces for which all examples are executed)",
          command::Parameter
          (
            "namespace",
            "",
            command::ParameterCheckAllowedNonConst( GetExamples().GetNamespaces())
          ),
          0,
          100
        )
      ),
      m_NamespaceExcludeList
      (
        new command::FlagDynamic
        (
          "exclude_namespace",
          "this is a list of namespaces for which all examples are excluded (App is always excluded except with -all or -namespace App)",
          command::Parameter
          (
            "namespace",
            "",
            command::ParameterCheckAllowedNonConst( GetExamples().GetNamespaces())
          ),
          0,
          100
        )
      ),
      m_ExampleList
      (
        new command::FlagDynamic
        (
          "exec",
          "list of examples to be executed",
          command::Parameter
          (
            "example",
            "",
            command::ParameterCheckEnumerate< ExampleClass>()
          ),
          0,
          1000
        )
      ),
      m_NamespacesExcludedDefault( "App")
    {
    }

    const ApplicationType Examples::Examples_Instance
    (
      GetAppGroups().AddAppToGroup( new Examples(), GetAppGroups().e_Bcl)
    );

  } // namespace app
} // namespace bcl
