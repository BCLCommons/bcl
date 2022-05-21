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
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_allowed_non_const.h"
#include "io/bcl_io_file.h"
#include "storage/bcl_storage_template_instantiations.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class WriteObject
    //! @brief TODO: add brief class comment
    //!
    //! @author woetzen
    //! @date Oct 31, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API WriteObject :
      public Interface
    {
    private:

    //////////
    // data //
    //////////

      //! output path
      util::ShPtr< command::FlagInterface> m_PrefixFlag;

      //! choice of objects to be written
      util::ShPtr< command::FlagInterface> m_DefaultObjectsToBeWritten;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! @brief default constructor
      WriteObject();

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      WriteObject *Clone() const
      {
        return new WriteObject( *this);
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
        return "Write a desired bcl object for use in other applications";
      }

      // This function is used to retrieve prior names for applications that are still valid (for backwards compatibility)
      // but which should not be displayed, e.g. for help
      storage::Vector< std::string> GetDeprecatedAppNames() const
      {
        return storage::Vector< std::string>( size_t( 1), "WriteBCLObject");
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

      static const ApplicationType WriteObject_Instance;

    }; // class WriteObject

    //! Main
    int WriteObject::Main() const
    {
      // generate a list of of all enums to be written
      const storage::Vector< std::string> enmus_commandline
      (
        m_DefaultObjectsToBeWritten->GetStringList()
      );

      // iterate over all given objects that are supposed to be written into files
      for
      (
        storage::Vector< std::string>::const_iterator
          itr( enmus_commandline.Begin()), itr_end( enmus_commandline.End());
        itr != itr_end;
        ++itr
      )
      {
        io::OFStream write;

        util::SiPtr< const util::ObjectInterface> object( GetObjectInstances().GetPtrToObjectFromIdentifier( *itr));
        io::File::MustOpenOFStream
        (
          write,
          io::File::ConvertClassIdentifierToFilename
          (
            m_PrefixFlag->GetFirstParameter()->GetValue() + object->GetClassIdentifier() + ".bcl"
          )
        );
        write << *object;
        io::File::CloseClearFStream( write);
      }

      // end
      return 0;
    }

    util::ShPtr< command::Command> WriteObject::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // input tag of sequence
      sp_cmd->AddFlag( m_PrefixFlag);

      //input density filename
      sp_cmd->AddFlag( m_DefaultObjectsToBeWritten);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief default constructor
    WriteObject::WriteObject() :
      m_PrefixFlag
      (
        new command::FlagStatic
        (
          "prefix",
          "prefix all files will get before they are written",
          command::Parameter
          (
            "output_prefix",
            "prefix that is prepended to the actual filename",
            ""
          )
        )
      ),
      m_DefaultObjectsToBeWritten
      (
        new command::FlagDynamic
        (
          "write_default_object",
          "give list of class names, and default objects will be written to \"bcl_namespace_classname.bcl\"",
          command::Parameter
          (
            "class_name",
            "name of bcl class for which the file will be written",
            command::ParameterCheckAllowedNonConst( GetObjectInstances().GetKnownObjectNames())
          ),
          0, 1000
        )
      )
    {
    }

    const ApplicationType WriteObject::WriteObject_Instance
    (
      GetAppGroups().AddAppToGroup( new WriteObject(), GetAppGroups().e_Utility)
    );

  } // namespace app
} // namespace bcl
