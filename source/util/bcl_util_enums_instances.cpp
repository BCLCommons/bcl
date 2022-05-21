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
#include "util/bcl_util_enums_instances.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_parameter.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

  //////////
  // data //
  //////////

    //! @brief GetFlagEnumsFiles gives command line flag to pass custom enums definitions
    //! @return ShPtr to a FlagInterface which is used to pass filenames for custom enum definitions
    ShPtr< command::FlagInterface> &EnumsInstances::GetFlagEnumsFiles()
    {
      static ShPtr< command::FlagInterface> s_flag_enums_files
      (
        new command::FlagDynamic
        (
          "enums_files",
          "files for enum data that adds enums or overrides data of existing enum data",
          command::Parameter
          (
            "enum_file",
            "file that is similar to a written Enums derived class"
          ),
          0,
          EnumsInstances::GetEnumsInstances().GetSize(),
          &EnumsInstances::AddEnumsFromCommandline
        )
      );

      return s_flag_enums_files;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    EnumsInstances::EnumsInstances() :
      m_EnumsInstances()
    {
    }

    //! @brief return the number of instances
    size_t EnumsInstances::GetSize() const
    {
      return m_EnumsInstances.size();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief add enums from files that are given in the commandline
    void EnumsInstances::AddEnumsFromCommandline()
    {
      // for each filename given in commandline
      for
      (
        ShPtrVector< command::ParameterInterface>::const_iterator
          itr_param( GetFlagEnumsFiles()->GetParameterList().Begin()),
          itr_param_end( GetFlagEnumsFiles()->GetParameterList().End());
        itr_param != itr_param_end;
        ++itr_param
      )
      {
        // create and open stream to given file
        io::IFStream read;
        io::File::MustOpenIFStream( read, ( *itr_param)->GetValue());

        const std::string enums_name( ObjectInterface::ExtractIdentifier( read));

        // find Enums derived class with given identifier
        iterator itr_enum( GetEnumsInstances().m_EnumsInstances.find( enums_name));

        if( itr_enum == GetEnumsInstances().m_EnumsInstances.end())
        {
          BCL_MessageCrt
          (
            "Enums instance with name \"" + enums_name + "\" given in file \"" + ( *itr_param)->GetValue() +
            " \" could not be found as instance. It might be that this Enums derived class is not instantiated since"
            " Singleton constructor was not called yet!"
          );
        }
        else
        {
          // read enums to that
          itr_enum->second->Read( read);
        }

        // close and clear stream
        io::File::CloseClearFStream( read);
      }
    }

    //! @brief the one and only instance of this class
    EnumsInstances &EnumsInstances::GetEnumsInstances()
    {
      static EnumsInstances s_enums_insances;

      return s_enums_insances;
    }

    //! @brief insert an instance of an Enums derived class
    void EnumsInstances::Insert( ObjectInterface &ENUMS_OBJECT, const std::string &ENUMS_NAME)
    {
      // check if Enums derived class only has one instance
      const_iterator itr( m_EnumsInstances.find( ENUMS_NAME));

      if( itr != m_EnumsInstances.end())
      {
        BCL_Exit( "more than one instance of Enums derived class: " + ENUMS_NAME + " is not supported!", -1);
      }

      // insert that instance
      m_EnumsInstances[ ENUMS_NAME] = &ENUMS_OBJECT;

      GetObjectInstances().AddInstanceWithName( &ENUMS_OBJECT, ENUMS_NAME);
    }

  } // namespace util
} // namespace bcl
