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
#include "biol/bcl_biol_rotamer_library.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_types.h"
#include "command/bcl_command_command_state.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> RotamerLibrary::s_Instance
    (
      GetObjectInstances().AddInstance( new RotamerLibrary())
    );

    //! map of loop libraries
    storage::HashMap< std::string, util::ShPtr< RotamerLibrary> > RotamerLibrary::s_RotamerLibraries =
      storage::HashMap< std::string, util::ShPtr< RotamerLibrary> >();

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RotamerLibrary::RotamerLibrary() :
      m_LibraryFileName()
    {
    }

    //! @brief construct from members
    //! @param LIBRARY_FILE_NAME path to the rotamer library
    RotamerLibrary::RotamerLibrary( const std::string &LIBRARY_FILE_NAME) :
      m_LibraryFileName( LIBRARY_FILE_NAME)
    {
      // read in the rotamer library at the given path
      std::ostringstream oss;
      ReadInitializerSuccessHook( util::ObjectDataLabel(), oss);
    }

    //! @brief returns a loop library from the given file name
    //! @param LIBRARY_FILE_NAME path to the rotamer library
    util::ShPtr< RotamerLibrary> RotamerLibrary::CreateRotamerLibrary( const std::string &LIBRARY_FILE_NAME)
    {
      if( command::CommandState::GetGlobalCommandState().IsInStaticInitialization())
      {
        return util::ShPtr< RotamerLibrary>();
      }

      // compute hash key for the given parameters
      const std::string key( LIBRARY_FILE_NAME);

      // if a library with this key does not exist, create one
      if( s_RotamerLibraries.Find( key) == s_RotamerLibraries.End())
      {
        util::ShPtr< RotamerLibrary> sp_library( new RotamerLibrary( LIBRARY_FILE_NAME));
        s_RotamerLibraries[ key] = sp_library;
        return sp_library;
      }
      return s_RotamerLibraries[ key];
    }

    //! @brief copy constructor
    //! @return pointer to a new RotamerLibrary
    RotamerLibrary *RotamerLibrary::Clone() const
    {
      return new RotamerLibrary( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &RotamerLibrary::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &RotamerLibrary::GetAlias() const
    {
      static const std::string s_alias( "RotamerLibrary");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RotamerLibrary::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Library containing rotamers and corresponding probabilities for amino acid side chains.");
      serializer.AddInitializer
      (
        "file path",
        "path to the rotamer library file",
        io::Serialization::GetAgent( &m_LibraryFileName)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    //! @return return code indicating success or failure
    bool RotamerLibrary::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      if( !command::CommandState::IsInStaticInitialization())
      {
        // read the rotamer library file
        io::IFStream lib_file;
        io::File::MustOpenIFStream( lib_file, m_LibraryFileName);
        ReadLibrary( lib_file);
        io::File::CloseClearFStream( lib_file);
      }

      return true;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief read a loop template library from the given input stream
    //! @param ISTREAM input stream from which to read the loop template library
    //! @param ANGLE_BIN_WIDTH bin width for chi angles
    void RotamerLibrary::ReadLibrary( std::istream &ISTREAM, double ANGLE_BIN_WIDTH)
    {
      // read in all lines in the library
      storage::Vector< std::string> lines( util::StringLineListFromIStream( ISTREAM));

      // iterate over all lines and add each valid rotamer to the library
      storage::HashMap< std::string, storage::HashMap< size_t, size_t> > rotamer_probability;
      for( auto line_itr( lines.Begin()), line_itr_end( lines.End()); line_itr != line_itr_end; ++line_itr)
      {
        // skip empty lines or comment lines indicated by a leading # or !
        if( line_itr->empty() || ( *line_itr)[ 0] == '!' || ( *line_itr)[ 0] == '#')
        {
          continue;
        }

        // read in the relevant data for this rotamer
        AAType aa_type;                     // type of this amino acid
        size_t count;                             // number of side chains
        storage::Vector< size_t> rotamers( 4, 0); // rotamer configuration
        storage::Vector< std::string> split_line( util::SplitString( *line_itr, " \t"));
        std::string aa_type_tmp( *split_line[ 0]);
        util::TryConvertFromString( count, *split_line[ 3], util::GetLogger());
        for( size_t rotamer_angles( 0); rotamer_angles < 4; ++rotamer_angles)
        {
          util::TryConvertFromString( rotamers( rotamer_angles), *split_line[ rotamer_angles + 4], util::GetLogger());
        }

        // in this implementation the backbone dihedral angles are ignored and the rotamer conformations are aggregated over
        // all backbone conformations
        size_t key( 0), radix( 1);
        for( size_t i( 0); i < 4; ++i, radix *= 10)
        {
          key += radix * rotamers( i);
        }
        rotamer_probability[ aa_type_tmp][ key] += count;
      }

    }

  } // namespace biol
} // namespace bcl
