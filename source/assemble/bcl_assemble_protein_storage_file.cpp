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
#include "assemble/bcl_assemble_protein_storage_file.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "command/bcl_command_command_state.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "opti/bcl_opti_tracker.h"
#include "pdb/bcl_pdb_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProteinStorageFile::s_Instance
    (
      GetObjectInstances().AddInstance( new ProteinStorageFile())
    );

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief InitializerType as string
    //! @param INIT_TYPE the InitializerType
    //! @return the string for the INIT_TYPE
    const std::string &ProteinStorageFile::GetInitializerTypeDescriptor( const InitializerType &INIT_TYPE)
    {
      static const std::string s_descriptors[] =
      {
        "Attach",
        "Overwrite",
        GetStaticClassName< InitializerType>()
      };
      return s_descriptors[ size_t( INIT_TYPE)];
    }

    //! @brief default storage as given in the storage flag used in the command line
    util::ShPtr< ProteinStorageFile> ProteinStorageFile::GetDefaultStorage()
    {
      return util::ShPtr< ProteinStorageFile>
             (
               new ProteinStorageFile
               (
                 GetDefaultStorageFlag()->GetParameterList()( 0)->GetValue(),
                 ProteinStorageFile::TypeEnum( GetDefaultStorageFlag()->GetParameterList()( 1)->GetValue()),
                 io::StreamBufferClass( GetDefaultStorageFlag()->GetParameterList()( 2)->GetValue())
               )
             );
    }

    //! @brief flag for default ProteinStorage
    //! @return ShPtr to flag that is used for default protein storage
    const util::ShPtr< command::FlagInterface> &ProteinStorageFile::GetDefaultStorageFlag()
    {
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "protein_storage",
          "choice of protein storage",
          util::ShPtrVector< command::ParameterInterface>::Create
          (
            util::ShPtr< command::ParameterInterface>
            (
              new command::Parameter
              (
                "storage_initializer",
                "initializer that is passed to the storage",
                ""
              )
            ),
            util::ShPtr< command::ParameterInterface>
            (
              new command::Parameter
              (
                "init_type",
                "type of initialization - attach to existing storage (no overwrite) or overwrite a storage (will overwrite existing files)",
                command::ParameterCheckSerializable( TypeEnum()),
                GetInitializerTypeDescriptor( e_Attach)
              )
            ),
            util::ShPtr< command::ParameterInterface>
            (
              new command::Parameter
              (
                "compression", "choice of compression",
                command::ParameterCheckEnumerate< io::StreamBufferClasses>(),
                io::GetStreamBufferClasses().e_Uncompressed.GetName()
              )
            )
          )
        )
      );

      // end
      return s_flag;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief initialize with stream buffer class type
    ProteinStorageFile::ProteinStorageFile
    (
      const std::string &INITIALIZER,
      const InitializerType &INIT_TYPE,
      const io::StreamBufferClass &COMPRESSION
    ) :
      m_InitializerType( INIT_TYPE),
      m_CompressionType( COMPRESSION)
    {
      if( !m_CompressionType.IsDefined() && !command::CommandState::IsInStaticInitialization())
      {
        m_CompressionType = io::GetStreamBufferClasses().e_Uncompressed;
      }
      if( !INITIALIZER.empty())
      {
        SetDirectory( INITIALIZER);
      }
    }

    //! @brief Clone function
    //! @return pointer to new ProteinStorageFile
    ProteinStorageFile *ProteinStorageFile::Clone() const
    {
      return new ProteinStorageFile( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinStorageFile::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief sets the pdb factory
    //! @param FACTORY PDB factory to use
    void ProteinStorageFile::SetFactory( const util::ShPtr< pdb::Factory> &FACTORY)
    {
      m_Factory = FACTORY;
    }

    //! @brief initialize the molecule storage
    //! @param INITIALIZER path (including prefix) in which proteins are stored
    //! @param INIT_FLAG flag for type of initialization
    //! @return true if initialize was successful
    bool ProteinStorageFile::SetDirectory( const std::string &INITIALIZER)
    {
      if( !m_CompressionType.IsDefined())
      {
        m_CompressionType = io::GetStreamBufferClasses().e_Uncompressed;
      }

      // set factory
      m_Factory = util::ShPtr< pdb::Factory>( new pdb::Factory());

      m_Directory = util::ShPtr< io::Directory>( new io::Directory( INITIALIZER));

      // check if directory exists
      const bool dir_exists( m_Directory->DoesExist());
      bool make_success( dir_exists);
      if( !dir_exists)
      {
        make_success = m_Directory->Make();
      }

      // end
      return make_success;
    }

    //! @brief get the initializer
    //! @return string used to initialize molecule storage
    const std::string &ProteinStorageFile::GetInitializer() const
    {
      if( m_Directory.IsDefined())
      {
        return m_Directory->GetPath();
      }
      // uninitialized path is the current directory
      static const std::string s_uninitialized;
      return s_uninitialized;
    }

    //! @brief get the initialization type
    //! @return true if existing files in the storage will be overwritten
    bool ProteinStorageFile::WillOverwrite() const
    {
      return m_InitializerType == e_Overwrite;
    }

    //! @brief number of molecules in source
    //! @param SOURCE source of molecule eg. given smallmolecule source like mGluR5 potentiators,
    //!        protein xray structures
    //! @return number of molecules
    size_t ProteinStorageFile::GetSize( const std::string &SOURCE) const
    {
      return GetAllKeys( SOURCE).GetSize();
    }

    //! @brief gets all sources in this storage
    //! @return all sources in this storage
    storage::Set< std::string> ProteinStorageFile::GetAllSources() const
    {
      // directory content
      const storage::List< io::DirectoryEntry> dir_entries
      (
        m_Directory->ListEntries( io::Directory::e_File, "", ( *m_CompressionType)->AddExtension( "." + pdb::GetDefaultFileExtension()))
      );

      // found sources
      storage::Set< std::string> sources;

      // iterate over all entries
      for
      (
        storage::List< io::DirectoryEntry>::const_iterator itr( dir_entries.Begin()), itr_end( dir_entries.End());
        itr != itr_end; ++itr
      )
      {
        // store the name
        const std::string name( io::File::RemoveLastExtension( io::File::RemoveCompressionExtension( itr->GetName())));

        // check for key sizes = s_SmallKeySize (final models)
        if( name.size() > s_SmallKeySize)
        {
          // get the substring
          const std::string key_short( name.substr( name.length() - s_SmallKeySize, s_SmallKeySize));

          // if this key is valid
          if( IsValidKey( key_short))
          {
            // add the rest of the name as the source
            sources.Insert( name.substr( 0, name.length() - s_SmallKeySize));
          }

          // check for key sizes = s_LargeKeySize
          if( name.size() > s_LargeKeySize)
          {
            // get the substring
            const std::string key_long( name.substr( name.length() - s_LargeKeySize, s_LargeKeySize));

            // if this key is valid
            if( IsValidKey( key_long))
            {
              // add the rest of the name as the source
              sources.Insert( name.substr( 0, name.length() - s_LargeKeySize));
            }
          }
        }
      }

      return sources;
    }

    //! @brief get all keys for given source and sort them
    //! @param SOURCE source of molecule eg. given smallmolecule source like mGluR5 potentiators,
    //!        protein xray structures
    //! @return all keys of given source
    storage::Vector< std::string> ProteinStorageFile::GetAllKeys( const std::string &SOURCE) const
    {
      // directory content
      const storage::List< io::DirectoryEntry> dir_entries( GetAllEntries( SOURCE));

      // keys
      storage::Vector< std::string> keys;

      // iterate over all entries
      for
      (
        storage::List< io::DirectoryEntry>::const_iterator itr( dir_entries.Begin()), itr_end( dir_entries.End());
        itr != itr_end; ++itr
      )
      {
        // filename without extension
        const std::string::size_type source_pos( itr->GetName().find( SOURCE));
        std::string key;
        if( m_CompressionType != io::GetStreamBufferClasses().e_Uncompressed)
        {
          key = io::File::RemoveLastExtension( io::File::RemoveCompressionExtension( itr->GetName()));
        }
        else
        {
          key = io::File::RemoveLastExtension( itr->GetName());
        }
        // extract the key
        key.erase( source_pos, SOURCE.length());

        key = opti::Tracker< ProteinModel, double>::RemoveTrackerTag( key);

        // key needs to be valid
        if( !IsValidKey( key))
        {
          continue;
        }

        keys.PushBack( key);
      }

      return keys;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get protein
    //! @param SOURCE prefix of protein eg. "run1"
    //! @param KEY key identifier for specific protein in given source like "00001"
    //! @return shptr to protein of interest, undefined if there is no such protein
    util::ShPtr< ProteinModel> ProteinStorageFile::Retrieve( const std::string &SOURCE, const std::string &KEY) const
    {
      // model
      util::ShPtr< ProteinModel> ptr_model;

      // check the key
      if( !IsValidKey( KEY))
      {
        return ptr_model;
      }

      // filename and check if exists
      const io::DirectoryEntry filename( Filename( SOURCE, KEY));

      // check if file with that name exists
      if( !filename.DoesExist())
      {
        return util::ShPtr< ProteinModel>();
      }

      // retrieve protein model
      ptr_model = util::ShPtr< ProteinModel>( m_Factory->ProteinModelFromPDBFilename( filename.GetFullName()).Clone());
      return ptr_model;
    }

    //! @brief get ensemble of proteins
    //! @param SOURCE corresponds to the prefix for the files e.g "run1"
    //! @return shptr list of proteins from given source
    util::ShPtrList< ProteinModel> ProteinStorageFile::RetrieveEnsemble( const std::string &SOURCE) const
    {
      // acquire all keys for that source and return ensemble
      return RetrieveEnsemble( SOURCE, GetAllKeys( SOURCE));
    }

    //! @brief get ensemble of proteins for given keys
    //! @param SOURCE corresponds to the prefix for the files e.g "run1"
    //! @param KEYS vector of identifiers for specific proteins in given source like "00001", "00002"
    //! @return shptr list of proteins from given source
    util::ShPtrList< ProteinModel> ProteinStorageFile::RetrieveEnsemble
    (
      const std::string &SOURCE,
      const storage::Vector< std::string> &KEYS
    ) const
    {
      util::ShPtrList< ProteinModel> ensemble;

      // iterate over all keys
      for( storage::Vector< std::string>::const_iterator itr( KEYS.Begin()), itr_end( KEYS.End()); itr != itr_end; ++itr)
      {
        util::ShPtr< ProteinModel> current( Retrieve( SOURCE, *itr));
        if( current.IsDefined())
        {
          ensemble.PushBack( current);
        }
        else
        {
          BCL_MessageCrt( "could not find protein " + Filename( SOURCE, *itr));
        }
      }

      // end
      return ensemble;
    }

    //! @brief get ensemble of proteins for given keys which does not map to keys, following plain numbering
    //! @param SOURCE corresponds to the prefix for the files e.g "run1"
    //! @param RANGE range of proteins
    //! @return shptr list of proteins from given source
    util::ShPtrList< ProteinModel> ProteinStorageFile::RetrieveEnsemble
    (
      const std::string &SOURCE,
      const math::Range< size_t> &RANGE
    ) const
    {
      // get all keys
      const storage::Vector< std::string> keys( GetAllKeys( SOURCE));

      // start position
      const size_t start_pos( RANGE.GetMin() + RANGE.GetLeftCondition());

      // width
      const size_t width( RANGE.GetWidth() + 1 - RANGE.GetLeftCondition() - RANGE.GetRightCondition());

      // subvector for the range
      storage::Vector< std::string> subkeys( keys, start_pos, width);

      // return sub keys ensemble
      return RetrieveEnsemble( SOURCE, subkeys);
    }

    //! @brief get any PDB remark lines
    //! @param SOURCE prefix of protein eg. "run1"
    //! @param KEY key identifier for specific protein in given source like "00001"
    //! @return PDB remark lines
    util::ShPtrList< pdb::Line> ProteinStorageFile::RetrieveRemarkLines
    (
      const std::string &SOURCE,
      const std::string &KEY
    ) const
    {
      // filename and check if exists
      const io::DirectoryEntry filename( Filename( SOURCE, KEY));

      // check if file with that name exists
      if( !filename.DoesExist())
      {
        return util::ShPtrList< pdb::Line>();
      }

      // read in pdb to handler
      io::IFStream read;
      io::File::MustOpenIFStream( read, filename.GetFullName());
      pdb::Handler handler( read);
      io::File::CloseClearFStream( read);

      // locate remark lines
      return handler.GetLines( pdb::GetLineTypes().REMARK);
    }

    //! @brief store protein
    //! @param MOLECULE protein to store
    //! @param SOURCE corresponds to the prefix for the files e.g "run1"
    //! @return key associated with protein
    std::string ProteinStorageFile::Store( const ProteinModel &MOLECULE, const std::string &SOURCE)
    {
      // if overwriting files and this source has not been written to yet
      if( m_InitializerType == e_Overwrite && !m_Sources.Contains( SOURCE))
      {
        // remove all entries w/ same source
        RemoveEntries( SOURCE);
        m_Sources.Insert( SOURCE);
      }

      // initialize current keys
      size_t current_key( 0);

      std::string key;
      while( true)
      {
        key = KeyToString( current_key);
        const io::DirectoryEntry filename( Filename( SOURCE, key));

        // check if a file with the name exists
        const bool file_exists( filename.DoesExist());

        // in the mean time, somebody might have written to the directory
        if( file_exists)
        {
          // advance to next key
          ++current_key;
          continue;
        }

        io::OFStream write;
        io::File::MustOpenOFStream( write, filename.GetFullName());
        m_Factory->WriteModelToPDB( MOLECULE, write);
        io::File::CloseClearFStream( write);
        break;
      }

      // increment key
      return key;
    }

    //! @brief store protein
    //! @param MOLECULE protein to store
    //! @param SOURCE corresponds to the prefix for the files e.g "run1"
    //! @param KEY key under which protein was stored
    //! @return true if store was successful
    bool ProteinStorageFile::Store( const ProteinModel &MOLECULE, const std::string &SOURCE, const std::string &KEY)
    {
      // check that key is valid
      if( !IsValidKey( KEY))
      {
        BCL_MessageStd( "KEY " + util::Format()( KEY) + " is not valid");
        return false;
      }

      // filename and check if exists
      io::DirectoryEntry filename( Filename( SOURCE, KEY));

      // if the file exists and should not be overwritten
      if( filename.DoesExist() && m_InitializerType != e_Overwrite)
      {
        std::string new_key( KEY);
        size_t i( s_SmallKeySize);
        do
        {
          i = s_SmallKeySize;
          while( --i)
          {
            if( new_key[ i] == '9')
            {
              new_key[ i] = '0';
            }
            else
            {
              new_key[ i] += 1;
              break;
            }
          }
          filename = Filename( SOURCE, new_key);
        } while( filename.DoesExist() && i);
        if( filename.DoesExist())
        {
          // only return false if the storage is full
          return false;
        }
      }

      // write
      io::OFStream write;
      io::File::MustOpenOFStream( write, filename.GetFullName());
      BCL_MessageDbg( "writing model to pdb file " + filename.GetFullName());
      m_Factory->WriteModelToPDB( MOLECULE, write);
      io::File::CloseClearFStream( write);

      // end
      return true;
    }

    //! @brief store ensemble of proteins
    //! @param ENSEMBLE shptr list of proteins
    //! @param SOURCE corresponds to the prefix for the files e.g "run1"
    //! @return vector of keys
    storage::Vector< std::string> ProteinStorageFile::Store
    (
      const util::ShPtrList< ProteinModel> &ENSEMBLE,
      const std::string &SOURCE
    )
    {
      // keys
      storage::Vector< std::string> keys;
      keys.AllocateMemory( ENSEMBLE.GetSize());

      // iterate over ensemble
      for
      (
        util::ShPtrList< ProteinModel>::const_iterator itr( ENSEMBLE.Begin()), itr_end( ENSEMBLE.End());
        itr != itr_end;
        ++itr
      )
      {
        // store molecule and insert key into keys
        keys.PushBack( Store( **itr, SOURCE));
      }

      return keys;
    }

    //! @brief deletes all entries in this storage
    void ProteinStorageFile::Delete()
    {
      // delete the storage
      if( m_Directory.IsDefined())
      {
        m_Directory->Remove( true);
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProteinStorageFile::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Directory, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProteinStorageFile::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Directory, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief construct complete filename from source and key
    //! @brief SOURCE the source for the model
    //! @brief KEY the key for that protein
    //! @brief filename of form {initializer}/{SOURCE}{KEY}.pdb
    std::string ProteinStorageFile::Filename( const std::string &SOURCE, const std::string &KEY) const
    {
      return m_Directory->AppendFilename
             (
               ( *m_CompressionType)->AddExtension
               (
                 SOURCE + KEY + io::File::GetExtensionDelimiter()
                 + pdb::GetDefaultFileExtension()
               )
             );
    }

    //! @brief converts a key into a string
    //! @param KEY key to be converted
    //! @return converted string
    std::string ProteinStorageFile::KeyToString( const size_t &KEY)
    {
      static const util::Format s_KeyToString
      (
        util::Format().W( 4).R().Fill( '0')
      );

      return "_" + s_KeyToString( KEY);
    }

    //! @brief check if key is valid string
    //! @param KEY the key to be checked
    //! @return true if the key is valid
    bool ProteinStorageFile::IsValidKey( const std::string &KEY)
    {
      // check that key is not empty
      if( KEY.empty())
      {
        return false;
      }

      // if this is a small string size
      if( KEY.size() == s_SmallKeySize)
      {
        // create iterators on the string
        std::string::const_iterator itr( KEY.begin()), itr_end( KEY.end());

        // first char should be "_"
        if( *itr != '_')
        {
          return false;
        }

        // other chars should be digits
        ++itr;
        for( ; itr != itr_end; ++itr)
        {
          if( !isdigit( *itr))
          {
            return false;
          }
        }
      }
      // if this is a large string size
      else if( KEY.size() == s_LargeKeySize)
      {
        // iterate over the string
        for
        (
          std::string::const_iterator itr( KEY.begin()), itr_end( KEY.end());
          itr != itr_end; ++itr
        )
        {
          // if the char is not a number
          if( !isdigit( *itr))
          {
            // check for '_' at the right positions
            if( *itr == '_' && ( itr == KEY.begin() || itr == KEY.end() - 6))
            {
              continue;
            }

            BCL_MessageStd
            (
              "non numeric \"" + util::Format()( *itr) + "\" in key \"" + KEY +
              "\" should be a \"_\" and either at the beginning or at end - 6"
            );

            return false;
          }
        }
      }
      // wrong size
      else
      {
        BCL_MessageStd
        (
          "key size does not match " + util::Format()( 11) + " or "
          + util::Format()( 5) + " as it has size " + util::Format()( KEY.size())
          + " candidate key was: " + KEY
        );
        return false;
      }

      // passed all checks
      return true;
    }

    //! @brief gets all entries that have the given source
    //! @param SOURCE source of entries
    //! @return list of found entries
    storage::List< io::DirectoryEntry> ProteinStorageFile::GetAllEntries( const std::string &SOURCE) const
    {
      // directory content
      return m_Directory->ListEntries
             (
               io::Directory::e_File,
               SOURCE,
               ( *m_CompressionType)->AddExtension( "." + pdb::GetDefaultFileExtension())
             );
    }

    //! @brief removes all entries that have the given source
    //! @param SOURCE source of entries to be removed
    void ProteinStorageFile::RemoveEntries( const std::string &SOURCE) const
    {
      // get the matching entries
      storage::List< io::DirectoryEntry> entries( GetAllEntries( SOURCE));

      // iterate over the list
      for
      (
        storage::List< io::DirectoryEntry>::iterator itr( entries.Begin()), itr_end( entries.End());
        itr != itr_end; ++itr
      )
      {
        // remove the entry
        if( !itr->Remove())
        {
          BCL_MessageStd( "Unable to remove file, " + itr->GetName());
        }
      }
    }

  } // namespace assemble
} // namespace bcl
