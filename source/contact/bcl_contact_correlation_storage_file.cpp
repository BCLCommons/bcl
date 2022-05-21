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
#include "align/bcl_align_alignment_node.h"
#include "align/bcl_align_handler_pir.h"
#include "contact/bcl_contact_correlation_matrix.h"
#include "contact/bcl_contact_correlation_storage_file.h"
#include "io/bcl_io_directory_entry.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {
  //////////
  // data //
  //////////

    //! format to convert Key to string                                                       // CHECK FORMAT
    const util::Format CorrelationStorageFile::s_KeyToString
    (
      util::Format().W( 6).R().Fill( '0')
    );

    //! default file extension for objects
    const std::string CorrelationStorageFile::s_FileExtension( "bcl");

    //! default prefix for filename of stored model information
    const std::string CorrelationStorageFile::s_FilePrefix( "ampair");

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> CorrelationStorageFile::s_Instance
    (
      GetObjectInstances().AddInstance( new CorrelationStorageFile())
    );

    //! Default delimiter to be used for reading/writing alignments to file
    const char CorrelationStorageFile::s_Delim( '$');

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  /////////////////
  // data access //
  /////////////////
    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CorrelationStorageFile::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief initialize the correlation storage
    //! @param INITIALIZER encodes where data is stored
    //! @param INIT_FLAG flag for type of initialization
    //! @return true if initialize was successful
    bool CorrelationStorageFile::Initialize( const std::string &INITIALIZER, const InitializerType INIT_FLAG)
    {
      BCL_MessageStd( "Initializer: " + INITIALIZER);
      util::ShPtr< io::Directory> sp_dir( new io::Directory( INITIALIZER));

      // check if directory exists
      const bool dir_exists( sp_dir->DoesExist());

      switch( INIT_FLAG)
      {
        // attach to existing storage - directory needs to exist
        case e_Attach:
        {
          if( dir_exists)
          {
            m_Directory = sp_dir;
          }
          return dir_exists;
        }
        // create storage, directory does need to not exist
        case e_Create:
        {
          // if it exists, it will not be overwritten
          if( dir_exists)
          {
            BCL_MessageCrt( "cannot create " + INITIALIZER + " since the directory already exists!");
            return false;
          }
          const bool make_success( sp_dir->Make());
          if( make_success)
          {
            m_Directory = sp_dir;
          }
          return make_success;
        }
        // overwrite directory, when it exist, will be overwritten
        case e_Overwrite:
        {
          if( dir_exists)
          {
            // clear out any existing model files; ignore other files
            m_Directory = sp_dir;

            // get the keys of all current models
            storage::Vector< std::string> keys( GetAllKeys());

            bool clear_success( true);

            // iterate over the keys and delete them
            for
            (
              storage::Vector< std::string>::const_iterator itr( keys.Begin()), itr_end( keys.End());
              itr != itr_end && clear_success;
              ++itr
            )
            {
              io::DirectoryEntry current_entry( sp_dir, s_FilePrefix + *itr + io::File::GetExtensionDelimiter() + s_FileExtension);
              clear_success &= current_entry.Remove();
            }
            if( !clear_success)
            {
              BCL_MessageCrt( "cannot clear " + INITIALIZER + " to overwrite!");
              return clear_success;
            }
            m_Directory = sp_dir;

            // initalize as create
            return clear_success;
          }
          else // directory does not exist
          {
            const bool make_success( sp_dir->Make());
            if( make_success)
            {
              m_Directory = sp_dir;
            }
            return make_success;
          }
        }
        default:
        {
          BCL_Exit( "unknown InitializerType supplied", -1);
          break;
        }
      }

      // end
      return true;
    }

    //! @brief number of data items in source
    //! @param SOURCE source of data file
    //! @return number of MSA correlation matrix pairs in storage List at SOURCE
    size_t CorrelationStorageFile::GetSize() const
    {
      return GetAllKeys().GetSize();
    }

    //! @brief get all keys for given source
    //! @return all keys of given source
    storage::Vector< std::string> CorrelationStorageFile::GetAllKeys() const
    {
      // keys
      storage::Vector< std::string> keys;

      // Check that directory exists
      if( m_Directory.IsDefined())
      {
        // directory content
        const storage::List< io::DirectoryEntry> dir_entries
        (
          m_Directory->ListEntries( io::Directory::e_File, s_FilePrefix, "." + s_FileExtension)
        );

        // iterate over all entries
        storage::List< io::DirectoryEntry>::const_iterator itr( dir_entries.Begin()), itr_end( dir_entries.End());
        for( ; itr != itr_end; ++itr)
        {
          // filename without extension
          std::string key( io::File::RemoveLastExtension( itr->GetName()));

          // extract the key and add it to keys vector
          key.erase( 0, s_FilePrefix.length());
          keys.PushBack( key);
        }
      }

      return keys;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get AlignmentMatrixPair
    //! @param SOURCE source of molecule eg. given AlignmentMatrixPairs
    //! @param KEY key identifier for specific molecule in given source
    //! @return shptr to AlignmentMatrixPair of interest and shptr to NULL if not found
    util::ShPtr< CorrelationStorageFile::AlignmentMatrixPair> CorrelationStorageFile::Retrieve( const std::string &KEY) const
    {
      // Create filename using key
      std::string filename( StorageFilename( KEY));

      // Check for file
      if( !io::DirectoryEntry( filename).DoesExist())
      {
        // BCL Message for verbose
        BCL_MessageVrb( "File :" + filename + " does not exist");

        // Return ShPtr to NULL if failed
        return util::ShPtr< AlignmentMatrixPair>( 0);
      }

      // BCL Message for verbose
      BCL_MessageVrb( "File :" + filename + " found");
      // Open file, create object, and read in
      io::IFStream read;
      io::File::MustOpenIFStream( read, filename);
      std::string alignment;
      CorrelationMatrix correlation_matrix;
      BCL_MessageDbg( "Reading alignment portion " + filename);
      io::Serialize::Read( alignment, read);
      BCL_MessageDbg( "Reading alignment portion " + filename + " done");
      BCL_MessageDbg( "Reading correlation matrix " + filename);
      io::Serialize::Read( correlation_matrix, read);
      BCL_MessageDbg( "Reading correlation matrix " + filename + " done");

      // Create alignment handler to read in string portion into alignment object
      align::HandlerPIR< biol::AABase> pir_handler( s_Delim);
      std::stringstream stream( alignment);
      BCL_MessageDbg( "Reading alignment " + filename);
      BCL_MessageDbg( "Read alignment :" + alignment);
      util::ShPtr< align::AlignmentNode< biol::AABase> > alignment_ptr
      (
        pir_handler.ReadAlignment( stream, biol::AASequence())
      );
      BCL_MessageDbg( "Reading alignment " + filename + " done");

      // Translate input into AlignmentMatrixPair object
      util::ShPtr< CorrelationStorageFile::AlignmentMatrixPair> am_pair_ptr( new AlignmentMatrixPair( *alignment_ptr, correlation_matrix));

      // Close file
      io::File::CloseClearFStream( read);

      // Return ShPtr to AMPair
      return am_pair_ptr;
    }

    //! @brief get ensemble of AlignmentMatrixPair
    //! @return shptr list of AlignmentMatrixPair from given source directory
    util::ShPtrList< CorrelationStorageFile::AlignmentMatrixPair> CorrelationStorageFile::RetrieveEnsemble() const
    {
      // Get keys
      storage::Vector< std::string> key_vector( GetAllKeys());

      // Return retrieve for all keys
      return RetrieveEnsemble( key_vector);
    }

    //! @brief get ensemble of AlignmentMatrixPair for given keys
    //! @param KEYS vector of keys
    //! @return shptr list of molecules from given source
    util::ShPtrList< CorrelationStorageFile::AlignmentMatrixPair> CorrelationStorageFile::RetrieveEnsemble
    (
      const storage::Vector< std::string> &KEYS
    ) const
    {
      // Retrieve given keys and save them in the list
      util::ShPtrList< AlignmentMatrixPair> am_pair_ptr_list;
      storage::Vector< std::string>::const_iterator iter( KEYS.Begin());
      storage::Vector< std::string>::const_iterator iter_end( KEYS.End());
      for( ; iter != iter_end; ++iter)
      {
        am_pair_ptr_list.PushBack( Retrieve( *iter));
      }

      // Return list of results
      return am_pair_ptr_list;
    }

//    //! @brief get ensemble of AlignmentMatrixPair for given keys which does not map to keys, following plain numbering
//    //! @param RANGE range of AlignmentMatrixPair
//    //! @return shptr list of AlignmentMatrixPair from given source
//    util::ShPtrList< AlignmentMatrixPair> CorrelationStorageFile::RetrieveEnsemble
//    (
//      const math::Range< size_t> &RANGE
//    ) const
//    {
//        // IMPLEMENT
//    }

    //! @brief store AlignmentMatrixPairs
    //! @param ALIGNMENT_MATRIX_PAIR to store
    //! @return key string associated with target sequence - pdb sequence ID, otherwise empty string
    std::string CorrelationStorageFile::Store( const AlignmentMatrixPair &ALIGNMENT_MATRIX_PAIR)
    {
      // key of stored AlignmentMatrixPair
      const std::string key( ALIGNMENT_MATRIX_PAIR.First().GetSequenceIds().FirstElement());

      // Store with key and return key
      if( CorrelationStorageFile::Store( ALIGNMENT_MATRIX_PAIR, key))
      {
        return key;
      }
      else
      {
        // else return empty string
        return std::string();
      }
    }

    //! @brief store AlignmentMatrixPairs as pairs of alignment string in wrapper and correlation matrix
    //! @param ALIGNMENT_MATRIX_PAIR to store
    //! @param KEY key under which AlignmentMatrixPair was stored
    //! @return true if store was successful false if unsuccessful or file already exists
    bool CorrelationStorageFile::Store( const AlignmentMatrixPair &ALIGNMENT_MATRIX_PAIR, const std::string &KEY)
    {
      // OFstream for writing into bcl file
      io::OFStream write;

      // filename for writing into
      std::string filename( StorageFilename( KEY));

      // Write out filename if verbose
      BCL_MessageVrb( "File path for storage: " + filename);

      // Check for file already of that pdb_id and prefix
      if( io::DirectoryEntry( filename).DoesExist())
      {
        return false;
      }

      // open file for reading
      io::File::MustOpenOFStream( write, filename);

      // translate Alignment given by node to string stream for storage
      align::HandlerPIR< biol::AABase> pir_handler( s_Delim);
      std::stringstream alignment_string;
      pir_handler.WriteAlignment( alignment_string, ALIGNMENT_MATRIX_PAIR.First());

      // write AlignmentMatrixPair
      io::Serialize::Write( alignment_string.str(), write) << '\n';
      io::Serialize::Write( ALIGNMENT_MATRIX_PAIR.Second(), write);
      // close file stream
      io::File::CloseClearFStream( write);

      return true;
    }

    //! @brief store ensemble of AlignmentMatrixPair
    //! @param ENSEMBLE shptr list of AlignmentMatrixPair
    //! @return vector of keys
    storage::Vector< std::string> CorrelationStorageFile::Store( const util::ShPtrList< AlignmentMatrixPair> &ENSEMBLE)
    {
      // Create vector of return keys
      storage::Vector< std::string> key_vector;

      // Iterate through ShPtrList of AlignmentMatrixPairs and store each
      util::ShPtrList< AlignmentMatrixPair>::const_iterator itr( ENSEMBLE.Begin());
      util::ShPtrList< AlignmentMatrixPair>::const_iterator itr_end( ENSEMBLE.End());
      for( ; itr != itr_end; ++itr)
      {
        // Store AlignmentMatrixPair and store returning key in vector to be returned
        key_vector.PushBack( Store( **itr));
      }

      // return vector of keys some of which may be empty strings
      return key_vector;
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &CorrelationStorageFile::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Directory, ISTREAM);
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &CorrelationStorageFile::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Directory, OSTREAM, INDENT);
      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace contact
} // namespace bcl
