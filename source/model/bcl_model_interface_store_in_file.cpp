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
#include "model/bcl_model_interface_store_in_file.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "model/bcl_model_interface_retrieve_from_file.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    //! format to convert Key to string
    const util::Format InterfaceStoreInFile::s_KeyToString
    (
      util::Format().W( 6).R().Fill( '0')
    );

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> InterfaceStoreInFile::s_Instance
    (
      util::Enumerated< StoreInterface>::AddInstance( new InterfaceStoreInFile())
    );

    //! @brief ExtensionType as string
    //! @param INIT_TYPE the ExtensionType
    //! @return the string for the INIT_TYPE
    const std::string &InterfaceStoreInFile::GetExtensionTypeDescriptor( const ExtensionType &TYPE)
    {
      static const std::string s_extensions[] =
      {
        "model",
        "descriptor",
        "result",
        "info",
        "independent",
        GetStaticClassName< ExtensionType>()
      };
      return s_extensions[ size_t( TYPE)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    InterfaceStoreInFile::InterfaceStoreInFile() :
      m_DirectoryName(),
      m_WriteDescriptors( true),
      m_WriteExperimentalPredicted( false),
      m_FilePrefix( ""),
      m_ModelKey( util::GetUndefinedSize_t()),
      m_HaveOpenLockFile( false),
      m_Compression( io::GetStreamBufferClasses().e_Uncompressed)
    {
      if( io::GetStreamBufferClasses().HaveEnumWithName( "BZ2"))
      {
        m_Compression = io::StreamBufferClass( "BZ2");
      }
      else if( io::GetStreamBufferClasses().HaveEnumWithName( "GZ"))
      {
        m_Compression = io::StreamBufferClass( "GZ");
      }
    }

    //! @brief default constructor
    InterfaceStoreInFile::InterfaceStoreInFile( const std::string &DIRECTORY_NAME) :
      m_DirectoryName(),
      m_WriteDescriptors( true),
      m_WriteExperimentalPredicted( false),
      m_FilePrefix( ""),
      m_ModelKey( util::GetUndefinedSize_t()),
      m_HaveOpenLockFile( false),
      m_Compression( io::GetStreamBufferClasses().e_Uncompressed)
    {
      m_Directory = io::Directory( DIRECTORY_NAME);

      // check if directory exists or could be created
      if( !m_Directory.DoesExist())
      {
        m_Directory.Make();
      }

      if( io::GetStreamBufferClasses().HaveEnumWithName( "BZ2"))
      {
        m_Compression = io::StreamBufferClass( "BZ2");
      }
      else if( io::GetStreamBufferClasses().HaveEnumWithName( "GZ"))
      {
        m_Compression = io::StreamBufferClass( "GZ");
      }
    }

    //! @brief destructor
    InterfaceStoreInFile::~InterfaceStoreInFile()
    {
      CleanUp();
    }

    //! @brief Clone function
    //! @return pointer to new InterfaceStoreInFile
    InterfaceStoreInFile *InterfaceStoreInFile::Clone() const
    {
      return new InterfaceStoreInFile( *this);
    }

    //! @brief cleanup function, to remove lock file if it exists
    void InterfaceStoreInFile::CleanUp()
    {
      // if the lock file is still open at exit
      // (this can happen due to the user terminating the program precisely when the result file is written out)
      // then remove the lock file
      if( m_HaveOpenLockFile)
      {
        io::DirectoryEntry lock_file( Filename( "", e_Result) + ".lock");
        if( lock_file.DoesExist())
        {
          lock_file.Remove();
        }
      }
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &InterfaceStoreInFile::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &InterfaceStoreInFile::GetAlias() const
    {
      static const std::string s_Name( "File");
      return s_Name;
    }

    //! @brief number of objects in source
    //! @return number of objects with that source prefix
    size_t InterfaceStoreInFile::GetSize() const
    {
      return GetAllKeys().GetSize();
    }

    //! @brief get all keys for given range
    //! @param RANGE range of key values
    //! @return all keys of given prefix
    storage::Vector< std::string> InterfaceStoreInFile::GetKeys( const math::Range< size_t> &RANGE) const
    {
      // keys
      storage::Vector< std::string> keys( GetAllKeys());
      storage::Vector< std::string> valid_keys;

      for
      (
        storage::Vector< std::string>::const_iterator itr( keys.Begin()), itr_end( keys.End());
        itr != itr_end;
        ++itr
      )
      {
        const size_t key_id( util::ConvertStringToNumericalValue< size_t>( *itr));

        // check for valid id
        if( RANGE.IsWithin( key_id))
        {
          valid_keys.PushBack( *itr);
        }
      }

      // return valid keys
      return valid_keys;
    }

    //! @brief get all keys for given source
    //! @return all keys of given prefix
    storage::Vector< std::string> InterfaceStoreInFile::GetAllKeys() const
    {
      // keys
      storage::Vector< std::string> keys;

      // directory content
      storage::List< io::DirectoryEntry> dir_entries;
      for
      (
        io::StreamBufferClasses::const_iterator
          itr_compression( io::GetStreamBufferClasses().Begin()),
          itr_compression_end( io::GetStreamBufferClasses().End());
        itr_compression != itr_compression_end;
        ++itr_compression
      )
      {
        dir_entries.Append
        (
          m_Directory.ListEntries
          (
            io::Directory::e_File,
            m_FilePrefix,
            ( **itr_compression)->AddExtension( "." + GetExtensionTypeDescriptor( e_Model))
          )
        );
      }

      // iterate over all entries
      for( storage::List< io::DirectoryEntry>::const_iterator itr( dir_entries.Begin()), itr_end( dir_entries.End()); itr != itr_end; ++itr)
      {
        // filename without extension
        std::string key( io::File::RemoveLastExtension( io::File::RemoveCompressionExtension( itr->GetName())));

        // extract the key
        key.erase( 0, std::string( m_FilePrefix).length());

        // key needs to be valid
        if( IsValidKey( key))
        {
          keys.PushBack( key);
        }
      }

      return keys;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief try to store the result descriptor
    //! should assert if the labels are already stored and differ from those given
    //! @param RESULT_DESCRIPTOR result object data label for all models in this storage
    void InterfaceStoreInFile::StoreResultDescriptor( const util::ObjectDataLabel &RESULT_DESCRIPTOR)
    {
      // filename and check if exists
      const std::string filename( Filename( "", e_Result));
      io::DirectoryEntry file( filename);
      const std::string lock_filename( filename + ".lock");
      io::DirectoryEntry lock_file( lock_filename);

      io::OFStream output;

      // check if the lock file already exists, if so, wait for it to vanish
      if( file.DoesExist() || lock_file.DoesExist())
      {
        // ensure that the already-stored descriptors match with those given
        const util::ObjectDataLabel retrieved_result
        (
          InterfaceRetrieveFromFile::RetrieveResultDescriptor( filename)
        );
        BCL_Assert
        (
          retrieved_result == RESULT_DESCRIPTOR,
          "Tried to store result descriptor\n" + RESULT_DESCRIPTOR.ToString()
          + " but models in this storage already have a different result descriptor:\n"
          + retrieved_result.ToString()
        );
        return;
      }

      // create the lock file
      BCL_Assert( io::File::TryOpenOFStream( output, lock_filename), "Could not create lock file!");
      m_HaveOpenLockFile = true;
      io::File::CloseClearFStream( output);

      // create the descriptor file
      io::File::MustOpenOFStream( output, filename);
      output << RESULT_DESCRIPTOR.ToString();
      io::File::CloseClearFStream( output);

      // remove the lock file
      lock_file.Remove();
      m_HaveOpenLockFile = false;
    }

    //! @brief store model with additional information in filesystem
    //! @param MODEL model::Interface that will be stored
    //! @param RESULT result value that was evaluated by objective function with current model
    //! @param DESCRIPTORS descriptors used to create dataset the current model was trained on
    //! @param METHOD_NAME name of iterate label that defines the used machine learning algorithm
    //! @param OBJ_FUNCTION objective function used in calculating the result value
    //! @param CV_INFO cross validation info
    //! @return key of stored model
    std::string InterfaceStoreInFile::Store
    (
      const util::ShPtr< Interface> &MODEL,
      const float RESULT,
      const util::ObjectDataLabel &DESCRIPTORS,
      const util::ObjectDataLabel &METHOD_NAME,
      const util::ObjectDataLabel &OBJ_FUNCTION,
      const CrossValidationInfo &CV_INFO
    )
    {
      // initialize current keys
      size_t current_key( 0);

      // get all keys
      storage::Vector< std::string> keys( GetAllKeys());

      // if there are keys existent for that source
      if( util::IsDefined( m_ModelKey))
      {
        current_key = m_ModelKey;
      }
      else if( !keys.IsEmpty())
      {
        // kind largest key
        storage::Vector< std::string>::const_iterator itr( std::max_element( keys.Begin(), keys.End()));

        // convert string to number
        current_key = util::ConvertStringToNumericalValue< size_t>( *itr);
      }

      std::string key;
      while( true)
      {
        key = s_KeyToString( current_key);
        io::DirectoryEntry filename( Filename( key, e_Model));

        // check if a file with the name exists
        const bool file_exists( filename.DoesExist());

        // in the mean time, somebody might have written to the directory
        if( file_exists)
        {
          // advance to next key
          ++current_key;
          continue;
        }

        Store( MODEL, RESULT, DESCRIPTORS, METHOD_NAME, OBJ_FUNCTION, CV_INFO, util::Format()( current_key));

        break;
      }

      // increment key
      return key;
    }

    //! @brief store model with additional information in filesystem
    //! @param MODEL model::Interface that will be stored
    //! @param RESULT result value that was evaluated by objective function with current model
    //! @param DESCRIPTORS descriptors used to create dataset the current model was trained on
    //! @param KEY prefered key for model to store
    //! @param METHOD_NAME name of iterate label that defines the used machine learning algorithm
    //! @param OBJ_FUNCTION objective function used in calculating the result value
    //! @param CV_INFO cross validation info
    //! @return key of stored model
    std::string InterfaceStoreInFile::Store
    (
      const util::ShPtr< Interface> &MODEL,
      const float RESULT,
      const util::ObjectDataLabel &DESCRIPTORS,
      const util::ObjectDataLabel &METHOD_NAME,
      const util::ObjectDataLabel &OBJ_FUNCTION,
      const CrossValidationInfo &CV_INFO,
      const std::string &KEY
    )
    {
      // convert key into right format
      const std::string key( s_KeyToString( KEY));

      // check that key is valid
      if( !IsValidKey( key))
      {
        BCL_MessageStd( "Model file with invalid key: " + key + " . Omitting storage!");
        return std::string();
      }

      WriteToFile( MODEL, key, e_Model);
      WriteToFile( CV_INFO, key, e_Info);

      if( m_WriteDescriptors)
      {
        WriteToFile( DESCRIPTORS, key, e_Descriptor);
      }

      // end
      return key;
    }

    //! @brief store model with additional information in filesystem
    //! @param MODEL model::Interface that will be stored
    //! @param DESCRIPTORS descriptors used to create dataset the current model was trained on
    //! @return key of stored model
    std::string InterfaceStoreInFile::Store
    (
      const util::ShPtr< Interface> &MODEL,
      const util::ObjectDataLabel &DESCRIPTORS
    )
    {
      // convert key into right format
      storage::Vector< std::string> keys( GetAllKeys());
      const std::string key
      (
        s_KeyToString
        (
          keys.GetSize() > 0
          ? util::ConvertStringToNumericalValue< size_t>( keys.LastElement()) + 1
          : size_t( 0)
        )
      );

      // check that key is valid
      if( !IsValidKey( key))
      {
        BCL_MessageStd( "Model file with invalid key: " + key + " . Omitting storage!");
        return std::string();
      }

      WriteToFile( MODEL, key, e_Model);

      if( m_WriteDescriptors)
      {
        WriteToFile( DESCRIPTORS, key, e_Descriptor);
      }

      // end
      return key;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief construct complete filename from source and key
    //! @brief KEY the key for that protein
    //! @brief filename of form {initializer}/m_FilePrefix{KEY}.s_ModelFileExtension
    std::string InterfaceStoreInFile::Filename( const std::string &KEY, const ExtensionType &EXTENSION_TYPE) const
    {
      return m_Directory.AppendFilename( ( *m_Compression)->AddExtension( m_FilePrefix + KEY + io::File::GetExtensionDelimiter() + GetExtensionTypeDescriptor( EXTENSION_TYPE)));
    }

    //! @brief write util::ObjectInterface from source and key
    //! @param OBJECT current object that should be written in file
    //! @param KEY current key
    //! @param EXTENSION_TYPE extention type of file the object will be written into
    void InterfaceStoreInFile::WriteToFile( const util::ObjectInterface &OBJECT, const std::string &KEY, const ExtensionType &EXTENSION_TYPE) const
    {
      // Filename(  and check if exists
      io::DirectoryEntry filename( Filename( KEY, EXTENSION_TYPE));

      if( filename.DoesExist())
      {
        BCL_MessageStd( "file with key: " + KEY + " already exists. Omitting storage!");
        return;
      }

      // initialize stream
      io::OFStream write;
      io::File::MustOpenOFStream( write, filename.GetFullName());
      // write result
      io::Serialize::Write( OBJECT, write);
      io::File::CloseClearFStream( write);
    }

    //! @brief write util::ObjectDataLabel from source and key
    //! @param OBJECT current object that should be written in file
    //! @param KEY current key
    //! @param EXTENSION_TYPE extention type of file the object will be written into
    void InterfaceStoreInFile::WriteToFile( const util::ObjectDataLabel &LABEL, const std::string &KEY, const ExtensionType &EXTENSION_TYPE) const
    {
      // Filename(  and check if exists
      io::DirectoryEntry filename( Filename( KEY, EXTENSION_TYPE));

      if( filename.DoesExist())
      {
        BCL_MessageStd( "file with key: " + KEY + " already exists. Omitting storage!");
        return;
      }

      // initialize stream
      io::OFStream write;
      io::File::MustOpenOFStream( write, filename.GetFullName());
      write << LABEL.ToString();
      io::File::CloseClearFStream( write);
    }

    //! @brief check if key is valid string
    //! @param KEY the key to be checked
    //! @return true if the key is valid
    bool InterfaceStoreInFile::IsValidKey( const std::string &KEY)
    {
      // check that the key is of size 6 and that it is an unsigned integer
      return KEY.length() == 6 && util::LengthOfUnsignedIntegerType( KEY) == 6;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer InterfaceStoreInFile::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Storage of trained model on the file system"
      );
      parameters.AddInitializer
      (
        "directory",
        "directory to stored models (for particular session)",
        io::Serialization::GetAgent( &m_DirectoryName)
      );
      parameters.AddInitializer
      (
        "prefix",
        "file prefix to model name eg. MYPREFIX000000." + GetExtensionTypeDescriptor( e_Model),
        io::Serialization::GetAgent( &m_FilePrefix),
        ""
      );
      parameters.AddInitializer
      (
        "write_descriptors",
        "write out descriptors to file",
        io::Serialization::GetAgent( &m_WriteDescriptors),
        "1"
      );
      parameters.AddOptionalInitializer
      (
        "key",
        "model key if not given it will be automatically set! in case of key duplicates the model will be overwritten!",
        io::Serialization::GetAgent( &m_ModelKey)
      );

      return parameters;
    }

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    bool InterfaceStoreInFile::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      m_Directory = io::Directory( m_DirectoryName);

      if( io::GetStreamBufferClasses().HaveEnumWithName( "BZ2"))
      {
        m_Compression = io::StreamBufferClass( "BZ2");
      }
      else if( io::GetStreamBufferClasses().HaveEnumWithName( "GZ"))
      {
        m_Compression = io::StreamBufferClass( "GZ");
      }
      else
      {
        m_Compression = io::GetStreamBufferClasses().e_Uncompressed;
      }
      // check if directory exists or could be created
      if( !m_Directory.DoesExist())
      {
        return m_Directory.Make();
      }
      return true;
    }

  } // namespace model
} // namespace bcl
