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
#include "model/bcl_model_meta_data_storage_file.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_dataset.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
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
    const util::Format MetaDataStorageFile::s_KeyToString
    (
      util::Format().W( 6).R().Fill( '0')
    );

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MetaDataStorageFile::s_Instance
    (
      util::Enumerated< MetaDataStorageInterface>::AddInstance( new MetaDataStorageFile())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new MetaDataStorageFile
    MetaDataStorageFile *MetaDataStorageFile::Clone() const
    {
      return new MetaDataStorageFile( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MetaDataStorageFile::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MetaDataStorageFile::GetAlias() const
    {
      static const std::string s_Name( "File");
      return s_Name;
    }

    //! @brief initialize the smallmolecule storage
    //! @return true if initialize was successful
    bool MetaDataStorageFile::Initialize()
    {
      util::ShPtr< io::Directory> sp_dir( new io::Directory( m_DirectoryName));

      // check if directory exists or could be created
      if( sp_dir->DoesExist() || sp_dir->Make())
      {
        m_Directory = sp_dir;
        return true;
      }

      BCL_MessageCrt( "Could not create directory: " + m_DirectoryName);
      // end
      return false;
    }

    //! @brief number of objects in source
    //! @return number of objects with that source prefix
    size_t MetaDataStorageFile::GetSize() const
    {
      return GetAllKeys().GetSize();
    }

    //! @brief get all keys for given source
    //! @return all keys of given prefix
    storage::Vector< std::string> MetaDataStorageFile::GetAllKeys() const
    {
      // directory content
      const storage::List< io::DirectoryEntry> dir_entries
      (
        m_Directory->ListEntries
        (
          io::Directory::e_File,
          m_FilePrefix,
          "." + GetExtensionTypeDescriptor( e_MetaData)
        )
      );

      // keys
      storage::Vector< std::string> keys;

      // iterate over all entries
      for( storage::List< io::DirectoryEntry>::const_iterator itr( dir_entries.Begin()), itr_end( dir_entries.End()); itr != itr_end; ++itr)
      {
        // filename without extension
        std::string key( io::File::RemoveLastExtension( itr->GetName()));

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

//    //! @brief get smallmolecule
//    //! @param SOURCE prefix of smallmolecule eg. "prefix1"
//    //! @param KEY key identifier for specific smallmolecule in given source like "00001"
//    //! @return shptr to smallmolecule of interest, undefined if there is no such smallmolecule
//    util::ShPtr< Interface> MetaDataStorageFile::Retrieve( const std::string &KEY) const
//    {
//      // model
//      util::ShPtr< Interface> sp_interface;
//
//      const std::string key( s_KeyToString( KEY));
//
//      // check the key
//      if( !IsValidKey( key))
//      {
//        return sp_interface;
//      }
//
//      // filename and check if exists
//      const std::string filename( Filename( key, e_MetaData));
//
//      // check if file with that name exists
//      if( !io::DirectoryEntry( filename).DoesExist())
//      {
//        return sp_interface;
//      }
//
//      // result stored with model
//      float result;
//
//      // create factory and retrieve protein model
//      io::IFStream input;
//      io::File::MustOpenIFStream( input, filename);
//      io::Serialize::Read( result, input);
//      util::ObjectDataLabel descriptors( input);
//      io::Serialize::Read( sp_interface, input);
//
//      io::File::CloseClearFStream( input);
//
//      // if smaller result is better
//      if( m_BestResultSmallest)
//      {
//        // get key to best result so far
//        m_KeyToBestModelByResult = result < m_ResultToBestModel ? key : m_KeyToBestModelByResult;
//      }
//      else
//      {
//        // get key to best result so far
//        m_KeyToBestModelByResult = result > m_ResultToBestModel ? key : m_KeyToBestModelByResult;
//      }
//
//      // return retrieved smallmolecule
//      return sp_interface;
//    }

    //! @brief get descriptor set of model by key
    //! @param KEY key identifier for specific model/descriptor set in given source
    //! @return shptr to descriptorset of interest
    util::ObjectDataLabel MetaDataStorageFile::RetrieveDescriptorSet( const std::string &KEY) const
    {
      // check the key
      if( !IsValidKey( KEY))
      {
        return util::ObjectDataLabel();
      }

      // filename and check if exists
      const std::string filename( Filename( KEY, e_Descriptor));

      // check if file with that name exists
      if( !io::DirectoryEntry( filename).DoesExist())
      {
        return util::ObjectDataLabel();
      }

      // result stored with model
      float result;

      // create factory and retrieve protein model
      io::IFStream input;
      io::File::MustOpenIFStream( input, filename);
      io::Serialize::Read( result, input);
      util::ObjectDataLabel descriptors( input);
      io::File::CloseClearFStream( input);

      // return retrieved descriptor set
      return descriptors;
    }

    //! @brief retrieve best descriptor set by result
    //! @return shptr to best descriptor set of model with best score/result
    util::ObjectDataLabel MetaDataStorageFile::RetrieveBestDescriptorSetByResult() const
    {
      // if there is no key to best model determined yet
      if( m_KeyToBestModelByResult.empty())
      {
        // acquire all keys for that source and return ensemble
        storage::List< util::ObjectDataLabel> models = RetrieveAllDescriptors();
      }

      // best model
      return RetrieveDescriptorSet( m_KeyToBestModelByResult);
    }

    //! @brief get ensemble of descriptors associated each associated with a model
    //! @return shptr to list of descriptors of interest
    storage::List< util::ObjectDataLabel> MetaDataStorageFile::RetrieveAllDescriptors() const
    {
      // get all keys
      storage::Vector< std::string> keys( GetAllKeys());

      // list for retrieved descriptors
      storage::List< util::ObjectDataLabel> descriptors;

      for
      (
        storage::Vector< std::string>::const_iterator itr_keys( keys.Begin()), itr_keys_end( keys.End());
        itr_keys != itr_keys_end;
        ++itr_keys
      )
      {
        // insert descriptor set according to key
        descriptors.PushBack( RetrieveDescriptorSet( *itr_keys));
      }

      // return list with descriptors
      return descriptors;
    }

    //! @brief store model with additional information in filesystem
    //! @param RESULT result value that was evaluated by objective function with current model
    //! @param DESCRIPTORS descriptors used to create dataset the current model was trained on
    //! @param METHOD_NAME name of iterate label that defines the used machine learning algorithm
    //! @param OBJECTIVE_FUNCTION name of objective function label that defines the objective function
    //! @return key of stored model
    std::string MetaDataStorageFile::Store
    (
      const float RESULT,
      const util::ObjectDataLabel &DESCRIPTORS,
      const util::ObjectDataLabel &METHOD_NAME,
      const util::ObjectDataLabel &OBJECTIVE_FUNCTION
    )
    {
      // initialize current keys
      size_t current_key( 0);

      // get all keys
      storage::Vector< std::string> keys( GetAllKeys());

      // if there are keys existent for that source
      if( !keys.IsEmpty())
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
        io::DirectoryEntry filename( Filename( key, e_MetaData));

        // check if a file with the name exists
        const bool file_exists( filename.DoesExist());

        // in the mean time, somebody might have written to the directory
        if( file_exists)
        {
          // advance to next key
          ++current_key;
          continue;
        }

        // initialize stream
        io::OFStream write;

        // write out meta data
        io::File::MustOpenOFStream( write, filename.GetFullName());
        io::Serialize::Write( RESULT, write) << std::endl;
        io::Serialize::Write( m_Round, write) << std::endl;
        io::Serialize::Write( m_Iteration, write) << std::endl;
        io::Serialize::Write( m_CrossValidationIdIndependent, write) << std::endl;
        io::Serialize::Write( m_CrossValidationIdMonitoring, write) << std::endl;
        io::Serialize::Write( m_CrossValidationTotalChunks, write) << std::endl;
        io::File::CloseClearFStream( write);

        // write out descriptors
        filename = io::DirectoryEntry( Filename( key, e_Descriptor));
        io::File::MustOpenOFStream( write, filename.GetFullName());
        io::Serialize::Write( DESCRIPTORS, write);
        io::File::CloseClearFStream( write);

        // stop while loop since meta data was written
        break;
      }

      // invalidate best key
      m_KeyToBestModelByResult = "";

      // increment key
      return key;
    }

    //! @brief store model with additional information in filesystem
    //! @param RESULT result value that was evaluated by objective function with current model
    //! @param DESCRIPTORS descriptors used to create dataset the current model was trained on
    //! @param KEY prefered key for model to store
    //! @param EXPERIMENTAL_PREDICTED_VALUES serialized object of a list with pairs of
    //!        experimental and predicted values
    //! @param METHOD_NAME name of iterate label that defines the used machine learning algorithm
    //! @param OBJECTIVE_FUNCTION name of objective function label that defines the objective function
    //! @return key of stored model
    std::string MetaDataStorageFile::Store
    (
      const float RESULT,
      const util::ObjectDataLabel &DESCRIPTORS,
      const util::ObjectDataLabel &METHOD_NAME,
      const util::ObjectDataLabel &OBJECTIVE_FUNCTION,
      const std::string &KEY
    )
    {
      // check that key is valid
      if( !IsValidKey( KEY))
      {
        return std::string();
      }

      // Filename(  and check if exists
      io::DirectoryEntry filename( Filename( KEY, e_MetaData));

      if( filename.DoesExist())
      {
        return std::string();
      }

      // initialize stream
      io::OFStream write;

      // write out meta data
      io::File::MustOpenOFStream( write, filename.GetFullName());
      io::Serialize::Write( RESULT, write) << std::endl;
      io::Serialize::Write( m_Round, write) << std::endl;
      io::Serialize::Write( m_Iteration, write) << std::endl;
      io::Serialize::Write( m_CrossValidationIdIndependent, write) << std::endl;
      io::Serialize::Write( m_CrossValidationIdMonitoring, write) << std::endl;
      io::Serialize::Write( m_CrossValidationTotalChunks, write) << std::endl;
      io::File::CloseClearFStream( write);

      // write out descriptors
      filename = io::DirectoryEntry( Filename( KEY, e_Descriptor));
      io::File::MustOpenOFStream( write, filename.GetFullName());
      io::Serialize::Write( DESCRIPTORS, write);
      io::File::CloseClearFStream( write);

      // invalidate best key
      m_KeyToBestModelByResult = "";

      // end
      return KEY;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MetaDataStorageFile::Read( std::istream &ISTREAM)
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
    std::ostream &MetaDataStorageFile::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Directory, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief construct complete filename from key
    //! @brief KEY the key for that protein
    //! @brief filename of form {path}/FilePrefix{KEY}.FileExtension
    std::string MetaDataStorageFile::Filename( const std::string &KEY, const ExtensionType &TYPE) const
    {
      return m_Directory->AppendFilename( m_FilePrefix + KEY + io::File::GetExtensionDelimiter() + GetExtensionTypeDescriptor( TYPE));
    }

    //! @brief check if key is valid string
    //! @param KEY the key to be checked
    //! @return true if the key is valid
    bool MetaDataStorageFile::IsValidKey( const std::string &KEY)
    {
      // check that the key is of size 6 and that it is an unsigned integer
      return KEY.length() == 6 && util::LengthOfUnsignedIntegerType( KEY) == 6;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MetaDataStorageFile::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Storage of meta data on file system"
      );
      parameters.AddInitializer
      (
        "directory",
        "directory containing meta data information about cross-validation or descriptor selection",
        io::Serialization::GetAgent( &m_DirectoryName),
        "dir_storage_meta_data"
      );
      parameters.AddInitializer
      (
        "prefix",
        "file prefix to model name eg. MYPREFIX000000." + GetExtensionTypeDescriptor( e_MetaData),
        io::Serialization::GetAgent( &m_FilePrefix),
        ""
      );
      parameters.AddInitializer
      (
        "round",
        "round of descriptor selection",
        io::Serialization::GetAgent( &m_Round),
        "0"
      );
      parameters.AddInitializer
      (
        "iteration",
        "iteration number of particular assembly of descriptor groups",
        io::Serialization::GetAgent( &m_Iteration),
        "0"
      );
      parameters.AddInitializer
      (
        "cv_ind_id",
        "cross-validation id of chunk assigned as independent dataset",
        io::Serialization::GetAgent( &m_CrossValidationIdIndependent),
        "0"
      );
      parameters.AddInitializer
      (
        "cv_mon_id",
        "cross-validation id of chunk assigned as monitoring dataset",
        io::Serialization::GetAgent( &m_CrossValidationIdMonitoring),
        "1"
      );
      parameters.AddInitializer
      (
        "cv_total",
        "number of total cross-validation chunks",
        io::Serialization::GetAgent( &m_CrossValidationTotalChunks),
        "10"
      );

      return parameters;
    }

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    bool MetaDataStorageFile::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      // invoke method initialize
      return Initialize();
    }

  } // namespace model
} // namespace bcl
