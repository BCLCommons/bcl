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
#include "model/bcl_model_interface_retrieve_from_file.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_interface.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "model/bcl_model.h"
#include "sched/bcl_sched_mutex.h"
#include "util/bcl_util_sh_ptr_list.h"
#include "util/bcl_util_stopwatch.h"
#include "util/bcl_util_time.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    //! @brief ExtensionType as string
    //! @param INIT_TYPE the ExtensionType
    //! @return the string for the INIT_TYPE
    const std::string &InterfaceRetrieveFromFile::GetExtensionTypeDescriptor( const ExtensionType &TYPE)
    {
      static const std::string s_extensions[] =
      {
        "model",
        "descriptor",
        "result",
        "info",
        GetStaticClassName< ExtensionType>()
      };
      return s_extensions[ size_t( TYPE)];
    }

    //! format to convert Key to string
    const util::Format InterfaceRetrieveFromFile::s_KeyToString
    (
      util::Format().W( 6).R().Fill( '0')
    );

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> InterfaceRetrieveFromFile::s_Instance
    (
      util::Enumerated< RetrieveInterface>::AddInstance
      (
        new InterfaceRetrieveFromFile()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    InterfaceRetrieveFromFile::InterfaceRetrieveFromFile() :
      m_DirectoryName(),
      m_FilePrefix( ""),
      m_SelectBestModels( false)
    {
    }

    //! @brief construct from directory name and file prefix
    //! @param DIRECTORY_NAME name of directory, e.g. /home/user/models/
    //! @param FILE_PREFIX prefix to all files for this model
    InterfaceRetrieveFromFile::InterfaceRetrieveFromFile
    (
      const std::string &DIRECTORY_NAME,
      const std::string &FILE_PREFIX
    ) :
      m_DirectoryName( DIRECTORY_NAME),
      m_Directory( m_DirectoryName),
      m_FilePrefix( FILE_PREFIX),
      m_SelectBestModels( false)
    {
      ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief Clone function
    //! @return pointer to new InterfaceRetrieveFromFile
    InterfaceRetrieveFromFile *InterfaceRetrieveFromFile::Clone() const
    {
      return new InterfaceRetrieveFromFile( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &InterfaceRetrieveFromFile::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &InterfaceRetrieveFromFile::GetAlias() const
    {
      static const std::string s_Name( "File");
      return s_Name;
    }

    //! @brief number of objects in source
    //! @return number of objects with that source prefix
    size_t InterfaceRetrieveFromFile::GetSize() const
    {
      return GetAllKeys().GetSize();
    }

    //! @brief get all keys for given range
    //! @param RANGE range of key values
    //! @return all keys of given prefix
    storage::Vector< std::string> InterfaceRetrieveFromFile::GetKeys( const math::Range< size_t> &RANGE) const
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
    storage::Vector< std::string> InterfaceRetrieveFromFile::GetAllKeys() const
    {
      if( !m_Key.empty())
      {
        // singular key, just return it
        return storage::Vector< std::string>( size_t( 1), m_Key);
      }

      // cache
      static sched::Mutex s_mutex;
      static storage::Map< util::ObjectDataLabel, storage::Vector< std::string> > s_keys;

      s_mutex.Lock();
      storage::Vector< std::string> &keys( s_keys[ GetLabel()]);

      if( !keys.IsEmpty())
      {
        s_mutex.Unlock();
        return keys;
      }

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
      keys.Sort( std::less< std::string>());

      if( m_SelectBestModels)
      {
        keys = FindBestModelForEachIndependentSet( keys);
      }
      s_mutex.Unlock();

      return keys;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get model
    //! @param KEY key identifier for specific smallmolecule in given source like "000001"
    //! @return shptr to smallmolecule of interest, undefined if there is no such smallmolecule
    RetrieveInterface::t_ModelPtr InterfaceRetrieveFromFile::Retrieve( const std::string &KEY) const
    {
      // model
      util::ShPtr< Interface> sp_interface;

      const std::string key( s_KeyToString( KEY));

      // check the key
      if( !IsValidKey( key))
      {
        BCL_MessageCrt( KEY + " is not a valid model key should be 6 integers, e.g. 123456");
        return RetrieveInterface::t_ModelPtr();
      }

      const std::string basename
      (
        m_AbsolutePath + "/" + m_FilePrefix + key + "." + GetExtensionTypeDescriptor( e_Model)
      );

      // if storing models ever becomes a memory issue, move this into another class and add a reset function
      // As of this writing, memory usage for models is typically trivial compared with dataset storage
      static sched::Mutex s_mutex;
      static storage::Map< std::string, RetrieveInterface::t_ModelPtr> s_results;

      s_mutex.Lock();
      RetrieveInterface::t_ModelPtr &interface( s_results[ basename]);
      if( interface.IsDefined())
      {
        s_mutex.Unlock();
        return t_ModelPtr( interface.operator ->(), false);
      }
      // filename and check if exists
      const std::string filename( Filename( key, e_Model));

      // check if file with that name exists
      if( !io::DirectoryEntry( filename).DoesExist())
      {
        BCL_MessageCrt( "could not find model " + filename);
        s_mutex.Unlock();
        return RetrieveInterface::t_ModelPtr();
      }

      BCL_MessageVrb( "Reading model with key = " + key + " at " + filename);
      // retrieve stored model
      io::IFStream input;

      static util::Stopwatch s_timer( "Reading models", util::Time( 0, 500), util::Message::e_Verbose, true, false);
      s_timer.Start();
      io::File::MustOpenIFStream( input, filename);
      io::Serialize::Read( sp_interface, input);
      io::File::CloseClearFStream( input);
      s_timer.Stop();

      // return retrieved small molecule
      interface = RetrieveInterface::t_ModelPtr( sp_interface->Clone());
      s_mutex.Unlock();
      return t_ModelPtr( interface.operator ->(), false);
    }

    //! @brief get result descriptor for the storage
    //! @return ObjectDataLabel for the results in this storage
    util::ObjectDataLabel InterfaceRetrieveFromFile::RetrieveResultDescriptor() const
    {
      static sched::Mutex s_mutex;
      static storage::Map< std::string, util::ObjectDataLabel> s_results;

      s_mutex.Lock();
      util::ObjectDataLabel &result( s_results[ m_AbsolutePath + "/" + m_FilePrefix]);
      if( !result.IsEmpty())
      {
        s_mutex.Unlock();
        return result;
      }
      result = this->RetrieveResultDescriptor( Filename( "", e_Result));
      s_mutex.Unlock();
      return result;
    }

    //! @brief get result descriptor located in the given filename
    //! @param FILENAME filename in which the results descriptor is located
    //! @return ObjectDataLabel for the results in this storage
    util::ObjectDataLabel InterfaceRetrieveFromFile::RetrieveResultDescriptor
    (
      const std::string &FILENAME
    )
    {
      // check if file with that name exists
      if( !io::DirectoryEntry( FILENAME).DoesExist())
      {
        return util::ObjectDataLabel();
      }

      // create factory and retrieve protein model
      io::IFStream input;

      // check that the writer's lock file isn't currently locked
      const io::DirectoryEntry lock_file( FILENAME + ".lock");
      while( lock_file.DoesExist())
      {
        // sleep a millisecond and check again
        util::Time::Delay( util::Time( 0, 1));
      }

      // clean up the path name for purposes of creating a semaphore
      io::File::MustOpenIFStream( input, FILENAME);
      util::ObjectDataLabel descriptors( input);
      io::File::CloseClearFStream( input);

      // return retrieved descriptor set
      return descriptors;
    }

    //! @brief get descriptor set of model by key
    //! @param KEY key identifier for specific model/descriptor set in given source
    //! @return shptr to descriptorset of interest
    util::ObjectDataLabel
    InterfaceRetrieveFromFile::RetrieveDescriptorSet( const std::string &KEY) const
    {
      // check the key
      if( !IsValidKey( KEY))
      {
        return util::ObjectDataLabel();
      }

      // filename and check if exists
      std::string filename( Filename( KEY, e_Descriptor));

      // check if file with that name exists
      if( !io::DirectoryEntry( filename).DoesExist())
      {
        // Often the descriptors are identical for all models in a given directory. In this case, allow there to be a
        // default descriptors file in the folder to avoid duplication
        filename = Filename( "", e_Descriptor);
        if( !io::DirectoryEntry( filename).DoesExist())
        {
          return util::ObjectDataLabel();
        }
      }

      // create factory and retrieve protein model
      io::IFStream input;
      io::File::MustOpenIFStream( input, filename);
      util::ObjectDataLabel descriptors( input);
      io::File::CloseClearFStream( input);

      // return retrieved descriptor set
      return descriptors;
    }

    //! @brief get cv info of model by key
    //! @param KEY key identifier for specific model/descriptor set in given source
    //! @return cross validation information for the model
    CrossValidationInfo InterfaceRetrieveFromFile::RetrieveCVInfo( const std::string &KEY) const
    {
      // check the key
      if( !IsValidKey( KEY))
      {
        return CrossValidationInfo();
      }

      // filename and check if exists
      const std::string filename( Filename( KEY, e_Info));

      // check if file with that name exists
      if( !io::DirectoryEntry( filename).DoesExist())
      {
        return CrossValidationInfo();
      }

      // return retrieved cv info
      return CrossValidationInfo( filename);
    }

    //! @brief get ensemble of descriptors associated each associated with a model
    //! @return shptr to list of descriptors of interest
    storage::List< util::ObjectDataLabel> InterfaceRetrieveFromFile::RetrieveEnsembleDescriptors() const
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

    //! @brief get ensemble of descriptors associated each associated with a model
    //! @return shptr to list of descriptors of interest
    storage::List< CrossValidationInfo> InterfaceRetrieveFromFile::RetrieveEnsembleCVInfo() const
    {
      // get all keys
      storage::Vector< std::string> keys( GetAllKeys());

      // list for retrieved descriptors
      storage::List< CrossValidationInfo> cv_infos;

      for
      (
        storage::Vector< std::string>::const_iterator itr_keys( keys.Begin()), itr_keys_end( keys.End());
        itr_keys != itr_keys_end;
        ++itr_keys
      )
      {
        // insert descriptor set according to key
        cv_infos.PushBack( RetrieveCVInfo( *itr_keys));
      }

      // return list with descriptors
      return cv_infos;
    }

    //! @brief get ensemble of objects
    //! @return shptr list of smallmolececules from given source
    RetrieveInterface::t_Container InterfaceRetrieveFromFile::RetrieveEnsemble() const
    {
      // acquire all keys for that source and return ensemble
      return RetrieveEnsemble( GetAllKeys());
    }

    //! @brief get ensemble of objects for given keys
    //! @param KEYS vector of identifiers for specific objects in given source like "000001", "000002"
    //! @return shptr list of objects from given source
    RetrieveInterface::t_Container InterfaceRetrieveFromFile::RetrieveEnsemble
    (
      const storage::Vector< std::string> &KEYS
    ) const
    {

      RetrieveInterface::t_Container ensemble;

      // iterate over all keys
      for( storage::Vector< std::string>::const_iterator itr( KEYS.Begin()), itr_end( KEYS.End()); itr != itr_end; ++itr)
      {
        RetrieveInterface::t_ModelPtr current( Retrieve( *itr));
        if( current.IsDefined())
        {
          ensemble.PushBack( current);
        }
      }
      // end
      return ensemble;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief construct complete filename from source and key
    //! @brief KEY the key for that protein
    //! @brief filename of form {initializer}/m_FilePrefix{KEY}.s_ModelFileExtension
    std::string InterfaceRetrieveFromFile::Filename( const std::string &KEY, const ExtensionType &EXTENSION_TYPE) const
    {
      const std::string basename
      (
        m_FilePrefix + KEY + io::File::GetExtensionDelimiter() + GetExtensionTypeDescriptor( EXTENSION_TYPE)
      );
      util::ShPtr< io::Directory> directory( m_Directory.Clone());
      // test compressed extensions
      for
      (
        io::StreamBufferClasses::const_iterator
          itr( io::GetStreamBufferClasses().Begin()), itr_end( io::GetStreamBufferClasses().End());
        itr != itr_end;
        ++itr
      )
      {
        const io::DirectoryEntry compressed_file( directory, ( **itr)->AddExtension( basename));
        if( compressed_file.DoesExist())
        {
          return compressed_file.GetFullName();
        }
      }
      // file does not exist with any compression type; return filename sans extension
      return m_Directory.AppendFilename( basename);
    }

    //! @brief helper function to select the best models for each unique independent set
    //! @param KEYS the set of keys to begin with
    //! @return a set of keys; one key for each unique independent set
    storage::Vector< std::string> InterfaceRetrieveFromFile::FindBestModelForEachIndependentSet
    (
      const storage::Vector< std::string> &KEYS
    ) const
    {
      if( KEYS.GetSize() < size_t( 2))
      {
        return KEYS;
      }

      // list for retrieved descriptors
      storage::Map< util::ObjectDataLabel, storage::Pair< float, std::string> > independent_set_to_best_result_and_key;

      // get the 1st cross-validation info
      const CrossValidationInfo first_cross_validation_info( RetrieveCVInfo( KEYS.FirstElement()));

      // determine the improvement type
      const opti::ImprovementTypeEnum improvement_type( first_cross_validation_info.GetImprovementType());

      // insert this value into the map
      independent_set_to_best_result_and_key[ first_cross_validation_info.GetIndependentDatasetRetriever()] =
        storage::Pair< float, std::string>( first_cross_validation_info.GetResult(), KEYS.FirstElement());

      for
      (
        storage::Vector< std::string>::const_iterator itr_keys( KEYS.Begin() + 1), itr_keys_end( KEYS.End());
        itr_keys != itr_keys_end;
        ++itr_keys
      )
      {
        const CrossValidationInfo cross_validation_info( RetrieveCVInfo( *itr_keys));
        BCL_Assert
        (
          improvement_type == cross_validation_info.GetImprovementType(),
          "Models in a given storage should all have the same improvement type but did not in " + this->GetString()
        );
        // find the current independent set in the map
        storage::Map< util::ObjectDataLabel, storage::Pair< float, std::string> >::iterator
          itr_independent
          (
            independent_set_to_best_result_and_key.Find( cross_validation_info.GetIndependentDatasetRetriever())
          );

        // if the independent set is not current in the map, insert this result and continue
        if( itr_independent == independent_set_to_best_result_and_key.End())
        {
          independent_set_to_best_result_and_key[ cross_validation_info.GetIndependentDatasetRetriever()] =
            storage::Pair< float, std::string>( cross_validation_info.GetResult(), *itr_keys);
          continue;
        }

        // independent set is already in the map.  Check whether this model performed better on it
        const float prior_result( itr_independent->second.First());
        if( opti::DoesImprove( cross_validation_info.GetResult(), prior_result, improvement_type))
        {
          // update with this key, since the corresponding model was better
          itr_independent->second.First() = cross_validation_info.GetResult();
          itr_independent->second.Second() = *itr_keys;
        }
      }

      // copy the best key for each cross validation into a vector and return it
      storage::Vector< std::string> best_keys;
      best_keys.AllocateMemory( independent_set_to_best_result_and_key.GetSize());
      for
      (
        storage::Map< util::ObjectDataLabel, storage::Pair< float, std::string> >::const_iterator
          itr( independent_set_to_best_result_and_key.Begin()), itr_end( independent_set_to_best_result_and_key.End());
        itr != itr_end;
        ++itr
      )
      {
        best_keys.PushBack( itr->second.Second());
      }
      best_keys.Sort( std::less< std::string>());
      return best_keys;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer InterfaceRetrieveFromFile::GetSerializer() const
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
        "pick best",
        "If set, only the best model for each unique independent set will be retrieved. "
        "The \"best\" model is the one with the best final objective function value on its independent set",
        io::Serialization::GetAgent( &m_SelectBestModels),
        "False"
      );
      parameters.AddOptionalInitializer
      (
        "key",
        "Single model key to retrieve; if omitted, all models from the directory with the given prefix are used",
        io::Serialization::GetAgent( &m_Key)
      );
      return parameters;
    }

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    bool InterfaceRetrieveFromFile::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      // create directory
      m_Directory = io::Directory( m_DirectoryName);
      if( !m_Directory.DoesExist() && !io::File::IsAbsolutePath( m_DirectoryName) && Model::GetModelPathFlag()->GetFlag())
      {
        m_Directory = io::Directory( Model::AddModelPath( m_DirectoryName));
      }
      m_AbsolutePath = io::File::MakeAbsolutePath( m_Directory.GetPath());
      m_Directory = io::Directory( m_AbsolutePath);

      // check whether storage directory exists
      if( !m_Directory.DoesExist())
      {
        ERROR_STREAM << "storage directory does not exist: " << m_Directory.GetPath();
        return false;
      }

      if( !m_Key.empty())
      {
        m_Key = s_KeyToString( m_Key);
        if( !io::DirectoryEntry( Filename( m_Key, e_Model)).DoesExist())
        {
          ERROR_STREAM
            << "given key did not exist in the filename at: "
            << Filename( m_Key, e_Model);
          return false;
        }
      }

      return true;
    }

  } // namespace model
} // namespace bcl
