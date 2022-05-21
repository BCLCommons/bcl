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

#ifndef BCL_MODEL_INTERFACE_RETRIEVE_FROM_FILE_H_
#define BCL_MODEL_INTERFACE_RETRIEVE_FROM_FILE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_retrieve_interface.h"
#include "io/bcl_io_directory.h"
#include "io/bcl_io_serialization_interface.h"
#include "util/bcl_util_object_data_label.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class InterfaceRetrieveFromFile
    //! @brief model::Interfaces or their associated descriptor sets are retrieved from a directory structure that
    //! follows that logic:
    //! {your directory}/{prefix}{KEY}.{model|descriptor|exp_pred}
    //! @details KEY is always 6 chars long "0000001", "0000002"
    //!
    //! the key is auto incremented an is only numeric
    //!
    //! @see @link example_model_interface_retrieve_from_file.cpp @endlink
    //! @author butkiem1
    //! @date Feb 24, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API InterfaceRetrieveFromFile :
      public RetrieveInterface
    {

    private:

    //////////
    // data //
    //////////

      //! @enum ExtensionType
      //! @brief enumerator for extension types
      enum ExtensionType
      {
        e_Model,       //!< file extension for model serialization
        e_Descriptor,  //!< file extension for descriptor label
        e_Result,      //!< file extension for result labels
        e_Info,        //!< file extension for cross validation information
        s_MaxExtensionType //!< max number of extension types
      };

      //! @brief ExtensionType as string
      //! @param INIT_TYPE the ExtensionType
      //! @return the string for the INIT_TYPE
      static const std::string &GetExtensionTypeDescriptor( const ExtensionType &TYPE);

      // string representation of directory
      std::string m_DirectoryName;

      //! directory e.g. /home/user/run5/
      io::Directory m_Directory;

      //! prefix for filename of stored model information
      std::string m_FilePrefix;

      //! Single key of model in the storage, if specified by the user
      std::string m_Key;

      //! absolute path of the directory; cached to reduce calls to getcwd
      std::string m_AbsolutePath;

      //! Bool that allows for choosing the best model for each unique independent set in the storage
      bool m_SelectBestModels;

    public:

    //////////
    // data //
    //////////

      //! format to convert Key to string
      static const util::Format s_KeyToString;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      InterfaceRetrieveFromFile();

      //! @brief construct from directory name and file prefix
      //! @param DIRECTORY_NAME name of directory, e.g. /home/user/models/
      //! @param FILE_PREFIX prefix to all files for this model
      InterfaceRetrieveFromFile
      (
        const std::string &DIRECTORY_NAME,
        const std::string &FILE_PREFIX
      );

      //! @brief Clone function
      //! @return pointer to new InterfaceRetrieveFromFile
      InterfaceRetrieveFromFile *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get member variable directory
      //! @return shptr on current directory
      const io::Directory &GetDirectory() const
      {
        return m_Directory;
      }

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief number of models in source
      //! @return number of models with that source prefix
      size_t GetSize() const;

      //! @brief get all keys for given source
      //! @return all keys of given prefix
      storage::Vector< std::string> GetAllKeys() const;

      //! @brief get all keys for given range
      //! @param RANGE range of key values
      //! @return all keys of given prefix
      storage::Vector< std::string> GetKeys( const math::Range< size_t> &RANGE) const;

    ////////////////
    // operations //
    ////////////////

      //! @brief get model
      //! @param KEY key identifier for specific model in given source like "000001"
      //! @return shptr to model of interest, undefined if there is no such model
      RetrieveInterface::t_ModelPtr Retrieve( const std::string &KEY) const;

      //! @brief get result descriptor for the storage
      //! @return ObjectDataLabel for the results in this storage
      util::ObjectDataLabel RetrieveResultDescriptor() const;

      //! @brief get result descriptor located in the given filename
      //! @param FILENAME filename in which the results descriptor is located
      //! @return ObjectDataLabel for the results in this storage
      static util::ObjectDataLabel RetrieveResultDescriptor( const std::string &FILENAME);

      //! @brief get descriptor set of model by key
      //! @param KEY key identifier for specific model/descriptor set in given source
      //! @return shptr to descriptorset of interest
      util::ObjectDataLabel RetrieveDescriptorSet( const std::string &KEY) const;

      //! @brief get cv info of model by key
      //! @param KEY key identifier for specific model/descriptor set in given source
      //! @return cross validation information for the model
      CrossValidationInfo RetrieveCVInfo( const std::string &KEY) const;

      //! @brief get ensemble of descriptors associated each associated with a model
      //! @return shptr to list of descriptors of interest
      storage::List< util::ObjectDataLabel> RetrieveEnsembleDescriptors() const;

      //! @brief get ensemble of descriptors associated each associated with a model
      //! @return shptr to list of descriptors of interest
      storage::List< CrossValidationInfo> RetrieveEnsembleCVInfo() const;

      //! @brief get list of models
      //! @return shptr list of models from given source
      RetrieveInterface::t_Container RetrieveEnsemble() const;

      //! @brief get list of models for given keys
      //! @param KEYS vector of identifiers for specific models in given source like "000001", "000002"
      //! @return shptr list of models from given source
      RetrieveInterface::t_Container RetrieveEnsemble( const storage::Vector< std::string> &KEYS) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief construct complete filename from key
      //! @brief KEY the key for that protein
      //! @brief filename of form {initializer}/s_FilePrefix{KEY}.s_ModelFileExtension
      std::string Filename( const std::string &KEY, const ExtensionType &EXTENSION_TYPE) const;

      //! @brief helper function to select the best models for each unique independent set
      //! @param KEYS the set of keys to begin with
      //! @return a set of keys; one key for each unique independent set
      storage::Vector< std::string> FindBestModelForEachIndependentSet
      (
        const storage::Vector< std::string> &KEYS
      ) const;

    public:

      //! @brief Set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read of this class
      //! @param ERROR_STREAM stream with which to write errors
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    }; // class InterfaceRetrieveFromFile

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_INTERFACE_RETRIEVE_FROM_FILE_H_
