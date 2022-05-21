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

#ifndef BCL_MODEL_META_DATA_STORAGE_FILE_H_
#define BCL_MODEL_META_DATA_STORAGE_FILE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_meta_data_storage_interface.h"
#include "io/bcl_io_directory.h"
#include "util/bcl_util_object_data_label.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MetaDataStorageFile
    //! @brief model::Interfaces are stored in a directory structure that follows that logic:
    //! {INITIALIZER}/{SOURCE}{KEY}.model
    //! @details KEY is always 6 chars long "000001", "000002"
    //!
    //! the key is auto incremented an is only numeric
    //! when initialized as attached, the largest key is located
    //! when initialized as created or overwrite, the key starts with 0
    //!
    //! @see @link example_model_meta_data_storage_file.cpp @endlink
    //! @author butkiem1
    //! @date Aug 19, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MetaDataStorageFile :
      public MetaDataStorageInterface
    {

    private:

    //////////
    // data //
    //////////

      //! @enum ExtensionType
      //! @brief enumerator for extension types
      enum ExtensionType
      {
        e_MetaData,       //!< file extension for meta data stored for descriptor selection
        e_Descriptor,  //!< file extension for descriptor label
        e_Predictions, //!< experimental and predicted data
        s_MaxExtensionType //!< max number of extension types
      };

      //! @brief ExtensionType as string
      //! @param INIT_TYPE the ExtensionType
      //! @return the string for the INIT_TYPE
      static const std::string &GetExtensionTypeDescriptor( const ExtensionType &TYPE)
      {
        static const std::string s_extensions[] =
        {
          "meta",
          "descriptor",
          "exp_pred",
          GetStaticClassName< ExtensionType>()
        };
        return s_extensions[ size_t( TYPE)];
      }

      //! directory name
      std::string m_DirectoryName;

      //! directory e.g. /home/user/run5/
      util::ShPtr< io::Directory> m_Directory;

      //! key to best model
      mutable std::string m_KeyToBestModelByResult;

      //! result of best model
      float m_ResultToBestModel;

      //! flag that indicates whether smallest result (true) or largest result (false) is best result
      bool m_BestResultSmallest;

      //! round number
      size_t m_Round;
      //! iteration number denotes a particular descriptor group chose at a particular iteration
      size_t m_Iteration;
      //! id of the independent cross-validation chunk
      size_t m_CrossValidationIdIndependent;
      //! id of the monitoring cross-validation chunk
      size_t m_CrossValidationIdMonitoring;
      //! number of total cross-validation chunks
      size_t m_CrossValidationTotalChunks;

      //! prefix for filename of stored model information
      std::string m_FilePrefix;

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

      //! @brief Clone function
      //! @return pointer to new MetaDataStorageFile
      MetaDataStorageFile *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get member variable directory
      //! @return shptr on current directory
      util::ShPtr< io::Directory> GetDirectory() const
      {
        return m_Directory;
      }

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief initialize the model storage
      //! @return true if initialize was successful
      bool Initialize();

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief number of models in source
      //! @return number of models with that source prefix
      size_t GetSize() const;

      //! @brief get all keys for given source
      //! @return all keys of given prefix
      storage::Vector< std::string> GetAllKeys() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief get descriptor set of model by key
      //! @param KEY key identifier for specific model/descriptor set in given source
      //! @return shptr to descriptorset of interest
      util::ObjectDataLabel RetrieveDescriptorSet( const std::string &KEY) const;

      //! @brief retrieve best descriptor set by result
      //! @return shptr to best descriptor set of model with best score/result
      util::ObjectDataLabel RetrieveBestDescriptorSetByResult() const;

      //! @brief get ensemble of descriptors associated each associated with a model
      //! @return shptr to list of descriptors of interest
      storage::List< util::ObjectDataLabel> RetrieveAllDescriptors() const;

      //! @brief store model with additional information in filesystem
      //! @param RESULT result value that was evaluated by objective function with current model
      //! @param DESCRIPTORS descriptors used to create dataset the current model was trained on
      //! @param METHOD_NAME name of iterate label that defines the used machine learning algorithm
      //! @param OBJECTIVE_FUNCTION name of objective function label that defines the objective function
      //! @return key of stored model
      std::string Store
      (
        const float RESULT,
        const util::ObjectDataLabel &DESCRIPTORS,
        const util::ObjectDataLabel &METHOD_NAME,
        const util::ObjectDataLabel &OBJECTIVE_FUNCTION
      );

      //! @brief store model with additional information in filesystem
      //! @param RESULT result value that was evaluated by objective function with current model
      //! @param DESCRIPTORS descriptors used to create dataset the current model was trained on
      //! @param METHOD_NAME name of iterate label that defines the used machine learning algorithm
      //! @param OBJECTIVE_FUNCTION name of objective function label that defines the objective function
      //! @param KEY preferred key for model to store
      //! @return key of stored model
      std::string Store
      (
        const float RESULT,
        const util::ObjectDataLabel &DESCRIPTORS,
        const util::ObjectDataLabel &METHOD_NAME,
        const util::ObjectDataLabel &OBJECTIVE_FUNCTION,
        const std::string &KEY
      );

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief construct complete filename from key
      //! @brief KEY the key for that protein
      //! @brief filename of form {path}/FilePrefix{KEY}.FileExtension
      std::string Filename( const std::string &KEY, const ExtensionType &TYPE) const;

      //! @brief Set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read of this class
      //! @param ERROR_STREAM stream with which to write errors
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    public:

      //! @brief check if key is valid string
      //! @param KEY the key to be checked
      //! @return true if the key is valid
      static bool IsValidKey( const std::string &KEY);

    }; // class MetaDataStorageFile

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_META_DATA_STORAGE_FILE_H_
