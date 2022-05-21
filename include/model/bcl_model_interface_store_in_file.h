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

#ifndef BCL_MODEL_INTERFACE_STORE_IN_FILE_H_
#define BCL_MODEL_INTERFACE_STORE_IN_FILE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_store_interface.h"
#include "io/bcl_io_directory.h"
#include "io/bcl_io_serialization_interface.h"
#include "util/bcl_util_cleanable_interface.h"
#include "util/bcl_util_object_data_label.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class InterfaceStoreInFile
    //! @brief model::Interfaces are stored in a directory structure that follows that logic:
    //! {INITIALIZER}/{SOURCE}{KEY}.model
    //! @details KEY is always 6 chars long "0000001", "0000002"
    //!
    //! the key is auto incremented an is only numeric
    //! when initialized as attached, the largest key is located
    //! when initialized as created or overwrite, the key starts with 0
    //!
    //! @see @link example_model_interface_store_in_file.cpp @endlink
    //! @author butkiem1
    //! @date Feb 25, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API InterfaceStoreInFile :
      public StoreInterface,
      public util::CleanableInterface
    {

    public:

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
        e_Info,        //!< file extension for cross validation info
        e_Predictions, //!< experimental and predicted data
        s_MaxExtensionType //!< max number of extension types
      };

      //! @brief ExtensionType as string
      //! @param INIT_TYPE the ExtensionType
      //! @return the string for the INIT_TYPE
      static const std::string &GetExtensionTypeDescriptor( const ExtensionType &TYPE);

    private:

      // string representation of directory
      std::string m_DirectoryName;

      //! directory e.g. /home/user/run5/
      io::Directory m_Directory;

      //! flag for serializing descriptors
      bool m_WriteDescriptors;

      //! flag for serializing experimental and predicted data
      bool m_WriteExperimentalPredicted;

      //! prefix for filename of stored model information
      std::string m_FilePrefix;

      //! model key if given
      size_t m_ModelKey;

      //! True if this storage currently has a lock file open
      bool m_HaveOpenLockFile;

      //! Compression type
      io::StreamBufferClass m_Compression;

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
      InterfaceStoreInFile();

      //! @brief default constructor
      InterfaceStoreInFile( const std::string &DIRECTORY_NAME);

      //! @brief destructor
      ~InterfaceStoreInFile();

      //! @brief Clone function
      //! @return pointer to new InterfaceStoreInFile
      InterfaceStoreInFile *Clone() const;

      //! @brief cleanup function, to remove lock file if it exists
      void CleanUp();

    /////////////////
    // data access //
    /////////////////

      //! @brief get member variable directory
      //! @return current directory
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

      //! @brief try to store the result descriptor
      //! should assert if the labels are already stored and differ from those given
      //! @param RESULT_DESCRIPTOR result object data label for all models in this storage
      void StoreResultDescriptor( const util::ObjectDataLabel &RESULT_DESCRIPTOR);

      //! @brief store model with additional information
      //! @param MODEL model::Interface that will be stored
      //! @param RESULT result value that was evaluated by objective function with current model
      //! @param DESCRIPTORS descriptors used to create dataset the current model was trained on
      //! @param EXPERIMENTAL_PREDICTED_VALUES serialized object of a list with pairs of
      //!        experimental and predicted values
      //! @param METHOD_NAME name of iterate label that defines the used machine learning algorithm
      //! @param OBJ_FUNCTION objective function used in calculating the result value
      //! @param CV_INFO cross validation info
      //! @return key of stored model
      std::string Store
      (
        const util::ShPtr< Interface> &MODEL,
        const float RESULT,
        const util::ObjectDataLabel &DESCRIPTORS,
        const util::ObjectDataLabel &METHOD_NAME,
        const util::ObjectDataLabel &OBJ_FUNCTION,
        const CrossValidationInfo &CV_INFO
      );

      //! @brief store model with additional information
      //! @param MODEL model::Interface that will be stored
      //! @param RESULT result value that was evaluated by objective function with current model
      //! @param DESCRIPTORS descriptors used to create dataset the current model was trained on
      //! @param EXPERIMENTAL_PREDICTED_VALUES serialized object of a list with pairs of
      //!        experimental and predicted values
      //! @param KEY preferred key for model to store
      //! @param METHOD_NAME name of iterate label that defines the used machine learning algorithm
      //! @param OBJ_FUNCTION objective function used in calculating the result value
      //! @param CV_INFO cross validation info
      //! @return key of stored model
      std::string Store
      (
        const util::ShPtr< Interface> &MODEL,
        const float RESULT,
        const util::ObjectDataLabel &DESCRIPTORS,
        const util::ObjectDataLabel &METHOD_NAME,
        const util::ObjectDataLabel &OBJ_FUNCTION,
        const CrossValidationInfo &CV_INFO,
        const std::string &KEY
      );

      //! @brief store model with additional information
      //! @param MODEL model::Interface that will be stored
      //! @param RESULT result value that was evaluated by objective function with current model
      //! @param DESCRIPTORS descriptors used to create dataset the current model was trained on
      //! @return key of stored model
      std::string Store
      (
        const util::ShPtr< Interface> &MODEL,
        const util::ObjectDataLabel &DESCRIPTORS
      );

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

      //! @brief write util::ObjectInterface from source and key
      //! @param OBJECT current object that should be written in file
      //! @param KEY current key
      //! @param EXTENSION_TYPE extention type of file the object will be written into
      void WriteToFile( const util::ObjectInterface &OBJECT, const std::string &KEY, const ExtensionType &EXTENSION_TYPE) const;

      //! @brief write util::ObjectDataLabel from source and key
      //! @param OBJECT current object that should be written in file
      //! @param KEY current key
      //! @param EXTENSION_TYPE extention type of file the object will be written into
      void WriteToFile( const util::ObjectDataLabel &LABEL, const std::string &KEY, const ExtensionType &EXTENSION_TYPE) const;

    public:

      //! @brief Set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read of this class
      //! @param ERROR_STREAM stream with which to write errors
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    private:

      //! @brief check if key is valid string
      //! @param KEY the key to be checked
      //! @return true if the key is valid
      static bool IsValidKey( const std::string &KEY);

    }; // class InterfaceStoreInFile

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_INTERFACE_STORE_IN_FILE_H_
