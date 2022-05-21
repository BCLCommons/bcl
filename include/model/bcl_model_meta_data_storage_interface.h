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

#ifndef BCL_MODEL_META_DATA_STORAGE_INTERFACE_H_
#define BCL_MODEL_META_DATA_STORAGE_INTERFACE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_model_interface.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MetaDataStorageInterface
    //! @brief interface to define model::Interface storage
    //! @details interface that defines methods for storing and retrieving model::Interfaces in different data sources
    //!
    //! @see @link example_model_interface_storage_interface.cpp @endlink
    //! @author butkiem1
    //! @date Sep 3, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MetaDataStorageInterface :
      public util::SerializableInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new MetaDataStorageInterface
      virtual MetaDataStorageInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief initialize the model storage
      //! @return true if initialize was successful
      virtual bool Initialize() = 0;

      //! @brief number of models in source
      //! @return number of models
      virtual size_t GetSize() const = 0;

      //! @brief get all keys for given source
      //! @return all keys of given source
      virtual storage::Vector< std::string> GetAllKeys() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief get descriptor set of model by key
      //! @param KEY key identifier for specific model/descriptor set in given source
      //! @return shptr to descriptorset of interest
      virtual util::ObjectDataLabel RetrieveDescriptorSet( const std::string &KEY) const = 0;

      //! @brief retrieve best descriptor set by result
      //! @param SMALLEST_RESULT flag wether the best model has smallest result
      //! @return shptr to best descriptor set of model with best score/result
      virtual util::ObjectDataLabel RetrieveBestDescriptorSetByResult() const = 0;

      //! @brief get ensemble of descriptors associated each associated with a model
      //! @return shptr to list of descriptors of interest
      virtual storage::List< util::ObjectDataLabel> RetrieveAllDescriptors() const = 0;

      //! @brief store model with additional information
      //! @param RESULT result value that was evaluated by objective function with current model
      //! @param DESCRIPTORS descriptors used to create dataset the current model was trained on
      //! @param METHOD_NAME name of iterate label that defines the used machine learning algorithm
      //! @param OBJECTIVE_FUNCTION name of objective function label that defines the objective function
      //! @return key of stored model
      virtual std::string Store
      (
        const float RESULT,
        const util::ObjectDataLabel &DESCRIPTORS,
        const util::ObjectDataLabel &METHOD_NAME,
        const util::ObjectDataLabel &OBJECTIVE_FUNCTION
      ) = 0;

      //! @brief store model with additional information
      //! @param RESULT result value that was evaluated by objective function with current model
      //! @param DESCRIPTORS descriptors used to create dataset the current model was trained on
      //! @param KEY preferred key for model to store
      //! @param METHOD_NAME name of iterate label that defines the used machine learning algorithm
      //! @param OBJECTIVE_FUNCTION name of objective function label that defines the objective function
      //! @return key of stored model
      virtual std::string Store
      (
        const float RESULT,
        const util::ObjectDataLabel &DESCRIPTORS,
        const util::ObjectDataLabel &METHOD_NAME,
        const util::ObjectDataLabel &OBJECTIVE_FUNCTION,
        const std::string &KEY
      ) = 0;

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class MetaDataStorageInterface

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_META_DATA_STORAGE_INTERFACE_H_
