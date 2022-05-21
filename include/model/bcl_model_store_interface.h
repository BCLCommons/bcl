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

#ifndef BCL_MODEL_STORE_INTERFACE_H_
#define BCL_MODEL_STORE_INTERFACE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_model_cross_validation_info.h"
#include "bcl_model_interface.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class StoreInterface
    //! @brief interface to define methods model::Interfaces to store
    //! @details interface that defines methods for storing model::Interfaces in different data sources
    //!
    //! @see @link example_model_store_interface.cpp @endlink
    //! @author butkiem1
    //! @date Feb 24, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API StoreInterface :
      public util::SerializableInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new StoreInterface
      virtual StoreInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief number of models in source
      //! @return number of models
      virtual size_t GetSize() const = 0;

      //! @brief get all keys for given source
      //! @return all keys of given source
      virtual storage::Vector< std::string> GetAllKeys() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief try to store the result descriptor
      //! should assert if the labels are already stored and differ from those given
      //! @param RESULT_DESCRIPTOR result object data label for all models in this storage
      virtual void StoreResultDescriptor( const util::ObjectDataLabel &RESULT_DESCRIPTOR) = 0;

      //! @brief store model with additional information
      //! @param MODEL model::Interface that will be stored
      //! @param RESULT result value that was evaluated by objective function with current model
      //! @param DESCRIPTORS descriptors used to create dataset the current model was trained on
      //! @param METHOD_NAME name of iterate label that defines the used machine learning algorithm
      //! @param OBJ_FUNCTION objective function used in calculating the result value
      //! @param CV_INFO cross validation info
      //! @return key of stored model
      virtual std::string Store
      (
        const util::ShPtr< Interface> &MODEL,
        const float RESULT,
        const util::ObjectDataLabel &DESCRIPTORS,
        const util::ObjectDataLabel &METHOD_NAME,
        const util::ObjectDataLabel &OBJ_FUNCTION,
        const CrossValidationInfo &CV_INFO
      ) = 0;

      //! @brief store model with additional information
      //! @param MODEL model::Interface that will be stored
      //! @param RESULT result value that was evaluated by objective function with current model
      //! @param DESCRIPTORS descriptors used to create dataset the current model was trained on
      //! @param KEY preferred key for model to store
      //! @param METHOD_NAME name of iterate label that defines the used machine learning algorithm
      //! @param OBJ_FUNCTION objective function used in calculating the result value
      //! @param CV_INFO cross validation info
      //! @return key of stored model
      virtual std::string Store
      (
        const util::ShPtr< Interface> &MODEL,
        const float RESULT,
        const util::ObjectDataLabel &DESCRIPTORS,
        const util::ObjectDataLabel &METHOD_NAME,
        const util::ObjectDataLabel &OBJ_FUNCTION,
        const CrossValidationInfo &CV_INFO,
        const std::string &KEY
      ) = 0;

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class StoreInterface

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_STORE_INTERFACE_H_
