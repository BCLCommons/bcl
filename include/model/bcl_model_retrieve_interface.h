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

#ifndef BCL_MODEL_RETRIEVE_INTERFACE_H_
#define BCL_MODEL_RETRIEVE_INTERFACE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically
#include "descriptor/bcl_descriptor.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_model_cross_validation_info.h"
#include "bcl_model_interface.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_serializable_interface.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RetrieveInterface
    //! @brief interface to define model::Interface retrieve functionality
    //! @details interface that defines methods for retrieving model::Interfaces in different data sources
    //!
    //! @see @link example_model_retrieve_interface.cpp @endlink
    //! @author butkiem1
    //! @date Feb 24, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RetrieveInterface :
      public util::SerializableInterface
    {

    public:

      //! Typedef for the container the models are stored in. OwnPtr allows control over whether models are cached
      //! or not
      typedef util::OwnPtr< Interface>     t_ModelPtr;
      typedef storage::Vector< t_ModelPtr> t_Container;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new RetrieveInterface
      virtual RetrieveInterface *Clone() const = 0;

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

      //! @brief get model by key
      //! @param KEY key identifier for specific model in given source
      //! @return shptr to model of interest
      virtual t_ModelPtr Retrieve( const std::string &KEY) const = 0;

      //! @brief get result descriptor for the storage
      //! @return ObjectDataLabel for the results in this storage
      virtual util::ObjectDataLabel RetrieveResultDescriptor() const = 0;

      //! @brief get descriptor set of model by key
      //! @param KEY key identifier for specific model/descriptor set in given source
      //! @return shptr to descriptorset of interest
      virtual util::ObjectDataLabel RetrieveDescriptorSet( const std::string &KEY) const = 0;

      //! @brief get ensemble of descriptors associated each associated with a model
      //! @return shptr to list of descriptors of interest
      virtual storage::List< util::ObjectDataLabel> RetrieveEnsembleDescriptors() const = 0;

      //! @brief Retrieves all descriptors, asserts if they are not all identical
      //! @return the unique descriptor
      virtual util::ObjectDataLabel RetrieveUniqueDescriptor() const;

      //! @brief get ensemble of descriptors associated each associated with a model
      //! @return shptr to list of descriptors of interest
      virtual storage::List< CrossValidationInfo> RetrieveEnsembleCVInfo() const = 0;

      //! @brief get the cv info that is held in common across all CV in the dataset
      //! @return the commonly-held cross-validation info
      virtual CrossValidationInfo RetrieveCommonCVInfo() const;

      //! @brief read the merged independent dataset predictions into a dataset object, whose features are the predicted values
      //!        results are the same result values, and ids are the given ids
      virtual util::ShPtr< descriptor::Dataset> ReadMergedIndependentPredictions();

      //! @brief read the merged independent dataset predictions into a dataset object, whose features are the predicted values
      //!        results are the same result values, and ids are the given ids
      storage::Vector< math::ROCCurve> GetMergedIndependentROCCurves();

      //! @brief transform an experimental/predicted matrix into a series of ROC curves
      //! @param EXP experimental results
      //! @param PRED predicted results
      //! @param THRESHOLD threshold for defining Positive/Negative class split
      //! @param POS_ABOVE_THRESH true if positives are those above the threshold
      static storage::Vector< math::ROCCurve> ROCCurvesFromDataset
      (
        const linal::MatrixConstInterface< float> &EXP,
        const linal::MatrixConstInterface< float> &PRED,
        const float &THRESHOLD,
        const bool &POS_ABOVE_THRESH
      );

      //! @brief get model ensemble with all models
      //! @return shptr to list of model of interest
      virtual t_Container RetrieveEnsemble() const = 0;

      //! @brief get model ensemble by keys
      //! @param KEYS vector with all key ids of models
      //! @return shptr to list of model of interest
      virtual t_Container RetrieveEnsemble( const storage::Vector< std::string> &KEYS) const = 0;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief check if key is valid string
      //! @param KEY the key to be checked
      //! @return true if the key is valid
      bool IsValidKey( const std::string &KEY) const
      {
        // check that the key is of size 6 and that it is an unsigned integer
        return KEY.length() == 6 && util::LengthOfUnsignedIntegerType( KEY) == 6;
      }

    }; // class RetrieveInterface

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_RETRIEVE_INTERFACE_H_
