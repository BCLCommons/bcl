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

#ifndef BCL_DESCRIPTOR_ITERATIVE_PREDICTION_H_
#define BCL_DESCRIPTOR_ITERATIVE_PREDICTION_H_

// include the namespace header
#include "bcl_descriptor.h"

// other forward includes
#include "model/bcl_model.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base.h"
#include "bcl_descriptor_prediction.h"
#include "model/bcl_model_retrieve_interface.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class IterativePrediction
    //! @brief Performs iterative evaluation of a predictive model that uses predictions from other models internally
    //! @detail Models can be created that are trained with output of other models that were trained to predict the same
    //!         output measure.  For example, models may be created to predict secondary structure based on one dataset,
    //!         while a second layer is shown various descriptors based off the outputs of the first layer.  In these
    //!         cases, it can be useful at each layer beyond the first to predict using the first layers predictions,
    //!         then recompute all prediction-based descriptors using the output of the second layer in an iterative
    //!         fashion or until convergence.
    //!
    //! @see @link example_descriptor_iterative_prediction.cpp @endlink
    //! @author mendenjl
    //! @date Mar 13, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class IterativePrediction :
      public Prediction< t_DataType>
    {
    private:

    //////////
    // data //
    //////////

      //! The tertiary-level predictor; uses lower level predictors, as needed
      Prediction< t_DataType> m_TertiaryPredictor;

      //! storage for the models used earlier in the training process
      util::Implementation< model::RetrieveInterface> m_InitialModelsStorage;

      //! Total iterations to perform
      size_t m_Iterations;

      //! Matrix of return for the current object
      linal::Matrix< float> m_LastRoundsPredictions;

      //! prediction label for the tertiary-level's models
      util::ObjectDataLabel m_TertiaryPredictionLabel;

      //! prediction mean label for the tertiary-level's models
      util::ObjectDataLabel m_TertiaryPredictionMeanLabel;

      //! prediction label for the lower-level's models
      util::ObjectDataLabel m_LowLevelPredictionLabel;

      //! prediction label for the lower-level's models
      util::ObjectDataLabel m_LowLevelPredictionMeanLabel;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;
      static const util::SiPtr< const util::ObjectInterface> s_MeanInstance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor, accepts bool of whether to auto-compute mean
      IterativePrediction( const bool &MEAN = false);

      //! @brief virtual copy constructor
      IterativePrediction *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

    protected:

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
      //! implementations
      iterate::Generic< Base< t_DataType, float> > GetInternalDescriptors();

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      void SetObjectHook();

      //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
      //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest; will be for the native dimension of the class
      //! @param STORAGE storage for the descriptor
      //! The inheriting class can assume that ITR is of native dimension and that STORAGE is zeroed and of the proper
      //! dimension
      void RecalculateImpl
      (
        const Iterator< t_DataType> &ITR,
        linal::VectorReference< float> &STORAGE
      );

      //! @brief Iteratively perform the prediction
      //! @param ITR iterator to use
      void Iterate( const Iterator< t_DataType> &ITR);

    }; // class IterativePrediction

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API IterativePrediction< chemistry::AtomConformationalInterface>;
    BCL_EXPIMP_TEMPLATE template class BCL_API IterativePrediction< biol::AABase>;
    BCL_EXPIMP_TEMPLATE template class BCL_API IterativePrediction< biol::Mutation>;
    BCL_EXPIMP_TEMPLATE template class BCL_API IterativePrediction< char>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_ITERATIVE_PREDICTION_H_
