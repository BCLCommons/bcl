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

#ifndef BCL_MODEL_OBJECTIVE_FUNCTION_RMSD_H_
#define BCL_MODEL_OBJECTIVE_FUNCTION_RMSD_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_feature_data_set.h"
#include "bcl_model_objective_function_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ObjectiveFunctionRmsd
    //! @brief The objective function for training a NeuralNetwork
    //! @details This class determines the prediction error of the NeuralNetwork model on the monitoring data set.
    //! The NeuralNetwork improves as long as the error goes down.
    //!
    //! @see @link example_model_objective_function_rmsd.cpp @endlink
    //! @author mueller
    //! @date Jul 6, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ObjectiveFunctionRmsd :
      public ObjectiveFunctionInterface
    {

    public:

      //! Different types of normalization often used for compute RMSD
      enum Normalization
      {
        e_None,   //!< no normalization applied
        e_MinMax, //!< normalize each output by the max-min seen in that column
        e_Std,    //!< Normalize by the standard deviation seen in the column
        s_NumberNormalization
      };

    private:

    //////////
    // data //
    //////////

      //! Normalization applied to rmsd
      Normalization m_Normalization;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_RawInstance;
      static const util::SiPtr< const util::ObjectInterface> s_MinMaxInstance;
      static const util::SiPtr< const util::ObjectInterface> s_StdInstance;

    public:

    //////////
    // data //
    //////////

      //! @brief Clone function
      //! @return pointer to new ObjectiveFunctionRmsd
      ObjectiveFunctionRmsd *Clone() const;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor from type of normalization
      ObjectiveFunctionRmsd( const Normalization &NORM = e_None);

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief determine what sign of the derivative of this objective function indicates improvement
      //! @return the sign of the derivative of this objective function indicates improvement
      opti::ImprovementType GetImprovementType() const
      {
        return opti::e_SmallerEqualIsBetter;
      }

      //! @brief get the overall goal of the objective function
      //! @return the goal of the objective function
      Goal GetGoalType() const
      {
        return e_Regression;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief evaluate for a given datasets with experimental and predicted values the implemented objective function
      //! @param EXPERIMENTAL feature dataset with experimental values
      //! @param PREDICTED feature dataset with predicted values
      //! @return objective function value based on the given data of experimental and predicted values
      float operator()
      (
        const FeatureDataSetInterface< float> &EXPERIMENTAL,
        const FeatureDataSetInterface< float> &PREDICTED
      ) const;

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class ObjectiveFunctionRmsd

  } // namespace model

} // namespace bcl

#endif // BCL_MODEL_OBJECTIVE_FUNCTION_RMSD_H_
