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

#ifndef BCL_MODEL_APPROXIMATOR_LEVERAGE_MATRIX_H_
#define BCL_MODEL_APPROXIMATOR_LEVERAGE_MATRIX_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_approximator_base.h"
#include "bcl_model_multiple_linear_regression.h"
#include "bcl_model_objective_function_wrapper.h"
#include "linal/bcl_linal_matrix_inversion_interface.h"
#include "math/bcl_math_range.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ApproximatorLeverageMatrix
    //! @brief trains a weight matrix
    //!
    //! @see @link example_model_approximator_leverage_matrix.cpp @endlink
    //! @author mendenjl
    //! @date May 06, 2016
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ApproximatorLeverageMatrix :
      public ApproximatorBase
    {

    //////////
    // data //
    //////////

    public:

      //! static instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    private:

      float                   m_LastObjectiveFunctionResult; //!< last objective function result
      util::ShPtr< Interface> m_LeverageMatrix;              //!< current model

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ApproximatorLeverageMatrix();

      //! @brief clone function
      //! @return pointer to new ApproximatorLeverageMatrix
      ApproximatorLeverageMatrix *Clone() const
      {
        return new ApproximatorLeverageMatrix( *this);
      }

    /////////////////
    // data access //
    /////////////////

    public:

      //! @brief set training data set for a specific iterate in approximator framework
      //! @param DATA training data set
      void SetTrainingData( util::ShPtr< descriptor::Dataset> &DATA);

      //! @brief get objective function to evaluate a monitoring dataset
      //! @return dataset interface with training data
      const util::ShPtr< ObjectiveFunctionWrapper> &GetObjectiveFunction() const
      {
        return m_ObjectiveFunction;
      }

      //! @brief set objective function to evaluate a monitoring dataset
      //! @param OBJ objective function of interest
      void SetObjectiveFunction( const util::ShPtr< ObjectiveFunctionWrapper> &OBJ)
      {
        m_ObjectiveFunction = OBJ;
      }

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief construct a model from the current iterate
      //! @return shptr to the new model interface
      util::ShPtr< Interface> GetCurrentModel() const;

      //! @brief returns the current approximation
      //! @return current argument result pair
      const util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> > GetCurrentApproximation() const;

      //! @brief conducts the next approximation step and stores the approximation
      void Next();

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief evaluates whether the approximation can continue
      //! @return true, if the approximation can continue - otherwise false
      bool CanContinue() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class ApproximatorLeverageMatrix

  } // namespace model

} // namespace bcl

#endif // BCL_MODEL_APPROXIMATOR_LEVERAGE_MATRIX_H_

