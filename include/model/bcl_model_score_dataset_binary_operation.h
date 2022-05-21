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

#ifndef BCL_MODEL_SCORE_DATASET_BINARY_OPERATION_H_
#define BCL_MODEL_SCORE_DATASET_BINARY_OPERATION_H_

// include the namespace header
#include "bcl_model.h"

// other forward includes

// includes from bcl - sorted alphabetically
#include "bcl_model_score_dataset_interface.h"
#include "math/bcl_math_assignment_operation_interface.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ScoreDatasetBinaryOperation
    //! @brief perform a math::Assignment on several dataset scores
    //!
    //! @see @link example_model_score_dataset_binary_operation.cpp @endlink
    //! @author mendenjl
    //! @date Apr 06, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ScoreDatasetBinaryOperation :
      public ScoreDatasetInterface
    {
    private:

    //////////
    // data //
    //////////

      //! 2+ scores to perform operation on
      storage::Vector< util::Implementation< ScoreDatasetInterface> > m_Scores;

      //! the operation that is performed to merge the properties (e.g. +=,-=,*=,/=)
      util::Implementation< math::AssignmentOperationInterface< float> > m_Op;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Add;          //!< add two scores together
      static const util::SiPtr< const util::ObjectInterface> s_Subtract;     //!< subtract two scores
      static const util::SiPtr< const util::ObjectInterface> s_Divide;       //!< divide scores
      static const util::SiPtr< const util::ObjectInterface> s_Multiply;     //!< multiply scores
      static const util::SiPtr< const util::ObjectInterface> s_Exponentiate; //!< exponentiate a score

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from an operation and an alias
      explicit ScoreDatasetBinaryOperation( const math::AssignmentOperationInterface< float> &OPERATION);

      //! @brief virtual copy constructor
      ScoreDatasetBinaryOperation *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief score a given dataset
      //! @param DATASET dataset of interest
      //! @return scores of the dataset
      linal::Vector< float> Score( const descriptor::Dataset &DATASET) const;

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class ScoreDatasetBinaryOperation

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_SCORE_DATASET_BINARY_OPERATION_H_
