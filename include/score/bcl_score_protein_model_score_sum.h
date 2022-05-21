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

#ifndef BCL_SCORE_PROTEIN_MODEL_SCORE_SUM_H_
#define BCL_SCORE_PROTEIN_MODEL_SCORE_SUM_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "fold/bcl_fold.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_protein_model.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "math/bcl_math_sum_function.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinModelScoreSum
    //! @brief This is a class for collecting scoring functions for protein models
    //! @details This class is a specialized sum function class for scoring protein models. In addition to SumFunction
    //! functionalities, it allows initialization from score weight sets
    //!
    //! @see @link example_score_protein_model_score_sum.cpp @endlink
    //! @author woetzen, karakam
    //! @date 05.10.2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinModelScoreSum :
      public math::SumFunctionMixin< ProteinModel>
    {
    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ProteinModelScoreSum( const std::string &SCHEME = "");

      //! @brief constructor from a map of functions and weights
      //! @param SCORE_WEIGHT_MAP map of scores and corresponding weights
      ProteinModelScoreSum
      (
        const storage::Map< util::ShPtr< ProteinModel>, double> &SCORE_WEIGHT_MAP
      );

      //! @brief constructor from a map of Score enums and weights
      //! @param SCORE_WEIGHT_MAP map of scores and corresponding weights
      ProteinModelScoreSum
      (
        const storage::Map< fold::Score, double> &SCORE_WEIGHT_MAP
      );

      //! @brief constructor from a map of functions and from a weight set
      //! @param SCORE_MAP map of scoring functions to be used
      //! @param WEIGHT_SET map of function schemes and corresponding weights
      ProteinModelScoreSum
      (
        const storage::Map< std::string, util::ShPtr< ProteinModel> > &SCORE_MAP,
        const storage::Map< std::string, double> &WEIGHT_SET
      );

      //! @brief virtual copy constructor
      ProteinModelScoreSum *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns a vector of string that has readable schemes of the individual functions
      //! @return the vector of strings that has readable schemes of the individual functions
      storage::Vector< std::string> GetReadableFunctionSchemes() const;

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief creates a table with individual function values, weights and weighted values
      //! the table has scores as the columns
      //! @param PROTEIN_MODEL ProteinModel to be used for calculating the function value
      //! @return a table with individual functions and their weighted sums
      storage::Table< double> CreateValueTableHorizontal( const assemble::ProteinModel &PROTEIN_MODEL) const;

      //! @brief creates a table with individual function values and their weighted sums for the given argument
      //! the table has scores as rows
      //! @param PROTEIN_MODEL ProteinModel to be used for calculating the function value
      //! @return a table with individual functions and their weighted sums
      storage::Table< double> CreateValueTableVertical( const assemble::ProteinModel &PROTEIN_MODEL) const;

      //! @brief creates a table with individual function values and their weighted sums for the given argument
      //! the table has scores as rows, using readable names instead of scheme
      //! @param PROTEIN_MODEL ProteinModel to be used for calculating the function value
      //! @return a table with individual functions and their weighted sums
      storage::Table< double> CreateSortedReadableTable( const assemble::ProteinModel &PROTEIN_MODEL) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief calculate the score for the given ProteinModel
      //! @brief PROTEIN_MODEL ProteinModel to be evaluated
      //! @return the score for the given model
      double operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write detailed scheme and values to OSTREAM
      //! @param PROTEIN_MODEL ProteinModel to be used to evaluate the function
      //! @param OSTREAM std::ostream to be written to
      //! @return std::ostream which was written to
      std::ostream &WriteDetailedSchemeAndValues
      (
        const assemble::ProteinModel &PROTEIN_MODEL,
        std::ostream &OSTREAM
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief initialize the score sum using the given score weight map
      //! @param SCORE_WEIGHT_MAP map of scores and corresponding weights
      void AddScoresWithWeights( const storage::Map< util::ShPtr< ProteinModel>, double> &SCORE_WEIGHT_MAP);

      //! @brief initialize the score sum using the given score map and weight map
      //! @param SCORE_MAP map of scoring functions to be used
      //! @param WEIGHT_SET map of function schemes and corresponding weights
      void AddScoresWithWeights
      (
        const storage::Map< std::string, util::ShPtr< ProteinModel> > &SCORE_MAP,
        const storage::Map< std::string, double> &WEIGHT_SET
      );

      //! @brief combines scores with the same readable scheme
      //! @return scores with the same readable scheme
      storage::Map< ProteinModel::TypeEnum, storage::Map< std::string, double> > CombineSimilarScores() const;

    }; //class ProteinModelScoreSum

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_PROTEIN_MODEL_SCORE_SUM_H_
