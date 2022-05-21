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

#ifndef BCL_FOLD_SCORES_H_
#define BCL_FOLD_SCORES_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "score/bcl_score_protein_model.h"
#include "score/bcl_score_protein_model_score_sum.h"
#include "util/bcl_util_enumerate.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Scores
    //! @brief enumerator for all scores used in folding
    //! @details util::Enumerate derived class that allows enumeration of all ProteinModel scoring functions used in
    //! different folding protocols
    //!
    //! @remarks example unnecessary
    //! @author karakam
    //! @date Mar 27, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Scores :
      public util::Enumerate< util::ShPtr< score::ProteinModel>, Scores>
    {
      friend class util::Enumerate< util::ShPtr< score::ProteinModel>, Scores>;

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Scores();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief adds a score to the enumerated scores
      //! @param SCORE Score to add
      EnumType AddScore( const util::ShPtr< score::ProteinModel> &SCORE);

      //! @brief adds a vector of scores to the enumerated scores
      //! @param SCORES Scores to add
      storage::Vector< Scores::EnumType> AddScoreVector( const util::ShPtrVector< score::ProteinModel> &SCORES);

    private:

      //! @brief function for adding a new enum
      //! @param NAME name of the current enum
      //! @param OBJECT object to be enumerated
      EnumType &AddEnum
      (
        const std::string &NAME,
        const util::ShPtr< score::ProteinModel> &OBJECT
      );

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief wrap a given unary function for ProteinModel into a FunctionCached object
      //! @param SP_FUNCTION ShPtr to function be wrapped
      //! @return a FunctionCached object
      static util::ShPtr< score::ProteinModel>
      WrapCacheProteinModelScore
      (
        const util::ShPtr< score::ProteinModel> &SP_FUNCTION
      );

      //! @brief wrap a given unary function for sse into a FunctionCached object
      //! @param SP_FUNCTION ShPtr to function be wrapped
      //! @return a FunctionCached object
      static util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> > >
      WrapCacheSSEScore
      (
        const util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> > > &SP_FUNCTION
      );

      //! @brief wrap a given binary function for two sse into a BinaryFunctionCached object
      //! @param SP_FUNCTION ShPtr to binary function be wrapped
      //! @param SYMMETRIC is function symmetric (a,b is the same as b,a)
      //! @return a FunctionCached object
      static util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, assemble::SSE, double> >
      WrapCacheSSEPairScore
      (
        const util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, assemble::SSE, double> > &SP_FUNCTION,
        const bool SYMMETRIC
      );

    }; // class Scores

    //! @brief construct on access function for all Scores
    //! @return reference to only instances of Scores
    BCL_API
    Scores &GetScores();

  } // namespace fold

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< score::ProteinModel>, fold::Scores>;

  } // namespace util
} // namespace bcl

#endif // BCL_FOLD_SCORES_H_
