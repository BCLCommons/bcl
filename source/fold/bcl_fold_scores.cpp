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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "fold/bcl_fold_scores.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_membrane.h"
#include "fold/bcl_fold_default_flags.h"
#include "math/bcl_math_binary_function_cached.h"
#include "math/bcl_math_function_cached.h"
#include "score/bcl_score_protein_model_wrapper.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Scores::Scores()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Scores::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief adds a score to the enumerated scores
    //! @param SCORE Score to add
    Scores::EnumType Scores::AddScore( const util::ShPtr< score::ProteinModel> &SCORE)
    {
      return AddEnum( SCORE->GetScheme(), SCORE);
    }

    //! @brief adds a vector of scores to the enumerated scores
    //! @param SCORES Scores to add
    storage::Vector< Scores::EnumType>
    Scores::AddScoreVector( const util::ShPtrVector< score::ProteinModel> &SCORES)
    {
      storage::Vector< Scores::EnumType> score_enums;

      // iterate through the scores to add them
      for
      (
        util::ShPtrVector< score::ProteinModel>::const_iterator
          score_itr( SCORES.Begin()), score_itr_end( SCORES.End());
        score_itr != score_itr_end;
        ++score_itr
      )
      {
        score_enums.PushBack( AddEnum( ( *score_itr)->GetScheme(), *score_itr));
      }

      return score_enums;
    }

    //! @brief function for adding a new enum
    //! @param NAME name of the current enum
    //! @param OBJECT object to be enumerated
    Scores::EnumType &Scores::AddEnum
    (
      const std::string &NAME,
      const util::ShPtr< score::ProteinModel> &OBJECT
    )
    {
      // make sure an enum with the given name does not exist already
      if( HaveEnumWithName( NAME))
      {
        BCL_MessageTop( "A score enum with the given name already exists: " + NAME);
      }

      // call the add enum
      return util::Enumerate< util::ShPtr< score::ProteinModel>, Scores>::AddEnum( NAME, OBJECT);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief wrap a given unary function for ProteinModel into a FunctionCached object
    //! @param SP_FUNCTION ShPtr to function be wrapped
    //! @return a FunctionCached object
    util::ShPtr< score::ProteinModel>
    Scores::WrapCacheProteinModelScore
    (
      const util::ShPtr< score::ProteinModel> &SP_FUNCTION
    )
    {
      // wrap function into cache function
      util::ShPtr< math::FunctionCached< assemble::ProteinModel, double> > sp_cache_function
      (
        new math::FunctionCached< assemble::ProteinModel, double>
        (
          SP_FUNCTION,
          &assemble::ProteinModel::GetDestructorSignal,
          &assemble::ProteinModel::GetChangeSignal
        )
      );

      // create score
      util::ShPtr< score::ProteinModel> sp_score
      (
        new score::ProteinModelWrapper( sp_cache_function, SP_FUNCTION->GetType(), SP_FUNCTION->GetReadableScheme())
      );

      // end
      return sp_score;
    }

    //! @brief wrap a given unary function for sse into a FunctionCached object
    //! @param SP_FUNCTION ShPtr to function be wrapped
    //! @return a FunctionCached object
    util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> > >
    Scores::WrapCacheSSEScore
    (
      const util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> > > &SP_FUNCTION
    )
    {
      // wrap function into cache function
      util::ShPtr< math::BinaryFunctionCached< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> > > sp_cache_function
      (
        new math::BinaryFunctionCached< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> >
        (
          SP_FUNCTION,
          &assemble::SSE::GetDestructorSignal,
          &biol::Membrane::GetDestructorSignal
        )
      );

      // add signal handler for coordinate changes
      sp_cache_function->AddSignalHandlerForArgument( &assemble::SSE::GetCoordinateChangeSignal);

      // end
      return sp_cache_function;
    }

    //! @brief wrap a given binary function for two sse into a BinaryFunctionCached object
    //! @param SP_FUNCTION ShPtr to binary function be wrapped
    //! @param SYMMETRIC is function symmetric (a,b is the same as b,a)
    //! @return a FunctionCached object
    util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, assemble::SSE, double> >
    Scores::WrapCacheSSEPairScore
    (
      const util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, assemble::SSE, double> > &SP_FUNCTION,
      const bool SYMMETRIC
    )
    {
      // wrap function into cache function
      util::ShPtr< math::BinaryFunctionCached< assemble::SSE, assemble::SSE, double> > sp_cache_function
      (
        new math::BinaryFunctionCached< assemble::SSE, assemble::SSE, double>
        (
          SP_FUNCTION,
          &assemble::SSE::GetDestructorSignal,
          SYMMETRIC
        )
      );

      // add signal handler for coordinate changes
      sp_cache_function->AddSignalHandlerForArgument( &assemble::SSE::GetCoordinateChangeSignal);

      // end
      return sp_cache_function;
    }

    //! @brief construct on access function for all Scores
    //! @return reference to only instances of Scores
    Scores &GetScores()
    {
      return Scores::GetEnums();
    }

  } // namespace fold

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< score::ProteinModel>, fold::Scores>;

  } // namespace util
} // namespace bcl
