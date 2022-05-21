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

// include forward header of this class
#include "model/bcl_model_collect_features_top.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> CollectFeaturesTop::s_Instance
    (
      util::Enumerated< CollectFeaturesInterface>::AddInstance( new CollectFeaturesTop( util::GetUndefined< size_t>()))
    );

    //! @brief constructor from parameter
    //! @param NUMBER_TO_KEEP number of features to keep
    CollectFeaturesTop::CollectFeaturesTop( const size_t &NUMBER_TO_KEEP) :
      m_NumberToKeep( NUMBER_TO_KEEP)
    {
    }

    //! @brief Clone function
    //! @return pointer to new CollectFeaturesTop
    CollectFeaturesTop *CollectFeaturesTop::Clone() const
    {
      return new CollectFeaturesTop( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &CollectFeaturesTop::GetClassIdentifier() const
    {
      return GetStaticClassName< CollectFeaturesTop>();
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &CollectFeaturesTop::GetAlias() const
    {
      static const std::string s_Name( "Top");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief collect the top m_NumberToKeep value indices by score
    //! @param SCORES any values
    //! @return top m_NumberToKeep value indices
    storage::Vector< size_t> CollectFeaturesTop::Collect( const linal::Vector< float> &SCORES) const
    {
      // check if it is desired to keep more scores than were even given
      if( m_NumberToKeep >= SCORES.GetSize())
      {
        // keep all the scores
        storage::Vector< size_t> to_keep( SCORES.GetSize());
        for( size_t i( 0), i_max( to_keep.GetSize()); i < i_max; ++i)
        {
          to_keep( i) = i;
        }
        return to_keep;
      }

      // initialize a vector to hold the scores and their indices
      storage::Vector< storage::Pair< float, size_t> > score_and_index;

      // allocate enough memory to hold all the columns
      score_and_index.AllocateMemory( SCORES.GetSize());

      // add all the scores and column indices to the vector
      for( size_t feature_index( 0), feature_size( SCORES.GetSize()); feature_index < feature_size; ++feature_index)
      {
        score_and_index.PushBack( storage::Pair< float, size_t>( SCORES( feature_index), feature_index));
      }

      // sort the vector by scores
      score_and_index.Sort( std::less< storage::Pair< float, size_t> >());

      // make a vector that will hold the top m_NumberToKeep indices
      storage::Vector< size_t> chosen_scores;

      // allocate enough memory to hold the columns that will be added
      chosen_scores.AllocateMemory( m_NumberToKeep);

      size_t best_score_index( 0);
      // add the best (e.g. the last) N_HIGHEST columns to the vector
      for
      (
        storage::Vector< storage::Pair< float, size_t> >::const_reverse_iterator
          rev_itr_score( score_and_index.ReverseBegin()),
          rev_itr_score_end( score_and_index.ReverseEnd());
        best_score_index < m_NumberToKeep;
        ++rev_itr_score, ++best_score_index
      )
      {
        chosen_scores.PushBack( rev_itr_score->Second());
      }

      // sort the vector of scores by index
      chosen_scores.Sort( std::less< size_t>());

      return chosen_scores;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer CollectFeaturesTop::GetSerializer() const
    {
      io::Serializer member_data;
      member_data.SetClassDescription( "Selects the indices of the highest N values");

      member_data.AddInitializer
      (
        "",
        "number of highest valued indices to collect",
        io::Serialization::GetAgentWithMin( &m_NumberToKeep, size_t( 1)),
        "1"
      );

      return member_data;
    }

  } // namespace model
} // namespace bcl
