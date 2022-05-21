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
#include "assemble/bcl_assemble_collector_protein_model_conformation_by_score.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_binary_function_stl_wrapper.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> CollectorProteinModelConformationByScore::s_Instance
    (
      util::Enumerated< find::CollectorInterface< util::SiPtrList< const ProteinModel>, ProteinModel> >::AddInstance
      (
        new CollectorProteinModelConformationByScore()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    CollectorProteinModelConformationByScore::CollectorProteinModelConformationByScore() :
      m_ScoreFunction(),
      m_NumberToCollect(),
      m_Sort( util::BinaryFunctionSTLWrapper< std::less< double> >()),
      m_ConsiderCurrentConformation()
    {
    }

    //! @brief constructor taking parameters
    //! @param SCORE the scoring function to be used
    //! @param NUM_TO_COLLECT the number of conformations to collect
    //! @param SORT the method for sorting scores to determine whether the best or worst scored models are collected
    //! @param CONSIDER_CURRENT if true the current conformation will be considered with the other conformations
    CollectorProteinModelConformationByScore::CollectorProteinModelConformationByScore
    (
      const util::ShPtr< math::FunctionInterfaceSerializable< ProteinModel, double> > &SCORE,
      const size_t NUM_TO_COLLECT,
      const util::ShPtr< util::BinaryFunctionInterfaceSerializable< double, double, bool> > &SORT,
      const bool CONSIDER_CURRENT
    ) :
      m_ScoreFunction( *SCORE),
      m_NumberToCollect( NUM_TO_COLLECT),
      m_Sort( *SORT),
      m_ConsiderCurrentConformation( CONSIDER_CURRENT)
    {
    }

    //! @brief constructor taking parameters
    //! @param NUM_TO_COLLECT the number of conformations to collect
    //! @param SORT the method for sorting scores to determine whether the best or worst scored models are collected
    //! @param CONSIDER_CURRENT if true the current conformation will be considered with the other conformations
    CollectorProteinModelConformationByScore::CollectorProteinModelConformationByScore
    (
      const size_t NUM_TO_COLLECT,
      const util::ShPtr< util::BinaryFunctionInterfaceSerializable< double, double, bool> > &SORT,
      const bool CONSIDER_CURRENT
    ) :
      m_ScoreFunction(),
      m_NumberToCollect( NUM_TO_COLLECT),
      m_Sort( *SORT),
      m_ConsiderCurrentConformation( CONSIDER_CURRENT)
    {
    }

    //! @brief Clone function
    //! @return pointer to new CollectorProteinModelConformationByScore
    CollectorProteinModelConformationByScore *CollectorProteinModelConformationByScore::Clone() const
    {
      return new CollectorProteinModelConformationByScore( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CollectorProteinModelConformationByScore::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &CollectorProteinModelConformationByScore::GetAlias() const
    {
      static const std::string s_alias( "CollectorProteinModelConformationByScore");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer CollectorProteinModelConformationByScore::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Optimization implementation for Monte Carlo Metropolis algorithms.");
      serializer.AddInitializer
      (
        "score function",
        "score function to evaluate the sampled protein models",
        io::Serialization::GetAgent( &m_ScoreFunction)
      );
      serializer.AddInitializer
      (
        "number",
        "number of conformations to collect",
        io::Serialization::GetAgent( &m_NumberToCollect)
      );
      serializer.AddInitializer
      (
        "sort method",
        "method for sorting scores",
        io::Serialization::GetAgent( &m_Sort)
      );
      serializer.AddInitializer
      (
        "consider current",
        "consider the current conformation",
        io::Serialization::GetAgent( &m_ConsiderCurrentConformation)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! Collect the t_ReturnType objects in t_ArgumentType
    //! @param PROTEIN_MODEL entity that contains a t_ReturnType
    //! @return returns Group of the collected t_ReturnType objects
    util::SiPtrList< const ProteinModel>
    CollectorProteinModelConformationByScore::Collect( const ProteinModel &PROTEIN_MODEL) const
    {
      // get the conformations in the protein model
      const ProteinEnsemble &conformation_ensemble( PROTEIN_MODEL.GetConformationalEnsemble());

      storage::List< storage::Pair< double, util::SiPtr< const ProteinModel> > > score_model;

      // iterate through the conformational ensemble to score the proteins and insert them into score list
      for
      (
        ProteinEnsemble::const_iterator
          model_itr( conformation_ensemble.Begin()), model_itr_end( conformation_ensemble.End());
        model_itr != model_itr_end;
        ++model_itr
      )
      {
        // get siptr to current model
        const util::SiPtr< const ProteinModel> current( **model_itr);

        // score the current model
        const double score( m_ScoreFunction->operator()( *current));

        // add the current model and its score into the score list
        score_model.PushBack( storage::Pair< double, util::SiPtr< const ProteinModel> >( score, current));
      }

      // true if the current conformation of the protein model is also being considered
      if( m_ConsiderCurrentConformation)
      {
        // get siptr to current model
        const util::SiPtr< const ProteinModel> current( PROTEIN_MODEL);

        // score the current model
        const double score( m_ScoreFunction->operator()( *current));

        // add the current model and its score into the score list
        score_model.PushBack( storage::Pair< double, util::SiPtr< const ProteinModel> >( score, current));
      }

      // sort the score list
      score_model.Sort
      (
        storage::PairBinaryPredicateFirst< double, util::SiPtr< const ProteinModel> >( *m_Sort)
      );

      // get the desired number of models starting at the beginning of the sorted list
      storage::List< storage::Pair< double, util::SiPtr< const ProteinModel> > > extracted_scores
      (
        score_model.ExtractElements( 0, m_NumberToCollect)
      );

      // to hold the models
      util::SiPtrList< const ProteinModel> desired_models;

      // iterate over the list of extracted scores to build up the list of desired models
      for
      (
        storage::List< storage::Pair< double, util::SiPtr< const ProteinModel> > >::const_iterator
          itr( extracted_scores.Begin()), itr_end( extracted_scores.End());
        itr != itr_end;
        ++itr
      )
      {
        desired_models.PushBack( itr->Second());
      }

      return desired_models;
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace assemble
} // namespace bcl
