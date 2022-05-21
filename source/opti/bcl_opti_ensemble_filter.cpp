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
#include "opti/bcl_opti_ensemble_filter.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> EnsembleFilter::s_Instance
    (
      util::Enumerated< OptimizationInterface< assemble::Ensemble< assemble::ProteinModel> > >::AddInstance( new EnsembleFilter())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    EnsembleFilter::EnsembleFilter() :
      m_ScoreFunction(),
      m_KeepPercentage()
    {
    }

    //! @brief construct from members
    //! @param SCORE_FUNCTION scoring function to evaluate the sampled models
    //! @param KEEP percentage of model to keep
    EnsembleFilter::EnsembleFilter( const score::ProteinModel &SCORE_FUNCTION, double KEEP) :
      m_ScoreFunction( SCORE_FUNCTION),
      m_KeepPercentage( KEEP)
    {
    }

    //! @brief copy constructor
    //! @return pointer to a new EnsembleFilter
    EnsembleFilter *EnsembleFilter::Clone() const
    {
      return new EnsembleFilter( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &EnsembleFilter::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &EnsembleFilter::GetAlias() const
    {
      static const std::string s_alias( "EnsembleFilter");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer EnsembleFilter::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Filters an ensemble using a provided scoring function.");
      serializer.AddInitializer
      (
        "score function",
        "score function to evaluate the sampled protein models",
        io::Serialization::GetAgent( &m_ScoreFunction)
      );
      serializer.AddInitializer
      (
        "keep",
        "percentage of models to keep",
        io::Serialization::GetAgent( &m_KeepPercentage)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief filters the provided ensemble
    //! @param ENSEMBLE ensemble to be filtered
    void EnsembleFilter::Optimize( assemble::Ensemble< assemble::ProteinModel> &ENSEMBLE) const
    {
      // compute the score for each element of the ensemble
      typedef assemble::Ensemble< assemble::ProteinModel>::Element Element;
      storage::Vector< storage::Pair< double, Element> > elements;
      storage::Vector< double> scores;
      for( auto el_it( ENSEMBLE.Begin()), el_it_end( ENSEMBLE.End()); el_it != el_it_end; ++el_it)
      {
        const assemble::ProteinModel &model( el_it->GetElement());
        const double score( ( *m_ScoreFunction)( model));
        const storage::Pair< double, Element> element_score( score, *el_it);
        elements.PushBack( element_score);
        scores.PushBack( score);
      }

      // sort the scores and determine the cutoff
      scores.Sort( std::less< double>());
      const size_t number_models( ENSEMBLE.GetSize());
      const size_t number_keep( m_KeepPercentage * number_models);
      const double cutoff( scores( number_keep - 1));

      // only keep elements with a score below the cutoff value
      assemble::Ensemble< assemble::ProteinModel> ensemble;
      for( auto el_it( elements.Begin()), el_it_end( elements.End()); el_it != el_it_end; ++el_it)
      {
        const double score( el_it->First());
        if( score <= cutoff)
        {
          ensemble.AddElement( el_it->Second());
        }
      }
      ENSEMBLE = ensemble;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace opti
} // namespace bcl
