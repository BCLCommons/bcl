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
#include "mc/bcl_mc_optimization_mcm.h"

// includes from bcl - sorted alphabetically
#include "mc/bcl_mc_approximator.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> OptimizationMCM::s_Instance
    (
      util::Enumerated< opti::OptimizationInterface< assemble::ProteinModel> >::AddInstance( new OptimizationMCM())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    OptimizationMCM::OptimizationMCM() :
      m_ScoreFunction(),
      m_Mutates(),
      m_Criterion(),
      m_Metropolis()
    {
    }

    //! @param SCORE_FUNCTION scoring function to evaluate the sampled models
    //! @param MUTATES mutates to sample models
    //! @param CRITERION termination criterion
    //! @param METROPOLIS metropolis criterion to decide whether mutates are accepted
    OptimizationMCM::OptimizationMCM
    (
      const score::ProteinModelScoreSum &SCORE_FUNCTION,
      const math::MutateInterface< assemble::ProteinModel> &MUTATES,
      const opti::CriterionInterface< assemble::ProteinModel, double> &CRITERION,
      const Metropolis< double> &METROPOLIS
    ) :
      m_ScoreFunction( SCORE_FUNCTION),
      m_Mutates( MUTATES),
      m_Criterion( CRITERION),
      m_Metropolis( METROPOLIS)
    {
    }

    //! @brief copy constructor
    //! @return pointer to a new OptimizationMCM
    OptimizationMCM *OptimizationMCM::Clone() const
    {
      return new OptimizationMCM( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &OptimizationMCM::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &OptimizationMCM::GetAlias() const
    {
      static const std::string s_alias( "MCMOptimizer");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer OptimizationMCM::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Optimization implementation for Monte Carlo Metropolis algorithms.");
      // serializer.Merge( opti::OptimizationInterface< assemble::ProteinModel>::GetSerializer());
      serializer.AddInitializer
      (
        "score function",
        "score function to evaluate the sampled protein models",
        io::Serialization::GetAgent( &m_ScoreFunction)
      );
      serializer.AddInitializer
      (
        "mutates",
        "mutates to sample the protein models",
        io::Serialization::GetAgent( &m_Mutates)
      );
      serializer.AddInitializer
      (
        "termination criterion",
        "criterion when the optimization will be terminated",
        io::Serialization::GetAgent( &m_Criterion)
      );
      serializer.AddInitializer
      (
        "metropolis",
        "metropolis criterion to decide which mutates are accepted",
        io::Serialization::GetAgent( &m_Metropolis)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief optimizes the provided protein model
    //! @param MODEL protein model to be optimized
    void OptimizationMCM::Optimize( assemble::ProteinModel &MODEL) const
    {
      // create the approximator
      Approximator< assemble::ProteinModel, double> approximator
      (
        m_ScoreFunction, *m_Mutates, m_Metropolis, *m_Criterion, MODEL
      );

      // optimize the given protein model
      approximator.Approximate();

      // get the result of the optimization
      MODEL = *approximator.GetTracker().GetBest()->First().HardCopy();
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief returns the default temperature control for Metropolis
    //! @return default temperature control for Metropolis
    util::ShPtr< TemperatureInterface> OptimizationMCM::GetDefaultTemperature()
    {
      return util::ShPtr< TemperatureInterface>( new TemperatureDefault);
    }

  } // namespace mc
} // namespace bcl
