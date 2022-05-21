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

#ifndef BCL_MC_OPTIMIZATION_MCM_H_
#define BCL_MC_OPTIMIZATION_MCM_H_

// include the namespace header
#include "bcl_mc.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_mc_metropolis.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_serializer.h"
#include "math/bcl_math_mutate_interface.h"
#include "opti/bcl_opti_criterion_interface.h"
#include "opti/bcl_opti_optimization_interface.h"
#include "score/bcl_score_protein_model_score_sum.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class OptimizationMCM
    //! @brief Performs an MCM optimization of a given argument.
    //! @detail An MCM approximation with the provided mutates, scoring function, and termination criteria is
    //! performed on the given argument.
    //!
    //! @see @link example_mc_optimization_mcm.cpp @endlink
    //! @author fischea
    //! @date Aug 08, 2016
    //!
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API OptimizationMCM :
      public opti::OptimizationInterface< assemble::ProteinModel>
    {

    //////////
    // data //
    //////////

    public:

      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    private:

      //! scoring function used to evaluate the sampled structures
      score::ProteinModelScoreSum m_ScoreFunction;

      //! mutates used to sample structures
      util::Implementation< math::MutateInterface< assemble::ProteinModel> > m_Mutates;

      //! termination criterion
      util::Implementation< opti::CriterionInterface< assemble::ProteinModel, double> > m_Criterion;

      //! Metropolis criterion
      Metropolis< double> m_Metropolis;

    ////////////////////////////////////
    //  construction and destruction  //
    ////////////////////////////////////

    public:

      //! @brief default constructor
      OptimizationMCM();

      //! @brief construct from members
      //! @param SCORE_FUNCTION scoring function to evaluate the sampled models
      //! @param MUTATES mutates to sample models
      //! @param CRITERION termination criterion
      //! @param METROPOLIS metropolis criterion to decide whether mutates are accepted
      OptimizationMCM
      (
        const score::ProteinModelScoreSum &SCORE_FUNCTION,
        const math::MutateInterface< assemble::ProteinModel> &MUTATES,
        const opti::CriterionInterface< assemble::ProteinModel, double> &CRITERION,
        const Metropolis< double> &METROPOLIS
      );

      //! @brief clone function
      //! @return pointer to a new OptimizationMCM
      OptimizationMCM *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief optimizes the provided protein model
      //! @param MODEL protein model to be optimized
      void Optimize( assemble::ProteinModel &MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief returns the default temperature control for Metropolis
      //! @return default temperature control for Metropolis
      static util::ShPtr< TemperatureInterface> GetDefaultTemperature();

    }; // class OptimizationMCM

  } // namespace mc
} // namespace bcl

#endif // BCL_MC_OPTIMIZATION_MCM_H_
