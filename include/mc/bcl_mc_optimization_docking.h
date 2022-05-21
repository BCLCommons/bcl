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

#ifndef BCL_MC_OPTIMIZATION_DOCKING_H_
#define BCL_MC_OPTIMIZATION_DOCKING_H_

// include header of this class

// includes from bcl - sorted alphabetically
#include "bcl_mc_optimization_mcm.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_serializer.h"
#include "math/bcl_math_mutate_interface.h"
#include "opti/bcl_opti_optimization_interface.h"
#include "quality/bcl_quality_measures.h"
#include "score/bcl_score_protein_model_score_sum.h"
#include "storage/bcl_storage_set.h"
#include "util/bcl_util_si_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
    namespace mc
    {

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class OptimizationDocking
      //! @brief Docking monte-carlo algorithm for proteins
      //!
      //! @see @link example_mc_optimization_docking.cpp @endlink
      //! @author lib14
      //! @date Nov 15, 2017
      //!
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      class BCL_API OptimizationDocking :
        public util::SerializableInterface
      {
    //////////
    // data //
    //////////

    public:

      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_instance;

    private:

      //! scoring function used to evaluate the sampled structures
      score::ProteinModelScoreSum m_ScoreFunction;

      //! mutates used to sample structures
      util::Implementation< math::MutateInterface< assemble::ProteinModel> > m_Mutates;

      //! termination criterion
      util::Implementation< opti::CriterionInterface< assemble::ProteinModel, double> > m_Criterion;

      //! Metropolis criterion
      Metropolis< double> m_Metropolis;

      //! quality measures used to evaluate how close the optimized model is to the native structure
      storage::Set< quality::Measure> m_QualityMeasures;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief default constructor
      //!@return an object of OptimizationMCMDocking
      OptimizationDocking();

      //! @brief constructor that takes in initializers
      //! @param SCORE_FUNCTION
      //! @param MUTATES
      //! @param CRITERION
      //! @param METROPOLIS
      //! @return
      OptimizationDocking
      (
        const score::ProteinModelScoreSum &SCORE_FUNCTION,
        const math::MutateInterface< assemble::ProteinModel> &MUTATES,
        const opti::CriterionInterface< assemble::ProteinModel, double> &CRITERION,
        const Metropolis< double> &METROPOLIS,
        const storage::Set< quality::Measure> &MEASURES
      );

      //! @brief clone function
      //! @return a pointer to a new Optimization object
      OptimizationDocking *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief gets the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

      //! @brief gets the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief get the scoring function
      //! @return the scoring function
      const score::ProteinModelScoreSum &GetScoreFunction() const;

      //! @brief get quality measures
      //! @return quality measures
      const storage::Set< quality::Measure> &GetQualityMeasures() const;

      //! @brief gets parameters for data members that are set up from data labels
      //!@return parameters for data members that are set up from data labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief optimizes the provided protein model
      //! @param MODEL protein model to be optimized
      void Optimize( assemble::ProteinModel &MODEL) const;

      }; // end of class OptimizationDocking
  } // namespace mc
} // namespace bcl

#endif // BCL_MC_OPTIMIZATION_DOCKING_H_
