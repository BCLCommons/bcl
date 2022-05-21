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

#include "biol/bcl_biol_atom.h"
#include "biol/bcl_biol_membrane.h"
#include "fold/bcl_fold_mutate_membrane_chain_move.h"
#include "fold/bcl_fold_mutate_protein_model_chain_move.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ofstream.h"
#include "mc/bcl_mc_approximator.h"
#include "mc/bcl_mc_optimization_docking.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_printer_quality_docking.h"
#include "pdb/bcl_pdb_printer_score.h"
#include "quality/bcl_quality_rmsd.h"
namespace bcl
{
    namespace mc
    {
  //////////
  // data //
  //////////

        //! single instance of this class
        const util::SiPtr< const util::ObjectInterface> OptimizationDocking::s_instance
        (
          GetObjectInstances().AddInstance( new OptimizationDocking())
        );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

        //! @brief default constructor
        OptimizationDocking::OptimizationDocking()
        {
        }

        //! @brief constructor that takes in initializers
        //! @param SCORE_FUNCTION scoring function for ranking docked models
        //! @param MUTATES Monte Carlo moves for sampling conformational space
        //! @param CRITERION termination criterion
        //! @param METROPOLIS Metropolis parameters
        //! @param MEASURES what measures to use to evaluate quality of models
        OptimizationDocking::OptimizationDocking
        (
          const score::ProteinModelScoreSum &SCORE_FUNCTION,
          const math::MutateInterface< assemble::ProteinModel> &MUTATES,
          const opti::CriterionInterface< assemble::ProteinModel, double> &CRITERION,
          const Metropolis< double> &METROPOLIS,
          const storage::Set< quality::Measure> &MEASURES
        ) :
          m_ScoreFunction( SCORE_FUNCTION),
          m_Mutates( MUTATES),
          m_Criterion( CRITERION),
          m_Metropolis( METROPOLIS),
          m_QualityMeasures( MEASURES)
        {
        }

        //! @brief clone function
        //! @return a pointer to a new OptimizationMCMDocking object
        OptimizationDocking *OptimizationDocking::Clone() const
        {
            return new OptimizationDocking( *this);
        }

  /////////////////
  // data access //
  /////////////////

        //! @brief gets the name of this class
        //! @return the name of this class
        const std::string &OptimizationDocking::GetClassIdentifier() const
        {
            return GetStaticClassName( *this);
        }

        //! @brief gets the name of the object when used in a dynamic context
        //! @return the name of the object when used in a dynamic context
        const std::string &OptimizationDocking::GetAlias() const
        {
            static const std::string s_alias( "DockingOptimizer");
            return s_alias;
        }

        //! @brief get the scoring function
        //! @return the scoring function
        const score::ProteinModelScoreSum &OptimizationDocking::GetScoreFunction() const
        {
          return m_ScoreFunction;
        }

        //! @brief get quality measures
        //! @return quality measures
        const storage::Set< quality::Measure> &OptimizationDocking::GetQualityMeasures() const
        {
          return m_QualityMeasures;
        }

        //! @brief gets parameters for data members that are set up from data labels
        //! @return parameters for data members that are set up from data labels
        io::Serializer OptimizationDocking::GetSerializer() const
        {
            // gets the serializer from the base class
            io::Serializer serializer;

            // sets class description
            serializer.SetClassDescription( "Monte Carlo sampling for protein-protein docking");

            // add initializers
            serializer.AddInitializer
            (
              "score function",
              "scoring function to be used along with Monte Carlo sampling",
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
            serializer.AddInitializer
            (
              "quality measures",
              "measures for evaluating how close the model is to the native structure",
              io::Serialization::GetAgent( &m_QualityMeasures),
              "(RMSD)"
            );

            return serializer;
        }

  ////////////////
  // operations //
  ////////////////

        //! @brief optimizes the provided protein model
        //! @param MODEL protein model to be optimized
        void OptimizationDocking::Optimize( assemble::ProteinModel &MODEL) const
        {
          // create the approximator
          Approximator< assemble::ProteinModel, double> approximator
          (
            m_ScoreFunction, *m_Mutates, m_Metropolis, *m_Criterion, MODEL
          );

          // optimize the given protein model
          approximator.Approximate();

          // get the result of the optimization
          MODEL = *util::ShPtr< assemble::ProteinModel>( approximator.GetTracker().GetBest()->First().HardCopy());
        }
  } // namespace mc
} // namespace bcl
