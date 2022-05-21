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

#ifndef BCL_MC_STAGE_H_
#define BCL_MC_STAGE_H_

// include the namespace header
#include "bcl_mc.h"

// include other forward headers - sorted alphabetically
#include "pdb/bcl_pdb.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_mc_approximator.h"
#include "fold/bcl_fold_mutate_tree.h"
#include "fold/bcl_fold_protocols.h"
#include "fold/bcl_fold_score_weight_set.h"
#include "fold/bcl_fold_setup.h"
#include "fold/bcl_fold_stage_interface.h"
#include "storage/bcl_storage_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Stage
    //! @brief This class encapsulates a specific part of MonteCarlo/Metropolis based protein folding
    //! @details This class is designed to allow multi-stage fold minimizations. The stage class comes with its own
    //! mutates, scores, mutate tree, score weightset and start model. In a multi-stage minimization, the start
    //! model will be the best model from the previous stage, unless this is the first stage, then the start model
    //! will be the start model as specified by fold::Setup class. Each stage also knows the maximum number of
    //! iterations it has as well as the protocols that define its scores and mutates. It also has its own temperature
    //! control.
    //!
    //! @remarks example unnecessary
    //! @author karakam, pinojc, fischea
    //! @date Feb 14, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Stage :
      public fold::StageInterface
    {

    //////////
    // data //
    //////////

    private:

      //! name
      std::string m_Name;

      //! number of the stage
      size_t m_Number;

      //! total number of stages
      bool m_IsLastStage;

      //! protocols
      storage::List< fold::Protocol> m_FoldProtocols;

      //! score protocols
      storage::List< fold::Protocol> m_ScoreProtocols;

      //! mutate protocols
      storage::List< fold::Protocol> m_MutateProtocols;

      //! max number of iterations
      size_t m_MaxNumberIterations;

      //! max number of unimproved iterations
      size_t m_MaxNumberUnimprovedIterations;

      //! scoring function
      util::ShPtr< score::ProteinModelScoreSum> m_ScoreFunction;

      //! score weight set
      util::ShPtr< fold::ScoreWeightSet> m_ScoreWeightSet;

      //! mutate tree
      util::ShPtr< fold::MutateTree> m_MutateTree;

      //! mutate
      util::ShPtr< math::MutateInterface< assemble::ProteinModel> > m_Mutate;

      //! temperature
      util::ShPtr< TemperatureInterface> m_Temperature;

      //! number of score to drop each round
      size_t m_NumberOfScoresDropped;

      //! indicates whether to modify the start model
      bool m_ModifyStartModel;

      //! indicates whether to print the start model
      bool m_PrintStartModel;

      //! indicates whether to print models every iteration
      bool m_PrintIterationModels;

      //! indicates whether to print the end model
      bool m_PrintEndModel;

      //! indicates whether to print the tracker history
      bool m_PrintTrackerHistory;

      //! pool filename postfix for use with this stage
      std::string m_PoolPostfix;

      //! number of the current model
      size_t m_RoundNumber;

      //! prefix
      std::string m_Prefix;

      //! output path
      std::string m_Path;

      //! quality measures to be computed for the generated models
      storage::Set< quality::Measure> m_QualityMeasures;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Stage();

      //! @brief constructor from members
      //! @param FOLD_PROTOCOLS Fold protocols
      //! @param SCORE_PROTOCOLS Score protocols
      //! @param MUTATE_PROTOCOLS Mutate protocols
      Stage
      (
        const storage::List< fold::Protocol> &FOLD_PROTOCOLS,
        const storage::List< fold::Protocol> &SCORE_PROTOCOLS,
        const storage::List< fold::Protocol> &MUTATE_PROTOCOLS
      );

      //! @brief Clone function
      //! @return pointer to new Stage
      Stage *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns name
      //! @return name
      const std::string &GetName() const
      {
        return m_Name;
      }

      //! @brief sets the name
      //! @param NAME new name
      void SetName( const std::string &NAME)
      {
        m_Name = NAME;
      }

      //! @brief returns the list of fold protocols used
      //! @return the list of fold protocols used
      const storage::List< fold::Protocol> &GetFoldProtocols() const;

      //! @brief sets the list of fold protocols used
      //! @param FOLD_PROTOCOLS list of fold protocols
      void SetFoldProtocols( const storage::List< fold::Protocol> &FOLD_PROTOCOLS);

      //! @brief returns the list of score protocols used
      //! @return the list of score protocols used
      const storage::List< fold::Protocol> &GetScoreProtocols() const;

      //! @brief sets the list of score protocols used
      //! @param SCORE_PROTOCOLS list of score protocols
      void SetScoreProtocols( const storage::List< fold::Protocol> &SCORE_PROTOCOLS);

      //! @brief returns the list of mutate protocols used
      //! @return the list of mutate protocols used
      const storage::List< fold::Protocol> &GetMutateProtocols() const;

      //! @brief sets the list of mutate protocols used
      //! @param MUTATE_PROTOCOLS list of mutate protocols
      void SetMutateProtocols( const storage::List< fold::Protocol> &MUTATE_PROTOCOLS);

      //! @brief returns maximum number of iterations
      //! @return maximum number of iterations
      size_t GetMaxNumberIterations() const
      {
        return m_MaxNumberIterations;
      }

      //! @brief set maximum number of iterations
      //! @param MAX_NUMBER_ITERATIONS new maximum number of iterations
      void SetMaxNumberIterations( const size_t MAX_NUMBER_ITERATIONS)
      {
        m_MaxNumberIterations = MAX_NUMBER_ITERATIONS;
      }

      //! @brief returns maximum number of unimproved iterations
      //! @return maximum number of unimproved iterations
      size_t GetMaxNumberUnimprovedIterations() const
      {
        return m_MaxNumberUnimprovedIterations;
      }

      //! @brief set maximum number of unimproved iterations
      //! @param MAX_NUMBER_UNIMPROVED_ITERATIONS new maximum number of unimproved iterations
      void SetMaxNumberUnimprovedIterations( const size_t MAX_NUMBER_UNIMPROVED_ITERATIONS)
      {
        m_MaxNumberUnimprovedIterations = MAX_NUMBER_UNIMPROVED_ITERATIONS;
      }

      //! @brief returns score weightset
      //! @return scoring weightset
      const util::ShPtr< fold::ScoreWeightSet> &GetScoreWeightSet() const;

      //! @brief returns score weightset
      //! @return scoring weightset
      void SetScoreWeightSet( const util::ShPtr< fold::ScoreWeightSet> &SCORE_WEIGHT_SET);

      //! @brief returns score function
      //! @return score function
      const util::ShPtr< score::ProteinModelScoreSum> &GetScoreFunction() const;

      //! @brief returns score function
      //! @return score function
      void SetScoreFunction( const util::ShPtr< score::ProteinModelScoreSum> &SCORE_FUNCTION);

      //! @brief get mutate tree
      //! @return mutate tree
      const util::ShPtr< fold::MutateTree> &GetMutateTree() const;

      //! @brief get mutate tree
      //! @return mutate tree
      void SetMutateTree( const util::ShPtr< fold::MutateTree> &MUTATE_TREE);

      //! @brief get mutate functions for protein models
      //! @return mutate functions for protein models
      const util::ShPtr< math::MutateInterface< assemble::ProteinModel> > &GetMutate() const;

      //! @brief get mutate functions for protein models
      //! @return mutate functions for protein models
      void SetMutate( const util::ShPtr< math::MutateInterface< assemble::ProteinModel> > &MUTATE);

      //! @brief return the temperature
      //! @return the temperature
      const util::ShPtr< TemperatureInterface> &GetTemperature() const;

      //! @brief return the tracker that is used to track the history of the MC process
      //! @param SP_TEMPERATURE new temperature to be used
      void SetTemperature( const util::ShPtr< TemperatureInterface> &SP_TEMPERATURE);

      //! @brief sets the boolean indicating if the start model should be modified or not
      //! @param MODIFY_START_MODEL indicates true if the start model should be modified - false if not
      void SetModifyStartModel( const bool MODIFY_START_MODEL);

      //! @brief sets the boolean indicating if the start model should be printed or not
      //! @param MODIFY_START_MODEL indicates true if the start model should be printed - false if not
      void SetPrintStartModel( const bool PRINT_START_MODEL);

      //! @brief sets the boolean indicating if every model should be printed or not
      //! @param PRINT_MODEL indicates true if the all models should be printed - false if not
      void SetPrintIterationModels( const bool PRINT_MODEL);

      //! @brief sets the boolean indicating if the end model should be printed or not
      //! @param MODIFY_END_MODEL indicates true if the end model should be printed - false if not
      void SetPrintEndModel( const bool PRINT_END_MODEL);

      //! @brief gets the filename of the pool for use with this stage
      //! @return string which is the filename of the pool for use with this stage
      const std::string &GetPoolPostfix() const;

      //! @brief sets the filename postfix of the pool for use with this stage
      //! @param POOL_POSTFIX the filename postfix of the pool for use with this stage
      void SetPoolPostfix( const std::string &POOL_POSTFIX);

      //! @brief sets the boolean indicating if the tracker history should be printed
      //! @param TRACKER_HISTORY true if the tracker history should be printer
      void SetPrintTrackerHistory( const bool &TRACKER_HISTORY);

      //! @brief initializes the stage with the current round number and stage number
      //! @param ROUND_NUMBER number of the model currently created
      //! @param STAGE_NUMBER number of the stage
      //! @param IS_LAST_STAGE indicates if it is the last stage
      void Initialize( const size_t &ROUND_NUMBER, const size_t &STAGE_NUMBER, const bool &IS_LAST_STAGE)
      {
        m_RoundNumber = ROUND_NUMBER;
        m_Number = STAGE_NUMBER;
        m_IsLastStage = IS_LAST_STAGE;
      }

      //! @brief sets the round number of the stage
      //! @param ROUND_NUMBER current round number
      void SetRoundNumber( const size_t ROUND_NUMBER)
      {
        m_RoundNumber = ROUND_NUMBER;
      }

      //! @brief sets the number of the stage
      //! @param STAGE_NUMBER number of this stage
      void SetStageNumber( const size_t STAGE_NUMBER)
      {
        m_Number = STAGE_NUMBER;
      }

      //! @brief sets the number of scores to drop for each protein model folding run
      //! @param NUMBER_SCORES the number of scores to drop
      void SetNumberScoresToDrop( const size_t &NUMBER_SCORES)
      {
        m_NumberOfScoresDropped = NUMBER_SCORES;
      }

      //! @brief return quality measures to be calculated
      //! @return quality measures to be calculated
      void SetQualityMeasures( const storage::Set< quality::Measure>);

      //! @brief gives the prefix object
      //! @return the prefix that is prepended to output files
      void SetPrefix( const std::string &PREFIX);

      //! @brief sets the path for the output
      //! @param PATH path for the output
      void SetPath( const std::string &PATH);

      //! @brief reset certain members using the round number and stage number passed
      //! @param ROUND_NUMBER round of optimization
      //! @param STAGE_NUMBER stage of optimization
      void Reset( const size_t ROUND_NUMBER, const size_t STAGE_NUMBER);

    ////////////////
    // operations //
    ////////////////

      //! @brief initiates the approximation and returns a shared pointer to the approximation result
      //! @param SP_PROTEIN_MODEL shared pointer to the starting model for the approximation
      //! @return shared pointer to the approximation result
      util::ShPtr< assemble::ProteinModel> Approximate( util::ShPtr< assemble::ProteinModel> &SP_PROTEIN_MODEL) const;

      //! @brief initialize all score
      void InitializeScores();

      //! @brief initialize all mutates
      void InitializeMutates();

      //! @brief informs the stage that it is the last stage in the approximation progress
      void SetIsLastStage();

    private:

      void InitializePool( util::ShPtr< assemble::ProteinModel> &SP_PROTEIN_MODEL) const;

      //! @brief modify the starting model
      //! @param START_MODEL starting model to be used
      void ModifyStartModel( assemble::ProteinModel &START_MODEL) const;

      //! @brief modify the terminate object
      //! @param CRITERION which will be modified by protocols
      void ModifyCriterion( opti::CriterionCombine< assemble::ProteinModel, double> &CRITERION) const;

      //! @brief modify the printer object
      //! @param PRINTER which will be modified by protocols
      void ModifyPrinter( PrinterCombined< assemble::ProteinModel, double> &PRINTER) const;

      //! @brief modify the pdb factory object
      //! @param FACTORY pdb factory to modify
      void ModifyFactory( util::ShPtr< pdb::Factory> &FACTORY) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class Stage

  } // namespace mc
} // namespace bcl

#endif // BCL_MC_STAGE_H_
