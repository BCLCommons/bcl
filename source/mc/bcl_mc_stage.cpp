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
#include "mc/bcl_mc_stage.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_printer_protein_model.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "fold/bcl_fold_default_flags.h"
#include "fold/bcl_fold_mutates.h"
#include "io/bcl_io_file.h"
#include "mc/bcl_mc_printer_combined.h"
#include "mc/bcl_mc_printer_with_criterion.h"
#include "opti/bcl_opti_criterion_all.h"
#include "opti/bcl_opti_criterion_phase.h"
#include "opti/bcl_opti_criterion_result_changed.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> Stage::s_Instance
    (
      GetObjectInstances().AddInstance( new Stage())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Stage::Stage() :
      m_Name( ""),
      m_Number( util::GetUndefined< size_t>()),
      m_FoldProtocols(),
      m_ScoreProtocols(),
      m_MutateProtocols(),
      m_MaxNumberIterations( util::GetUndefined< size_t>()),
      m_MaxNumberUnimprovedIterations( util::GetUndefined< size_t>()),
      m_ScoreFunction(),
      m_ScoreWeightSet(),
      m_MutateTree(),
      m_Mutate(),
      m_Temperature(),
      m_NumberOfScoresDropped( 0),
      m_ModifyStartModel( true),
      m_PrintStartModel( false),
      m_PrintIterationModels( false),
      m_PrintEndModel( false),
      m_PrintTrackerHistory( false),
      m_PoolPostfix( ""),
      m_RoundNumber(util::GetUndefined< size_t>()),
      m_Prefix( ""),
      m_Path(),
      m_QualityMeasures()
    {
    }

    //! @brief constructor from members
    //! @param FOLD_PROTOCOLS Fold protocols
    //! @param SCORE_PROTOCOLS Score protocols
    //! @param MUTATE_PROTOCOLS Mutate protocols
    Stage::Stage
    (
      const storage::List< fold::Protocol> &FOLD_PROTOCOLS,
      const storage::List< fold::Protocol> &SCORE_PROTOCOLS,
      const storage::List< fold::Protocol> &MUTATE_PROTOCOLS
    ) :
      m_Name( ""),
      m_Number( util::GetUndefined< size_t>()),
      m_FoldProtocols( FOLD_PROTOCOLS),
      m_ScoreProtocols( SCORE_PROTOCOLS),
      m_MutateProtocols( MUTATE_PROTOCOLS),
      m_MaxNumberIterations( util::GetUndefined< size_t>()),
      m_MaxNumberUnimprovedIterations( util::GetUndefined< size_t>()),
      m_ScoreFunction(),
      m_ScoreWeightSet(),
      m_MutateTree(),
      m_Mutate(),
      m_Temperature(),
      m_NumberOfScoresDropped( 0),
      m_ModifyStartModel( true),
      m_PrintStartModel( false),
      m_PrintIterationModels( false),
      m_PrintEndModel( false),
      m_PrintTrackerHistory( false),
      m_PoolPostfix( ""),
      m_RoundNumber(util::GetUndefined< size_t>()),
      m_Prefix( ""),
      m_Path(),
      m_QualityMeasures()
    {
    }

    //! @brief Clone function
    //! @return pointer to new Stage
    Stage *Stage::Clone() const
    {
      return new Stage( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Stage::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the list of fold protocols used
    //! @return the list of fold protocols used
    const storage::List< fold::Protocol> &Stage::GetFoldProtocols() const
    {
      return m_FoldProtocols;
    }

    //! @brief sets the list of fold protocols used
    //! @param FOLD_PROTOCOLS list of fold protocols
    void Stage::SetFoldProtocols( const storage::List< fold::Protocol> &FOLD_PROTOCOLS)
    {
      m_FoldProtocols = FOLD_PROTOCOLS;
    }

    //! @brief returns the list of score protocols used
    //! @return the list of score protocols used
    const storage::List< fold::Protocol> &Stage::GetScoreProtocols() const
    {
      return m_ScoreProtocols;
    }

    //! @brief sets the list of score protocols used
    //! @param SCORE_PROTOCOLS list of score protocols
    void Stage::SetScoreProtocols( const storage::List< fold::Protocol> &SCORE_PROTOCOLS)
    {
      m_ScoreProtocols = SCORE_PROTOCOLS;
    }

    //! @brief returns the list of mutate protocols used
    //! @return the list of mutate protocols used
    const storage::List< fold::Protocol> &Stage::GetMutateProtocols() const
    {
      return m_MutateProtocols;
    }

    //! @brief sets the list of mutate protocols used
    //! @param MUTATE_PROTOCOLS list of mutate protocols
    void Stage::SetMutateProtocols( const storage::List< fold::Protocol> &MUTATE_PROTOCOLS)
    {
      m_MutateProtocols = MUTATE_PROTOCOLS;
    }

    //! @brief returns score weightset
    //! @return scoring weightset
    const util::ShPtr< fold::ScoreWeightSet> &Stage::GetScoreWeightSet() const
    {
      return m_ScoreWeightSet;
    }

    //! @brief returns score weightset
    //! @return scoring weightset
    void Stage::SetScoreWeightSet( const util::ShPtr< fold::ScoreWeightSet> &SCORE_WEIGHT_SET)
    {
      m_ScoreWeightSet = SCORE_WEIGHT_SET;
    }

    //! @brief returns score function
    //! @return score function
    const util::ShPtr< score::ProteinModelScoreSum> &Stage::GetScoreFunction() const
    {
      return m_ScoreFunction;
    }

    //! @brief returns score function
    //! @return score function
    void Stage::SetScoreFunction( const util::ShPtr< score::ProteinModelScoreSum> &SCORE_FUNCTION)
    {
      m_ScoreFunction = SCORE_FUNCTION;
    }

    //! @brief get mutate tree
    //! @return mutate tree
    const util::ShPtr< fold::MutateTree> &Stage::GetMutateTree() const
    {
      return m_MutateTree;
    }

    //! @brief get mutate tree
    //! @return mutate tree
    void Stage::SetMutateTree( const util::ShPtr< fold::MutateTree> &MUTATE_TREE)
    {
      m_MutateTree = MUTATE_TREE;
    }

    //! @brief get mutate functions for protein models
    //! @return mutate functions for protein models
    const util::ShPtr< math::MutateInterface< assemble::ProteinModel> > &Stage::GetMutate() const
    {
      return m_Mutate;
    }

    //! @brief get mutate functions for protein models
    //! @return mutate functions for protein models
    void Stage::SetMutate( const util::ShPtr< math::MutateInterface< assemble::ProteinModel> > &MUTATE)
    {
      m_Mutate = MUTATE;
    }

    //! @brief return the temperature
    //! @return the temperature
    const util::ShPtr< TemperatureInterface> &Stage::GetTemperature() const
    {
      return m_Temperature;
    }

    //! @brief sets the temperature
    //! @param SP_TEMPERATURE new temperature to be usedif( is_last_stage)
    void Stage::SetTemperature( const util::ShPtr< TemperatureInterface> &SP_TEMPERATURE)
    {
      m_Temperature = SP_TEMPERATURE;
    }

    //! @brief sets the boolean indicating if the start model should be modified or not
    //! @param MODIFY_START_MODEL indicates true if the start model should be modified - false if not
    void Stage::SetModifyStartModel( const bool MODIFY_START_MODEL)
    {
      m_ModifyStartModel = MODIFY_START_MODEL;
    }

    //! @brief sets the boolean indicating if the start model should be printed or not
    //! @param MODIFY_START_MODEL indicates true if the start model should be printed - false if not
    void Stage::SetPrintStartModel( const bool PRINT_START_MODEL)
    {
      m_PrintStartModel = PRINT_START_MODEL;
    }

    //! @brief sets the boolean indicating if every model should be printed or not
    //! @param PRINT_MODEL indicates true if the all models should be printed - false if not
    void Stage::SetPrintIterationModels( const bool PRINT_MODEL)
    {
      m_PrintIterationModels = PRINT_MODEL;
    }

    //! @brief sets the boolean indicating if the end model should be printed or not
    //! @param MODIFY_END_MODEL indicates true if the end model should be printed - false if not
    void Stage::SetPrintEndModel( const bool PRINT_END_MODEL)
    {
      m_PrintEndModel = PRINT_END_MODEL;
    }

    //! @brief gets the filename postfix of the pool for use with this stage
    //! @return string which is the filename postfix of the pool for use with this stage
    const std::string &Stage::GetPoolPostfix() const
    {
      return m_PoolPostfix;
    }

    //! @brief sets the filename postfix of the pool for use with this stage
    //! @param POOL_POSTFIX the filename postfix of the pool for use with this stage
    void Stage::SetPoolPostfix( const std::string &POOL_POSTFIX)
    {
      m_PoolPostfix = POOL_POSTFIX;
    }
    //! @brief return quality measures to be calculated
    //! @return quality measures to be calculated
    void Stage::SetQualityMeasures( const storage::Set< quality::Measure> QUALITY_MEASURES)
    {
       m_QualityMeasures = QUALITY_MEASURES;
    }
    //! @brief gives the prefix object
    //! @return the prefix that is prepended to output files
    void Stage::SetPrefix( const std::string &PREFIX)
    {
       m_Prefix = PREFIX;
    }

    //! @brief sets the boolean indicating if the tracker history should be printed
    //! @param TRACKER_HISTORY true if the tracker history should be printer
    void Stage::SetPrintTrackerHistory( const bool &TRACKER_HISTORY)
    {
      m_PrintTrackerHistory = TRACKER_HISTORY;
    }

    //! @brief sets the path for the output
    //! @param PATH path for the output
    void Stage::SetPath( const std::string &PATH)
    {
       m_Prefix = PATH;
    }

    //! @brief initiates the approximation and returns a shared pointer to the approximation result
    //! @param SP_PROTEIN_MODEL shared pointer to the starting model for the approximation
    //! @return shared pointer to the approximation result
    util::ShPtr< assemble::ProteinModel> Stage::Approximate( util::ShPtr< assemble::ProteinModel> &SP_PROTEIN_MODEL) const
    {
      // if a new sse pool is needed for this stage
      if( !m_PoolPostfix.empty())
      {
        InitializePool( SP_PROTEIN_MODEL);
      }

      // modify the start model
      if( m_ModifyStartModel)
      {
        BCL_MessageStd( "modifying start model for stage " + util::Format()( m_Number + 1));
        ModifyStartModel( *SP_PROTEIN_MODEL);
      }

    /////////////////////////////////////
    // setting up printers and storage //
    /////////////////////////////////////

      // create the combined printer
      PrinterCombined< assemble::ProteinModel, double> printer_combine;

      // create the criterion for the protein model printer
      storage::Set< opti::PhaseEnum> phases;
      if( m_PrintStartModel)
      {
        phases.Insert( opti::e_Start);
      }
      if( m_PrintIterationModels)
      {
        phases.Insert( opti::e_Iteration);
      }
      if( m_PrintEndModel)
      {
        phases.Insert( opti::e_End);
      }

      util::ShPtr< opti::CriterionInterface< assemble::ProteinModel, double> > sp_print_criterion
      (
        new opti::CriterionPhase< assemble::ProteinModel, double>( phases)
      );
      if( m_PrintIterationModels)
      {
        opti::CriterionAll< assemble::ProteinModel, double> comb;
        comb.InsertCriteria( *sp_print_criterion);
        comb.InsertCriteria( opti::CriterionResultChanged< assemble::ProteinModel, double>());
        sp_print_criterion = util::ShPtr< opti::CriterionInterface< assemble::ProteinModel, double> >( comb.Clone());
      }
      // create the protein model printer and add it to the combined printer
      util::ShPtr< PrintInterface< assemble::ProteinModel, double> > sp_printer
      (
        new assemble::PrinterProteinModel
        (
          fold::GetSetup().GetPrefix() + ( !m_IsLastStage ? m_Name : ""),
          fold::GetSetup().GetStorage(),
          fold::GetSetup().GetSuperimposeMeasure()
        )
      );
      util::ShPtr< PrintInterface< assemble::ProteinModel, double> > sp_model_printer
      (
        new PrinterWithCriterion< assemble::ProteinModel, double>( sp_printer, sp_print_criterion)
      );
      printer_combine.Insert( sp_model_printer);
      printer_combine.Initialize( m_RoundNumber, m_Number);
      printer_combine.SetPrefix( m_Prefix);

      // create the pdb factory
      util::ShPtr< pdb::Factory> factory( new pdb::Factory());
      ModifyFactory( factory);

      // modify the factory so that it adds the scoring table to the pdb when printed
      util::ShPtr< assemble::ProteinStorageFile> sp_storage( fold::GetSetup().GetStorage());
      sp_storage->SetFactory( factory);

    ////////////////////////////////////////////
    // preparing and performing approximation //
    ////////////////////////////////////////////

      // create the metropolis criterion
      const double metropolis_min_change( 0.0001);
      const util::ShPtr< Metropolis< double> > sp_metropolis
      (
        new Metropolis< double>( m_Temperature, true, metropolis_min_change)
      );

      // create the termination criteria for the approximation
      opti::CriterionCombine< assemble::ProteinModel, double> criterion_combine;
      ModifyCriterion( criterion_combine);

      storage::Vector< std::string> dropped_functions( m_NumberOfScoresDropped);
      storage::Vector< double> dropped_values( m_NumberOfScoresDropped);
      util::ShPtr< score::ProteinModelScoreSum> sum( m_ScoreFunction);
      static storage::Set< std::string> s_prev_dropped_functions;
      static bool s_score_dropped_continuity( false);
      if( m_NumberOfScoresDropped)
      {
        storage::Vector< std::string> functions;
        if( s_score_dropped_continuity)
        {
          storage::Vector< std::string> prev_dropped_functions;
          functions = m_ScoreFunction->GetFunctionSchemes();
          functions.PopBack();
          functions.Shuffle();
          size_t n_not_dropped( 0);
          for( auto itr_a( functions.Begin()), itr_place( functions.Begin()), itr_end( functions.End()); itr_a != itr_end; ++itr_a)
          {
            if( !s_prev_dropped_functions.Contains( *itr_a))
            {
              if( itr_a != itr_place)
              {
                *itr_place = *itr_a;
              }
              ++n_not_dropped;
              ++itr_place;
            }
            else
            {
              prev_dropped_functions.PushBack( *itr_a);
            }
          }
          functions.Resize( n_not_dropped);
          functions.InsertElements( 0, prev_dropped_functions);
          s_prev_dropped_functions.Reset();
          s_prev_dropped_functions.InsertElements( functions.Begin(), functions.Begin() + m_NumberOfScoresDropped);
        }
        else
        {
          functions = m_ScoreFunction->GetFunctionSchemes();
          functions.PopBack();
          functions.Shuffle();
        }
        for( size_t i( 0); i < m_NumberOfScoresDropped; ++i)
        {
          dropped_functions( i) = functions( i);
          dropped_values( i) = sum->GetTerm( functions( i)).First();
          sum->SetCoefficient( functions( i), 0.0);
        }
      }

      // create the approximator
      Approximator< assemble::ProteinModel, double> approximator
      (
        *m_ScoreFunction,
        *m_Mutate,
        *sp_metropolis,
        criterion_combine,
        *SP_PROTEIN_MODEL
      );

      // add the combined printer to the approximator
      approximator.SetPrinter( printer_combine);

      // start the approximation
      approximator.Approximate();

      // get the approximation result
      util::ShPtr< storage::Pair< assemble::ProteinModel, double> > sp_final_model_and_score_pair
      (
        approximator.GetTracker().GetBest()
      );
      util::ShPtr< assemble::ProteinModel> sp_model_return
      (
        util::ShPtr< assemble::ProteinModel>( approximator.GetTracker().GetBest()->First().HardCopy())
      );

      // output information regarding the final model of this stage
      BCL_MessageCrt( "#SSEs: " + util::Format()( sp_model_return->GetNumberSSEs()));

      // if this stage is the last stage print the final model
      if( m_IsLastStage)
      {
        // initialize with round number
        sp_printer->Initialize( m_RoundNumber);

        // set the prefix
        sp_printer->SetPrefix( m_Prefix);

        // print the final structure
        sp_printer->Print( approximator.GetTracker());

        s_prev_dropped_functions.Reset();
      }

      if( m_NumberOfScoresDropped)
      {
        for( size_t i( 0); i < m_NumberOfScoresDropped; ++i)
        {
          sum->SetCoefficient( dropped_functions( i), dropped_values( i));
        }
      }

      return sp_model_return;
    }

    //! @brief informs the stage that it is the last stage in the approximation progress
    void Stage::SetIsLastStage()
    {
      m_IsLastStage = true;
    }

    //! @brief reset certain members using the round number and stage number passed
    //! @param ROUND_NUMBER round of optimization
    //! @param STAGE_NUMBER stage of optimization
    void Stage::Stage::Reset( const size_t ROUND_NUMBER, const size_t STAGE_NUMBER)
    {
      // reset the temperature
      m_Temperature->Reset();

      // iterate over the protocols
      for
      (
        storage::List< fold::Protocol>::iterator
          protocol_itr( m_FoldProtocols.Begin()), protocol_itr_end( m_FoldProtocols.End());
        protocol_itr != protocol_itr_end; ++protocol_itr
      )
      {
        // call reset
        ( **protocol_itr)->Reset();
      }
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initialize all score
    void Stage::InitializeScores()
    {
      // iterate over the list
      for
      (
        storage::List< fold::Protocol>::iterator itr( m_ScoreProtocols.Begin()), itr_end( m_ScoreProtocols.End());
        itr != itr_end; ++itr
      )
      {
        // initialize scores for this protocol
        ( **itr)->InitializeScores();
      }
    }

    //! @brief initialize all mutates
    void Stage::InitializeMutates()
    {
      // iterate over the list
      for
      (
        storage::List< fold::Protocol>::iterator itr( m_MutateProtocols.Begin()), itr_end( m_MutateProtocols.End());
        itr != itr_end; ++itr
      )
      {
        // initialize mutates for this protocol
        ( **itr)->InitializeMutates();
      }
    }

    void Stage::InitializePool( util::ShPtr< assemble::ProteinModel> &SP_PROTEIN_MODEL) const
    {
      // initialize empty pool
      util::ShPtr< assemble::SSEPool> sp_pool( new assemble::SSEPool());

      // initialize map to hold the min pool lengths
      storage::Map< biol::SSType, size_t> min_pool_sse_lengths( assemble::SSEPool::GetCommandLineMinSSELengths());

      // open pool file
      io::IFStream read;
      const std::string pool_file
      (
        assemble::SSEPool::GetFlagPoolPrefix()->GetFirstParameter()->GetValue() + "." + GetPoolPostfix()
      );
      BCL_MessageStd( "Reading pool from file for stage setting " + pool_file);
      io::File::MustOpenIFStream( read, pool_file);

      // read pool
      sp_pool->ReadSSEPool
      (
        read,
        *SP_PROTEIN_MODEL,
        min_pool_sse_lengths[ biol::GetSSTypes().HELIX],
        min_pool_sse_lengths[ biol::GetSSTypes().STRAND]
      );

      // if separate pool flag was provided and pool is not overlapping
      if( fold::DefaultFlags::GetFlagPoolSeparate()->GetFlag() && !sp_pool->IsOverlapping())
      {
        // separate pools
        BCL_MessageStd( "separating adjoining SSEs in the pool");
        const bool success
        (
          sp_pool->Separate
          (
            assemble::SSEPool::GetCommandLineMinSSELengths(),
            fold::DefaultFlags::GetFlagPoolSeparate()->GetFirstParameter()->GetNumericalValue< size_t>()
          )
        );
        // make sure to prune it to remove loops
        sp_pool->Prune( assemble::SSEPool::GetCommandLineMinSSELengths());

        BCL_MessageStd( "separating success: " + util::Format()( success));
      }

      BCL_MessageStd( "pool set to ");
      sp_pool->WriteSSEPool( util::GetLogger());

      SP_PROTEIN_MODEL->SetSSEPoolData( sp_pool);
    }

    //! @brief modify the starting model
    //! @param START_MODEL starting model to be used
    void Stage::ModifyStartModel( assemble::ProteinModel &START_MODEL) const
    {
      // iterate over the protocols
      for
      (
        storage::List< fold::Protocol>::const_iterator
          protocol_itr( m_FoldProtocols.Begin()), protocol_itr_end( m_FoldProtocols.End());
        protocol_itr != protocol_itr_end; ++protocol_itr
      )
      {
        // modify the start model
        ( **protocol_itr)->ModifyStartModel( START_MODEL);
      }
    }

    //! @brief modify the terminate object
    //! @param CRITERION which will be modified by protocols
    void Stage::ModifyCriterion( opti::CriterionCombine< assemble::ProteinModel, double> &CRITERION) const
    {
      // iterate over the protocols
      for
      (
        storage::List< fold::Protocol>::const_iterator
          protocol_itr( m_FoldProtocols.Begin()), protocol_itr_end( m_FoldProtocols.End());
        protocol_itr != protocol_itr_end; ++protocol_itr
      )
      {
        // modify the terminate
        ( **protocol_itr)->ModifyCriterion( CRITERION, *this);
      }
    }

    //! @brief modify the printer object
    //! @param PRINTER which will be modified by protocols
    void Stage::ModifyPrinter( PrinterCombined< assemble::ProteinModel, double> &PRINTER) const
    {
      // iterate over the protocols
      for
      (
        storage::List< fold::Protocol>::const_iterator
          protocol_itr( m_FoldProtocols.Begin()), protocol_itr_end( m_FoldProtocols.End());
        protocol_itr != protocol_itr_end; ++protocol_itr
      )
      {
        // modify the printer
        ( **protocol_itr)->ModifyPrinter( PRINTER, *this);
      }
    }

    //! @brief modify the pdb factory object
    //! @param FACTORY pdb factory to modify
    void Stage::ModifyFactory( util::ShPtr< pdb::Factory> &FACTORY) const
    {
      // reset factory printers
      FACTORY->ResetPrinters();

      // iterate over the protocols
      for
      (
        storage::List< fold::Protocol>::const_iterator
          protocol_itr( m_FoldProtocols.Begin()), protocol_itr_end( m_FoldProtocols.End());
        protocol_itr != protocol_itr_end; ++protocol_itr
      )
      {
        // modify the factory
        ( **protocol_itr)->ModifyFactory( FACTORY, *this);
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Stage::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Stage::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return the stream
      return OSTREAM;
    }

  } // namespace mc
} // namespace bcl
