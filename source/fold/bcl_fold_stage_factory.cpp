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
#include "fold/bcl_fold_stage_factory.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_default_flags.h"
#include "fold/bcl_fold_mutates.h"
#include "io/bcl_io_file.h"
#include "mc/bcl_mc_temperature_accepted.h"
#include "storage/bcl_storage_hash_map.h"
#include "util/bcl_util_logger_interface.h"
#include "util/bcl_util_message.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    const std::string StageFactory::s_LineTypeNames[ s_NumberLineTypes] =
    {
      "NUMBER_CYCLES",         // Number of cycles
      "STAGE",                 // Start of a stage definition, followed by an optional name
      "TYPE",                  // type of the approximation used for this stage
      "FOLD_PROTOCOLS",        // fold protocols to use
      "SCORE_PROTOCOLS",       // score protocols to use (optional)
      "MUTATE_PROTOCOLS",      // mutate protocols to use (optional)
      "SCORE_WEIGHTSET",       // score weight set table
      "SCORE_WEIGHTSET_FILE",  // file containing score weight set
      "SCORE_DROPOUT_RATE",    // Fraction of scores to drop randomly for each protein model
      "MUTATE_WEIGHTSET",      // mutate weight set table
      "MUTATE_WEIGHTSET_FILE", // file containing mutate weight set
      "NUMBER_ITERATIONS",     // followed by total number iterations and max number unimproved
      "MODIFY_START_MODEL",    // true if the stage should modify the start model, false otherwise
      "PRINT_START_MODEL",     // whether to print the start model for this stage
      "PRINT_ITERATION_MODELS",// whether to print the all models for this stage
      "PRINT_END_MODEL",       // whether to print the end model for this stage
      "PRINT_TRACKER_HISTORY", // whether to print the tracker history for this stage
      "POOL_POSTFIX",          // file postfix containing pool that should be used for this stage
      "STAGE_END"              // end of a stage definition
    };

    //! @brief finds the LineType enum that corresponds to given string
    //! @param LINE_NAME Line name of interest
    //! @return the LineType enum that corresponds to given string
    StageFactory::LineType StageFactory::LineTypeFromString( const std::string &LINE_NAME)
    {
      // iterate over line types
      for( size_t i( 0); i < s_NumberLineTypes; ++i)
      {
        // if strings match then return corresponding enum
        if( LINE_NAME == s_LineTypeNames[ i])
        {
          return LineType( i);
        }
      }
      // otherwise exit
      BCL_Exit( "The provided stage file line name is not valid " + LINE_NAME, -1);
      // this line needs to be added to prevent the compiler warning
      return e_NumberCycles;
    }

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> StageFactory::s_Instance
    (
      GetObjectInstances().AddInstance( new StageFactory())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    StageFactory::StageFactory()
    {
    }

    //! @brief Clone function
    //! @return pointer to new StageFactory
    StageFactory *StageFactory::Clone() const
    {
      return new StageFactory( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &StageFactory::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief construct stages
    //! @return vector of stages
    util::ShPtrVector< StageInterface> StageFactory::CreateStages()
    {
      // initialize stages
      util::ShPtrVector< StageInterface> stages;

      // if StageFile was specified
      if( DefaultFlags::GetFlagStagesFileRead()->GetFlag())
      {
        BCL_MessageStd( "Constructing multi-stage from stages file");
        BCL_MessageStd( "All related flags specified in commandline will be ignored");

        // get the number of cycles from command line
        size_t number_cycles
        (
          DefaultFlags::GetFlagStagesNumberCycles()->GetFirstParameter()->GetNumericalValue< size_t>()
        );

        // construct the stages from the stage file
        stages = CreateStagesFromFile( number_cycles);
      }
      else // construct from command line options
      {
        stages.PushBack( CreateStage());
      }

      return stages;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &StageFactory::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &StageFactory::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief sets the temperature for the given STAGE
    void StageFactory::SetTemperature( mc::Stage &STAGE)
    {
      // create the temperature control
      util::ShPtr< mc::TemperatureInterface> sp_temperature
      (
        new mc::TemperatureAccepted
        (
          // start fraction
          mc::TemperatureAccepted::GetParameterStartFraction()->GetNumericalValue< double>(),
          // end fraction
          mc::TemperatureAccepted::GetParameterEndFraction()->GetNumericalValue< double>(),
          // total number of steps
          STAGE.GetMaxNumberIterations(),
          // start temperature

          mc::TemperatureAccepted::GetParameterStartTemperature()->GetNumericalValue< double>(),
          // nr steps between each update
          mc::TemperatureAccepted::GetParameterUpdateInterval()->GetNumericalValue< size_t>()
        )
      );

      // set the temperature
      STAGE.SetTemperature( sp_temperature);
    }

    //! @brief construct and return stages from stage file
    //! @param NUMBER_CYCLES number of cycles which can be updated from command line
    //! @return stages constructed from stage file
    util::ShPtrVector< StageInterface> StageFactory::CreateStagesFromFile( size_t &NUMBER_CYCLES)
    {
      // stores the stages
      util::ShPtrVector< StageInterface> stages;

      // open a stream to the specified stage file
      io::IFStream read;
      const std::string stage_file( DefaultFlags::GetFlagStagesFileRead()->GetFirstParameter()->GetValue());
      io::File::MustOpenIFStream( read, stage_file);
      BCL_MessageStd( "Reading stage file " + stage_file)

      // stores the total number of stages defined in the stage file
      size_t number_stages( 0);

      // stores the number of cycles
      size_t number_cycles
      (
        DefaultFlags::GetFlagStagesNumberCycles()->GetFirstParameter()->GetNumericalValue< size_t>()
      );

      // parse the stage file
      std::string line;
      while( !read.eof() && std::getline( read, line))
      {
        // if the line is empty, skip it
        if( util::TrimString( line).empty())
        {
          continue;
        }

        // determine the line type
        storage::Vector< std::string> strings( util::SplitString( line));
        LineType line_type( LineTypeFromString( strings.FirstElement()));

        // if this line sets the number of cycles
        if( line_type == e_NumberCycles)
        {
          number_cycles = util::ConvertStringToNumericalValue< size_t>( strings( 1));
        }
        else if( line_type == e_StageStart) // This line marks the begin of a stage definition
        {
          // determine the name for the stage
          ++number_stages;
          std::string name;
          if( strings.GetSize() > 1)
          {
            name = strings( 1);
          }
          else
          {
            name = "stage_" + util::Format()( number_stages);
          }

          // construct the stage form the stage file
          stages.PushBack( CreateStage( read, name, number_stages));
        }
        else
        {
          BCL_Exit( strings.FirstElement() + " is unknown or cannot be used outside of a stage", -1);
        }
      }

      // repeat the stages based on the number of cycles provided
      util::ShPtrVector< StageInterface> original_stages( stages);
      for( size_t i( 1); i < number_cycles; ++i)
      {
        for
        (
          util::ShPtrVector< StageInterface>::const_iterator it( original_stages.Begin());
          it != original_stages.End();
          ++it
        )
        {
          stages.PushBack( *it);
        }
      }

      return stages;
    }

    //! @brief constructs one stage from an IFStream
    //! @param READ IFStream to create the stage from
    //! @param NAME the name of the stage
    //! @param NUMBER the number of the stage
    //! @return shared pointer to the constructed stage
    util::ShPtr< StageInterface> StageFactory::CreateStage( io::IFStream &READ, const std::string &NAME, const size_t &NUMBER)
    {
      // create a map to store line types and their respective values
      storage::HashMap< size_t, storage::Vector< std::string> > linetypes_values;

      // add the name of the stage
      linetypes_values.Insert
      (
        storage::Pair< size_t, storage::Vector< std::string> >
        (
          e_StageStart,
          storage::Vector< std::string>( 1, NAME)
        )
      );

      // parse the stage section until the end of the stage or the end of the file is reached
      std::string line;
      bool stage_closed( false);
      while( !READ.eof() && std::getline( READ, line))
      {
        // determine the line type
        const storage::Vector< std::string> strings( util::SplitString( line));
        const LineType line_type( LineTypeFromString( strings.FirstElement()));

        // if the end of the stage definition is reached
        if( line_type == e_StageEnd)
        {
          stage_closed = true;
          break;
        }

        // check if the same line type has already been read for this stage
        BCL_Assert( linetypes_values.Count( line_type) == 0, "Redefinition of " + s_LineTypeNames[ line_type]);

        // store the line type and its values for further processing
        storage::Vector< std::string> values;
        if( strings.GetSize() > 1)
        {
          values = storage::Vector< std::string>( strings.Begin() + 1, strings.Last() + 1);
        }
        linetypes_values.Insert( storage::Pair< size_t, storage::Vector< std::string> >( line_type, values));
      }

      // make sure that the definition for this stage was complete
      BCL_Assert( stage_closed, "Stage definition is not complete");

      // make sure that the stage type was given
      BCL_Assert( linetypes_values.Count( e_Type) == 1, "Stage type was not provided");

      // construct the stage according to the given type
      util::ShPtr< StageInterface> sp_stage;
      const std::string &type( linetypes_values[ e_Type].FirstElement());
      if( type.compare( "MCM") == 0)
      {
        sp_stage = util::ShPtr< StageInterface>( CreateMcStage( linetypes_values, NUMBER));
      }
      else if( type.compare( "GRADMIN") == 0)
      {

      }
      else
      {
        BCL_Exit( type + " is an unknown stage type", -1);
      }

      return sp_stage;
    }

    //! @brief constructs a stage from the command line
    //! @return shared pointer to the constructed stage
    util::ShPtr< StageInterface> StageFactory::CreateStage()
    {
      BCL_MessageStd( "Constructing single stage from command line");

      // constructing a single stage from the specified protocols
      util::ShPtr< mc::Stage> sp_stage
      (
        new mc::Stage
        (
          GetProtocols().GetCommandLineProtocolList(),
          GetProtocols().GetCommandLineScoreProtocolList(),
          GetProtocols().GetCommandLineMutateProtocolList()
        )
      );

      // initialize the stage
      sp_stage->InitializeScores();
      sp_stage->InitializeMutates();
      SetSingleStageScoreFunction( *sp_stage);
      SetSingleStageMutate( *sp_stage);

      sp_stage->SetMaxNumberIterations
      (
        DefaultFlags::GetFlagMCNumberIterations()->GetParameterList()( 0)->GetNumericalValue< size_t>()
      );
      sp_stage->SetMaxNumberUnimprovedIterations
      (
        DefaultFlags::GetFlagMCNumberIterations()->GetParameterList()( 1)->GetNumericalValue< size_t>()
      );
      sp_stage->SetName( "Stage_0");
      SetTemperature( *sp_stage);
      sp_stage->SetPrefix( GetSetup().GetPrefix());

      return sp_stage;
    }

    //! @brief creates a monte carlo metropolis stage based on the given parameters
    //! @param SETTINGS contains the parameters to construct the stage from
    //! @return shared pointer to the constructed stage
    util::ShPtr< StageInterface> StageFactory::CreateMcStage
    (
      storage::HashMap< size_t, storage::Vector< std::string> > &SETTINGS,
      const size_t &NUMBER
    )
    {
      // create the fold, score and mutate protocols from the arguments given in the stage file
      storage::List< Protocol> fold_protocols( CreateUniqueProtocolList( ConstructProtocolList( SETTINGS[ e_FoldProtocols])));
      storage::List< Protocol> score_protocols( CreateUniqueProtocolList( ConstructProtocolList( SETTINGS[ e_ScoreProtocols])));
      storage::List< Protocol> mutate_protocols( CreateUniqueProtocolList( ConstructProtocolList( SETTINGS[ e_MutateProtocols])));

      // the final unique protocols
      storage::List< Protocol> fold_protocols_final;
      storage::List< Protocol> score_protocols_final;
      storage::List< Protocol> mutate_protocols_final;

      // if the fold protocols have not been specified, they are constructed from the score and mutate protocols
      if( fold_protocols.IsEmpty())
      {
        BCL_MessageStd( "No fold protocols specified, constructing from score and mutate protocols");
        BCL_Assert
        (
          !score_protocols.IsEmpty() && !mutate_protocols.IsEmpty(),
          "Both score and mutate protocols need to be set if fold protocols are not provided for stage"
        );
        fold_protocols = storage::List< Protocol>( score_protocols);
        fold_protocols.Append( mutate_protocols);
        fold_protocols_final = CreateUniqueProtocolList( fold_protocols);
      }
      else // if the fold protocols have been specified
      {
        // if the score protocols have not been specified, they are constructed from the fold protocols
        if( score_protocols.IsEmpty())
        {
          BCL_MessageStd( "No score protocols specified, constructing from fold protocols");
          score_protocols = fold_protocols;
        }

        // if the mutate protocols have not been specified, they are constructed from the fold protocols
        if( mutate_protocols.IsEmpty())
        {
          BCL_MessageStd( "No mutate protocols specified, constructing from fold protocols");
          mutate_protocols = fold_protocols;
        }

        // use the provided fold protocols as the final fold protocols
        fold_protocols_final = fold_protocols;
      }

      // initialize the score protocols
      for( storage::List< Protocol>::iterator it( score_protocols.Begin()); it != score_protocols.End(); ++it)
      {
        // only initialize if the have not been initialized yet
        if( std::find( score_protocols_final.Begin(), score_protocols_final.End(), *it) == score_protocols_final.End())
        {
          BCL_MessageStd( "\t\tInitializing scores from protocol " + it->GetName());
          ( **it)->InitializeScores();
        }
      }

      // initialize the mutate protocols
      for( storage::List< Protocol>::iterator it( mutate_protocols.Begin()); it != mutate_protocols.End(); ++it)
      {
        // only initialize if the have not been initialized yet
        if( std::find( mutate_protocols_final.Begin(), mutate_protocols_final.End(), *it) == mutate_protocols_final.End())
        {
          BCL_MessageStd( "\t\tInitializing mutates from protocol " + it->GetName());
          // initialize mutate
          ( **it)->InitializeMutates();
          // insert into initialized mutates
          mutate_protocols_final.PushBack( *it);
        }
      }

      // create the score weight set
      util::ShPtr< ScoreWeightSet> sp_score_weight;
      // if a score weight set file was provided create the score weight set from it
      if( SETTINGS.Count( e_ScoreWeightSetFile) == 1)
      {
        const std::string &score_weightset_file( SETTINGS[ e_ScoreWeightSetFile]( 0));

        storage::Table< double> score_weightset = ReadTableFromFile( score_weightset_file);

        sp_score_weight = util::ShPtr< ScoreWeightSet>( new ScoreWeightSet( score_weightset));
      }
      else // otherwise construct the weight set from the score protocols
      {
        BCL_MessageStd( "Score weight set not specified, constructing from score protocols");
        sp_score_weight = util::ShPtr< ScoreWeightSet>( new ScoreWeightSet);
        for( storage::List< Protocol>::const_iterator it( score_protocols.Begin()); it != score_protocols.End(); ++it)
        {
          ( **it)->ModifyScoreWeightSet( *sp_score_weight);
        }
      }

      // create the scoring function
      util::ShPtr< score::ProteinModelScoreSum> sp_scoring_function( sp_score_weight->ConstructScoreSum());

      // create the mutate weight set
      util::ShPtr< MutateTree> sp_mutate_tree;
      // if a mutate weight set file was provided create the score weight set from it
      if( SETTINGS.Count( e_MutateWeightSetFile) == 1)
      {
        const std::string &mutate_weightset_file( SETTINGS[ e_MutateWeightSetFile]( 0));
        storage::Table< double> mutate_weightset = ReadTableFromFile( mutate_weightset_file);
        sp_mutate_tree = util::ShPtr< MutateTree>( new MutateTree( mutate_weightset));
      }
      else // otherwise construct the weight set from the mutate protocols
      {
        BCL_MessageStd( "Mutate weight set not specified, constructing from mutate protocols");
        sp_mutate_tree = util::ShPtr< MutateTree>( new MutateTree());
        for( storage::List< Protocol>::const_iterator it( mutate_protocols_final.Begin()); it != mutate_protocols_final.End(); ++it)
        {
          sp_mutate_tree = ( **it)->GetMutateTree();
          ( **it)->MergeAndModifyMutateTree( *sp_mutate_tree);
        }
      }

      // set the maximum number of steps / of unimproved steps in a row
      const size_t max_num_steps
      (
        util::ConvertStringToNumericalValue< size_t>( SETTINGS[ e_NumberIterations]( 0))
      );
      const size_t max_num_unimproved_steps
      (
        util::ConvertStringToNumericalValue< size_t>( SETTINGS[ e_NumberIterations]( 1))
      );
      const size_t scores_to_drop
      (
        SETTINGS[ e_ScoreDropoutRate].GetSize()
        ? util::ConvertStringToNumericalValue< double>( SETTINGS[ e_ScoreDropoutRate]( 0))
          * sp_score_weight->GetWeightMap().GetSize()
        : 0
      );

      // set the printing options
      const bool print_start_model( SETTINGS.Count( e_PrintStartModel) == 0 ? false : true);
      const bool print_itr_model( SETTINGS.Count( e_PrintIterationModel) == 0 ? false : true);
      const bool print_end_model( SETTINGS.Count( e_PrintEndModel) == 0 ? false : true);
      const bool print_tracker_history( SETTINGS.Count( e_PrintTrackerHistory) == 0 ? false : true);

      // set the pool for this stage if specified
      const std::string &pool_postfix( SETTINGS.Count( e_PoolPostfix) == 1 ? SETTINGS[ e_PoolPostfix]( 0) : "");

      // create the mcm stage
      util::ShPtr< mc::Stage> sp_stage( new mc::Stage( fold_protocols_final, score_protocols_final, mutate_protocols_final));
      sp_stage->SetStageNumber( NUMBER);
      sp_stage->SetScoreWeightSet( sp_score_weight);
      sp_stage->SetScoreFunction( sp_scoring_function);
      sp_stage->SetNumberScoresToDrop( scores_to_drop);
      sp_stage->SetMutateTree( sp_mutate_tree);
      sp_stage->SetMutate( sp_mutate_tree->ConstructMutates());
      sp_stage->SetMaxNumberIterations( max_num_steps);
      sp_stage->SetMaxNumberUnimprovedIterations( max_num_unimproved_steps);
      sp_stage->SetPrintStartModel( print_start_model);
      sp_stage->SetPrintIterationModels( print_itr_model);
      sp_stage->SetPrintEndModel( print_end_model);
      sp_stage->SetPrintTrackerHistory( print_tracker_history);
      sp_stage->SetPoolPostfix( pool_postfix);
      sp_stage->SetPrefix( GetSetup().GetPrefix());

      // set the temperature for this stage
      SetTemperature( *sp_stage);

      if( util::GetMessenger().GetCurrentMessageLevel() >= util::Message::e_Verbose)
      {
        WriteStage( *sp_stage, util::GetLogger());
      }

      return sp_stage;
    }

    //! @brief write stage to stream
    //! @param STAGE Stage to be written
    //! @param OSTREAM ostream to be written to
    //! @return ostream which was written to
    std::ostream &StageFactory::WriteStage( const mc::Stage &STAGE, std::ostream &OSTREAM)
    {
      // initialize indent
      static const std::string s_indent( "  ");

      // write name
      OSTREAM << s_LineTypeNames[ e_StageStart] << ' ' << STAGE.GetName() << '\n';

      // write protocols
      OSTREAM << s_indent << s_LineTypeNames[ e_FoldProtocols];
      WriteProtocolList( STAGE.GetFoldProtocols(), OSTREAM);
      OSTREAM << s_indent << s_LineTypeNames[ e_ScoreProtocols];
      WriteProtocolList( STAGE.GetScoreProtocols(), OSTREAM);
      OSTREAM << s_indent << s_LineTypeNames[ e_MutateProtocols];
      WriteProtocolList( STAGE.GetMutateProtocols(), OSTREAM);

      // write the number of iterations
      OSTREAM << s_indent << s_LineTypeNames[ e_NumberIterations] << ' '
             << STAGE.GetMaxNumberIterations() << ' '
             << STAGE.GetMaxNumberUnimprovedIterations() << '\n';

      // write score weightset
      OSTREAM << s_indent << s_LineTypeNames[ e_ScoreWeightSet] << '\n';
      STAGE.GetScoreWeightSet()->CreateTable().WriteFormatted( OSTREAM) << '\n';

      // write mutate weightset
      OSTREAM << s_indent << s_LineTypeNames[ e_MutateWeightSet] << '\n';
      STAGE.GetMutateTree()->CreateTable().WriteFormatted( OSTREAM) << '\n';

      // write end
      OSTREAM << s_LineTypeNames[ e_StageEnd] << '\n';

      // end
      return OSTREAM;
    }

    //! @brief write given protocol list to stream
    //! @param PROTOCOL_LIST Protocol list of interest
    //! @param OSTREAM ostream to be written to
    //! @return ostream which was written to
    std::ostream &StageFactory::WriteProtocolList
    (
      const storage::List< Protocol> &PROTOCOL_LIST,
      std::ostream &OSTREAM
    )
    {
      // iterate over the list
      for
      (
        storage::List< Protocol>::const_iterator itr( PROTOCOL_LIST.Begin()), itr_end( PROTOCOL_LIST.End());
        itr != itr_end; ++itr
      )
      {
        OSTREAM << ' ' << itr->GetName();
      }
      OSTREAM << '\n';

      // end
      return OSTREAM;
    }

    //! @brief removes duplicates if any from the given protocol list
    //! @param PROTOCOL_LIST Protocol list of interest
    //! @return list of unique protocols
    storage::List< Protocol> StageFactory::CreateUniqueProtocolList( const storage::List< Protocol> &PROTOCOL_LIST)
    {
      // initialize list
      storage::List< Protocol> list;

      // iterate over the given list
      for
      (
        storage::List< Protocol>::const_iterator itr( PROTOCOL_LIST.Begin()), itr_end( PROTOCOL_LIST.End());
        itr != itr_end; ++itr
      )
      {
        // insert into list only if not a duplicate
        if( std::find( list.Begin(), list.End(), *itr) == list.End())
        {
          list.PushBack( *itr);
        }
      }

      // end
      return list;
    }

    //! @brief construct a protocol list from the given vector of strings making sure there are no duplicates
    //! @param STRING_VECTOR Vector of strings
    //! @return list of unique protocols
    storage::List< Protocol> StageFactory::ConstructProtocolList( const storage::Vector< std::string> &STRING_VECTOR)
    {
      // construct list
      storage::List< Protocol> list;

      // iterate over string vector
      for
      (
        storage::Vector< std::string>::const_iterator itr( STRING_VECTOR.Begin()), itr_end( STRING_VECTOR.End());
        itr != itr_end; ++itr
      )
      {
        // get the corresponding protocol
        Protocol this_protocol( GetProtocols().GetEnumFromName( *itr));
        BCL_Assert( this_protocol.IsDefined(), "No protocol found with given name " + *itr);

        // insert into list only if not a duplicate
        if( std::find( list.Begin(), list.End(), this_protocol) == list.End())
        {
          list.PushBack( this_protocol);
        }
      }

      // end
      return list;
    }

    //! @brief read table from file
    //! @param FILENAME name for the file that contains table
    //! @return Table read from the file
    storage::Table< double> StageFactory::ReadTableFromFile( const std::string &FILENAME)
    {
      // initialize stream and tie it to filename
      io::IFStream read;
      io::File::MustOpenIFStream( read, FILENAME);

      // construct table and read from stream
      storage::Table< double> table;
      table.ReadFormatted( read, false);
      io::File::CloseClearFStream( read);

      // end
      return table;
    }

    //! @brief write table to file
    //! @param TABLE Table to write
    //! @param FILENAME name fo the file to which the table will be written to
    void StageFactory::WriteTableToFile( const storage::Table< double> &TABLE, const std::string &FILENAME)
    {
      // initialize stream and tie it to filename
      io::OFStream write;
      io::File::MustOpenOFStream( write, FILENAME);

      // write table
      TABLE.WriteFormatted( write);
      io::File::CloseClearFStream( write);
    }

    //! @brief function to set the scoring function for the given single stage
    //! @param SINGLE_STAGE Stage of interest
    void StageFactory::SetSingleStageScoreFunction( mc::Stage &SINGLE_STAGE)
    {
      // if a score file was given using the corresponding flag
      if( DefaultFlags::GetFlagScoreRead()->GetFlag())
      {
        // open the score file
        BCL_MessageStd
        (
          "reading scoring function from " + DefaultFlags::GetFlagScoreRead()->GetFirstParameter()->GetValue()
        );
        io::IFStream read;
        io::File::MustOpenIFStream( read, DefaultFlags::GetFlagScoreRead()->GetFirstParameter()->GetValue());

        // read the scoring function
        util::ShPtr< score::ProteinModelScoreSum> sp_score_function;
        read >> sp_score_function;
        io::File::CloseClearFStream( read);

        // set it
        SINGLE_STAGE.SetScoreFunction( sp_score_function);
      }
      // otherwise use defined scores
      else
      {
        BCL_MessageStd( "initializing scores from protocols chosen");

        // set the score weightset
        SetSingleStageScoreWeightSet( SINGLE_STAGE);

        // construct the score sum
        util::ShPtr< score::ProteinModelScoreSum> sp_score_function
        (
          SINGLE_STAGE.GetScoreWeightSet()->ConstructScoreSum()
        );
        SINGLE_STAGE.SetScoreFunction( sp_score_function);
      }

      // if score write flag was set
      if( DefaultFlags::GetFlagScoreWrite()->GetFlag())
      {
        BCL_MessageStd
        (
          "writing the scoring function to file " + DefaultFlags::GetFlagScoreWrite()->GetFirstParameter()->GetValue()
        );
        // initialize stream
        io::OFStream write;
        io::File::MustOpenOFStream( write, DefaultFlags::GetFlagScoreWrite()->GetFirstParameter()->GetValue());
        // output the stream
        write << SINGLE_STAGE.GetScoreFunction();
        io::File::CloseClearFStream( write);
      }
    }

    //! @brief function to set the scoring weightset for the given single stage
    //! @param SINGLE_STAGE Stage of interest
    void StageFactory::SetSingleStageScoreWeightSet( mc::Stage &SINGLE_STAGE)
    {
      // initialize scoring weight_set
      if( DefaultFlags::GetFlagScoreWeightSetRead()->GetFlag())
      {
        // store filename
        const std::string filename( DefaultFlags::GetFlagScoreWeightSetRead()->GetFirstParameter()->GetValue());
        // initialize read
        BCL_MessageStd( "reading score weight set from file " + filename);

        // read table from given file and construct weightset from it
        SINGLE_STAGE.SetScoreWeightSet
        (
          util::ShPtr< ScoreWeightSet>( new ScoreWeightSet( ReadTableFromFile( filename)))
        );

      }
      // if the flag was not given use the default scoring function
      else
      {
        // initialize a empty ScoreWeightSet
        util::ShPtr< ScoreWeightSet> sp_weightset( new ScoreWeightSet());

        BCL_MessageStd( "initializing score weights");
        // iterate over the score protocols again
        for
        (
          storage::List< Protocol>::const_iterator
            protocol_itr( SINGLE_STAGE.GetScoreProtocols().Begin()),
            protocol_itr_end( SINGLE_STAGE.GetScoreProtocols().End());
          protocol_itr != protocol_itr_end; ++protocol_itr
        )
        {
          BCL_MessageStd( "\t" + protocol_itr->GetName());
          // this time initialize the score weightsets
          ( **protocol_itr)->ModifyScoreWeightSet( *sp_weightset);
        }

        // set the weight set
        SINGLE_STAGE.SetScoreWeightSet( sp_weightset);
      }
    }

    //! @brief function to set the mutate for the given single stage
    //! @param SINGLE_STAGE Stage of interest
    void StageFactory::SetSingleStageMutate( mc::Stage &SINGLE_STAGE)
    {
      // if read mutate flag was set
      if( DefaultFlags::GetFlagMutateRead()->GetFlag())
      {
        // initialize empty mutate
        util::ShPtr< math::MutateInterface< assemble::ProteinModel> > sp_mutate;

        // initialize read
        BCL_MessageStd
        (
          "Reading mutates from " + DefaultFlags::GetFlagMutateRead()->GetFirstParameter()->GetValue()
        );
        io::IFStream read;
        io::File::MustOpenIFStream( read, DefaultFlags::GetFlagMutateRead()->GetFirstParameter()->GetValue()),
        // read the mutate
        read >> sp_mutate;
        io::File::CloseClearFStream( read);

        // set the mutate
        SINGLE_STAGE.SetMutate( sp_mutate);
      }
      // otherwise use the mutates already defined
      else
      {
        BCL_MessageStd( "initializing mutate tree from protocols chosen");

        // set the mutate tree
        SetSingleStageMutateTree( SINGLE_STAGE);

        // now construct the mutate and set it
        SINGLE_STAGE.SetMutate( SINGLE_STAGE.GetMutateTree()->ConstructMutates());
      }
      // if mutate write flag was provided
      if( DefaultFlags::GetFlagMutateWrite()->GetFlag())
      {
        // initialize write
        BCL_MessageStd
        (
          "writing mutates to file" + DefaultFlags::GetFlagMutateWrite()->GetFirstParameter()->GetValue()
        );
        io::OFStream write;
        io::File::MustOpenOFStream( write, DefaultFlags::GetFlagMutateWrite()->GetFirstParameter()->GetValue());

        // write the mutate
        write << SINGLE_STAGE.GetMutate();
        io::File::CloseClearFStream( write);
      }
    }

    //! @brief function to set the mutate tree for the given single stage
    //! @param SINGLE_STAGE Stage of interest
    void StageFactory::SetSingleStageMutateTree( mc::Stage &SINGLE_STAGE)
    {
      // if a weights file is given
      if( DefaultFlags::GetFlagMutateWeightSetRead()->GetFlag())
      {
        // initialize filename
        const std::string filename( DefaultFlags::GetFlagMutateWeightSetRead()->GetFirstParameter()->GetValue());

        // read table and construct mutate tree with it
        BCL_MessageStd( "reading mutate weight set from file " + filename);
        SINGLE_STAGE.SetMutateTree( util::ShPtr< MutateTree>( new MutateTree( ReadTableFromFile( filename))));
      }

      // if no weights file
      else
      {
        // initialize empty mutate weightset
        util::ShPtr< MutateTree> sp_mutate_tree( new MutateTree());

        BCL_MessageStd( "initializing mutate weightset from protocol chosen");
        // iterate over the mutate protocols again
        for
        (
          storage::List< Protocol>::const_iterator
            protocol_itr( SINGLE_STAGE.GetMutateProtocols().Begin()),
            protocol_itr_end( SINGLE_STAGE.GetMutateProtocols().End());
          protocol_itr != protocol_itr_end; ++protocol_itr
        )
        {
          BCL_MessageStd( "\t" + protocol_itr->GetName());
          // this time initialize the mutate tree
          ( **protocol_itr)->MergeAndModifyMutateTree( *sp_mutate_tree);
        }

        // set the tree
        SINGLE_STAGE.SetMutateTree( sp_mutate_tree);
      }
    }

  } // namespace fold
} // namespace bcl
