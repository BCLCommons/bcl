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

// include headers from the bcl - sorted alphabetically
#include "app/bcl_app_apps.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "chemistry/bcl_chemistry_aa_fragment_complete.h"
#include "chemistry/bcl_chemistry_conformation_comparison_psi_flex_field.h"
#include "chemistry/bcl_chemistry_fragment_evolve_base.h"
#include "chemistry/bcl_chemistry_fragment_evolve_implementations.h"
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "chemistry/bcl_chemistry_fragment_mutate_interface.h"
#include "chemistry/bcl_chemistry_fragment_react.h"
#include "chemistry/bcl_chemistry_fragment_split_gadd_fragments.h"
#include "chemistry/bcl_chemistry_molecule_evolution_info.h"
#include "chemistry/bcl_chemistry_molecule_evolutionary_optimizer.h"
#include "chemistry/bcl_chemistry_reaction_search.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "descriptor/bcl_descriptor_molecule_druglike.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_running_average.h"
#include "math/bcl_math_template_instantiations.h"
#include "pdb/bcl_pdb_factory.h"
#include "random/bcl_random_uniform_distribution.h"
#include "storage/bcl_storage_triplet.h"
#include "util/bcl_util_stopwatch.h"

// External includes - sorted alphabetically
#include <cstdio>
#include <iomanip>

namespace bcl
{
  namespace app
  {
    // typedef for convenience
    typedef chemistry::MoleculeEvolutionInfo MolInfo;
    typedef util::ShPtr< chemistry::MoleculeEvolutionInfo> MolInfo_sp;

    // Forward declarations of comparison functions
    bool operator <( const MolInfo_sp &FIRST, const MolInfo_sp &SECOND);
    bool operator >( const MolInfo_sp &FIRST, const MolInfo_sp &SECOND);
    bool operator <( const std::pair< const MolInfo *, size_t> &FIRST, const std::pair< const MolInfo *, size_t> &SECOND);
    bool operator >( const std::pair< const MolInfo *, size_t> &FIRST, const std::pair< const MolInfo *, size_t> &SECOND);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EvoGen
    //! @brief A de-novo design application for small molecules using a stochastic search algorithm
    //!
    //! @author geanesar
    //! @date 05/12/2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class EvoGen :
      public Interface
    {

    protected:

      //! native target ligands
      mutable chemistry::FragmentEnsemble m_NativeLigands;

      //! native target pockets
      mutable chemistry::FragmentEnsemble m_BindingPocket;

    private:

    //////////
    // data //
    //////////

      util::ShPtr< command::FlagInterface> m_StartingMolsFlag;
      util::ShPtr< command::FlagInterface> m_InsertMolsFlag;
      util::ShPtr< command::FlagInterface> m_ReactionsFlag;
      util::ShPtr< command::FlagInterface> m_ReagentsFlag;
      util::ShPtr< command::FlagInterface> m_ModelDirectoryFlag;
      util::ShPtr< command::FlagInterface> m_ModelPrefixFlag;
      util::ShPtr< command::FlagInterface> m_ExtScoringCommandFlag;
      util::ShPtr< command::FlagInterface> m_DruglikenessModelFlag;
      util::ShPtr< command::FlagInterface> m_WeightScoreFlag;
      util::ShPtr< command::FlagInterface> m_OutputDirFlag;
      util::ShPtr< command::FlagInterface> m_OutputPrefixFlag;
      util::ShPtr< command::FlagInterface> m_PopulationSizeFlag;
      util::ShPtr< command::FlagInterface> m_OversampleFlag;
      util::ShPtr< command::FlagInterface> m_ShuffleFlag;
      util::ShPtr< command::FlagInterface> m_IterationsFlag;
      util::ShPtr< command::FlagInterface> m_ActiveCutoffFlag;
      util::ShPtr< command::FlagInterface> m_SelectionTypeFlag;
      util::ShPtr< command::FlagInterface> m_EvolutionBalanceTypeFlag;
      util::ShPtr< command::FlagInterface> m_ReplaceTournSizeFlag;
      util::ShPtr< command::FlagInterface> m_ModifyTournSizeFlag;
      util::ShPtr< command::FlagInterface> m_RetirementTypeFlag;
      util::ShPtr< command::FlagInterface> m_ImplementationFlag;

    public:

    //////////////////////////////////
    // Construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      EvoGen();

      //! @brief clone function
      //! @return a pointer to a copy of this class
      EvoGen *Clone() const
      {
        return new EvoGen( *this);
      }

      //! @brief Get the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief initializes the command object for that executable
      //! @return a ShPtr containing the commands
      util::ShPtr< command::Command> InitializeCommand() const
      {
        // initialize command to be returned
        util::ShPtr< command::Command> sp_cmd( new command::Command());

        // insert all the flags and params
        sp_cmd->AddFlag( m_StartingMolsFlag);
        sp_cmd->AddFlag( m_InsertMolsFlag);
        sp_cmd->AddFlag( m_ReactionsFlag);
        sp_cmd->AddFlag( m_ReagentsFlag);
        sp_cmd->AddFlag( m_ModelDirectoryFlag);
        sp_cmd->AddFlag( m_ModelPrefixFlag);
        sp_cmd->AddFlag( m_ExtScoringCommandFlag);
        sp_cmd->AddFlag( m_DruglikenessModelFlag);
        sp_cmd->AddFlag( m_WeightScoreFlag);
        sp_cmd->AddFlag( m_OutputDirFlag);
        sp_cmd->AddFlag( m_OutputPrefixFlag);
        sp_cmd->AddFlag( m_PopulationSizeFlag);
        sp_cmd->AddFlag( m_OversampleFlag);
        sp_cmd->AddFlag( m_ShuffleFlag);
        sp_cmd->AddFlag( m_IterationsFlag);
        sp_cmd->AddFlag( m_ActiveCutoffFlag);
        sp_cmd->AddFlag( m_SelectionTypeFlag);
        sp_cmd->AddFlag( m_EvolutionBalanceTypeFlag);
        sp_cmd->AddFlag( m_ReplaceTournSizeFlag);
        sp_cmd->AddFlag( m_ModifyTournSizeFlag);
        sp_cmd->AddFlag( m_RetirementTypeFlag);
        sp_cmd->AddFlag( m_ImplementationFlag);

        // default flags
        command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

        // return assembled Command object
        return sp_cmd;
      }

      std::string SetupDescriptorString() const
      {

        std::string descriptor_str( "Combine(Multiply(");

        std::string model_prefix( m_ModelPrefixFlag->GetFlag() ? m_ModelPrefixFlag->GetFirstParameter()->GetValue() : "model");
        descriptor_str += "PredictionMean(storage=File(directory="
                          + m_ModelDirectoryFlag->GetFirstParameter()->GetValue()
                          + ",prefix=" + model_prefix
                          +"))";
        if( m_WeightScoreFlag->GetFlag())
        {
          descriptor_str += ",WithinRangeSmooth(descriptor=Weight,begin=150,end=550,left width=50,right width=50)";
        }
        descriptor_str += "))";
        return descriptor_str;
      }

      void SetupInitialPopulation( chemistry::FragmentEnsemble &INIT_POP, const storage::Vector< chemistry::FragmentComplete> &MOLECULES, size_t POP_SIZE) const
      {
        for( size_t mol_num( 0), added( 0); mol_num < MOLECULES.GetSize() && added < POP_SIZE; ++mol_num)
        {
          if( MOLECULES( mol_num).GetNumberAtoms() > 0)
          {
            INIT_POP.PushBack( MOLECULES( mol_num));
            ++added;
          }
        }
      }

      //! @brief the main part of the application
      //! @return 0 on success
      int Main() const
      {

      /////////////////////////
      // Initial population  //
      /////////////////////////

        BCL_MessageStd( "Reading in molecules for starting population");

        //! Molecules that will be put into initial population
        storage::Vector< chemistry::FragmentComplete> molecules;

        //! @brief the population size as requested by user
        size_t requested_pop_size( m_PopulationSizeFlag->GetFirstParameter()->GetNumericalValue< size_t>());

        // If shuffle is specified, read in everything.  Otherwise only read in however many the user wants in a population
        size_t n_to_read( m_ShuffleFlag->GetFlag() ? std::numeric_limits< size_t>::max() : requested_pop_size);

        // add indices so that the (alex stopped his comment here, indicating that there is no good reason to add indices - Ben)
        chemistry::FragmentFeed input_mols( m_StartingMolsFlag->GetStringList(), sdf::e_Saturate, n_to_read);
        for( size_t mol_no( 0); input_mols.NotAtEnd(); ++input_mols, ++mol_no)
        {
          if( input_mols->GetNumberAtoms() > 0)
          {
            chemistry::FragmentComplete next_mol( *input_mols);
            if( next_mol.IsPropertyStored( "EvoGenIdentifier"))
            {
              next_mol.RemoveProperty( "EvoGenIdentifier");
            }
            next_mol.StoreProperty( "EvoGenIdentifier", "P" + util::Format()( 0) + "M" + util::Format()( mol_no));
            molecules.PushBack( next_mol);
          }
        }

        BCL_MessageStd( "Read in a total of " + util::Format()( molecules.GetSize()) + " molecules");

        // If shuffling is specified, randomize the molecules
        if( m_ShuffleFlag->GetFlag())
        {
          BCL_MessageStd( "Shuffling molecules");
          molecules.Shuffle();
        }

        //! @brief the population size that will be used
        size_t set_pop_size = requested_pop_size;

        // Size of the initial population
        // Warn if there are not enough molecules to fill the initial population to user's request
        if( molecules.GetSize() < set_pop_size)
        {
          BCL_MessageCrt
          (
            "WARNING: population size specified as " + util::Format()( requested_pop_size)
            + " but input files only contain " + util::Format()( molecules.GetSize()) + "."
          );
        }

        chemistry::FragmentEnsemble init_pop;

        SetupInitialPopulation( init_pop, molecules, set_pop_size);

        BCL_MessageStd( "Initial population contains " + util::Format()( init_pop.GetSize()) + " members");

        // Write the status file to wherever the output filename will be written
        std::string output_dir( m_OutputDirFlag->GetFlag() ? m_OutputDirFlag->GetFirstParameter()->GetValue() : ".");

        std::string prefix( output_dir + "/");
        std::string out_prefix( m_OutputPrefixFlag->GetFirstParameter()->GetValue());
        std::string final_output_fn( prefix + out_prefix + "_final.sdf.gz");
        std::string checkpoint_filename( prefix + out_prefix + "_checkpoint.sdf.gz");

        std::string details_prefix( prefix + out_prefix + "_details");
        std::string member_history_filename = details_prefix + "_histories.txt";
        std::string member_details_filename = details_prefix + "_member_details.csv";
        std::string population_stats_filename = details_prefix + "_pop_stats.csv";
        std::string active_molecules_filename = details_prefix + "_active.sdf";

        std::string pop_sdf_prefix = details_prefix + "_sdf";

        util::ShPtr< chemistry::FragmentEnsemble> active_mols( new chemistry::FragmentEnsemble());

        BCL_MessageStd( "");
        BCL_MessageStd( "=============================");
        BCL_MessageStd( "RUN DETAILS:");
        BCL_MessageStd( "  Output file: --------------------- " + final_output_fn)
        BCL_MessageStd( "  Checkpoint file: ----------------- " + checkpoint_filename);

        if( !pop_sdf_prefix.empty())
        {
          BCL_MessageStd( "    (details file) Population SDFs:  " + pop_sdf_prefix + "_pop_X.sdf.gz");
        }

        BCL_MessageStd( "=============================\n");

        util::Stopwatch timer
        (
          "Molecule generation",
          util::Message::e_Standard,
          false
        );

        timer.Start();

        // Execute the evolutionary algorithm
        chemistry::MoleculeEvolutionaryOptimizer optimizer;

        // set up the population filter/selector
        std::string sel_type( m_SelectionTypeFlag->GetFirstParameter()->GetValue());
        if( sel_type == "top")
        {
          optimizer.SetSelectionTop();
        }
        else if( sel_type == "tournament")
        {
          optimizer.SetReplacementTypeTournament( m_ReplaceTournSizeFlag->GetFirstParameter()->GetNumericalValue< float>());
          optimizer.SetModifyTypeTournament( m_ModifyTournSizeFlag->GetFirstParameter()->GetNumericalValue< float>());
        }
        else
        {
          BCL_MessageStd( "Unknown selection type \"" + m_SelectionTypeFlag->GetFirstParameter()->GetValue() + "\" specified");
          return -1;
        }

        // set up the evolution balance type
        std::string evbal_type( m_EvolutionBalanceTypeFlag->GetFirstParameter()->GetValue());
        if( evbal_type == "alchemical_mutate")
        {
          optimizer.SetEvolutionBalanceAlchemicalMutate();
        }
        else if( evbal_type == "reaction_dominant")
        {
          optimizer.SetEvolutionBalanceReactionDominant();
        }
        else if( evbal_type == "reaction_insertion_only")
        {
          optimizer.SetEvolutionBalanceReactionInsertionOnly();
        }
        else if( evbal_type == "recombination_dominant")
        {
          optimizer.SetEvolutionBalanceRecombinationDominant();
        }
        else if( evbal_type == "recombination_insertion_only")
        {
          optimizer.SetEvolutionBalanceRecombinationInsertionOnly();
        }
        else if( evbal_type == "balanced")
        {
          optimizer.SetEvolutionBalancedBalanced();
        }
        else
        {
            BCL_MessageStd( "Unknown evolution balance type \"" + m_EvolutionBalanceTypeFlag->GetFirstParameter()->GetValue() + "\" specified");
            return -1;
        }

        // determine retirement type
        std::string ret_type( m_RetirementTypeFlag->GetFirstParameter()->GetValue());
        if( ret_type == "all")
        {
          optimizer.SetRetirementTypeAll();
        }
        else if( ret_type == "probabilistic")
        {
          optimizer.SetRetirementTypeProbabilistic();
        }
        else if( ret_type == "none")
        {
          optimizer.SetRetirementTypeNone();
        }
        else
        {
          BCL_Exit( "Retirement type not recognized", -1);
        }

        // determine population sizes and sampling factors
        float oversample_factor( m_OversampleFlag->GetFirstParameter()->GetNumericalValue< float>());
        optimizer.SetFinalPopSize( set_pop_size);
        optimizer.SetMaxToGenerate( set_pop_size * oversample_factor);

        // set up scoring functions
        BCL_Assert
        (
          ( m_ExtScoringCommandFlag->GetFlag() || m_ModelDirectoryFlag->GetFlag())
            && !( m_ExtScoringCommandFlag->GetFlag() && m_ModelDirectoryFlag->GetFlag()),
          "Exactly one of " + m_ExtScoringCommandFlag->GetName() + " or " + m_ModelDirectoryFlag->GetName() + " must be specified"
        );

        // get the model directory
        if( m_ModelDirectoryFlag->GetFlag())
        {
          // set up fitness function
          std::string descriptor_str = SetupDescriptorString();
          BCL_MessageStd( "Score descriptor: " + descriptor_str);
          optimizer.SetModelTypeInternal();
          optimizer.SetModelDescriptor( descriptor_str);
        }
        else if( m_ExtScoringCommandFlag->GetFlag())
        {
          optimizer.SetModelTypeExternal();
          optimizer.SetModelCmd( m_ExtScoringCommandFlag->GetFirstParameter()->GetValue());
        }

        // set up reaction/insertion operations
        optimizer.SetupReactOperation( m_ReagentsFlag->GetFirstParameter()->GetValue(), m_ReactionsFlag->GetFirstParameter()->GetValue());
        optimizer.SetupInsertOperation( m_InsertMolsFlag->GetFirstParameter()->GetValue());
        optimizer.SetupAlchemicalMutate( m_ImplementationFlag->GetFirstParameter()->GetValue());

        // consider adding flag later to allow recombination
        // with a special subset of fragments

        // set the initial population
        optimizer.SetInitialPopulation( init_pop);

        // determine the number of iterations to run for
        size_t n_iters( m_IterationsFlag->GetFirstParameter()->GetNumericalValue< size_t>());

        // write out initial population scores before we begin
        io::OFStream out;
        io::OFStream out_csv;
        io::File::MustOpenOFStream( out_csv, pop_sdf_prefix + "_all_scores.csv");

        out_csv << "pop_no,fitness\n";

        { // this is here on purpose so that we can re-use variable names later
          std::string out_filename( pop_sdf_prefix + "_0.sdf.gz");
          const std::vector< MolInfo> &mols( optimizer.GetMoleculeEvolutionInfos().back());

          io::File::MustOpenOFStream( out, out_filename);
          for( size_t i( 0), end_i( mols.size()); i < end_i; ++i)
          {
            chemistry::FragmentComplete mol( mols[ i].GetMolecule());
            if( mol.IsPropertyStored( "EvoGenFitness"))
            {
              mol.RemoveProperty( "EvoGenFitness");
            }
            if( mol.IsPropertyStored( "EvoGenHistory"))
            {
              mol.RemoveProperty( "EvoGenHistory");
            }
            if( mol.IsPropertyStored( "EvoGenIdentifier"))
            {
              mol.RemoveProperty( "EvoGenIdentifier");
            }
            mol.StoreProperty( "EvoGenFitness", util::Format()( mols[ i].GetMoleculeFitness()));
            mol.StoreProperty( "EvoGenHistory", util::Format()( mols[ i].GetMoleculeHistory()));
            mol.StoreProperty( "EvoGenIdentifier", util::Format()( mols[ i].GetMoleculeIdentifier()));
            mol.WriteMDL( out);
            out_csv << "0," << mols[ i].GetMoleculeFitness() << "\n";
          }
          out_csv.flush();
          io::File::CloseClearFStream( out);
        }

        // open the activity log
        std::string activity_log( pop_sdf_prefix + "_activity_log.json");
        optimizer.SetLogFile( activity_log);
//        optimizer.StartLogging();

        // execute the optimization until we have iterated enough times
        chemistry::FragmentEnsemble mols_across_generations;
        chemistry::ConstitutionSet unique_mols;
        while( optimizer.GetMoleculeEvolutionInfos().size() < n_iters)
        {
          // population ID number
          size_t pop_no( optimizer.GetMoleculeEvolutionInfos().size());

          // check for errors
          if( optimizer.Next() == -1)
          {
            BCL_MessageStd( "Could not optimize population " + util::Format()( optimizer.GetMoleculeEvolutionInfos().size()) + ", optimizer returned error");
            break;
          }

          // open up an sdf file to write molecule info to
          std::string out_filename( pop_sdf_prefix + "_" + util::Format()( pop_no) + ".sdf.gz");
          const std::vector< MolInfo> &mols( optimizer.GetMoleculeEvolutionInfos().back());

          // write all molecules, storing necessary properties
          io::File::MustOpenOFStream( out, out_filename);
          for( size_t i( 0), end_i( mols.size()); i < end_i; ++i)
          {
            chemistry::FragmentComplete mol( mols[ i].GetMolecule());
            if( mol.IsPropertyStored( "EvoGenFitness"))
            {
              mol.RemoveProperty( "EvoGenFitness");
            }
            if( mol.IsPropertyStored( "EvoGenHistory"))
            {
              mol.RemoveProperty( "EvoGenHistory");
            }
            if( mol.IsPropertyStored( "EvoGenIdentifier"))
            {
              mol.RemoveProperty( "EvoGenIdentifier");
            }
            mol.StoreProperty( "EvoGenFitness", util::Format()( mols[ i].GetMoleculeFitness()));
            mol.StoreProperty( "EvoGenHistory", util::Format()( mols[ i].GetMoleculeHistory()));
            mol.StoreProperty( "EvoGenIdentifier", util::Format()( mols[ i].GetMoleculeIdentifier()));
            mol.WriteMDL( out);
            out_csv << pop_no << "," << mols[ i].GetMoleculeFitness() << "\n";

            // Save current generation to all, but make sure the final selection does not include duplicates carried across generations
            if( unique_mols.Insert( chemistry::FragmentConstitutionShared( mols[ i].GetMolecule())).second)
            {
              mols_across_generations.PushBack( mol);
            }
          }

          // write info to the all_scores csv without closing the stream
          out_csv.flush();

          // close the molecule output file
          io::File::CloseClearFStream( out);

          // display message to user
          BCL_MessageStd
          (
            " Opti has generated " + util::Format()( optimizer.GetMoleculeEvolutionInfos().size()) + " populations with "
            + util::Format()( optimizer.GetMoleculeEvolutionInfos().back().size())
          );
        }

        // close down the json log file
//        optimizer.StopLogging();
        io::File::CloseClearFStream( out_csv);

        //Sort by fitness and save the top N
        size_t n_best( m_PopulationSizeFlag->GetFirstParameter()->GetNumericalValue< size_t>()), itr_index( 0);
        mols_across_generations.Sort( "EvoGenFitness"); // I cannot figure out how to reverse iterate with a FragmentEnsemble
        storage::List< chemistry::FragmentComplete> sorted_mols( mols_across_generations.GetMolecules());
        io::File::MustOpenOFStream( out, final_output_fn);
        BCL_MessageStd( "Collecting best molecules...");
        for( storage::List< chemistry::FragmentComplete>::reverse_iterator itr( sorted_mols.ReverseBegin()); itr != sorted_mols.ReverseEnd(); ++itr, ++itr_index)
        {
          if( itr_index < n_best)
          {
            itr->WriteMDL( out);
          }
        }

        timer.Stop();
        timer.WriteMessage();

        BCL_MessageStd( "Finished evolving molecules");

        return 0;
      } // Main()

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    private:

      // a static instance of this class
      static const ApplicationType EvoGen_Instance;

    }; // EvoGen

    //! @brief standard constructor
    EvoGen::EvoGen() :
      m_StartingMolsFlag
      (
        new command::FlagStatic
        (
          "starting_mols", "molecules to begin the run",
          command::Parameter
          (
            "filename", "sdf file of starting molecules",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_InsertMolsFlag
      (
        new command::FlagStatic
        (
          "insert_mols", "molecules to insert into the population during the run",
          command::Parameter
          (
            "filename", "sdf file of insertion molecules",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_ReactionsFlag
      (
        new command::FlagStatic
        (
          "reactions", "reactions to use for structure modification",
          command::Parameter
          (
            "filename", "rxn file containing desired reactions",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_ReagentsFlag
      (
        new command::FlagStatic
        (
          "reagents", "reagents to use in reactions",
          command::Parameter
          (
            "filename", "sdf file of reagnet molecules",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_ModelDirectoryFlag
      (
        new command::FlagDynamic
        (
          "model_dir", "QSAR model directory",
          command::Parameter
          (
            "directory", "directory containing models that could be used in a Prediction() descriptor",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_ModelPrefixFlag
      (
        new command::FlagDynamic
        (
          "model_prefix", "model prefix",
          command::Parameter
          (
            "prefix", "the prefix for models in the directory specified by -model_dir"
          )
        )
      ),
      m_ExtScoringCommandFlag
      (
        new command::FlagDynamic
        (
          "score_command", "the command to use for external scoring",
          command::Parameter
          (
            "command",
            "the command to use; should return 0 on success, with "
              "usage: program <input sdf> <output sdf>; output sdf should have a field \"EvoGenFitness\" where scores can be read"
          )
        )
      ),
      m_DruglikenessModelFlag
      (
        new command::FlagDynamic
        (
          "druglikeness_descriptor", "a descriptor that provides binary 1/0 output if a molecule is druglike or not",
          command::Parameter
          (
            "descriptor", "A descriptor that returns 1 if a molecule is druglike and 0 otherwise"
          )
        )
      ),
      m_WeightScoreFlag
      (
        new command::FlagDynamic
        (
          "use_weight_term", "whether to use a molecular weight scoring term to restrict molecule size",
          command::Parameter
          (
            "use term", "set this flag to use a weight term"
          )
        )
      ),
      m_OutputDirFlag
      (
        new command::FlagStatic
        (
          "output_dir", "the output directory for the run",
          command::Parameter
          (
            "directory", "the directory for run output",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_OutputPrefixFlag
      (
        new command::FlagStatic
        (
          "output_prefix",
          "what output files should be prefixed with",
          command::Parameter
          (
            "prefix", "identifier for the run",
            "evogen"
          )
        )
      ),
      m_PopulationSizeFlag
      (
        new command::FlagStatic
        (
          "population_size", "size of the population to use",
          command::Parameter
          (
            "size", "number of molecules present in each iteration of the population",
            command::ParameterCheckRanged< size_t>( 0, std::numeric_limits< size_t>::max()),
            "50"
          )
        )
      ),
      m_OversampleFlag
      (
        new command::FlagStatic
        (
          "oversample_factor",
          "the oversampling factor during the run.  the run will generate factor*pop size members during the run before performing replacement",
          command::Parameter
          (
            "factor", "oversampling factor",
            command::ParameterCheckRanged< float>( 1.0, std::numeric_limits< float>::max()),
            "10"
          )
        )
      ),
      m_ShuffleFlag
      (
        new command::FlagDynamic
        (
          "shuffle",
          "whether to randomize initial population",
          command::Parameter
          (
            "shuffle",
            "if set, the initial population will be chosen randomly from the input molecules"
          )
        )
      ),
      m_IterationsFlag
      (
        new command::FlagStatic
        (
          "iterations",
          "the criterion to use to stop the algorithm",
          command::Parameter
          (
            "iterations",
            "the number of iterations to run",
            command::ParameterCheckRanged< size_t>( 1, std::numeric_limits< size_t>::max()),
            "25"
          )
        )
      ),
      m_ActiveCutoffFlag
      (
        new command::FlagDynamic
        (
          "active_cutoff", "the model score value that differentiates actives from inactives",
          command::Parameter
          (
            "cutoff", "the model score cutoff",
            command::ParameterCheckRanged< float>( std::numeric_limits< float>::min(), std::numeric_limits< float>::max())
          )
        )
      ),
      m_SelectionTypeFlag
      (
        new command::FlagDynamic
        (
          "selection_type", "the type of molecule selection to perform",
          command::Parameter
          (
            "selection_type", "selection type",
            command::ParameterCheckAllowed( storage::Vector< std::string>::Create( "tournament", "top")),
            "tournament"
          )
        )
      ),
      m_EvolutionBalanceTypeFlag
      (
        new command::FlagDynamic
        (
          "evolution_balance_type", "the type of balance applied to molecule evolution processes",
          command::Parameter
          (
            "evolution_balance_type", "evolution balance type",
            command::ParameterCheckAllowed( storage::Vector< std::string>::Create( "alchemical_mutate", "reaction_dominant", "reaction_insertion_only", "recombination_dominant", "recombination_insertion_only", "balanced")),
            "reaction_insertion_only"
          )
        )
      ),
      m_ReplaceTournSizeFlag
      (
        new command::FlagDynamic
        (
          "rep_tourn_factor", "the fraction of molecules to consider for replacement tournament selection steps",
          command::Parameter
          (
            "size_factor", "the size factor to use",
            command::ParameterCheckRanged< float>( 0.0, 1.0),
            "0.2"
          )
        )
      ),
      m_ModifyTournSizeFlag
      (
        new command::FlagDynamic
        (
          "mod_tourn_factor", "the fraction of molecules to consider for modification tournament selection steps",
          command::Parameter
          (
            "size_factor", "the size factor to use",
            command::ParameterCheckRanged< float>( 0.0, 1.0),
            "0.2"
          )
        )
      ),
      m_RetirementTypeFlag
      (
        new command::FlagDynamic
        (
          "retirement_type", "the fraction of molecules to consider for modification tournament selection steps",
          command::Parameter
          (
            "type", "the type of retirement to use",
            command::ParameterCheckAllowed( storage::Vector< std::string>::Create( "all", "none", "probabilistic")),
            "none"
          )
        )
      ),
      m_ImplementationFlag
      (
        new command::FlagStatic
        (
          "implementation",
          "method to mutate molecules",
          command::Parameter
          (
            "mutate",
            "",
            command::ParameterCheckSerializable
            (
              util::Implementation< chemistry::FragmentMutateInterface>()
            )
          )
        )
      )
    {
    }

    //! @brief less-than operator for ShPtrs to MolInfos
    //! @return true if left operand is less than right operand
    bool operator <( const MolInfo_sp &FIRST, const MolInfo_sp &SECOND)
    {
      BCL_Assert( FIRST.IsDefined() && SECOND.IsDefined(), "comparison would have dereferenced null pointer");
      return *FIRST < *SECOND;
    }

    //! @brief greater-than operator for ShPtrs to MolInfos
    //! @return true if left operand is greater than right operand
    bool operator >( const MolInfo_sp &FIRST, const MolInfo_sp &SECOND)
    {
      BCL_Assert( FIRST.IsDefined() && SECOND.IsDefined(), "comparison would have dereferenced null pointer");
      return *FIRST > *SECOND;
    }

    //! @brief less-than operator for pair of MolInfo* and size_t
    //! @return true if left operand is less than right operand
    bool operator <( const std::pair< const MolInfo *, size_t> &FIRST, const std::pair< const MolInfo *, size_t> &SECOND)
    {
      return FIRST.first && SECOND.first ? *FIRST.first < *SECOND.first : 0;
    }

    //! @brief greater-than operator for pair of MolInfo* and size_t
    //! @return true if left operand is greater than right operand
    bool operator >( const std::pair< const MolInfo *, size_t> &FIRST, const std::pair< const MolInfo *, size_t> &SECOND)
    {
      return FIRST.first && SECOND.first ? *FIRST.first > *SECOND.first : 0;
    }

    // static instance of this class used for adding to the command line
    const ApplicationType EvoGen::EvoGen_Instance
    (
      GetAppGroups().AddAppToGroup( new EvoGen(), GetAppGroups().e_ChemInfo)
    );

  } // namespace app
} // namespace bcl
