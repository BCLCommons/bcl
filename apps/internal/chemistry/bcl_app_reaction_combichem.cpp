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
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "chemistry/bcl_chemistry_configurational_bond_types.h"
#include "chemistry/bcl_chemistry_conformation_comparison_psi_field.h"
#include "chemistry/bcl_chemistry_fragment_map_conformer.h"
#include "chemistry/bcl_chemistry_fragment_react.h"
#include "chemistry/bcl_chemistry_fragment_stochastic_pose_optimizer.h"
#include "chemistry/bcl_chemistry_merge_fragment_complete.h"
#include "chemistry/bcl_chemistry_reaction_search.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "chemistry/bcl_chemistry_sample_conformations.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "command/bcl_command_parameter_check_default.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_serialize.h"
#include "sched/bcl_sched_job_interface.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_unary_function_job_with_data.h"
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "bcl_app_reaction_combichem.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_quality.h"
#include "assemble/bcl_assemble_voxel_grid_atom.h"
#include "chemistry/bcl_chemistry_configuration_interface.h"
#include "chemistry/bcl_chemistry_configuration_set.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_configuration_shared.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_hydrogens_handler.h"
#include "chemistry/bcl_chemistry_reaction_complete.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "command/bcl_command_parameter_check_extensions_file_existence.h"
#include "graph/bcl_graph_common_subgraph_isomorphism.h"
#include "graph/bcl_graph_subgraph_isomorphism.h"
#include "io/bcl_io_file.h"
#include "iterate/bcl_iterate_generic.h"
#include "math/bcl_math_mutate_result.h"
#include "storage/bcl_storage_template_instantiations.h"
namespace bcl
{
  namespace app
  {

    // Static instance initialization
    const ApplicationType ReactionCombichem::ReactionCombichem_Instance
    (
      GetAppGroups().AddAppToGroup( new ReactionCombichem(), GetAppGroups().e_ChemInfo)
    );
    chemistry::ConformationComparisonPsiField ReactionCombichem::s_Aligner = chemistry::ConformationComparisonPsiField();

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief standard constructor
    ReactionCombichem::ReactionCombichem() :
      m_StartingMoleculesFlag
      (
        new command::FlagStatic
        (
          "starting_molecules",
          "the ensemble of molecules that will be reacted with the reagents",
          command::Parameter
          (
            "input_file",
            "pairwise reaction of each input molecule with reagent molecules",
            ""
          )
        )
      ),
      m_OutputFilenameFlag
      (
        new command::FlagStatic
        (
          "output_filename",
          "product molecules",
          command::Parameter
          (
            "output_filename",
            "SDF containing linked molecules",
            "reaction_products.sdf.gz"
          )
        )
      ),
      m_LigandLocalDockFlag
      (
        new command::FlagDynamic
        (
          "ligand_dock_score",
          "the protein-ligand interface scoring function to use for pose refinement",
          command::Parameter
          (
            "function",
            "the scoring function implementation to use"
          ),
          0,
          1
        )
      ),
      m_ReactionsFlag
      (
        new command::FlagDynamic
        (
          "reactions", "reactions to use for structure modification",
          command::Parameter
          (
            "directory",
            "path to directory containing rxn files with desired reactions"
          ),
          0,
          1
        )
      ),
      m_ReagentsFlag
      (
        new command::FlagStatic
        (
          "reagents", "reagents to use in reactions",
          command::Parameter
          (
            "filename", "sdf file of reagent molecules",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_NRoundsFlag
      (
        new command::FlagStatic
        (
          "rounds", "the number of rounds of combinatorial chemistry to perform",
          command::Parameter
          (
            "number of rounds",
            "for every round after round 1, products from the previous round will be "
            "used to begin a new round of reactions with the specified reagents/reactions",
            command::ParameterCheckRanged< size_t>( 1, math::GetHighestBoundedValue< size_t>()),
            "1"
          )
        )
      ),
      m_SaveAllRoundsFlag
      (
        new command::FlagStatic
        (
          "save_products_all_rounds",
          "write all products generated across all rounds to the output SDF; "
          "default behavior is to only save the final round of products (the intent is "
          "to mimic a reaction reaching equilibrium)"
        )
      ),
      m_LimitOneRxnPerRoundFlag
      (
        new command::FlagStatic
        (
          "limit_one_reaction_per_round",
          "limit the number of reactions that can occur each round; "
          "if multiple rounds are used then a unique reaction is used each round; "
          "if the number of rounds exceeds the number of unique reactions then the run is terminated"
        )
      ),
      m_MDLStringFlag
      (
        new command::FlagDynamic
        (
          "MDL_property",
          "the MDL property specifying the protein binding pocket used for pose-dependent scoring",
          command::Parameter
          (
            "property_name",
            "the name of the property"
          ),
          0,
          1
        )
      ),
      m_DrugLikenessTypeFlag
      (
        new command::FlagStatic
        (
          "druglikeness_type", "the type of druglikeness filter to use to determine when a molecule is skipped by the Monte Carlo algorithm",
          command::Parameter
          (
            "type", "the type of druglikenes to use",
            command::ParameterCheckAllowed( storage::Vector< std::string>::Create( "IsConstitutionDruglike", "IsConstitutionDruglikeAndHitlike", "None")),
            "IsConstitutionDruglike"
          )
        )
      ),
      m_Start
      (
        new command::FlagStatic
        (
          "starting_molecules_start",
          "flag for indicating which mol from ensemble a to load in first",
          command::Parameter
          (
            "index",
            "index of ensemble a molecules to start with",
            command::ParameterCheckRanged< size_t>( 0, math::GetHighestBoundedValue< size_t>()),
            "0"
          )
        )
      ),
      m_MaxMols
      (
        new command::FlagStatic
        (
          "starting_molecules_max",
          "flag for indicating maximum number of molecules to take from ensemble a",
          command::Parameter
          (
            "max",
            "max number of molecules from ensemble a",
            command::ParameterCheckRanged< size_t>( 0, math::GetHighestBoundedValue< size_t>()),
            util::Format()( math::GetHighestBoundedValue< size_t>())
          )
        )
      ),
      m_CorinaFlag
      (
        new command::FlagStatic
        (
          "corina",
          "make a system call to Corina to make the starting conformer during molecule cleaning;"
          "this means that if only 1 conformer is desired (i.e. pose-independent scoring) it will be the corina default conformer,"
          "while if multiple conformers are desired (i.e. refinement phase of pose-dependent scoring) there will be no effect; "
          "this option is meant primarily to allow backward compatibility for QSAR models generated with Corina conformers"
        )
      )
    {
    }

    //! @brief Clone function
    //! @return pointer to new FoldProtein
    ReactionCombichem *ReactionCombichem::Clone() const
    {
      return new ReactionCombichem( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ReactionCombichem::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string ReactionCombichem::GetDescription() const
    {
      return "Reaction-based de novo molecule design";
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ReactionCombichem::GetReadMe() const
    {
      static std::string s_read_me =
        "Performs random or enumerative reaction-based design on molecules with provided reactants and reagents. "
        "Supports reactions involving 1 - 4 reagents. Multicomponent enumerative reactions proceed to full completion "
        "unless it is specified that all products from each reaction stage should be saved. ";
      return s_read_me;
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> ReactionCombichem::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // molecule i/o
      sp_cmd->AddFlag( m_StartingMoleculesFlag);
      sp_cmd->AddFlag( m_Start);
      sp_cmd->AddFlag( m_MaxMols);
      sp_cmd->AddFlag( m_OutputFilenameFlag);
      sp_cmd->AddFlag( m_ReagentsFlag);
      sp_cmd->AddFlag( m_ReactionsFlag);
      sp_cmd->AddFlag( m_NRoundsFlag);
      sp_cmd->AddFlag( m_LimitOneRxnPerRoundFlag);
      sp_cmd->AddFlag( m_SaveAllRoundsFlag);
      sp_cmd->AddFlag( m_CorinaFlag);

      // protein-ligand interface scorer for pose refinement
      sp_cmd->AddFlag( m_LigandLocalDockFlag);
      sp_cmd->AddFlag( m_MDLStringFlag);

      // drug-likeness filter
      sp_cmd->AddFlag( m_DrugLikenessTypeFlag);

      // add defaults
      sdf::AddMoleculeIOPrefFlags( *sp_cmd);
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled command object
      return sp_cmd;
    }

  /////////////////
  // operations  //
  /////////////////

    //! @brief the Main function
    //! @return error code - 0 for success
    int ReactionCombichem::Main() const
    {

      // now load the reagents
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_ReagentsFlag->GetFirstParameter()->GetValue());
//      m_ReactantEnsemble = chemistry::FragmentEnsemble( input, sdf::e_Maintain);
      m_ReactantEnsemble = chemistry::FragmentEnsemble( input, sdf::GetCommandLineHydrogensPref());
      io::File::CloseClearFStream( input);

      // we must have reagents
      BCL_Assert( m_ReactantEnsemble.GetSize(), "No reagents specified! Exiting...");

      // read in starting molecules filename
      if( m_StartingMoleculesFlag->GetFlag())
      {
        io::File::MustOpenIFStream( input, m_StartingMoleculesFlag->GetFirstParameter()->GetValue());

        // get starting molecules adjusted for ensemble range
        math::Range< size_t> ens_a_load_rng( size_t( 0), math::GetHighestBoundedValue< size_t>());
        ens_a_load_rng.SetMin( util::ConvertStringToNumericalValue< size_t>( m_Start->GetFirstParameter()->GetValue()));

        // set max number of molecules to load for starting molecule ensemble
        const size_t n_to_load_a( util::ConvertStringToNumericalValue< size_t>( m_MaxMols->GetFirstParameter()->GetValue()));
        ens_a_load_rng.SetMax
        (
          math::GetHighestBoundedValue< size_t>() - n_to_load_a > ens_a_load_rng.GetMin() ?
              ens_a_load_rng.GetMin() + n_to_load_a - 1 :
              math::GetHighestBoundedValue< size_t>()
        );

        // get command-line preference for hydrogen atom handling
        chemistry::FragmentEnsemble start_ensemble( input, sdf::GetCommandLineHydrogensPref(), ens_a_load_rng);

        // close input stream
        io::File::CloseClearFStream( input);

        // finalize start_ensemble
        m_StartEnsemble = storage::Vector< chemistry::FragmentComplete>( start_ensemble.Begin(), start_ensemble.End());
      }
      else
      {
        BCL_MessageStd( "No starting materials specified! Reactions will be performed using only the reagent pool");
        m_StartEnsemble = storage::Vector< chemistry::FragmentComplete>( m_ReactantEnsemble.Begin(), m_ReactantEnsemble.End());
      }

      // find the reaction files directory
      if( m_ReactionsFlag->GetFlag())
      {
        m_ReactionsDirectory = m_ReactionsFlag->GetFirstParameter()->GetValue();
      }
      else
      {
        m_ReactionsDirectory = chemistry::RotamerLibraryFile::GetRotamerFinder().FindFile( "") + ( "functional_reactions");
      }

      // find suitable reactions with the provided reagents
      m_ReactionSearch = chemistry::ReactionSearch
      (
        m_ReagentsFlag->GetFirstParameter()->GetValue(),
        m_ReactionsDirectory
      );

      // running initialize somewhere somehow is absolutely mission critical
      m_ReactionSearch.Initialize();

      // may not even need this
      m_ReactionEnsemble = chemistry::ReactionEnsemble( storage::List< chemistry::ReactionComplete>( m_ReactionSearch.GetReactions()->Begin(), m_ReactionSearch.GetReactions()->End()));
      m_Reactions = *( m_ReactionSearch.GetReactions());
      if( m_LimitOneRxnPerRoundFlag->GetFlag())
      {
        for( size_t rxn_i( 0), rxn_sz( m_Reactions.GetSize()); rxn_i < rxn_sz; ++rxn_i)
        {
          m_UnusedReactions.Insert( std::make_pair( rxn_i, m_Reactions( rxn_i)));
        }
      }

      // setup pose-dependent flags
      if( m_MDLStringFlag->GetFlag())
      {
        m_MDLString = m_MDLStringFlag->GetFirstParameter()->GetValue();
        m_PocketFilename = m_StartEnsemble.Begin()->GetStoredPropertiesNonConst().GetMDLProperty( m_MDLString);
      }
      if( m_LigandLocalDockFlag->GetFlag())
      {
        BCL_Assert
        (
          m_PocketFilename.size(),
          "Either no MDL property was specified for the protein binding pocket filename, or "
          "the specified MDL property did not contain a filename. Check your input SDF and try again."
        );
      }

      // open output stream
      io::OFStream output, output_debug_mdl;
      io::File::MustOpenOFStream( output, m_OutputFilenameFlag->GetFirstParameter()->GetValue());
      io::File::MustOpenOFStream( output_debug_mdl, "Debug_" + m_OutputFilenameFlag->GetFirstParameter()->GetValue());

      // loop over rounds of combichem
      size_t n_total_rounds( size_t( 1));
      if( m_NRoundsFlag->GetFlag())
      {
        n_total_rounds = m_NRoundsFlag->GetFirstParameter()->GetNumericalValue< size_t>();
      }

      // collect products by round; outer vector indexes products by round index
      storage::Vector< storage::Vector< chemistry::FragmentComplete> > unique_final_products( n_total_rounds);

      // perform N rounds of reactions
      for( size_t round_index( 0); round_index < n_total_rounds; ++round_index)
      {
        // get a random reaction position
        size_t rand_rxn_key;
        if( m_LimitOneRxnPerRoundFlag->GetFlag())
        {
          storage::Vector< size_t> keys( m_UnusedReactions.GetKeysAsVector());
          keys.Shuffle();
          rand_rxn_key = keys.FirstElement();
        }
        if
        (
            round_index > 0 && //< if we are actually doing more than 1 round of combichem
            m_Products.GetSize() //< if we actually made products in the previous round
        )
        {
          // add old starting materials to the reagents list
          for
          (
              auto ens_itr( m_StartEnsemble.Begin()), ens_itr_end( m_StartEnsemble.End());
              ens_itr != ens_itr_end;
              ++ens_itr
          )
          {
//            ens_itr->SaturateWithH();
            m_ReactantEnsemble.PushBack( *ens_itr);
          }
          // set reactants to prior starting molecules
//          m_ReactantEnsemble = chemistry::FragmentEnsemble( storage::List< chemistry::FragmentComplete>( m_StartEnsemble.Begin(), m_StartEnsemble.End()));

          // set starting materials to products from previous round
          m_StartEnsemble = m_Products;

          // Reset the reaction search object with the new molecules
          m_ReactionSearch.Reset();
          m_ReactionSearch = chemistry::ReactionSearch
          (
            m_ReactantEnsemble,
            m_ReactionEnsemble
          );
          m_ReactionSearch.Initialize();

        }
        else if( round_index > 0 && !m_Products.GetSize())
        {
          // save products from previous round and exit
          for
          (
              auto mol_itr( unique_final_products( round_index - 1).Begin()),
              mol_itr_end( unique_final_products( round_index - 1).End());
              mol_itr != mol_itr_end;
              ++mol_itr
          )
          {
            mol_itr->WriteMDL( output);
          }
          break;
        }

        // create the reaction operator with any updates to the ReactionSearch object
        m_ReactOp = chemistry::FragmentReact( m_ReactionSearch);

        // we will want to track all unique products for a single round
        chemistry::ConfigurationSet unique_products_configs;

        // go over each start molecule
        size_t ens_index( 0);
        for
        (
            auto ens_itr( m_StartEnsemble.Begin()), ens_itr_end( m_StartEnsemble.End());
            ens_itr != ens_itr_end;
            ++ens_itr, ++ens_index
        )
        {
          // enumerate reactions for this molecule
          storage::Pair
          <
            storage::Vector< storage::Pair< util::SiPtr< const chemistry::ReactionComplete>, chemistry::FragmentEnsemble> >,
            storage::Vector< storage::Vector< std::string> >
          > rxn_exh;

          if( m_LimitOneRxnPerRoundFlag->GetFlag())
          {
            // perform reaction
            BCL_MessageStd( "React molecule " + util::Format()( ens_index) + " with reactants!");
            rxn_exh = m_ReactOp.ReactExhaustiveOneReaction( *ens_itr, m_UnusedReactions.Find( rand_rxn_key)->second);
            BCL_MessageStd( "Done reacting molecule " + util::Format()( ens_index) + " with reactants!");
          }
          else
          {
            BCL_MessageStd( "React molecule " + util::Format()( ens_index) + " with reactants!");
            rxn_exh = m_ReactOp.ReactExhaustive( *ens_itr);
            BCL_MessageStd( "Done reacting molecule " + util::Format()( ens_index) + " with reactants!");
          }

          // we only want to store the products from this round
          m_Products.Reset();

          // iterate over the outcomes of each reaction
          size_t rxn_info_v_index( 0);
          for
          (
              auto rxn_prod_itr( rxn_exh.First().Begin()), rxn_prod_itr_end( rxn_exh.First().End());
              rxn_prod_itr != rxn_prod_itr_end;
              ++rxn_prod_itr, ++rxn_info_v_index
          )
          {
            // a reaction may have generated more than one molecule; add each of them to the new population
            size_t rxn_info_index( 0);
            for
            (
                chemistry::FragmentEnsemble::const_iterator itr_mol( rxn_prod_itr->Second().Begin()), itr_mol_end( rxn_prod_itr->Second().End());
                itr_mol != itr_mol_end;
                ++itr_mol, ++rxn_info_index
            )
            {
              m_Products.PushBack( *itr_mol);
              m_Products.LastElement().StoreProperty( "Reaction", rxn_exh.Second()( rxn_info_v_index)( rxn_info_index));
            } // end product from a single reaction of one molecule
          } // end all reaction products for one start molecule

          // get only the unique products
          chemistry::ConfigurationSet unique_product_configs;
          storage::Vector< chemistry::FragmentComplete> unique_products;
          for
          (
              auto itr( m_Products.Begin()), itr_end( m_Products.End());
              itr != itr_end;
              ++itr
          )
          {
            if( unique_product_configs.Insert( chemistry::FragmentConfigurationShared( *itr)).second)
            {
              unique_products.PushBack( *itr);
            }
          }
          // save only unique products for this starting molecule
          m_Products = unique_products;

          // we will compare configurations again after we do all the reactions of all starting molecules
          unique_products.Reset();

          // refine docked pose if relevant
          chemistry::FragmentEnsemble mols;
          if( m_LigandLocalDockFlag->GetFlag())
          {
            // get protein-ligand interaction scorer
            m_LigandDockScorer = m_LigandLocalDockFlag->GetFirstParameter()->GetValue();

            // a little non-obvious, but doing this molecule-by-molecule instead of my ensemble
            // handles a bit better, especially with storing MDL properties
            for
            (
                auto prod_itr( m_Products.Begin()), prod_itr_end( m_Products.End());
                prod_itr != prod_itr_end;
                ++prod_itr
            )
            {
              // clean molecule
              util::ShPtr< chemistry::FragmentComplete> clean_mol( CallCleaner( *prod_itr, *ens_itr));
              if( clean_mol.IsDefined() && clean_mol->GetSize())
              {
                chemistry::FragmentEnsemble clean_mol_ens;
                clean_mol_ens.PushBack( *clean_mol);
                mols = chemistry::FragmentEnsemble
                    (
                      LigandLocalDock
                      (
                        clean_mol_ens,
                        m_LigandDockScorer,
                        *ens_itr
                      )
                    );
                // output all new molecules
                for
                (
                    auto mol_itr( mols.Begin()), mol_itr_end( mols.End());
                    mol_itr != mol_itr_end;
                    ++mol_itr
                )
                {
                  if( mol_itr->GetSize())
                  {
                    mol_itr->StoreProperty( "Reaction", prod_itr->GetMDLProperty( "Reaction"));
                    mol_itr->SetName( prod_itr->GetMDLProperty( "Reaction"));
                    mol_itr->WriteMDL( output);
                  }
                }
              }
            }
          }
          // otherwise just clean and be done
          else
          {
            for
            (
                auto prod_itr( m_Products.Begin()), prod_itr_end( m_Products.End());
                prod_itr != prod_itr_end;
                ++prod_itr
            )
            {
              prod_itr->WriteMDL( output_debug_mdl);
              chemistry::AtomVector< chemistry::AtomComplete> atoms( prod_itr->GetAtomVector());
              chemistry::HydrogensHandler::Remove( atoms);
              chemistry::FragmentMapConformer cleaner
              (
                m_DrugLikenessTypeFlag->GetFirstParameter()->GetValue(),
                m_CorinaFlag->GetFlag()
              );
              util::ShPtr< chemistry::FragmentComplete> clean_mol
              (
                cleaner.Clean
                (
                  atoms,
                  *ens_itr,
                  m_DrugLikenessTypeFlag->GetFirstParameter()->GetValue(),
                  false
                )
              );
              if( clean_mol.IsDefined() && clean_mol->GetSize())
              {
                if( unique_products_configs.Insert( chemistry::FragmentConfigurationShared( *prod_itr)).second)
                {
                  // save for next round
                  // TODO: fix these labels to track total history of molecule
                  clean_mol->StoreProperty( "Reaction", prod_itr->GetMDLProperty( "Reaction"));
                  clean_mol->SetName( prod_itr->GetMDLProperty( "Reaction"));
                  unique_final_products( round_index).PushBack( *clean_mol);

                  // write products if we are collecting all
                  if( m_SaveAllRoundsFlag->GetFlag())
                  {
                    clean_mol->WriteMDL( output);
                  }
                }
              }
            }
          }
        } // end reacting all start molecules

        // for next round
        m_Products = unique_final_products( round_index);

        if( m_LimitOneRxnPerRoundFlag->GetFlag())
        {
          size_t n_deleted_keys( m_UnusedReactions.Erase( rand_rxn_key));
        }

      } // end combichem rounds

      // take final round products
      if( !m_SaveAllRoundsFlag->GetFlag())
      {
        // iterate over product sets from each round starting from the last round
        size_t ri( m_NRoundsFlag->GetFirstParameter()->GetNumericalValue< size_t>() - 1);
        for
        (
            auto prod_itr_r( unique_final_products.ReverseBegin()), prod_itr_r_end( unique_final_products.ReverseEnd());
            prod_itr_r != prod_itr_r_end;
            ++prod_itr_r, --ri
        )
        {
          // the last round to have products in it are the products we want to keep
          if( prod_itr_r->GetSize())
          {
            for
            (
                auto mol_itr( prod_itr_r->Begin()), mol_itr_end( prod_itr_r->End());
                mol_itr != mol_itr_end;
                ++mol_itr
            )
            {
              mol_itr->WriteMDL( output);
            }
            break;
          }
        }
      }

      io::File::CloseClearFStream( output);
      return 0;
    }

  //////////////////////
  // helper functions //
  /////////////////////

    //! @brief remove whitespace (via isspace) from a string
    //! @param STR the string to remove whitespace from
    //! @return STR without any whitespace
    std::string ReactionCombichem::RemoveWhitespace( const std::string &STR) const
    {
      std::string str;
      str.reserve( STR.length());
      for
      (
        std::string::const_iterator itr( STR.begin()), itr_end( STR.end());
        itr != itr_end;
        ++itr
      )
      {
        if( !isspace( *itr))
        {
          str.push_back( *itr);
        }
      }
      return str;
    }

     //! @brief wrapper that calls the FragmentMapConformer::Clean function to generate a legitimate 3D conformer and fix bond lengths
     //! @param FRAGMENT the molecule that needs to be fixed
     //! @param REFERENCE scaffold molecule for substructure alignment reference
     //! @return pointer to cleaned molecule
     util::ShPtr< chemistry::FragmentComplete> ReactionCombichem::CallCleaner
     (
       const chemistry::FragmentComplete &FRAGMENT,
       const chemistry::FragmentComplete &REFERENCE
     ) const
     {
       // hold on to this
       static chemistry::HydrogensHandler hdyrogens_handler;

       // remove hydrogen atoms
       chemistry::AtomVector< chemistry::AtomComplete> new_frag_v( FRAGMENT.GetAtomVector());
       hdyrogens_handler.Remove( new_frag_v);

       // call cleaner
       static chemistry::FragmentMapConformer cleaner;
       if( m_CorinaFlag->GetFlag())
       {
         cleaner = chemistry::FragmentMapConformer( "None", true);
       }
       else
       {
         cleaner = chemistry::FragmentMapConformer( "None", false);
       }
       util::ShPtr< chemistry::FragmentComplete> clean_frag( cleaner.Clean
         (
           new_frag_v,
           REFERENCE,
           m_DrugLikenessTypeFlag->GetFirstParameter()->GetValue()
         ));
       return clean_frag;
     }

     //! @brief wrapper that calls ResolveClashes and OptimizePose with a combined global/local conf ensemble for each mol
     //! @param ENSEMBLE the molecule that needs to be fixed
     //! @param SCORER the property to be used to score the interface
     //! @param REFERENCE the scaffold for pre-alignment
     //! @return pointer to cleaned molecule
     chemistry::FragmentEnsemble ReactionCombichem::LigandLocalDock
     (
       const chemistry::FragmentEnsemble &ENSEMBLE,
       const descriptor::CheminfoProperty &SCORER,
       const chemistry::FragmentComplete &REFERENCE
     ) const
     {
       // so we can add the properties
       chemistry::FragmentEnsemble ensemble( ENSEMBLE), best_confs;

       // we only want to make once
       static chemistry::RotamerLibraryFile rotamer_library_file;
       static chemistry::SampleConformations sample_confs
       (
         rotamer_library_file,
         "",
         0.0,    // tolerance
         250,    // number of conformations
         2000,   // number of iterations
         false,  // change chirality
         0.0,    // random dihedral change weight
         true,   // generate 3d
         0.5     // clash tolerance
       );
       static chemistry::FragmentMapConformer conf_mapper
       (
         m_DrugLikenessTypeFlag->GetFirstParameter()->GetValue(),
         m_MDLString,
         m_PocketFilename,
         SCORER,
         true,
         storage::Vector< float>(),
         false
       );
       static storage::Vector< float> bfactors( conf_mapper.GetBFactors());

       // Sample conformers of each molecule
       for
       (
           auto mol_itr( ensemble.Begin()), mol_itr_end( ensemble.End());
           mol_itr != mol_itr_end;
           ++mol_itr
       )
       {
         // store MDL for interface scoring and generate confs
         mol_itr->GetStoredPropertiesNonConst().SetMDLProperty( conf_mapper.GetMDL(), conf_mapper.GetPocketFilename());
         chemistry::FragmentEnsemble confs( sample_confs( *mol_itr).First());

         // make sure we have conformers, then align them to common substructure
         if( confs.GetSize())
         {
           for
           (
               auto conf_itr( confs.Begin()), conf_itr_end( confs.End());
               conf_itr != conf_itr_end;
               ++conf_itr
           )
           {
             // perform substructure-based alignment to reference
             s_Aligner.ArbitraryScaffoldAlignment
             (
               *conf_itr,
               REFERENCE,
               chemistry::ConformationGraphConverter::AtomComparisonType::e_AtomType,
               chemistry::ConfigurationalBondTypeData::Data::e_BondOrderAmideOrAromaticWithRingness
             );
             //          conf_itr->WriteMDL( debug_mdl_out);
           }
           //        io::File::CloseClearFStream( debug_mdl_out);
         }

         // optimize pose
         util::ShPtr< chemistry::FragmentComplete> opti_pose;
         if( confs.GetSize())
         {
           BCL_MessageStd( "Getting best scoring non-clashing conformer!");
           static chemistry::FragmentStochasticPoseOptimizer pose_optimizer
           (
             SCORER,
             bfactors,
             conf_mapper.GetMDL(),
             conf_mapper.GetPocketFilename(),
             size_t( 20),
             size_t( 100),
             size_t( 100),
             float( 1.0),
             opti::e_LargerIsBetter,
             float( 5.0)
           );
           opti_pose = pose_optimizer.StochasticPoseOptimization( confs);
         }
         else
         {
           opti_pose = util::ShPtr< chemistry::FragmentComplete>( new chemistry::FragmentComplete( *mol_itr));
           continue;
         }

         // save best pose
         if( opti_pose->GetSize())
         {
           best_confs.PushBack( *opti_pose);
         }
       }
       return best_confs;
     }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ReactionCombichem::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &ReactionCombichem::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace app
} // namespace bcl
