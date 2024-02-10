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
#include "bcl_app_conformer_from_scaffold.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "chemistry/bcl_chemistry_conformation_comparison_interface.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_align_to_scaffold.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "descriptor/bcl_descriptor_molecule_similarity.h"
#include "graph/bcl_graph_csi_substructure.h"
#include "io/bcl_io_file.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_thunk_job.h"
#include "storage/bcl_storage_template_instantiations.h"
#include "util/bcl_util_implementation.h"

namespace bcl
{
  namespace app
  {

    const ApplicationType ConformerFromScaffold::ConformerFromScaffold_Instance
    (
      GetAppGroups().AddAppToGroup( new ConformerFromScaffold(), GetAppGroups().e_Molecule)
    );

  ///////////////////////////////////
  // construction and destruction //
  ///////////////////////////////////

    //! default constructor
    //! @brief standard constructor
    ConformerFromScaffold::ConformerFromScaffold() :
        m_InputFileFlag
        (
          new command::FlagStatic
          (
            "input_filenames",
            "molecules for which to generate new conformers",
            command::Parameter
            (
              "filename",
              "SDF filename for where to write out molecules",
              command::ParameterCheckFileExistence()
            )
          )
        ),
        m_InputRangeMinFlag
        (
          new command::FlagStatic
          (
            "input_range_min",
            "Operate on a subset of input molecules beginning with this molecule index",
            command::Parameter
            (
              "index",
              "",
              command::ParameterCheckRanged< size_t>( 0, std::numeric_limits< size_t>::max()),
              "0"
            )
          )
        ),
        m_InputRangeEndFlag
        (
          new command::FlagStatic
          (
            "input_range_max",
            "Operate on a subset of input molecules ending with this molecule index",
            command::Parameter
            (
              "index",
              "",
              command::ParameterCheckRanged< size_t>( 0, std::numeric_limits< size_t>::max()),
              util::Format()( std::numeric_limits< size_t>::max())
            )
          )
        ),
        m_ScaffoldFileFlag
        (
          new command::FlagStatic
          (
            "scaffold_filenames",
            "molecules against which input molecules will be compared for conformer generation",
            command::Parameter
            (
              "filename",
              "SDF filename for where to write out molecules",
              command::ParameterCheckFileExistence()
            )
          )
        ),
        m_ScaffoldRangeMinFlag
        (
          new command::FlagStatic
          (
            "scaffold_range_min",
            "Operate on a subset of scaffold molecules beginning with this molecule index",
            command::Parameter
            (
              "index",
              "",
              command::ParameterCheckRanged< size_t>( 0, std::numeric_limits< size_t>::max()),
              "0"
            )
          )
        ),
        m_ScaffoldRangeEndFlag
        (
          new command::FlagStatic
          (
            "scaffold_range_max",
            "Operate on a subset of scaffold molecules ending with this molecule index",
            command::Parameter
            (
              "index",
              "",
              command::ParameterCheckRanged< size_t>( 0, std::numeric_limits< size_t>::max()),
              util::Format()( std::numeric_limits< size_t>::max())
            )
          )
        ),
        m_OutputFileFlag
        (
          new command::FlagStatic
          (
            "output_filename",
            "new conformers of the input molecules will be written to this file; if unsuccessful, original input molecule can be written",
            command::Parameter
            (
              "filename",
              "SDF filename for where to write out molecules",
              ""
            )
          )
        ),
        m_OutputSimilarityFailureFileFlag
        (
          new command::FlagStatic
          (
            "output_similarity_failures_filename",
            "failed molecules; occurs due to insufficient similarity to a scaffold",
            command::Parameter
            (
              "filename",
              "SDF filename for where to write out molecules",
              ""
            )
          )
        ),
        m_OutputConfGenFailureFileFlag
        (
          new command::FlagStatic
          (
            "output_confgen_failures_filename",
            "failed molecules; occurs due to invalid conformer generation",
            command::Parameter
            (
              "filename",
              "SDF filename for where to write out molecules",
              ""
            )
          )
        ),
        m_AtomTypeFlag
        (
          new command::FlagStatic
          (
            "atom_type",
            "how to compare atom types",
            command::Parameter
            (
              "scheme",
              "used to compare atoms",
              command::ParameterCheckSerializable( chemistry::ConformationGraphConverter::AtomComparisonTypeEnum()),
              "ElementType"
            )
          )
        ),
        m_BondTypeFlag
        (
          new command::FlagStatic
          (
            "bond_type",
            "how to compare bond types",
            command::Parameter
            (
              "scheme",
              "used to compare bonds",
              command::ParameterCheckSerializable( chemistry::ConfigurationalBondTypeData::DataEnum()),
              "BondOrderAmideOrAromaticWithRingness"
            )
          )
        ),
        m_MinSizeFlag
        (
          new command::FlagStatic
          (
            "min_size",
            "minimum size of the fragments to keep it in the fragments list",
            command::Parameter
            (
              "min_size",
              "minimum size of the fragments to keep it in the fragments list",
              command::ParameterCheckRanged< size_t>(),
              "3"
            )
          )
        ),
        m_SolutionTypeFlag
        (
          new command::FlagStatic
          (
            "solution_type",
            "graph search solution type",
            command::Parameter
            (
              "",
              "",
              command::ParameterCheckAllowed
              (
                storage::Vector< std::string>::Create
                (
                  "Connected",
                  "Unconnected",
                  "GreedyUnconnected"
                )
              ),
              "Connected"
            )
          )
        ),
        m_SimilarityThresholdFlag
        (
          new command::FlagStatic
          (
            "similarity_threshold",
            "minimum similarity to a scaffold required to perform conformer generation",
            util::ShPtrVector< command::ParameterInterface>::Create
            (
              util::ShPtr< command::ParameterInterface>
              (
                new command::Parameter
                (
                  "minimum",
                  "the lowest allowable similarity that a molecule is allowed to have with a scaffold in order to be a candidate for conformer generation",
                  command::ParameterCheckRanged< float>(0.0, 1.0),
                  "0.0"
                )
              ),
              util::ShPtr< command::ParameterInterface>
              (
                new command::Parameter
                (
                  "maximum",
                  "the highest allowable similarity that a molecule is allowed to have with a scaffold in order to be a candidate for conformer generation",
                  command::ParameterCheckRanged< float>(0.0, 1.0),
                  "1.0"
                )
              )
            )
          )
        ),
        m_FindAllFlag
        (
          new command::FlagStatic
          (
            "find_all",
            "if this flag is provided, then follow the largest common substructure search "
            "with a search for all subgraph isomorphisms between the common substructures; "
            "this allows alternative poses to be discovered for symmetric molecules"
          )
        ),
        m_UniqueFlag
        (
          new command::FlagStatic
          (
            "unique",
            "distance between conformers required for a molecule to be unique; only applicable if 'find_all' is enabled",
            util::ShPtrVector< command::ParameterInterface>::Create
            (
              util::ShPtr< command::ParameterInterface>
              (
                new command::Parameter
                (
                  "rmsd_type",
                  "method to use to find rmsd between conformers/poses",
                  command::ParameterCheckSerializable( util::Implementation< chemistry::ConformationComparisonInterface>()),
                  "SymmetryRealSpaceRMSD"
                )
              ),
              util::ShPtr< command::ParameterInterface>
              (
                new command::Parameter
                (
                  "tolerance",
                  "all conformers/poses must be at least this rmsd from the previously accepted pose; "
                  "prior to evaluation, all alignment poses are scored with the MolAlign score; "
                  "if you want additional local conformational heterogeneity, then run ConformerGenerator on the output of this app.",
                  command::ParameterCheckRanged< float>( 0.0, std::numeric_limits< float>::max()),
                  "1.0"
                )
              )
            )
          )
        ),
        m_SaveEnsembleFlag
        (
          new command::FlagStatic
          (
            "save_ensemble",
            "keep the ensemble of unique generated conformers; "
            "only applicable if 'find_all' is enabled; "
            "occurs after application of 'unique'"
          )
        )
    {
    }

    //! copy constructor, only copy the flags
    ConformerFromScaffold::ConformerFromScaffold( const ConformerFromScaffold &PARENT) :
          m_InputFileFlag( PARENT.m_InputFileFlag),
          m_InputRangeMinFlag( PARENT.m_InputRangeMinFlag),
          m_InputRangeEndFlag( PARENT.m_InputRangeEndFlag),
          m_ScaffoldFileFlag( PARENT.m_ScaffoldFileFlag),
          m_ScaffoldRangeMinFlag( PARENT.m_ScaffoldRangeMinFlag),
          m_ScaffoldRangeEndFlag( PARENT.m_ScaffoldRangeEndFlag),
          m_OutputFileFlag( PARENT.m_OutputFileFlag),
          m_OutputSimilarityFailureFileFlag( PARENT.m_OutputSimilarityFailureFileFlag),
          m_OutputConfGenFailureFileFlag( PARENT.m_OutputConfGenFailureFileFlag),
          m_AtomTypeFlag( PARENT.m_AtomTypeFlag),
          m_BondTypeFlag( PARENT.m_BondTypeFlag),
          m_MinSizeFlag( PARENT.m_MinSizeFlag),
          m_SolutionTypeFlag( PARENT.m_SolutionTypeFlag),
          m_SimilarityThresholdFlag( PARENT.m_SimilarityThresholdFlag),
          m_FindAllFlag( PARENT.m_FindAllFlag),
          m_UniqueFlag( PARENT.m_UniqueFlag),
          m_SaveEnsembleFlag( PARENT.m_SaveEnsembleFlag)
    {
    }

    //! @brief Clone function
    //! @return pointer to new ConformerFromScaffold
    ConformerFromScaffold *ConformerFromScaffold::Clone() const
    {
      return new ConformerFromScaffold( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ConformerFromScaffold::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> ConformerFromScaffold::InitializeCommand() const
    {
      // initialize command to be returned
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // Input flags
      sp_cmd->AddFlag( m_InputFileFlag);
      sp_cmd->AddFlag( m_InputRangeMinFlag);
      sp_cmd->AddFlag( m_InputRangeEndFlag);
      sp_cmd->AddFlag( m_ScaffoldFileFlag);
      sp_cmd->AddFlag( m_ScaffoldRangeMinFlag);
      sp_cmd->AddFlag( m_ScaffoldRangeEndFlag);
      sp_cmd->AddFlag( sdf::GetNeutralizeChargesFlag());

      // Run flags
      sp_cmd->AddFlag( m_AtomTypeFlag);
      sp_cmd->AddFlag( m_BondTypeFlag);
      sp_cmd->AddFlag( m_MinSizeFlag);
      sp_cmd->AddFlag( m_SolutionTypeFlag);
      sp_cmd->AddFlag( m_SimilarityThresholdFlag);
      sp_cmd->AddFlag( m_FindAllFlag);
      sp_cmd->AddFlag( m_UniqueFlag);

      // Output flags
      sp_cmd->AddFlag( m_OutputFileFlag);
      sp_cmd->AddFlag( m_OutputSimilarityFailureFileFlag);
      sp_cmd->AddFlag( m_OutputConfGenFailureFileFlag);
      sp_cmd->AddFlag( m_SaveEnsembleFlag);
      sp_cmd->AddFlag( sdf::GetExplicitAromaticityFlag());

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd, storage::Set< command::FlagTypeEnum>( command::e_AppGeneric));

      // return assembled Command object
      return sp_cmd;
    }

  ////////////////
  //    main    //
  ////////////////

    //! @brief the Main function
    //! @return error code - 0 for success
    int ConformerFromScaffold::Main() const
    {
      BCL_MessageStd("Reading molecules...");

      // read in ensemble without hydrogen atoms to speedup substructure search
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_InputFileFlag->GetFirstParameter()->GetValue());
      chemistry::FragmentEnsemble input_ensemble
      (
        input,
        sdf::e_Remove,
        math::Range< size_t>( m_InputRangeMinFlag->GetFirstParameter()->GetNumericalValue< size_t>(), m_InputRangeEndFlag->GetFirstParameter()->GetNumericalValue< size_t>()),
        sdf::e_CmdLine
      );
      const storage::Vector< chemistry::FragmentComplete> input_molecules( input_ensemble.Begin(), input_ensemble.End());
      io::File::CloseClearFStream( input);
      const size_t ensemble_size( input_ensemble.GetSize());
      BCL_Assert( ensemble_size, "Must have at least one molecule in the input ensemble. Exiting...\n");

      // read in the scaffold ensemble without hydrogen atoms to speedup substructure search
      io::File::MustOpenIFStream( input, m_ScaffoldFileFlag->GetFirstParameter()->GetValue());
      const chemistry::FragmentEnsemble scaffold_ensemble
      (
        input,
        sdf::e_Remove,
        math::Range< size_t>( m_ScaffoldRangeMinFlag->GetFirstParameter()->GetNumericalValue< size_t>(), m_ScaffoldRangeEndFlag->GetFirstParameter()->GetNumericalValue< size_t>()),
        sdf::e_CmdLine
      );
      const storage::Vector< chemistry::FragmentComplete> scaffold_molecules( scaffold_ensemble.Begin(), scaffold_ensemble.End());
      io::File::CloseClearFStream( input);
      const size_t scaffold_ensemble_size( scaffold_molecules.GetSize());
      BCL_Assert( scaffold_ensemble_size, "Must have at least one molecule in the scaffold ensemble. Exiting...\n");

      // initialize output so that we can write as we go
      InitializeOutputFiles();

      // Get the atom and bond type resolution for substructure comparisons
      util::Implementation< chemistry::ConformationComparisonInterface> similarity_metric( InitializeSimilarityMetric());
      graph::CommonSubgraphIsomorphismBase::SolutionType solution_type( InitializeSolutionType());
      chemistry::FragmentAlignToScaffold ats( InitializeAlignmentObject( solution_type));

      // generate conformers for each input file based on template substructure
      size_t total_confs( 0);
      BCL_MessageStd("Generating new conformer ensemble for input molecules...");
      for
      (
          auto mol_itr( input_ensemble.Begin()), mol_itr_end( input_ensemble.End());
          mol_itr != mol_itr_end;
          ++mol_itr, ++m_MoleculeIndex
      )
      {
        // status update
        if( m_MoleculeIndex % 100 == 0)
        {
          util::GetLogger().LogStatus( "Completed " + std::to_string( m_MoleculeIndex) + "/" + std::to_string( ensemble_size) + " molecules...\n");
        }

        // find the most similar molecule
        storage::Pair< size_t, float> similarity_result( FindBestScaffoldBySimilarity( *mol_itr, scaffold_molecules, similarity_metric ) );
        if( !util::IsDefined( similarity_result.First()))
        {
          continue;
        }

        // add an alignment scorer
        descriptor::MoleculeSimilarity alignment_scorer
        (
          "PropertyFieldDistance",
          chemistry::FragmentEnsemble( storage::List< chemistry::FragmentComplete>( 1, scaffold_molecules( similarity_result.First())) )
        );

        // build candidate conformers based on substructure comparisons to the most similar molecule in scaffolds list
        chemistry::FragmentEnsemble confs( Run( *mol_itr, scaffold_molecules( similarity_result.First()), ats, alignment_scorer, m_FindAllFlag->GetFlag() ) );
        bool success( confs.GetSize());
        if( !success)
        {
          if( m_OutputConfGenFailureFileFlag->GetFlag())
          {
            mol_itr->StoreProperties( input_molecules( m_MoleculeIndex));
            mol_itr->WriteMDL( m_OutputConfgenFailures);
          }

          BCL_MessageVrb
          (
            "Failed to generate conformer for molecule " + std::to_string( m_MoleculeIndex) +
            " using scaffold " + std::to_string( similarity_result.First()) +
            " with a Tanimoto similarity of " + std::to_string( similarity_result.Second())
          );
          ++m_ConfgenFailureCount;
        }
        else
        {
          // assign perfect alignment score if it is an identical molecule
//          for( auto conf_itr( confs.Begin()), conf_itr_end( confs.End()); conf_itr != conf_itr_end; ++conf_itr)
//          {
//            if( similarity_result.Second() == 1.0)
//            {
//              conf_itr->GetStoredPropertiesNonConst().SetMDLProperty("alignment_scorer.GetAlias()", linal::Vector< float>(1, 0.0));
//            }
//          }

          // lowest alignment score first
          confs.Sort( alignment_scorer.GetAlias());

          // filter to get unique poses
          if ( m_UniqueFlag->GetFlag() && confs.GetSize() > 1)
          {
            util::Implementation<chemistry::ConformationComparisonInterface> unique( m_UniqueFlag->GetFirstParameter()->GetValue());
            chemistry::FragmentEnsemble uniq_confs;

            // add the first conformer
            auto conf_itr_i( confs.Begin());
            uniq_confs.PushBack( *conf_itr_i);
            for( ++conf_itr_i; conf_itr_i != confs.End(); ++conf_itr_i) {
              float min_rmsd( std::numeric_limits<float>::max());
              for( auto conf_itr_j( uniq_confs.Begin()); conf_itr_j != uniq_confs.End(); ++conf_itr_j)
              {
                float rmsd( (*unique)(*conf_itr_i, *conf_itr_j) );
                min_rmsd = std::min(min_rmsd, rmsd);
              }
              // Check if min_rmsd is >= the threshold with respect to all previously observed conformers
              bool meets_threshold( min_rmsd >= m_UniqueFlag->GetParameterList().LastElement()->GetNumericalValue<float>() );
              if( meets_threshold)
              {
                uniq_confs.PushBack(*conf_itr_i);
              }
            }
            confs = uniq_confs;
          }

          // save all remaining
          auto confs_itr( confs.Begin());
          if( m_SaveEnsembleFlag->GetFlag())
          {
            for( auto conf_itr_end( confs.End()); confs_itr != conf_itr_end; ++confs_itr)
            {
              // set all of the properties from the original molecule as long as they have not been re-cached on the new conformer
              const chemistry::SmallMoleculeMiscProperties &original_properties( input_molecules( m_MoleculeIndex).GetStoredProperties() );
              for
              (
                  auto prop_itr( original_properties.Begin()), prop_itr_end( original_properties.End());
                  prop_itr != prop_itr_end;
                  ++prop_itr
              )
              {
                if( confs_itr->GetMDLProperty( prop_itr->first).empty())
                {
                  confs_itr->GetStoredPropertiesNonConst().SetMDLProperty(prop_itr->first, prop_itr->second);
                }
              }
              confs_itr->GetStoredPropertiesNonConst().SetMDLProperty( "ConformerFromScaffold_scaffold_filename", m_ScaffoldFileFlag->GetFirstParameter()->GetValue() );
              confs_itr->GetStoredPropertiesNonConst().SetMDLProperty( "ConformerFromScaffold_scaffold_molecule_index", linal::Vector<float>(1, similarity_result.First()) );
              confs_itr->GetStoredPropertiesNonConst().SetMDLProperty( "ConformerFromScaffold_similarity_to_scaffold_molecule", linal::Vector< float>(1, similarity_result.Second()) );
              confs_itr->WriteMDL( m_Output);
              ++total_confs;
            }
          }
          else
          {
            const chemistry::SmallMoleculeMiscProperties &original_properties( input_molecules( m_MoleculeIndex).GetStoredProperties() );
            for
            (
                auto prop_itr( original_properties.Begin()), prop_itr_end( original_properties.End());
                prop_itr != prop_itr_end;
                ++prop_itr
            )
            {
              if( confs_itr->GetMDLProperty( prop_itr->first).empty())
              {
                confs_itr->GetStoredPropertiesNonConst().SetMDLProperty(prop_itr->first, prop_itr->second);
              }
            }
            confs_itr->GetStoredPropertiesNonConst().SetMDLProperty( "ConformerFromScaffold_scaffold_filename", m_ScaffoldFileFlag->GetFirstParameter()->GetValue() );
            confs_itr->GetStoredPropertiesNonConst().SetMDLProperty( "ConformerFromScaffold_scaffold_molecule_index", linal::Vector<float>(1, similarity_result.First()) );
            confs_itr->GetStoredPropertiesNonConst().SetMDLProperty( "ConformerFromScaffold_similarity_to_scaffold_molecule", linal::Vector< float>(1, similarity_result.Second()) );
            confs_itr->WriteMDL( m_Output);
            ++total_confs;
          }
          ++m_SuccessCount;
        }
      }

      // close the file stream
      io::File::CloseClearFStream( m_Output);
      io::File::CloseClearFStream( m_OutputSimilarityFailures);
      io::File::CloseClearFStream( m_OutputConfgenFailures);
      BCL_MessageStd( "Summary report: ");
      BCL_MessageStd( std::to_string( m_SuccessCount) + "/" + std::to_string( ensemble_size) + " molecules were successful.");
      BCL_MessageStd( std::to_string( m_SimilarityFailureCount) + "/" + std::to_string( ensemble_size) + " molecules failed due to out of range similarity scores.");
      BCL_MessageStd( std::to_string( m_ConfgenFailureCount) + "/" + std::to_string( ensemble_size) + " molecules failed during conformer generation.");
      BCL_MessageStd( std::to_string( total_confs) + " total conformers were saved.");
      return 0;
    }

  //////////////////////
  // helper functions //
  /////////////////////

    void ConformerFromScaffold::InitializeOutputFiles() const
    {
      BCL_MessageStd("Initializing output files...");
      if( m_OutputFileFlag->GetFlag())
      {
        BCL_MessageStd("Conformers will be written to "
          + util::Format()( m_OutputFileFlag->GetFirstParameter()->GetValue()));
        io::File::MustOpenOFStream( m_Output, m_OutputFileFlag->GetFirstParameter()->GetValue());
      }
      if( m_OutputSimilarityFailureFileFlag->GetFlag())
      {
        BCL_MessageStd("Molecules with insufficient similarity to scaffold(s) will be written to "
          + util::Format()( m_OutputSimilarityFailureFileFlag->GetFirstParameter()->GetValue()));
        io::File::MustOpenOFStream( m_OutputSimilarityFailures, m_OutputSimilarityFailureFileFlag->GetFirstParameter()->GetValue());
      }
      if( m_OutputConfGenFailureFileFlag->GetFlag())
      {
        BCL_MessageStd("Molecules that fail during conformer generation will be written to "
          + util::Format()( m_OutputConfGenFailureFileFlag->GetFirstParameter()->GetValue()));
        io::File::MustOpenOFStream( m_OutputConfgenFailures, m_OutputConfGenFailureFileFlag->GetFirstParameter()->GetValue());
      }
    }

    util::Implementation< chemistry::ConformationComparisonInterface> ConformerFromScaffold::InitializeSimilarityMetric() const
    {
      util::Implementation< chemistry::ConformationComparisonInterface> similarity_metric("LargestCommonSubstructureTanimoto");
      if( m_SolutionTypeFlag->GetFirstParameter()->GetValue() == "Unconnected")
      {
        similarity_metric = util::Implementation< chemistry::ConformationComparisonInterface>("LargestCommonDisconnectedSubstructureTanimoto");
      }
      else if( m_SolutionTypeFlag->GetFirstParameter()->GetValue() == "GreedyUnconnected")
      {
        similarity_metric = util::Implementation< chemistry::ConformationComparisonInterface>("LargestCommonDisconnectedSubstructureTanimoto");
      }
      return similarity_metric;
    }

    graph::CommonSubgraphIsomorphismBase::SolutionType ConformerFromScaffold::InitializeSolutionType() const
    {
      graph::CommonSubgraphIsomorphismBase::SolutionType solution_type( graph::CommonSubgraphIsomorphismBase::SolutionType::e_Connected);
      if( m_SolutionTypeFlag->GetFirstParameter()->GetValue() == "Unconnected")
      {
        solution_type = graph::CommonSubgraphIsomorphismBase::SolutionType::e_Unconnected;
      }
      else if( m_SolutionTypeFlag->GetFirstParameter()->GetValue() == "GreedyUnconnected")
      {
        solution_type = graph::CommonSubgraphIsomorphismBase::SolutionType::e_GreedyUnconnected;
      }
      return solution_type;
    }

    chemistry::FragmentAlignToScaffold ConformerFromScaffold::InitializeAlignmentObject
    (
      const graph::CommonSubgraphIsomorphismBase::SolutionType &SOLUTION_TYPE
    ) const
    {
      return chemistry::FragmentAlignToScaffold
      (
        chemistry::ConformationGraphConverter::AtomComparisonTypeEnum( m_AtomTypeFlag->GetFirstParameter()->GetValue()),
        chemistry::ConfigurationalBondTypeData::DataEnum( m_BondTypeFlag->GetFirstParameter()->GetValue()),
        m_MinSizeFlag->GetFirstParameter()->GetNumericalValue< size_t>(),
        SOLUTION_TYPE
      );
    }

    storage::Pair< size_t, float> ConformerFromScaffold::FindBestScaffoldBySimilarity
    (
      const chemistry::FragmentComplete &MOLECULE,
      const storage::Vector< chemistry::FragmentComplete> &SCAFFOLD_MOLECULES,
      const util::Implementation< chemistry::ConformationComparisonInterface> &SIMILARITY_METRIC
    ) const
    {
      // TODO: allow users to pass pre-computed similarity matrix to avoid computing similarity at this step
      // compute the largest common substructure tanimoto similarity of current molecule to each scaffold
      float best_similarity( 0.0);
      size_t best_similarity_index( 0), scaffold_index( 0);
      for
      (
          auto scaffold_itr( SCAFFOLD_MOLECULES.Begin() ), scaffold_itr_end( SCAFFOLD_MOLECULES.End());
          scaffold_itr != scaffold_itr_end;
          ++scaffold_itr, ++scaffold_index
      )
      {
        float current_similarity( ( *SIMILARITY_METRIC)( MOLECULE, *scaffold_itr ) );
        BCL_MessageVrb( "Current similarity for molecule "
          + std::to_string( m_MoleculeIndex) + " against scaffold "
          + std::to_string( scaffold_index) + ": "
          + std::to_string( current_similarity)
        );
        if( current_similarity > best_similarity)
        {
          best_similarity = current_similarity;
          best_similarity_index = scaffold_index;
          if( best_similarity == 1.0)
          {
            break;
          }
        }
      }
      BCL_MessageVrb
      (
        "Best similarity for molecule " + std::to_string( m_MoleculeIndex)
        + " is against scaffold " + std::to_string( best_similarity_index)
        + ": " + std::to_string( best_similarity )
      );
      if
      (
          best_similarity < m_SimilarityThresholdFlag->GetFirstParameter()->GetNumericalValue< float>() ||
          best_similarity > m_SimilarityThresholdFlag->GetParameterList().LastElement()->GetNumericalValue< float>()
      )
      {
        if( m_OutputSimilarityFailureFileFlag->GetFlag())
        {
          MOLECULE.WriteMDL( m_OutputSimilarityFailures);
        }

        BCL_MessageVrb
        (
          "Molecule " + std::to_string( m_MoleculeIndex) + " failed similarity threshold requirement! "
          "Maximum Tanimoto similarity to scaffold " + std::to_string( best_similarity_index)  +
          " with a value of " + std::to_string( best_similarity)
        );
        ++m_SimilarityFailureCount;
        return storage::Pair< size_t, float>( util::GetUndefinedSize_t(), util::GetUndefined< float>());
      }
      return storage::Pair< size_t, float>( best_similarity_index, best_similarity);
    }

    chemistry::FragmentEnsemble ConformerFromScaffold::Run
    (
      const chemistry::FragmentComplete &TARGET_MOLECULE,
      const chemistry::FragmentComplete &SCAFFOLD_MOLECULE,
      const chemistry::FragmentAlignToScaffold &ALIGNMENT_OBJECT,
      const descriptor::MoleculeSimilarity &ALIGNMENT_SCORER,
      const bool FIND_ALL
    ) const
    {
      chemistry::FragmentEnsemble confs;
      if( FIND_ALL)
      {
        confs = chemistry::FragmentEnsemble
        (
          ALIGNMENT_OBJECT.ConformersFromScaffoldIterative
          (
            TARGET_MOLECULE,
            SCAFFOLD_MOLECULE,
            storage::Vector< size_t>(),
            storage::Vector< size_t>(),
            ALIGNMENT_SCORER
          )
        );
      }
      else
      {
        // generate the new conformer based on the most similar scaffold
        chemistry::FragmentComplete target_mol( TARGET_MOLECULE);
        bool success(
          ALIGNMENT_OBJECT.ConformerFromScaffoldMCS
          (
            target_mol,
            SCAFFOLD_MOLECULE,
            storage::Vector< size_t>(),
            storage::Vector< size_t>(),
            ALIGNMENT_SCORER
          )
        );
        if( success)
        {
          confs = chemistry::FragmentEnsemble( storage::List< chemistry::FragmentComplete>(1, target_mol));
        }
      }
      return confs;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ConformerFromScaffold::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &ConformerFromScaffold::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace app
} // namespace bcl
