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
#include "chemistry/bcl_chemistry_fragment_map_conformer.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
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
              "input molecules",
              "SDF format",
              command::ParameterCheckFileExistence()
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
              "scaffold molecules",
              "SDF format",
              command::ParameterCheckFileExistence()
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
              "output molecules",
              "SDF format"
            )
          )
        ),
        m_OutputFailureFileFlag
        (
          new command::FlagDynamic
          (
            "output_failures_filename",
            "failed molecules; occurs due to invalid conformer generation or insufficient similarity to a scaffold",
            command::Parameter
            (
              "output failure molecules",
              "SDF format"
            )
          )
        ),
        m_ModeFlag
        (
          new command::FlagStatic
          (
            "mode",
            "mode to determine how input molecules will be compared to output",
            command::Parameter
            (
              "mode",
              "",
              command::ParameterCheckAllowed( storage::Vector< std::string>::Create( "all", "similarity" ) ),
              "all"
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
              command::ParameterCheckAllowed(
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
            command::Parameter
            (
              "",
              "",
              command::ParameterCheckRanged< float>(0.0, 1.0),
              "0.0"
            )
          )
        )
    {
    }

    //! copy constructor, only copy the flags
    ConformerFromScaffold::ConformerFromScaffold( const ConformerFromScaffold &PARENT) :
          m_InputFileFlag( PARENT.m_InputFileFlag),
          m_ScaffoldFileFlag( PARENT.m_ScaffoldFileFlag),
          m_OutputFileFlag( PARENT.m_OutputFileFlag),
          m_OutputFailureFileFlag( PARENT.m_OutputFailureFileFlag),
          m_ModeFlag( PARENT.m_ModeFlag),
          m_AtomTypeFlag( PARENT.m_AtomTypeFlag),
          m_BondTypeFlag( PARENT.m_BondTypeFlag),
          m_MinSizeFlag( PARENT.m_MinSizeFlag),
          m_SolutionTypeFlag( PARENT.m_SolutionTypeFlag),
          m_SimilarityThresholdFlag( PARENT.m_SimilarityThresholdFlag)
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

      // insert all the flags and params
      sp_cmd->AddFlag( m_InputFileFlag);
      sp_cmd->AddFlag( m_ScaffoldFileFlag);
      sp_cmd->AddFlag( m_OutputFileFlag);
      sp_cmd->AddFlag( m_OutputFailureFileFlag);
      sp_cmd->AddFlag( m_ModeFlag);
      sp_cmd->AddFlag( m_AtomTypeFlag);
      sp_cmd->AddFlag( m_BondTypeFlag);
      sp_cmd->AddFlag( m_MinSizeFlag);
      sp_cmd->AddFlag( m_SolutionTypeFlag);
      sp_cmd->AddFlag( m_SimilarityThresholdFlag);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int ConformerFromScaffold::Main() const
    {
      // read in ensemble without hydrogen atoms to speedup substructure search
      io::IFStream input;
      BCL_MessageStd("Reading input molecules...");
      io::File::MustOpenIFStream( input, m_InputFileFlag->GetFirstParameter()->GetValue());
      chemistry::FragmentEnsemble input_ensemble( input, sdf::e_Remove);
      io::File::CloseClearFStream( input);
      const size_t ensemble_size( input_ensemble.GetSize());
      BCL_Assert( ensemble_size, "Must have at least one molecule in the input ensemble. Exiting...\n");

      // read in the scaffold ensemble without hydrogen atoms to speedup substructure search
      BCL_MessageStd("Reading scaffold molecules...");
      io::File::MustOpenIFStream( input, m_ScaffoldFileFlag->GetFirstParameter()->GetValue());
      const chemistry::FragmentEnsemble scaffold_ensemble( input, sdf::e_Remove);
//      storage::Vector< chemistry::FragmentComplete> scaffold_molecules( scaffold_ensemble.Begin(), scaffold_ensemble.End());
      io::File::CloseClearFStream( input);
      const size_t scaffold_ensemble_size( scaffold_ensemble.GetSize());
      BCL_Assert( scaffold_ensemble_size, "Must have at least one molecule in the scaffold ensemble. Exiting...\n");

      // initialize output so that we can write as we go
      io::OFStream output, output_failures;
      io::File::MustOpenOFStream( output, m_OutputFileFlag->GetFirstParameter()->GetValue());
      io::File::MustOpenOFStream( output_failures, m_OutputFailureFileFlag->GetFirstParameter()->GetValue());

      // cleaning the scaffold molecules prior to comparison and conformer generation
      storage::Vector< chemistry::FragmentComplete> scaffold_molecules;
      for
      (
          auto scaffold_itr( scaffold_ensemble.Begin()), scaffold_itr_end( scaffold_ensemble.End());
          scaffold_itr != scaffold_itr_end;
          ++scaffold_itr
      )
      {
        chemistry::AtomVector< chemistry::AtomComplete> scaffold_atom_vector( scaffold_itr->GetAtomVector());
        scaffold_atom_vector = chemistry::FragmentMapConformer::CleanAtoms( scaffold_atom_vector, "None", true, true);
        chemistry::FragmentComplete clean_scaffold( scaffold_atom_vector, scaffold_itr->GetName());
        clean_scaffold.StoreProperties( *scaffold_itr);
        scaffold_molecules.PushBack( clean_scaffold);
      }

      // Get the atom and bond type resolution for substructure comparisons
      BCL_MessageStd("Initializing scaffold alignment object...");
      graph::CommonSubgraphIsomorphismBase::SolutionType solution_type( graph::CommonSubgraphIsomorphismBase::SolutionType::e_Connected);
      if( m_SolutionTypeFlag->GetFirstParameter()->GetValue() == "Unconnected")
      {
        solution_type = graph::CommonSubgraphIsomorphismBase::SolutionType::e_Unconnected;
      }
      else if( m_SolutionTypeFlag->GetFirstParameter()->GetValue() == "GreedyUnconnected")
      {
        solution_type = graph::CommonSubgraphIsomorphismBase::SolutionType::e_GreedyUnconnected;
      }
      chemistry::FragmentAlignToScaffold ats
      (
        chemistry::ConformationGraphConverter::AtomComparisonTypeEnum( m_AtomTypeFlag->GetFirstParameter()->GetValue()),
        chemistry::ConfigurationalBondTypeData::DataEnum( m_BondTypeFlag->GetFirstParameter()->GetValue()),
        m_MinSizeFlag->GetFirstParameter()->GetNumericalValue< size_t>(),
        solution_type
      );

      // generate conformers for each input file based on template substructure
      const util::Implementation< chemistry::ConformationComparisonInterface> similarity_metric("LargestCommonSubstructureTanimoto");
      BCL_MessageStd("Generating new conformer ensemble for input molecules...");
      size_t mol_index( 0);
      BCL_MessageStd( "Completed " + std::to_string( mol_index) + "/" + std::to_string( ensemble_size) + " molecules.");
      for
      (
          auto mol_itr( input_ensemble.Begin()), mol_itr_end( input_ensemble.End());
          mol_itr != mol_itr_end;
          ++mol_itr, ++mol_index
      )
      {
        // status update
        if( mol_index % 10)
        {
          BCL_MessageStd( "Completed " + std::to_string( mol_index) + "/" + std::to_string( ensemble_size) + " molecules.");
//          util::GetLogger().LogStatus( "Completed " + std::to_string( mol_index) + "/" + std::to_string( ensemble_size) + " molecules.");
        }

        // clean molecule before comparison
        chemistry::AtomVector< chemistry::AtomComplete> mol_atom_vector( mol_itr->GetAtomVector());
        mol_atom_vector = chemistry::FragmentMapConformer::CleanAtoms( mol_atom_vector, "None", true, true);
        chemistry::FragmentComplete clean_mol( mol_atom_vector, mol_itr->GetName());
        clean_mol.StoreProperties( *mol_itr);
        *mol_itr = clean_mol;

        // TODO: allow users to pass pre-computed similarity matrix to avoid computing similarity at this step
        // compute the largest common substructure tanimoto similarity of current molecule to each scaffold
        float best_similarity( 0.0);
        size_t best_similarity_index( 0), scaffold_index( 0);
        for
        (
            auto scaffold_itr( scaffold_molecules.Begin() ), scaffold_itr_end( scaffold_molecules.End());
            scaffold_itr != scaffold_itr_end;
            ++scaffold_itr, ++scaffold_index
        )
        {
          float current_similarity( ( *similarity_metric)( *mol_itr, *scaffold_itr ) );
          BCL_MessageVrb( "Current similarity for molecule " + std::to_string( mol_index) + " against scaffold " + std::to_string( scaffold_index) + ": " + std::to_string( current_similarity));
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

        BCL_MessageVrb( "Best similarity for molecule " + std::to_string( mol_index) + " is against scaffold " + std::to_string( best_similarity_index) + ": " + std::to_string( best_similarity ) );
        if( best_similarity < m_SimilarityThresholdFlag->GetFirstParameter()->GetNumericalValue< float>() )
        {
          if( m_OutputFailureFileFlag->GetFlag())
          {
            mol_itr->WriteMDL( output_failures);
          }

          BCL_MessageVrb
          (
            "Molecule " + std::to_string( mol_index) + " failed similarity threshold requirement! "
            "Maximum Tanimoto similarity to scaffold " + std::to_string( best_similarity_index)  +
            " with a value of " + std::to_string( best_similarity)
          );
          continue;
        }

        // generate the new conformer based on the most similar scaffold
        bool success( ats.GenerateConformerBasedOnScaffold( *mol_itr, scaffold_molecules( best_similarity_index) ) );
        if( !success)
        {
          if( m_OutputFailureFileFlag->GetFlag())
          {
            mol_itr->WriteMDL( output_failures);
          }

          BCL_MessageVrb
          (
            "Failed to generate conformer for molecule " + std::to_string( mol_index) +
            " using scaffold " + std::to_string( best_similarity_index) +
            " with a Tanimoto similarity of " + std::to_string( best_similarity)
          );
        }
        else
        {
          mol_itr->WriteMDL( output);
        }
      }

      // close the file stream
      io::File::CloseClearFStream( output);
      io::File::CloseClearFStream( output_failures);
      BCL_MessageStd( "Completed " + std::to_string( mol_index) + "/" + std::to_string( ensemble_size) + " molecules.");
      BCL_MessageStd("Done!");
      return 0;
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
