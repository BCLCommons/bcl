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
#include "bcl_app_align_binding_poses.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_quality.h"
#include "assemble/bcl_assemble_voxel_grid_atom.h"
#include "biol/bcl_biol_dssp.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "command/bcl_command_parameter_check_extensions_file_existence.h"
#include "graph/bcl_graph_common_subgraph_isomorphism.h"
#include "graph/bcl_graph_subgraph_isomorphism.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_mutate_result.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "quality/bcl_quality_rmsd.h"
#include "quality/bcl_quality_superimpose_measures.h"
#include "storage/bcl_storage_template_instantiations.h"
namespace bcl
{
  namespace app
  {

    // Static instance initialization
    const ApplicationType AlignBindingPoses::AlignBindingPoses_Instance
    (
      GetAppGroups().AddAppToGroup( new AlignBindingPoses(), GetAppGroups().e_Molecule)
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief standard constructor
    AlignBindingPoses::AlignBindingPoses() :
      m_InputMolInBindingPos
      (
        new command::Parameter
        (
          "template_mol",
          "molecule in the binding pose. Will be output as first molecule in the ensemble",
          command::ParameterCheckFileExistence()
        )
      ),
      m_TemplatePDB
      (
        new command::Parameter
        (
          "template_pdb",
          "protein model to which all others will be aligned",
          command::ParameterCheckFileExistence()
        )
      ),
      m_OutputFilename
      (
        new command::Parameter( "output", "filename for output sdf of aligned ensemble")
      ),
      m_Molecules
      (
        new command::FlagDynamic
        (
          "ligands",
          "ligands; each file may have multiple ligands in the same binding pocket that were already aligned or alternate "
          "ligands in the same binding pocket. each ligand should have a corresponding entry in -proteins for the associated PDB",
          command::Parameter
          (
            "ligand",
            "ligands of interest",
            command::ParameterCheckFileExistence()
          ),
          0,
          10000000
        )
      ),
      m_Pdbs
      (
        new command::FlagDynamic
        (
          "proteins",
          "proteins associated with each ligand",
          command::Parameter
          (
            "protein",
            "protein of interest",
            command::ParameterCheckFileExistence()
          ),
          0,
          10000000
        )
      ),
      m_SuperimposeMethodFlag
      (
        new command::FlagStatic
        (
          "superimpose_method",
          "metric to optimize with the superimposition",
          command::Parameter
          (
            "superimpose_measure",
            "\tquality measure to use for superimposition",
            command::ParameterCheckEnumerate< quality::SuperimposeMeasures>(),
            quality::GetSuperimposeMeasures().e_RMSD.GetName()
          )
        )
      )
    {
    }

    //! @brief Clone function
    //! @return pointer to new FoldProtein
    AlignBindingPoses *AlignBindingPoses::Clone() const
    {
      return new AlignBindingPoses( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AlignBindingPoses::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> AlignBindingPoses::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // scaffold sdf that ensemble will be aligned to
      sp_cmd->AddParameter( m_InputMolInBindingPos);
      sp_cmd->AddParameter( m_TemplatePDB);
      sp_cmd->AddParameter( m_OutputFilename);

      // hydrogen preferences
      sdf::AddMoleculeIOPrefFlags( *sp_cmd);

      // ensemble containing the molecules to be aligned to scaffold
      sp_cmd->AddFlag( m_Molecules);
      // ensemble that will be written out
      sp_cmd->AddFlag( m_Pdbs);
      // atoms in the scaffold to align
      sp_cmd->AddFlag( m_SuperimposeMethodFlag);
      sp_cmd->AddFlag( pdb::Factory::GetFlagAAClass());
      sp_cmd->AddFlag( pdb::Factory::GetFlagWritePDBAtomID());
      sp_cmd->AddFlag( pdb::Factory::GetFlagWritePDBResID());

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    void AlignBindingPoses::WriteContacts( const chemistry::FragmentComplete &FRAG, const assemble::ProteinModel &MODEL) const
    {
      assemble::VoxelGridAtom grid( 4.0);
      auto aas( *MODEL.GetSequences().FirstElement());

      storage::Set< std::string> pdb_ids;
      storage::Map< size_t, util::SiPtr< const biol::AABase> > atom_pdb_ids_to_aa;
      util::SiPtrVector< const biol::Atom> all_the_atoms;
      for( auto itr( aas.Begin()), itr_end( aas.End()); itr != itr_end; ++itr)
      {
        all_the_atoms.Append( ( *itr)->GetAtoms());
        for( auto itr_atom( ( *itr)->GetAtoms().Begin()), itr_atom_end( ( *itr)->GetAtoms().End()); itr_atom != itr_atom_end; ++itr_atom)
        {
          atom_pdb_ids_to_aa[ ( *itr_atom)->GetPdbID()] = *itr;
        }
      }
      grid.SetObjects( all_the_atoms);
      for( auto itr( FRAG.GetAtomsIterator()); itr.NotAtEnd(); ++itr)
      {
        auto neigh( grid.GetNeighbors( itr->GetPosition(), 4.0));
        for( auto itr_neigh( neigh.Begin()), itr_neigh_end( neigh.End()); itr_neigh != itr_neigh_end; ++itr_neigh)
        {
          pdb_ids.Insert
          (
            util::Format()( atom_pdb_ids_to_aa[ itr_neigh->First()->GetPdbID()]->GetPdbID()) + " "
            + util::Format()( atom_pdb_ids_to_aa[ itr_neigh->First()->GetPdbID()]->GetType()->GetOneLetterCode())
          );
        }
      }
      BCL_MessageStd( "Contacting Atom PDBIDs: " + util::Format()( util::Join( ",", storage::Vector< std::string>( pdb_ids.Begin(), pdb_ids.End()))));
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int AlignBindingPoses::Main() const
    {
      BCL_Assert
      (
        m_Pdbs->GetStringList().GetSize() == m_Molecules->GetStringList().GetSize(),
        "Must have the same number of entries for -pdbs and -ligands!"
      );

      // if the superimposition flag is set
      // read the native pdb in
      pdb::Factory factory;
      io::IFStream read;
      io::File::MustOpenIFStream( read, m_TemplatePDB->GetValue());
      pdb::Handler template_pdb_reader( read);
      io::File::CloseClearFStream( read);

      // instantiate protein model
      assemble::ProteinModel template_model( factory.ProteinModelFromPDB( template_pdb_reader));

      // process model with dssp
      biol::DSSP model_dssp;
//      math::MutateResult< assemble::ProteinModel> new_model( model_dssp( template_model));
//      if( new_model.GetArgument().IsDefined())
//      {
//        template_model = *new_model.GetArgument();
//      }

      // get the enum from name
      const quality::SuperimposeMeasure measure
      (
        m_SuperimposeMethodFlag->GetFirstParameter()->GetValue()
      );

      io::File::MustOpenIFStream( read, m_InputMolInBindingPos->GetValue());
      chemistry::FragmentEnsemble native_mols( read, sdf::GetCommandLineHydrogensPref());
      io::File::CloseClearFStream( read);

      io::OFStream output, output_pdb;

      io::File::MustOpenOFStream( output_pdb, m_TemplatePDB->GetValue() + ".transformed.pdb");
      factory.WriteModelToPDB( template_model, output_pdb);
      io::File::CloseClearFStream( output_pdb);

      io::File::MustOpenOFStream( output, m_OutputFilename->GetValue());
      native_mols.WriteMDL( output);
      for
      (
        auto itr_ligand( native_mols.Begin()), itr_ligand_end( native_mols.End());
        itr_ligand != itr_ligand_end;
        ++itr_ligand
      )
      {
        WriteContacts( *itr_ligand, template_model);
      }

      const size_t n_mols( m_Pdbs->GetStringList().GetSize());
      for( size_t mol_n( 0); mol_n < n_mols; ++mol_n)
      {
        // read in the pdb model
        io::File::MustOpenIFStream( read, m_Pdbs->GetStringList()( mol_n));
        pdb::Handler reader( read);
        io::File::CloseClearFStream( read);

        assemble::ProteinModel aligned_model( factory.ProteinModelFromPDB( reader));
//        math::MutateResult< assemble::ProteinModel> new_model( model_dssp( aligned_model));
//        if( new_model.GetArgument().IsDefined())
//        {
//          aligned_model = *new_model.GetArgument();
//        }
        io::File::MustOpenIFStream( read, m_Molecules->GetStringList()( mol_n));
        chemistry::FragmentEnsemble transformed_ligands( read, sdf::GetCommandLineHydrogensPref());
        io::File::CloseClearFStream( read);
        for
        (
          auto itr_ligand( transformed_ligands.Begin()), itr_ligand_end( transformed_ligands.End());
          itr_ligand != itr_ligand_end;
          ++itr_ligand
        )
        {
          WriteContacts( *itr_ligand, aligned_model);
        }
        // transform the model
        auto pair_transformation
        (
          assemble::Quality::SuperimposeModelWithAlignment
          (
            measure,
            aligned_model,
            template_model,
            biol::GetAtomTypes().GetBackBoneAtomTypes()
          )
        );
        if( !pair_transformation.First())
        {
          BCL_MessageCrt
          (
            "Could not align " + m_Pdbs->GetStringList()( mol_n)
            + " to the template pdb; associated molecules "
            " will not be output"
          );
          continue;
        }
        io::File::MustOpenOFStream( output_pdb, m_Pdbs->GetStringList()( mol_n)+".transformed.pdb");
        factory.WriteModelToPDB( aligned_model, output_pdb);
        io::File::CloseClearFStream( output_pdb);

        for
        (
          auto itr_ligand( transformed_ligands.Begin()), itr_ligand_end( transformed_ligands.End());
          itr_ligand != itr_ligand_end;
          ++itr_ligand
        )
        {
          chemistry::FragmentComplete frag( *itr_ligand);
          frag.Transform( pair_transformation.Second());
          frag.WriteMDL( output);
        }
      }
      io::File::CloseClearFStream( output);
      return 0;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AlignBindingPoses::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &AlignBindingPoses::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace app
} // namespace bcl
