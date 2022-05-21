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
#include "bcl_app_molecule_extract_protein_pocket.h"
#include "biol/bcl_biol_aa_complete.h"
#include "biol/bcl_biol_aa_side_chain_factory.h"
#include "chemistry/bcl_chemistry_aa_fragment_complete.h"
#include "chemistry/bcl_chemistry_conformation_comparison_interface.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_voxel_grid_atom.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_histogram.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_thunk_job.h"
#include "sched/bcl_sched_unary_function_job_with_data.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_template_instantiations.h"
#include "util/bcl_util_implementation.h"
namespace bcl
{
  namespace app
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief standard constructor
    MoleculeExtractProteinPocket::MoleculeExtractProteinPocket() :
      m_InputMoleculeSDF
      (
        new command::FlagStatic
        (
          "molecule_filename",
          "filename for input sdf of fragment ensemble",
          command::Parameter( "input_sdf", "a single SDF file of ligands", "")
        )
      ),
      m_InputPocketPDB
      (
        new command::FlagStatic
        (
          "receptor_filename",
          "filename for a single receptor",
          command::Parameter( "input_pdb", "a single PDB file of the protein receptor containing the ligand binding pocket", "")
        )
      ),
      m_OutputPrefix
      (
        new command::FlagStatic
        (
          "output_prefix",
          "output prefix for the resulting binding pocket pdb file",
          command::Parameter( "output prefix", "<output_prefix>.pdb", "output_bcl_pocket")
        )
      ),
      m_NeighborDistanceCutoff
      (
        new command::FlagStatic
        (
          "neighbor_distance",
          "distance within which two atoms are neighbors",
          command::Parameter( "neighbor_distance", "determines which atoms from the protein will be considered pocket atoms of the ligand", "8.0")
        )
      ),
      m_OutputAtomsOnly
      (
        new command::FlagStatic
        (
          "atoms_only",
          "only output pocket atoms and not full residues",
          command::Parameter( "atoms_only", "output pocket atoms in PDB format", "")
        )
      )
    {
    }

    //! @brief Clone function
    //! @return pointer to new FoldProtein
    MoleculeExtractProteinPocket *MoleculeExtractProteinPocket::Clone() const
    {
      return new MoleculeExtractProteinPocket( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MoleculeExtractProteinPocket::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string MoleculeExtractProteinPocket::GetDescription() const
    {
      return "extract binding pocket residues from a protein";
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> MoleculeExtractProteinPocket::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // ensemble containing the small molecules
      sp_cmd->AddFlag( m_InputMoleculeSDF);

      // PDB file for a single protein binding pocket
      sp_cmd->AddFlag( m_InputPocketPDB);

      // output prefix
      sp_cmd->AddFlag( m_OutputPrefix);

      // distance specifying neighbors
      sp_cmd->AddFlag( m_NeighborDistanceCutoff);

      // specify only output neighbor atoms
      sp_cmd->AddFlag( m_OutputAtomsOnly);

      // hydrogen flags
      sdf::AddMoleculeIOPrefFlags( *sp_cmd);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int MoleculeExtractProteinPocket::Main() const
    {
      chemistry::FragmentEnsemble ensemble_a;
      chemistry::FragmentComplete combined_mol, atom_pocket;
      chemistry::AAFragmentComplete aa_fragment, pocket;
      const pdb::Factory factory( biol::GetAAClasses().e_AAComplete);
      assemble::ProteinModel protein_model;

      BCL_Assert( m_InputMoleculeSDF->GetFlag() == m_InputPocketPDB->GetFlag(), "Both an SDF for molecules and a PDB for proteins are required");
      if( m_InputMoleculeSDF->GetFlag() && m_InputPocketPDB->GetFlag())
      {
        // read in ligand(s)
        io::IFStream input_sdf;
        io::File::MustOpenIFStream( input_sdf, m_InputMoleculeSDF->GetFirstParameter()->GetValue());
        ensemble_a = chemistry::FragmentEnsemble( input_sdf, sdf::e_Maintain);
        BCL_MessageStd( "Binding pocket will be determined based on the positions of " + util::Format()( ensemble_a.GetSize()) + " small molecules.");

        // If more than one molecule, connect without bonds
        // The goal here is to select residues for the binding pocket based on the entire ensemble, so we need one voxel grid
        if( ensemble_a.GetSize() > 1)
        {
          storage::Vector< sdf::BondInfo> empty( 0);
          chemistry::AtomVector< chemistry::AtomComplete> ensemble_atoms( ensemble_a.Begin()->GetAtomVector());
          for( auto ens_a_itr( ++ensemble_a.Begin()); ens_a_itr != ensemble_a.End(); ++ens_a_itr)
          {
            ensemble_atoms.AddAtomsWithConnectivity( ens_a_itr->GetAtomVector(), empty);
          }
          combined_mol = chemistry::FragmentComplete( ensemble_atoms, "");
        }
        else
        {
          combined_mol = *ensemble_a.Begin();
        }

        // read in protein
        protein_model = assemble::ProteinModel( factory.ProteinModelFromPDBFilename( m_InputPocketPDB->GetFirstParameter()->GetValue(), pdb::Factory::GetSSETypeMinSizes( 0, 0, 0)));

        //instantiate AASideChainFactory no hydrogens include backbone atoms for sidechain placement superimposition
        biol::AASideChainFactory side_chain_factory( false, true);

        // add side chains to model and set it to model
        //   protein_model = *side_chain_factory.ProteinModelWithSideChains( protein_model);
        aa_fragment = chemistry::AAFragmentComplete( protein_model.GetAminoAcids(), true);

        io::File::CloseClearFStream( input_sdf);
      }

      //Build voxel grids for the protein and the ligand(s)
      double neighbor_distance( m_NeighborDistanceCutoff->GetFirstParameter()->GetNumericalValue< double>());
      chemistry::VoxelGridAtom voxel_grid_pocket( neighbor_distance), voxel_grid_mol( neighbor_distance);

      voxel_grid_mol.SetObjects( util::SiPtrVector< const chemistry::AtomConformationalInterface>( combined_mol.GetAtomsIterator(), combined_mol.GetAtomsIterator().End()));
      voxel_grid_pocket.SetObjects( util::SiPtrVector< const chemistry::AtomConformationalInterface>( aa_fragment.GetAtomsIterator(), aa_fragment.GetAtomsIterator().End()));
      auto neighbors( voxel_grid_mol.GetNeighborsIn( voxel_grid_pocket, neighbor_distance));

      // default behavior returns protein binding-pocket residues
      if( !m_OutputAtomsOnly->GetFlag())
      {
        //Get residue IDs back from atoms found in the protein
        util::SiPtrVector< biol::AABase> pocket_residues;
        for( auto itr_neighbors( neighbors.Begin() + 1), itr_neighbors_end( neighbors.End()); itr_neighbors != itr_neighbors_end; ++itr_neighbors)
        {
          pocket_residues.PushBack( new biol::AAComplete( aa_fragment.GetAtomsResidue( aa_fragment.GetAtomIndex( *itr_neighbors->Second()))));
          BCL_MessageStd( aa_fragment.GetAtomsResidue( aa_fragment.GetAtomIndex( *itr_neighbors->Second())).GetIdentification());
        }

        //Construct a set to get only the unique residues
        storage::Set< util::SiPtr< biol::AABase>> pocket_set( pocket_residues.Begin(), pocket_residues.End());
        util::SiPtrVector< biol::AABase> unique_pocket_residues( pocket_set.Begin(), pocket_set.End());

        //Make new molecule with residues
        pocket = chemistry::AAFragmentComplete( unique_pocket_residues, true);

        io::OFStream output_pocket;
        io::File::MustOpenOFStream( output_pocket, m_OutputPrefix->GetFirstParameter()->GetValue() + ".pdb");
        pdb::Factory().WriteModelToPDB( pocket.ReconstructProteinModel( protein_model), output_pocket);
      }

      // alternatively, just return the atoms in SDF format (for you, Jeff)
      else
      {
        // Make a pseudomolecule from (1) atominfo in neighbor atoms, and (2) empty bond connections
        storage::Vector< sdf::AtomInfo> atom_info;
        storage::Vector< sdf::BondInfo> empty;
        for( auto itr_neighbors( neighbors.Begin() + 1), itr_neighbors_end( neighbors.End()); itr_neighbors != itr_neighbors_end; ++itr_neighbors)
        {
          atom_info.PushBack( itr_neighbors->Second()->GetAtomInfo());
        }
        chemistry::AtomVector< chemistry::AtomComplete> pocket_atoms( atom_info, empty);

        // Make new molecule with pocket atoms
        atom_pocket = chemistry::FragmentComplete( pocket_atoms, "");

        // Write out pocket atoms in SDF format
        io::OFStream output_pocket;
        io::File::MustOpenOFStream( output_pocket, m_OutputPrefix->GetFirstParameter()->GetValue() + ".sdf");
        atom_pocket.WriteMDL( output_pocket);
      }

      return 0;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MoleculeExtractProteinPocket::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &MoleculeExtractProteinPocket::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &MoleculeExtractProteinPocket::GetReadMe() const
    {
      static std::string s_read_me =
        "MoleculeExtractProteinPocket creates a PDB file containing a binding pocket provided the following:\n"
         "-a small molecule ligand or reference coordinate\n"
         "-a protein\n"
         "-a distance indicating how far from the reference coordinate(s) to check for protein atoms/residues\n";
      return s_read_me;
    }

    const ApplicationType MoleculeExtractProteinPocket::MoleculeExtractProteinPocket_Instance
    (
      GetAppGroups().AddAppToGroup( new MoleculeExtractProteinPocket(), GetAppGroups().e_Molecule)
    );

  } // namespace app
} // namespace bcl

