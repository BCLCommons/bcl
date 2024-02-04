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

// App header
#include "app/bcl_app_apps.h"

// include header for this application
#include "molecule/bcl_app_generate_partial_charges_file.h"

// bcl includes
#include "chemistry/bcl_chemistry_atom_types.h"
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_bond_isometry_handler.h"
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "chemistry/bcl_chemistry_configurational_bond_types.h"
#include "chemistry/bcl_chemistry_conformation_comparison_psi_field.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "chemistry/bcl_chemistry_fragment_make_conformers.h"
#include "chemistry/bcl_chemistry_fragment_split_interface.h"
#include "chemistry/bcl_chemistry_fragment_split_largest_component.h"
#include "chemistry/bcl_chemistry_merge_fragment_complete.h"
#include "chemistry/bcl_chemistry_molecule_feature_mapper.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "chemistry/bcl_chemistry_stereocenters_handler.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "graph/bcl_graph_common_subgraph_isomorphism.h"
#include "graph/bcl_graph_common_subgraph_isomorphism_base.h"
#include "graph/bcl_graph_connectivity.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_atom_info.h"
#include "sdf/bcl_sdf_bond_info.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_format.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"

namespace bcl
{
  namespace app
  {

    //! @brief initializes the command object for this application
    //! @return a ShPtr to a Command containing all of this applications parameters
    util::ShPtr< command::Command> GeneratePartialChargesFile::InitializeCommand() const
    {
      // initialize command to be returned
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // ensembles containing the molecules to be filtered
       chemistry::FragmentFeed::AddFlags( *sp_cmd);

      // add command line options to add/remove hydrogens
      sdf::AddMoleculeIOPrefFlags( *sp_cmd);

      // add AtomMdlLine to molecule
      sdf::MdlHandler::AddAtomMdlLineFlag( *sp_cmd);

      // default flags
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // unique flags
      sp_cmd->AddFlag( m_OutputFlag);
      sp_cmd->AddFlag( m_PartialChargeTypeFlag);
      return sp_cmd;

    } // InitializeCommand

    //! @brief the Main function
    //! @return 0 for success
    int GeneratePartialChargesFile::Main() const
    {
      // Set partial charge type
      m_PartialChargeType = PartialChargeType( m_PartialChargeTypeFlag->GetFirstParameter()->GetNumericalValue< size_t>());
      BCL_MessageVrb( "Partial charge type: " + util::Format()( m_PartialChargeType));

      // input molecules; default constructor will read commandline preferences for hydrogen atoms and neutralization
      chemistry::FragmentFeed feed;
      for( ; feed.NotAtEnd(); ++feed)
      {
        // keep track of position in input
        const size_t mol_index( feed.GetPosition());
        const chemistry::FragmentComplete &current_mol( *feed);

        // Generate partial charge file
        const linal::Vector< float> partial_charges( ComputeAtomicPartialCharges( current_mol));
        if( partial_charges.GetSize() == current_mol.GetSize())
        {
          // Output those partial charges with the element
          WritePartialChargeFile( mol_index, partial_charges, current_mol.GetAtomVector());
        }
        else
        {
          BCL_MessageStd( "Error: Fails to generate partial charge file for molecule # " + util::Format()( mol_index));
        }
      }
      return 0;
    } // Main

    //! @brief Compute partial charges of MOL
    //! @param MOL the molecule for which atomic partial charges will be computed
    //! @return atomic partial charges of MOL
    const linal::Vector< float>
    GeneratePartialChargesFile::ComputeAtomicPartialCharges
    (
      const chemistry::FragmentComplete &MOL
    ) const
    {
      linal::Vector< float> partial_charges;
      if( m_PartialChargeType == e_SigmaCharge)
      {
        partial_charges = descriptor::GetCheminfoProperties().calc_SigmaCharge->CollectValuesOnEachElementOfObject( MOL);
      }
      else if( m_PartialChargeType == e_PiCharge)
      {
        partial_charges = descriptor::GetCheminfoProperties().calc_PiCharge->CollectValuesOnEachElementOfObject( MOL);
      }
      else if( m_PartialChargeType == e_TotalCharge)
      {
        partial_charges = descriptor::GetCheminfoProperties().calc_TotalCharge->CollectValuesOnEachElementOfObject( MOL);
      }
      else if( m_PartialChargeType == e_VCharge)
      {
        partial_charges = descriptor::GetCheminfoProperties().calc_VCharge->CollectValuesOnEachElementOfObject( MOL);
      }
      else if( m_PartialChargeType == e_VCharge2)
      {
        partial_charges = descriptor::GetCheminfoProperties().calc_VChargeV2->CollectValuesOnEachElementOfObject( MOL);
      }
      else // default TotalCharge
      {
        partial_charges = descriptor::GetCheminfoProperties().calc_TotalCharge->CollectValuesOnEachElementOfObject( MOL);
      }
      return partial_charges;
    }

    //! @brief write the file that contains partial charge and element type for each atom in NCAA
    //! @parameter MOL_INDEX: index of the NCAA in the
    //! @parameter PARTIAL_CHARGES: current index of the C backbone atom
    //! @parameter ATOMS: current index of the N backbone atom
    void GeneratePartialChargesFile::WritePartialChargeFile
    (
      const size_t &MOL_INDEX,
      const linal::Vector< float> &PARTIAL_CHARGES,
      const chemistry::AtomVector< chemistry::AtomComplete> &ATOMS
    ) const
    {
      // open output file
      io::OFStream out;
      io::File::MustOpenOFStream
      (
        out,
        m_OutputFlag->GetFirstParameter()->GetValue() + "." + util::Format()( MOL_INDEX) + ".txt"
      );

      // Write out the element name and partial charge for each atom in the dipeptide form of the NCAA
      size_t index( 0);
      for
      (
          linal::Vector< float>::const_iterator iter( PARTIAL_CHARGES.Begin()), end_iter( PARTIAL_CHARGES.End());
          iter != end_iter;
          ++iter, ++index
      )
      {
        out <<
            util::Format()( index + 1) <<
            " " <<
            util::Format()( ATOMS( index).GetElementType()->GetChemicalSymbol()) <<
            " " <<
            util::Format()( *iter) <<
            std::endl; // the number of the N-terminal atom
      }

      // close output file
      io::File::CloseClearFStream( out);
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string GeneratePartialChargesFile::GetDescription() const
    {
      return "GeneratePartialChargesFile returns a partial charges file that can be used with "
          "Rosetta molfile_to_params scripts to generate small molecules and non-canonical amino acids "
          "with the specified atomic partial charges";
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &GeneratePartialChargesFile::GetReadMe() const
    {
      static std::string s_read_me =
        "GeneratePartialChargesFile computes atomic partial charges for molecules and outputs them "
        "in a format that is readable by Rosetta molfile_to_params scripts for param file creation.";
      return s_read_me;
    }

    //! @brief standard constructor
    GeneratePartialChargesFile::GeneratePartialChargesFile() :
      m_OutputFlag
      (
        new command::FlagStatic
        (
          "output_prefix",
          "prefix for files to which output will be written",
          command::Parameter
          (
            "output_prefix",
            ""
          )
        )
      ),
      m_PartialChargeTypeFlag
      (
        new command::FlagStatic
        (
          "partial_charge_type",
          "specify the type of partial charges to assign if 'generate_partial_charge_file' is set.",
          command::Parameter
          (
            "partial_charge_type",
            "0 - SigmaCharge \n"
            "1 - PiCharge \n"
            "2 - TotalCharge \n"
            "3 - VCharge \n"
            "4 - VCharge2 \n",
            command::ParameterCheckRanged< size_t>( 0, 4),
            "2"
          )
        )
      )
    {
    }

    // Construct the static instance of this application, and add it to the ChemInfo group
    const ApplicationType GeneratePartialChargesFile::GeneratePartialChargesFile_Instance
    (
      GetAppGroups().AddAppToGroup( new GeneratePartialChargesFile(), GetAppGroups().e_Molecule)
    );

  } // namespace app
} // namespace bcl
