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
#include "command/bcl_command_parameter_check_ranged.h"
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "app/bcl_app_apps.h"
#include "chemistry/bcl_chemistry_atom_vector.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "chemistry/bcl_chemistry_fragment_map_conformer.h"
#include "chemistry/bcl_chemistry_fragment_split_interface.h"
#include "chemistry/bcl_chemistry_hydrogens_handler.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "graph/bcl_graph_common_subgraph_isomorphism.h"
#include "graph/bcl_graph_const_graph.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ifstream.h"
#include "molecule/bcl_app_set_sample_by_parts_atoms.h"
#include "sdf/bcl_sdf_fragment_factory.h"
#include "storage/bcl_storage_template_instantiations.h"
#include "util/bcl_util_own_ptr.h"
namespace bcl
{
  namespace app
  {
    const ApplicationType SetSampleByPartsAtoms::SetSampleByPartsAtoms_Instance
    (
      GetAppGroups().AddAppToGroup( new SetSampleByPartsAtoms(), GetAppGroups().e_Molecule)
    );

    //! @brief standard constructor
    SetSampleByPartsAtoms::SetSampleByPartsAtoms() :
      m_OutputFilenameFlag
      (
        new command::FlagStatic
        (
          "output",
          "file to write split molecules into",
          command::Parameter
          (
            "output",
            "file to write split molecules into"
          )
        )
      ),
      m_AtomTypeDataFlag
      (
        new command::FlagStatic
        (
          "atom_type_comparison",
          "how to compare atom data",
          command::Parameter
          (
            "scheme",
            "used to compare atoms",
            command::ParameterCheckSerializable( chemistry::ConformationGraphConverter::AtomComparisonTypeEnum()),
            "ElementType"
          )
        )
      ),
      m_BondTypeDataFlag
      (
        new command::FlagStatic
        (
          "bond_type_comparison",
          "how to compare bond data",
          command::Parameter
          (
            "scheme",
            "used to compare bonds",
            command::ParameterCheckSerializable( chemistry::ConfigurationalBondTypeData::DataEnum()),
            "BondOrderAmideOrAromaticWithRingness"
          )
        )
      ),
      m_ReferenceFragmentFlag
      (
        new command::FlagStatic
        (
          "reference_mol",
          "reference molecule for the substructure comparison to input molecules",
          command::Parameter
          (
            "reference_mol",
            "used to determine which atom indices will be set in SampleByParts",
            command::ParameterCheckFileExistence(),
            ""
          )
        )
      ),
      m_DisableComplementFlag
      (
        new command::FlagStatic
        (
          "disable_complement_indices",
          "sets the SampleByParts indices to be the matched atom indices to the reference structure; "
          "by default this application will return the complement to the matched substructure atoms"
        )
      )
    {
    }

    //! copy constructor; skips i/o streams
    SetSampleByPartsAtoms::SetSampleByPartsAtoms( const SetSampleByPartsAtoms &PARENT) :
          m_OutputFilenameFlag( PARENT.m_OutputFilenameFlag),
          m_AtomTypeDataFlag( PARENT.m_AtomTypeDataFlag),
          m_BondTypeDataFlag( PARENT.m_BondTypeDataFlag),
          m_ReferenceFragmentFlag( PARENT.m_ReferenceFragmentFlag),
          m_DisableComplementFlag( PARENT.m_DisableComplementFlag)
    {
    }

    //! @brief Clone function
    //! @return pointer to new FoldProtein
    SetSampleByPartsAtoms *SetSampleByPartsAtoms::Clone() const
    {
      return new SetSampleByPartsAtoms( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &SetSampleByPartsAtoms::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> SetSampleByPartsAtoms::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

    ////////////////////
    // common options //
    ////////////////////

      // ensembles containing the molecules to be edited
      chemistry::FragmentFeed::AddFlags( *sp_cmd);

      // add command line options to add/remove hydrogens
      sdf::AddMoleculeIOPrefFlags( *sp_cmd);

      //! output filename
      sp_cmd->AddFlag( m_OutputFilenameFlag);
      sp_cmd->AddFlag( m_AtomTypeDataFlag);
      sp_cmd->AddFlag( m_BondTypeDataFlag);
      sp_cmd->AddFlag( m_ReferenceFragmentFlag);
      sp_cmd->AddFlag( m_DisableComplementFlag);

    ///////////////////
    // default flags //
    ///////////////////

      // default flags are unnecessary for this application, but message level and readme are useful
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags
        (
          *sp_cmd,
          storage::Set< command::FlagTypeEnum>( command::e_Pthread)
        );

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int SetSampleByPartsAtoms::Main() const
    {
      // open the output file
      io::File::MustOpenOFStream( m_OutputFile, m_OutputFilenameFlag->GetFirstParameter()->GetValue());

      // setup the reference molecule
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_ReferenceFragmentFlag->GetFirstParameter()->GetValue());
      m_ReferenceFragment = sdf::FragmentFactory::MakeFragment( input, sdf::e_Remove);
      io::File::CloseClearFStream( input);

      // setup data
      m_AtomComparison = chemistry::ConformationGraphConverter::AtomComparisonTypeEnum( m_AtomTypeDataFlag->GetFirstParameter()->GetValue());
      m_BondComparison = chemistry::ConfigurationalBondTypeData::DataEnum( m_BondTypeDataFlag->GetFirstParameter()->GetValue());
      chemistry::ConformationGraphConverter graph_converter( m_AtomComparison, m_BondComparison);

      // prepare graph of reference structure
      graph::ConstGraph< size_t, size_t> ref_mol_graph( graph_converter( m_ReferenceFragment));
      util::OwnPtr< graph::ConstGraph< size_t, size_t> > ref_mol_graph_ptr( &ref_mol_graph, false);

      // set the reference graph to be graph_A in the common subgraph isomorphism
      graph::CommonSubgraphIsomorphism< size_t, size_t> common_subgraph_iso;
      common_subgraph_iso.SetGraphA( ref_mol_graph_ptr);

      // iterate over all input molecules
      chemistry::FragmentFeed itr_fragments;
      bool has_h( false);
      if( m_ReferenceFragmentFlag->GetFlag())
      {
        // get matched atom indices
        for( ; itr_fragments.NotAtEnd(); ++itr_fragments)
        {
          // so we know how to setup output molecules
          if( itr_fragments->GetNumberHydrogens())
          {
            has_h = true;
          }

          // speed up isomorphism search by removing hydrogen atoms
          chemistry::AtomVector< chemistry::AtomComplete> new_mol_atoms_noh( itr_fragments->GetAtomVector());
          chemistry::HydrogensHandler::Remove( new_mol_atoms_noh);
          chemistry::FragmentComplete labeled_mol( new_mol_atoms_noh, "");

          // get the graph of our current molecule
          graph::ConstGraph< size_t, size_t> mol_graph( graph_converter( labeled_mol));
          util::OwnPtr< graph::ConstGraph< size_t, size_t> > mol_graph_ptr( &mol_graph, false);
          common_subgraph_iso.SetGraphB( mol_graph_ptr);

          // get the isomorphism betwene our current and reference molecules
          common_subgraph_iso.FindIsomorphism( common_subgraph_iso.EstimateUpperBounds(), 1);
          graph::Subgraph< size_t, size_t> subgraph_b
          (
            common_subgraph_iso.GetSubgraphIsomorphismsOfGraphB().FirstElement()
          );

          // get the moveable indices
          storage::Vector< size_t> moveable_indices;
          if( m_DisableComplementFlag->GetFlag())
          {
            // get the indices of the matched atoms
            moveable_indices = subgraph_b.GetVertexIndices();
          }
          else
          {
            // get the indices of the unmatched atoms
            moveable_indices = subgraph_b.GetComplement().GetVertexIndices();
          }

          // now set the MDL property on the current molecule
          BCL_MessageVrb( "Getting atom indices for conformer sampling...");
          labeled_mol.GetStoredPropertiesNonConst().SetMDLProperty( "SampleByParts", storage::Vector< size_t>( moveable_indices.Begin(), moveable_indices.End()));
          BCL_MessageVrb( "SampleByParts with the following atom indices: " + util::Format()( labeled_mol.GetMDLProperty( "SampleByParts")));

          // add back hydrogen atoms?
          if( sdf::GetCommandLineHydrogensPref() != sdf::e_Remove && has_h)
          {
            labeled_mol.SaturateWithH();
          }
          else if( sdf::GetCommandLineHydrogensPref() == sdf::e_Saturate)
          {
            labeled_mol.SaturateWithH();
          }

          // write out our molecule
          Write( labeled_mol);
        }
      }
      // if no reference specified, just output the input molecules
      else
      {
        for( ; itr_fragments.NotAtEnd(); ++itr_fragments)
        {
          chemistry::FragmentComplete fragment( *itr_fragments);
          Write( fragment);
        }
      }

      io::File::CloseClearFStream( m_OutputFile);

      // end
      return 0;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SetSampleByPartsAtoms::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &SetSampleByPartsAtoms::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief writes a single molecule out
    //! @param MOLECULE the molecule or fragment to write
    void SetSampleByPartsAtoms::Write( chemistry::FragmentComplete &MOLECULE) const
    {
      MOLECULE.WriteMDL( m_OutputFile);
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &SetSampleByPartsAtoms::GetReadMe() const
    {
      static std::string s_read_me =
        "SetSampleByPartsAtoms assigns atoms to the SampleByParts MDL property to be used with SampleConformations. "
        "Atoms to be sampled are determined by performing a substructure comparison of each input molecule to a reference molecule. "
        "The atom indices to be returned can either be the complement to the subgraph or the common subgraph atoms themselves.";
      return s_read_me;
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string SetSampleByPartsAtoms::GetDescription() const
    {
      return "Label molecules with atom indices to be sampled in SampleConformations based on common substructure to a reference molecule";
    }

  } // namespace app
} // namespace bcl
