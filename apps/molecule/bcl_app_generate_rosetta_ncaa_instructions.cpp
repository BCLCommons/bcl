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
#include "molecule/bcl_app_generate_rosetta_ncaa_instructions.h"
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
    util::ShPtr< command::Command> GenerateRosettaNCAAInstructions::InitializeCommand() const
    {
      // initialize command to be returned
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // ensembles containing the molecules to be filtered
       chemistry::FragmentFeed::AddFlags( *sp_cmd);

      // add command line options to add/remove hydrogens
      sdf::AddMoleculeIOPrefFlags( *sp_cmd);

      // add AtomMdlLine to molecule
      sdf::MdlHandler::AddAtomMdlLineFlag( *sp_cmd);

      // flags for input/output
      sp_cmd->AddFlag( m_OutputPrefixFlag);
      sp_cmd->AddFlag( m_Generate3DFlag);
      sp_cmd->AddFlag( m_SideChainSampleBypartsFlag);
      sp_cmd->AddFlag( m_ExtraPropertiesFlag);
      sp_cmd->AddFlag( m_ChiralityFlag);
      // Flag for output the partial charge file, which can be used to
      // assign partial charges to Rosetta params file
      sp_cmd->AddFlag( m_GeneratePartialChargeFileFlag);
      sp_cmd->AddFlag( m_CaAndChi1IndicesFlag);

      // default flags
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    } // InitializeCommand

    //! @brief the Main function
    //! @return 0 for success
    int GenerateRosettaNCAAInstructions::Main() const
    {
      // load ncaa base
      storage::Pair< bool, chemistry::FragmentComplete> glycine( ReadNCAABase());
      BCL_Assert( glycine.First(), "Someone altered the indices of the ncaa base neutral glycine reference structure file!!");

      // Input CA chirality
      std::string input_ca_chirality( m_ChiralityFlag->GetFirstParameter()->GetValue());
      // meh, may as well align by MCS
      static chemistry::ConformationComparisonPsiField ats;

      // generate graph of our ncaa base
      chemistry::ConformationGraphConverter graph_maker
      (
        chemistry::ConformationGraphConverter::e_ElementType,
        chemistry::ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness
      );
      graph::ConstGraph< size_t, size_t> base_graph( graph_maker( glycine.Second()));
      util::OwnPtr< graph::ConstGraph< size_t, size_t> > base_graph_ptr( &base_graph, false);

      // set the base molecule as the first graph for the isomorphism
      graph::CommonSubgraphIsomorphism< size_t, size_t> common_subgraph_iso( graph::CommonSubgraphIsomorphismBase::e_Connected);
      common_subgraph_iso.SetGraphA( base_graph_ptr);

      // input molecules
      chemistry::FragmentFeed feed;

      // go over all molecules from input
      for( ; feed.NotAtEnd(); ++feed)
      {
        // minimum atom count requirement
        if( feed->GetNumberAtoms() - feed->GetNumberHydrogens() < size_t( 6))
        {
          BCL_MessageStd( "Error: Not enough atoms for reliable substructure comparison to glycine backbone. Exiting...");
          continue;
        }
        size_t mol_index( feed.GetPosition());

        // initialize constant output data; these should remain the same so long as the ncaa base is not modified in the library
        // NOTE that these are based on 0-indexed atom IDs for coding purposes; in the output file we do +1 on all indices
        const size_t nter_index( 2), ca_index( 3), c_index( 1), o_index( 4),
            upper_n_index( 5), lower_c_index( 7);
        const storage::Vector< size_t> ignore_indices
        (
          storage::Vector< size_t>::Create
          (
            8,
            9,
            10,
            11,
            12,
            13,
            14,
            15,
            16,
            17
          )
        ); // consists of the capping atoms minus upper_n and lower_c

        // make molecule
        chemistry::FragmentComplete current_mol( *feed);
        // remove salts and whatnot
        chemistry::FragmentSplitLargestComponent splitter( size_t( 2));
        chemistry::FragmentEnsemble largest_component( splitter( current_mol));
        current_mol = largest_component.GetMolecules().FirstElement();

        // reachable will contains the indices of ncaa sc
        storage::Vector< size_t> reachable;
        size_t chi1_atom_index, mol_ca_index;
        // If user did not provide CA and Chi 1 index through command line
        // or user wants to determine the chirarity automatically

        if
        (
            !m_CaAndChi1IndicesFlag->GetFlag() ||
            m_ChiralityFlag->GetFirstParameter()->GetValue().compare( "auto") == 0
        )
        {
          // Need to keep Hydrogens to not mess up the indices of Ca and chi1
          // If the input of CA and Chi1 indices are provided
          if( !m_CaAndChi1IndicesFlag->GetFlag())
          {
            current_mol.RemoveH();
          }

          // align to ncaa base
          ats.ArbitraryScaffoldAlignment
          (
            current_mol,
            glycine.Second(),
            chemistry::ConformationGraphConverter::e_ElementType,
            chemistry::ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness
          );

          // make graph of molecule
          graph::ConstGraph< size_t, size_t> mol_graph( graph_maker( current_mol));
          util::OwnPtr< graph::ConstGraph< size_t, size_t> > mol_graph_ptr( &mol_graph, false);

          // set the base molecule as the second graph for the isomorphism
          common_subgraph_iso.SetGraphB( mol_graph_ptr);

          // find the largest subgraph isomorphism
          common_subgraph_iso.FindIsomorphism( common_subgraph_iso.EstimateUpperBounds(), size_t( 1));
          //common_subgraph_iso.FindIsomorphism( size_t( 6), size_t( 3));

          // obtain the mapped vertices between the ncaa base and the current molecule
          storage::Vector< storage::Map< size_t, size_t>> largest_iso( common_subgraph_iso.GetIsomorphisms());

          // get subgraph of our current molecule
          storage::Vector< graph::Subgraph< size_t, size_t> > mol_subgraph
          (
            common_subgraph_iso.GetSubgraphIsomorphismsOfGraphB()
          );

          // iterate over isomorphisms
          for
          (
            auto iso_itr( largest_iso.Begin()), iso_itr_end( largest_iso.End());
            iso_itr != iso_itr_end;
            ++iso_itr
          )
          {
            // we need the component off of the carbon alpha of our reference backbone
            // our CA is always indexed at 4
            size_t mol_ca_index( iso_itr->Find( size_t( 2))->second);

            //BCL_Debug( mol_ca_index);

            for
            (
                auto b_iso_itr( mol_subgraph.Begin()), b_iso_itr_end( mol_subgraph.End());
                b_iso_itr != b_iso_itr_end;
                ++b_iso_itr
            )
            {
              // get edges connecting subgraph to supergraph
              storage::List< storage::Pair< size_t, size_t> > subgraph_mol_edges( b_iso_itr->GetAdjacentEdgeIndices());
              //BCL_Debug( subgraph_mol_edges);
              for
              (
                  auto subgraph_mol_itr( subgraph_mol_edges.Begin()),
                  subgraph_mol_itr_end( subgraph_mol_edges.End());
                  subgraph_mol_itr != subgraph_mol_itr_end;
                  ++subgraph_mol_itr
              )
              {
                // only care about CA
                if( subgraph_mol_itr->First() != mol_ca_index)
                {
                  continue;
                }
                // find the part of the parent graph that contains the vertex subgraph_parent_itr->Second()
                // (subgraph vertex) but not subgraph_parent_itr->First(), or the sidechain subgraph
                chi1_atom_index = subgraph_mol_itr->Second();

                // If the user provided Ca and Chi 1 indices, check to see if it matches with the alignment output
                if( m_CaAndChi1IndicesFlag->GetFlag())
                {
                  // If the provided Ca index does not match with the alignment output, reject this isomorphism
                  if
                  (
                      m_CaAndChi1IndicesFlag->GetFirstParameter()->GetNumericalValue< size_t>() - 1 != mol_ca_index &&
                      m_CaAndChi1IndicesFlag->GetParameterList()( 1)->GetNumericalValue< size_t>() - 1 != chi1_atom_index
                  )
                  {
                    BCL_MessageStd
                    (
                      "Error: Molecule # " + util::Format()( mol_index) +
                      ": cannot automatically determine CA chirality using the input CA and Chi 1 indices." +
                      " Please try again with explicit input CA chirality of L_AA or D_AA. Skip this molecule..."
                    )
                    continue;
                  }
                }
                // If we need to determine the original Ca chirality automatically
                if( m_ChiralityFlag->GetFirstParameter()->GetValue().compare( "auto") == 0)
                {
                  size_t mol_n_index( iso_itr->Find( size_t( 1))->second);
                  //BCL_Debug( mol_n_index);
                  size_t mol_c_index( iso_itr->Find( size_t( 3))->second);
                  //BCL_Debug( mol_c_index);
                  std::ostringstream error_stream;
                  //chemistry::AtomVector< chemistry::AtomComplete> ncaa( current_mol.GetAtomVector());

                  input_ca_chirality =
                    FindCAChirarity
                    (
                      current_mol.GetAtomVector(),
                      mol_ca_index,
                      mol_c_index,
                      mol_n_index,
                      chi1_atom_index
                    );
                  BCL_MessageStd(
                    "Molecule # " + util::Format()( mol_index) +
                    ": Computed CA chirality from input structure is " + input_ca_chirality
                  );
                }

                util::ShPtr< storage::Vector< size_t> > reachable_vertices_to_keep
                (
                  graph::Connectivity::GetVerticesReachableFromDirectedEdge
                  (
                    mol_graph,
                    chi1_atom_index,
                    mol_ca_index
                  )
                );

                reachable = *reachable_vertices_to_keep;
                //BCL_Debug( reachable);
                break;
              }
              if( reachable.GetSize())
              {
                break;
              }
            } // edges
            if( reachable.GetSize())
            {
              break;
            }
          } // isomorphisms
        }
        // If the indices of CA and Chi 1 are provided by the user, and CA chirality is not auto
        else
        {
          mol_ca_index = m_CaAndChi1IndicesFlag->GetFirstParameter()->GetNumericalValue< size_t>() - 1;
          chi1_atom_index = m_CaAndChi1IndicesFlag->GetParameterList()( 1)->GetNumericalValue< size_t>() - 1;

          // make graph of molecule
          graph::ConstGraph< size_t, size_t> mol_graph( graph_maker( current_mol));
          util::ShPtr< storage::Vector< size_t> > reachable_vertices_to_keep
          (
            graph::Connectivity::GetVerticesReachableFromDirectedEdge
            (
              mol_graph,
              chi1_atom_index,
              mol_ca_index
            )
          );
          reachable = *reachable_vertices_to_keep;
        }
        if( !reachable.GetSize())
        {
           BCL_MessageStd
           (
             "Error: molecule # " + util::Format()( mol_index) + " No substructure match allows incorporation of non-canonical into backbone. "
             "May require manual addition of correct Chi 1 and CA atom indices for matching. Exiting..."
           );
           continue;
        }

        // New CA id = size of sidechain
        // New chi 1 index is the order of this atom in the sidechain
        size_t new_mol_ca_index( 0);
        size_t ncaa_sc_chi1_index( 0);
        if( util::IsDefined( chi1_atom_index))
        {
          for( size_t i( 0); i < reachable.GetSize(); ++i)
          {
            if( reachable( i) == chi1_atom_index)
            {
              ncaa_sc_chi1_index = i;
              break;
            }
            else
            {
              ++new_mol_ca_index;
            }
          }
        }
        else
        {
          BCL_MessageStd( "Error: Undefined atom attached to CA!, skip molecule # " + util::Format()( mol_index));
          continue;
        }

        // make sidechain fragment of ncaa from reachable atoms
        chemistry::AtomVector< chemistry::AtomComplete> ncaa_sc_v( current_mol.GetAtomVector());

        ncaa_sc_v.Reorder( reachable);
        chemistry::FragmentComplete ncaa_sc( ncaa_sc_v, "");

        // connect sidechain to backbone via AddMedChem
        m_AddMedChem.SetMutableAtomIndices( storage::Vector< size_t>( 1, new_mol_ca_index));
        storage::Triplet< bool, size_t, chemistry::FragmentComplete> glycine_dipeptide
        (
          ReadDipeptideBackbone( "ALPHA_AA")
        );
        if( !glycine_dipeptide.First())
        {
          BCL_MessageStd
          (
            " Cant read in ALPHA_AA backbone for molecule # " + util::Format()( mol_index)
          );
          continue;
        }
        storage::Triplet< chemistry::FragmentComplete, size_t, size_t> glycine_ov;
        if( input_ca_chirality.compare( "L_AA") == 0)
        {
          // open the 0-index hydrogen atom valence
          glycine_ov = storage::Triplet< chemistry::FragmentComplete, size_t, size_t>
          (
            m_AddMedChem.OpenValence
            (
              glycine_dipeptide.Third(),
              glycine_dipeptide.Second(),
              false,
              false
            )
          );
        }
        else if( input_ca_chirality.compare( "D_AA") == 0)
        {
          // open the 1-index hydrogen atom valence
          glycine_ov = storage::Triplet< chemistry::FragmentComplete, size_t, size_t>
          (
            m_AddMedChem.OpenValence
            (
              glycine_dipeptide.Third(),
              glycine_dipeptide.Second(),
              false,
              true
            )
          );
        }
        else
        {
          BCL_MessageStd
          (
            "Error: molecule " + util::Format()( mol_index) +
            " : only input CA chirality of L_AA and D_AA are supported for now."
          );
          continue;
        }
        // join the glycine fragment to our base molecule
        // the first atom of ncaa_sc should be the chi 1 atom
        storage::Pair< bool, chemistry::FragmentComplete> new_ncaa
        (
          chemistry::MergeFragmentComplete::MergeFragments
          (
            glycine_ov.First(),
            ncaa_sc,
            chemistry::GetConfigurationalBondTypes().e_NonConjugatedSingleBond,
            storage::Pair< size_t, size_t>( glycine_ov.Second(), ncaa_sc_chi1_index)
          )
        );
        if( !new_ncaa.First())
        {
          BCL_MessageStd( "Error: Unable to add sidechain from input molecule # " + util::Format()( mol_index) +
            " to backbone due to MergeFragments failure.");
          continue;
        }

        // clean at atom_vector level
        chemistry::AtomVector< chemistry::AtomComplete> new_mol_v( new_ncaa.Second().GetAtomVector());
        chemistry::AtomsCompleteStandardizer standardizer( new_mol_v, "", true);
        standardizer.SetConjugationOfBondTypes( new_mol_v);

        // add isometry info, which may be different for a subgraph
        chemistry::BondIsometryHandler::AddIsometryInformation( new_mol_v, true);
        // add stereocenter information
        chemistry::StereocentersHandler::AddChiralityFromConformation( new_mol_v);

        // add stereocenter information
        chemistry::FragmentComplete new_mol( new_mol_v, "");
        new_mol.SaturateWithH();

        // resample the movable indices
        if( m_SideChainSampleBypartsFlag->GetFlag())
        {
          // set movable indices; full backbone indices are the first 18 atoms
          storage::Set< size_t> movable_indices;
          movable_indices.InsertElement( size_t( 3));
          //movable_indices.InsertElement( new_mol.GetSize() - 1);
          for( size_t sc_i( 18); sc_i < new_mol.GetSize(); ++sc_i)
          {
            movable_indices.InsertElement( sc_i);
          }

          new_mol.GetStoredPropertiesNonConst().SetMDLProperty
              (
                "SampleByParts", storage::Vector< size_t>( movable_indices.Begin(), movable_indices.End())
              );
          BCL_MessageStd
          (
            "SampleByParts with the following atom indices: " + util::Format()( new_mol.GetMDLProperty( "SampleByParts"))
          );
        }

        // generate better 3D conf
        static chemistry::RotamerLibraryFile rotamer_library;
        static chemistry::FragmentMakeConformers make_one_conf
        (
          rotamer_library,
          2000,
          0.01,
          m_Generate3DFlag->GetFlag(),
          false
        );
        util::ShPtr< chemistry::FragmentComplete> gen_mol_3d_sp( make_one_conf.MakeConformer( new_mol));

        // get the chi1 index
        size_t chi1_index;
        if( gen_mol_3d_sp.IsDefined())
        {
          chi1_index = FindChi1Index
              (
                gen_mol_3d_sp->GetAtomVector(),
                ca_index,
                c_index,
                nter_index
              );
        }
        // get the formal charge of the final molecule
        const float formal_charge
        (
          descriptor::GetCheminfoProperties().calc_FormalCharge->SumOverObject( *gen_mol_3d_sp)( 0)
        );

        // compute the total numbers of aromatic rings
        const size_t aromatic_ring_num
        (
          descriptor::GetCheminfoProperties().calc_NAromaticRings->SumOverObject( ncaa_sc)( 0)
        );

        // compute the total numbers of rings
        const size_t ring_num
        (
          descriptor::GetCheminfoProperties().calc_NRings->SumOverObject( ncaa_sc)( 0)
        );

        // add chirality properties
        std::string chirality;
        if( util::IsDefined( chi1_index))
        {
          chirality = FindCAChirarity
              (
                gen_mol_3d_sp->GetAtomVector(),
                ca_index,
                c_index,
                nter_index,
                chi1_index
              );
        }
        if( chirality.compare( input_ca_chirality) != 0)
        {
          BCL_MessageStd
          (
            "Error: molecule " + util::Format()( mol_index) + " input CA chirarity " + input_ca_chirality +
            " is not the same with the output CA chirality: " + chirality
          );
          continue;
        };

        std::string properties_string
        (
          GetSidechainPropertiesList
          (
            formal_charge,
            ring_num,
            aromatic_ring_num,
            chirality,
            m_ExtraPropertiesFlag->GetFlag() ?
                m_ExtraPropertiesFlag->GetStringList() :
                storage::Vector< std::string>()
          )
        );

        // write the instructions file
        std::string instructions
        (
          WriteRosettaInstructions
          (
            nter_index,
            ca_index,
            c_index,
            o_index,
            chi1_index,
            ignore_indices,
            upper_n_index,
            lower_c_index,
            formal_charge,
            properties_string
          )
        );
        //gen_mol_3d_sp->RemoveProperty( "SampleByParts");
        gen_mol_3d_sp->GetStoredPropertiesNonConst().SetMDLProperty
            (
              "RosettaParamsInstructions",
              instructions
            );

        // output the final molecule
        io::OFStream mol_output;
        io::File::MustOpenOFStream
        (
          mol_output, m_OutputPrefixFlag->GetFirstParameter()->GetValue() + "_" + util::Format()( mol_index) + ".sdf"
        );
        if( gen_mol_3d_sp.IsDefined())
        {
          gen_mol_3d_sp->WriteMDL( mol_output);
        }
        io::File::CloseClearFStream( mol_output);

        // Generate partial charge file
        if( m_GeneratePartialChargeFileFlag->GetFlag())
        {
          // compute the total of sigma and pi charge per atom of the NCAA
          const linal::Vector< float> partial_charges
          (
            descriptor::GetCheminfoProperties().calc_TotalCharge->CollectValuesOnEachElementOfObject( *gen_mol_3d_sp)
          );
          if( partial_charges.GetSize() == gen_mol_3d_sp->GetSize())
          {
            // Output those partial charges with the element
            WritePartialChargeFile( mol_index, partial_charges, gen_mol_3d_sp->GetAtomVector());
          }
          else
          {
            BCL_MessageStd( "Error: Fails to generate partial charge file for molecule # " + util::Format()( mol_index));
          }
        }
      }
      return 0;
    } // Main

    //! @brief load neutral glycine residue from library
    //! @return the neutral glycine as the ncaa base
    const storage::Pair< bool, chemistry::FragmentComplete> GenerateRosettaNCAAInstructions::ReadNCAABase() const
    {
      // Begin
      chemistry::FragmentEnsemble glycine;

      // Read in neutral glycine file
      io::IFStream file;
      io::File::MustOpenIFStream
      (
        file,
        chemistry::RotamerLibraryFile::GetRotamerFinder().FindFile( "") + "ncaa_base/glycine_bb.sdf.gz"
        //chemistry::RotamerLibraryFile::GetRotamerFinder().FindFile( "") + "ncaa_base/glycine_neutral.sdf.gz"
      );
      glycine.ReadMoreFromMdl( file, sdf::e_Maintain);

      // return the glycine residue
      io::File::CloseClearFStream( file);

      // get the ncaa base (which will provide our reference to Rosetta for peptide backbone atoms)
      chemistry::FragmentComplete ncaa_base( glycine.GetMolecules().FirstElement());

      // make sure no one fucked with the reference file in a harmful way
      const chemistry::AtomVector< chemistry::AtomComplete> atom_v( ncaa_base.GetAtomVector());
      if
      (
          // central amide oxygen atom
          atom_v( 4).GetElementType() != chemistry::GetElementTypes().e_Oxygen ||

          // central amide nitrogen atom
          atom_v( 1).GetElementType() != chemistry::GetElementTypes().e_Nitrogen ||

          // alpha and beta carbon atoms
          atom_v( 3).GetElementType() != chemistry::GetElementTypes().e_Carbon ||
          atom_v( 2).GetElementType() != chemistry::GetElementTypes().e_Carbon ||

          // alpha and beta carbon hydrogen atoms
          atom_v( 6).GetElementType() != chemistry::GetElementTypes().e_Hydrogen ||
          atom_v( 7).GetElementType() != chemistry::GetElementTypes().e_Hydrogen
      )
      {
        return storage::Pair< bool, chemistry::FragmentComplete>
        (
          std::make_pair( bool( false), chemistry::FragmentComplete( ncaa_base))
        );
      }
      return storage::Pair< bool, chemistry::FragmentComplete>
      (
        std::make_pair( bool( true), chemistry::FragmentComplete( ncaa_base))
      );
    }

    //! @brief load neutral glycine dipeptide from library as backbone for alpha AA
    //! @return the neutral glycine dipeptide as the ncaa base
    const storage::Triplet< bool, size_t, chemistry::FragmentComplete> GenerateRosettaNCAAInstructions::ReadDipeptideBackbone
    (
      const std::string &BACKBONE_TYPE
    ) const
    {
      // Begin
      chemistry::FragmentEnsemble backbone;

      if( BACKBONE_TYPE.compare( "ALPHA_AA") == 0)
      {
        // Read in neutral glycine file
        io::IFStream file;
        io::File::MustOpenIFStream
        (
          file,
          chemistry::RotamerLibraryFile::GetRotamerFinder().FindFile( "") + "ncaa_base/glycine_neutral.sdf.gz"
        );
        backbone.ReadMoreFromMdl( file, sdf::e_Maintain);

        // return the glycine residue
        io::File::CloseClearFStream( file);

        // get the ncaa base (which will provide our reference to Rosetta for peptide backbone atoms)
        chemistry::FragmentComplete alpha_bb( backbone.GetMolecules().FirstElement());

        // make sure no one fucked with the reference file in a harmful way
        const chemistry::AtomVector< chemistry::AtomComplete> atom_v( alpha_bb.GetAtomVector());
        if
        (
            //central amide oxygen atom
            atom_v( 5).GetElementType() != chemistry::GetElementTypes().e_Oxygen ||

            // central amide nitrogen atom
            atom_v( 3).GetElementType() != chemistry::GetElementTypes().e_Nitrogen ||

            // alpha and beta carbon atoms
            atom_v( 4).GetElementType() != chemistry::GetElementTypes().e_Carbon ||
            atom_v( 2).GetElementType() != chemistry::GetElementTypes().e_Carbon ||

            // alpha and beta carbon hydrogen atoms
            atom_v( 0).GetElementType() != chemistry::GetElementTypes().e_Hydrogen ||
            atom_v( 1).GetElementType() != chemistry::GetElementTypes().e_Hydrogen
        )
        {
          return storage::Triplet< bool, size_t, chemistry::FragmentComplete>
          (
            bool( false), size_t(), chemistry::FragmentComplete()
          );
        }
        return storage::Triplet< bool, size_t, chemistry::FragmentComplete>
        (
          bool( true), size_t( 4), chemistry::FragmentComplete( alpha_bb)
        );
      }
      else
      {
        BCL_MessageStd( "Type of backbone are not supported. Allowed type: ALPHA_AA");
        return storage::Triplet< bool, size_t, chemistry::FragmentComplete>
        (
          bool( false), size_t(), chemistry::FragmentComplete()
        );
      }
    }

    //! @brief return the chi1 atom indices
    //! @param NCAA: the atom vector of NCAA
    //! @param C_INDEX: the index of backbone C atom
    //! @param N_INDEX: the index of the backbone N atom
    //! @return the chi1 atom index
    const size_t GenerateRosettaNCAAInstructions::FindChi1Index
    (
      const chemistry::AtomVector< chemistry::AtomComplete> &NCAA,
      const size_t &CA_INDEX,
      const size_t &C_INDEX,
      const size_t &N_INDEX
    ) const
    {
      // iterate over bonds of carbon alpha
      for
      (
          auto bond_itr( NCAA( CA_INDEX).GetBonds().Begin()),
          bond_itr_end( NCAA( CA_INDEX).GetBonds().End());
          bond_itr != bond_itr_end;
          ++bond_itr
      )
      {
        // skip backbone atoms
        const size_t &atom_index( NCAA.GetAtomIndex( bond_itr->GetTargetAtom()));
        if( atom_index == C_INDEX || atom_index == N_INDEX)
        {
          continue;
        }

        // if the next atom is a heavy atom then accept it as chi1;
        // i realize that if there are two heavy atoms then this is not perfect;
        // R/S will be assigned below a
        if( bond_itr->GetTargetAtom().GetElementType() != chemistry::GetElementTypes().e_Hydrogen)
        {
          return atom_index;
        }
      }
      return util::GetUndefinedSize_t();
    }

    //! brief return the correct type of NCAA backbone (so this can be added to the instruction file) later
    //! @parameter NCAA: the atom vectors of NCAA
    //! @parameter C_INDEX: the index of backbone C atom
    //! @parameter N_INDEX: the index of the backbone N atom
    //! @param CHI1_INDEX: the index of the chi1 angle atom
    const std::string GenerateRosettaNCAAInstructions::FindCAChirarity
    (
      const chemistry::AtomVector< chemistry::AtomComplete> &NCAA,
      const size_t &CA_INDEX,
      const size_t &C_INDEX,
      const size_t &N_INDEX,
      const size_t &CHI1_INDEX
    ) const
    {
      std::string ca_chirarity;
      // Purpose: ranking the priority among the sidechain branch, the C=O and the amine group

      //chemistry::StereocentersHandler::AddChiralityFromConformation( NCAA);

      //BCL_Debug( NCAA[ CA_INDEX]->GetAtomInfo());
      //BCL_Debug( NCAA[ CA_INDEX]->GetBonds());

      //BCL_Debug( NCAA[ CHI1_INDEX]->GetAtomInfo());
      //BCL_Debug( NCAA[ C_INDEX]->GetAtomInfo());
      //BCL_Debug( NCAA[ N_INDEX]->GetAtomInfo());
      if
      (
          GetChiralityName( NCAA[ CA_INDEX]->GetChirality()) == "R" ||
          GetChiralityName( NCAA[ CA_INDEX]->GetChirality()) == "S"
      )
      {
        // make a vector of N, C, chi1 atoms.
        // If this is L-AA, the order of those three atoms will be counter clockwise
        // To determine the order of those atoms in space, we can determine the sign of the determinant
        // of the matrix formed by the 3D coordinate of those three points
        // For the explanation, see https://math.stackexchange.com/questions/2400911/ordering-points-in-mathbbr3-using-the-sign-of-determinant
        storage::Vector< size_t> ca_branches( storage::Vector< size_t>::Create( N_INDEX, C_INDEX, CHI1_INDEX));
        // Create a matrix that will hold the x,y, and z coordinates for the three atoms of highest priority.
        linal::Matrix3x3< float> xyz_coordinates( 0.0); // make a matrix of size 3 X 3
        const linal::Vector3D &root_position( NCAA[ CA_INDEX]->GetPosition());

        // put the positions (relative to the root atom) of the 3 highest priority substituents into a matrix
        // the sign of the determinant of the matrix will give us clock-wise (D-AA) or counter-clockwise (L-AA)
        for( size_t index( 0), size( 3); index < size; ++index)
        {
          // get the atom position out of the vector, which has been sorted by priority
          const linal::Vector3D &position( NCAA[ ca_branches( index)]->GetPosition());

          // put the positions of this atom relative to the root atom into the rows of the matrix
          xyz_coordinates( index, 0) = position( 0) - root_position( 0);
          xyz_coordinates( index, 1) = position( 1) - root_position( 1);
          xyz_coordinates( index, 2) = position( 2) - root_position( 2);
        }

        // Calculate the determinant of a 3 X 3 matrix that has rows sorted in descending order of ca_branches.
        // Opposite orders will have opposite signs.
        const float determinant( xyz_coordinates.Determinant());
        //BCL_Debug( xyz_coordinates.Determinant());
        if( determinant > 0.0)
        {
          ca_chirarity = "L_AA";
        }
        else if( determinant < 0.0)
        {
          ca_chirarity = "D_AA";
        }
        else
        {
          BCL_MessageStd( "Warning: cannot determine whether the NCAA is L or D stereoisomer");
          ca_chirarity = "ACHIRAL_BACKBONE";
        }
      }
      else
      {
        BCL_MessageStd( "Warning: cannot determine whether the NCAA is L or D stereoisomer");
        ca_chirarity = "ACHIRAL_BACKBONE";
      }
      //BCL_Debug( ca_chirarity);
      return ca_chirarity;
    }

    //! brief write the final Rosetta instructions file
    const std::string GenerateRosettaNCAAInstructions::WriteRosettaInstructions
    (
      const size_t &NTER_INDEX,
      const size_t &CA_INDEX,
      const size_t &C_INDEX,
      const size_t &O_INDEX,
      const size_t &CHI1_INDEX,
      const storage::Vector< size_t> &IGNORE_INDICES,
      const size_t &UPPER_N_INDEX,
      const size_t &LOWER_C_INDEX,
      const float &FORMAL_CHARGE,
      const std::string &PROPERTIES
    ) const
    {
      // open output file
      std::stringstream out;

      // NOTE that the +1 is so that we correspond to SDF and Rosetta numbering
      out << "M  ROOT " << util::Format()( NTER_INDEX + 1) << std::endl; // the number of the N-terminal atom
      out << "M  POLY_N_BB " << util::Format()( NTER_INDEX + 1) << std::endl; // ditto, unless for some reason the most N terminal atom in your residue type is not N
      out << "M  POLY_CA_BB " << util::Format()( CA_INDEX + 1) << std::endl; // the CA atom number
      out << "M  POLY_C_BB " << util::Format()( C_INDEX + 1) << std::endl; // the C atom number
      out << "M  POLY_O_BB " << util::Format()( O_INDEX + 1) << std::endl; // the O atom number
      out << "M  POLY_IGNORE "; // all the atoms of the capping groups except UPPER and LOWER
      for( size_t i( 0), sz( IGNORE_INDICES.GetSize()); i < sz; ++i)
      {
        out << util::Format()( IGNORE_INDICES( i) + 1) << " ";
      }
      out << std::endl;
      out << "M  POLY_UPPER " << util::Format()( UPPER_N_INDEX + 1) << std::endl; // the nitrogen atom number of the C-terminal methyl amide
      out << "M  POLY_LOWER " << util::Format()( LOWER_C_INDEX + 1) << std::endl; // the carbonyl atom number of the N-terminal acetyl
      out << "M  POLY_CHG " << util::Format()( FORMAL_CHARGE) << std::endl; // the charge
      out << "M  POLY_PROPERTIES " << PROPERTIES << std::endl; // any properties, e.g. PROTEIN, CHARGED, L_AA, etc
      out << "M  END" << std::endl;
      // Done, we dont need to output cap atom file anymore
      return out.str();
    }

    //! brief output the final property list for the SIDECHAIN of NCAA
    const std::string GenerateRosettaNCAAInstructions::GetSidechainPropertiesList
    (
      const float &FORMAL_CHARGE,
      const size_t &RING_NUM,
      const size_t &AROMATIC_NUM,
      const std::string &CA_CHIRARITY,
      const storage::Vector< std::string> &EXTRAS
    ) const
    {
      storage::Map< std::string, std::string> poly_properties
      (
        storage::Map< std::string, std::string>::Create
        (
          // Let's just set the backbone to be alpha_AA for now.
          // Will add more option later
          std::pair< std::string, std::string>( "molecule", "PROTEIN"),
          std::pair< std::string, std::string>( "type", "POLYMER"),
          std::pair< std::string, std::string>( "backbone", "ALPHA_AA"),
          std::pair< std::string, std::string>( "cachirality", CA_CHIRARITY)
        )
      );

      // add poly properties based on molecule type
      if( FORMAL_CHARGE != float( 0.0))
      {
        poly_properties.Insert( std::pair< std::string, std::string>( "charge", "CHARGED"));
      }

      // add the sidechain ring info
      if( RING_NUM < 1)
      {
        poly_properties.Insert( std::pair< std::string, std::string>( "scring", "ALIPHATIC"));
      }
      else if( AROMATIC_NUM > 0)
      {
        poly_properties.Insert( std::pair< std::string, std::string>( "scring", "AROMATIC"));
      }
      else // empty string if has non-aromatic rings
      {
        poly_properties.Insert( std::pair< std::string, std::string>( "scring", ""));
      }

      // Combine EXTRA properties and computed properties
      if( !EXTRAS.IsEmpty())
      {
        for
        (
            auto prop_itr( EXTRAS.Begin()),
            prop_itr_end( EXTRAS.End());
            prop_itr != prop_itr_end;
            ++prop_itr
        )
        {
          std::string property( *prop_itr);
          // Only add the property if this is a new value
          if( poly_properties.GetMappedValues().Find( property) == poly_properties.GetMappedValues().GetSize())
          {
            // ignore input CA chirarity info if incorrect
            if( property == "L_AA" || property == "D_AA" || property == "ACHIRAL_BACKBONE")
            {
              BCL_MessageStd
              (
                "Input extra CA chirarity is " + property +
                ", but the computed CA chirarity is " + poly_properties.GetValue( "cachirality") +
                ". Ignore extra CA chirarity input"
              );
            }
            // ignore input backbone type info if incorrect
            else if
            (
                property == "PEPTOID" || property == "BETA_AA" || property == "GAMMA_AA" ||
                property == "R_PEPTOID" || property == "S_PEPTOID"
            )
            {
              BCL_MessageStd( "non alpha amino acid backbones are not supported. Ignore input property " + property);
            }
            // Check on number of rings and aromatic rings
            else if
            (
                ( property == "ALIPHATIC" || property == "AROMATIC") &&
                property != poly_properties.GetValue( "scring")
            )
            {
              std::string computed_scring;
              if( poly_properties.GetValue( "scring") == "")
              {
                computed_scring = "N/A";
              }
              else
              {
                computed_scring = poly_properties.GetValue( "scring");
              }
              BCL_MessageStd
              (
                "Input extra sidechain type is " + property +
                ", but the computed sidechain type is " + computed_scring +
                ". Ignore extra sidechain type input."
              );
            }
            else
            {
              poly_properties.Insert( std::pair< std::string, std::string>( "extra", property));
            }
          }
        }
      }
      bcl::storage::Vector< std::string> poly_properties_v
      (
        poly_properties.GetMappedValues()
      );
      // put all properties together into a string
      std::string ncaa_properties;
      for
      (
          auto itr( poly_properties_v.Begin()),
          itr_end( poly_properties_v.End());
          itr != itr_end;
          ++itr
      )
      {
        ncaa_properties.append( *itr);
        ncaa_properties.append( " ");
      }
      return ncaa_properties;
    }

    //! brief write the file that contains partial charge and element type for each atom in NCAA
    //! @parameter MOL_INDEX: index of the NCAA in the
    //! @parameter PARTIAL_CHARGES: current index of the C backbone atom
    //! @parameter ATOMS: current index of the N backbone atom
    void GenerateRosettaNCAAInstructions::WritePartialChargeFile
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
        m_OutputPrefixFlag->GetFirstParameter()->GetValue() + "_" + util::Format()( MOL_INDEX) + ".PartialCharges.txt"
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
    std::string GenerateRosettaNCAAInstructions::GetDescription() const
    {
      return "GenerateRosettaNCAAInstructions returns an instructions file from which to generate non-canonical amino acids in Rosetta";
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &GenerateRosettaNCAAInstructions::GetReadMe() const
    {
      static std::string s_read_me =
        "GenerateRosettaNCAAInstructions is a helper application to remove manual design of residue "
        "instruction files for generating non-canonical amino acids";
      return s_read_me;
    }

    //! @brief standard constructor
    GenerateRosettaNCAAInstructions::GenerateRosettaNCAAInstructions() :
      m_OutputPrefixFlag
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
      m_Generate3DFlag
      (
        new command::FlagStatic
        (
          "generate_3D",
          "de novo initialization coordinates of NCAA CA and sidechain atoms; "
          "may give better starting structure, but may not work for certain complex ring systems / macrocycles "
          "in the NCAA sidechain"
        )
      ),
      m_SideChainSampleBypartsFlag
      (
        new command::FlagStatic
        (
          "sidechain_sample_by_parts",
          "restrict sample by parts atoms to sidechain atoms plus carbon alpha; note that this will result in "
          "fewer rotamers if conformers are subsequently generated - specifically, the rotamers will be biased "
          "by the reference glycine backbone conformation."
        )
      ),
      m_ExtraPropertiesFlag
      (
        new command::FlagDynamic
        (
          "extra_properties",
          "additional properties to append to Rosetta params 'PROPERTIES' line",
          command::Parameter
          (
            "properties",
            "The application will override some basic incorrect extra properties with the correct ones."
          )
        )
      ),
      m_ChiralityFlag
      (
        new command::FlagStatic
        (
          "chirality",
          "specify whether the new amino acid should be 'L_AA' or 'D_AA' for alpha NCAAs",
          command::Parameter
          (
            "Sidechain chirality",
            "auto if the user want to determine the CA chirality automatically from input backbone or L_AA or D_AA.",
            command::ParameterCheckAllowed
            (
              storage::Vector< std::string>::Create( "L_AA", "D_AA", "auto")
            ),
            "auto"
          )
        )
      ),
      m_GeneratePartialChargeFileFlag
      (
        new command::FlagStatic
        (
          "generate_partial_charge_file",
          "Generate a file that includes the atomic partial charges; "
          "The partial charge is the total of sigma and pi charges for each atom in the dipeptide form of the NCAA"
        )
      ),
      m_CaAndChi1IndicesFlag
      (
        new command::FlagStatic
        (
          "ca_chi1_ids",
          "Specifies CA and Chi 1 indices to help identification of backbone and sidechain atoms of NCAAs",
          util::ShPtrVector< command::ParameterInterface>::Create
          (
            util::ShPtr< command::ParameterInterface>
            (
              new command::Parameter
              (
                "CA index",
                "Index of CA atom in the input NCAA",
                command::ParameterCheckRanged< size_t>( 0, std::numeric_limits< size_t>::max()),
                ""
              )
            ),
            util::ShPtr< command::ParameterInterface>
            (
              new command::Parameter
              (
                "Chi 1 index",
                "Index of CA atom in the input NCAA",
                command::ParameterCheckRanged< size_t>( 0, std::numeric_limits< size_t>::max()),
                ""
              )
            )
          )
        )
      )
    {
    }

    // Construct the static instance of this application, and add it to the ChemInfo group
    const ApplicationType GenerateRosettaNCAAInstructions::GenerateRosettaNCAAInstructions_Instance
    (
      GetAppGroups().AddAppToGroup( new GenerateRosettaNCAAInstructions(), GetAppGroups().e_Molecule)
    );

  } // namespace app
} // namespace bcl
