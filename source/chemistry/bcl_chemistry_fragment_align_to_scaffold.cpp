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
#include <chemistry/bcl_chemistry_fragment_split_rings.h>
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "chemistry/bcl_chemistry_fragment_align_to_scaffold.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_complete.h"
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_bond_isometry_handler.h"
#include "chemistry/bcl_chemistry_bond_lengths.h"
#include "chemistry/bcl_chemistry_conformation_comparison_by_rmsd.h"
#include "chemistry/bcl_chemistry_conformation_comparison_property_field_correlation.h"
#include "chemistry/bcl_chemistry_conformation_comparison_property_rmsd_x.h"
#include "chemistry/bcl_chemistry_conformation_comparison_psi_field.h"
#include "chemistry/bcl_chemistry_conformation_interface.h"
#include "chemistry/bcl_chemistry_fragment_evolve_base.h"
#include "chemistry/bcl_chemistry_fragment_map_conformer.h"
#include "chemistry/bcl_chemistry_fragment_make_conformers.h"
#include "chemistry/bcl_chemistry_fragment_split_largest_component.h"
#include "chemistry/bcl_chemistry_fragment_split_rings_with_unsaturated_substituents.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "chemistry/bcl_chemistry_sample_conformations.h"
#include "chemistry/bcl_chemistry_stereocenters_handler.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "graph/bcl_graph_common_subgraph_isomorphism_base.h"
#include "graph/bcl_graph_const_graph.h"
#include "graph/bcl_graph_common_subgraph_isomorphism.h"
#include "graph/bcl_graph_subgraph.h"
#include "linal/bcl_linal_vector_3d.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "quality/bcl_quality_rmsd.h"
#include "random/bcl_random_distribution_interface.h"
#include "random/bcl_random_uniform_distribution.h"
#include "sdf/bcl_sdf_atom_info.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_own_ptr.h"

// external includes - sorted alphabetically

#undef AddAtom  // for Windows compilation
#undef ATOMS    // for windows compilation

namespace bcl
{
  namespace chemistry
  {

    //////////
    // data //
    //////////


    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    //! @brief default constructor
    FragmentAlignToScaffold::FragmentAlignToScaffold() :
        m_AtomType( ConformationGraphConverter::AtomComparisonType::e_ElementType),
        m_BondType( ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness),
        m_MinIsoSize( size_t( 3)),
        m_SolutionType( graph::CommonSubgraphIsomorphismBase::e_GreedyUnconnected)
    {
    }

    //! @brief constructor with graph isomorphism settings
    //! @param ATOM_TYPE the atom comparison type to be used for substructure comparison
    //! @param BOND_TYPE the bond comparison type to be used for substructure comparison
    //! @param MIN_ISO_SIZE the minimum size a substructure can be and be a solution
    //! @param SOLUTION_TYPE the solution type for the graph isomorphism
    FragmentAlignToScaffold::FragmentAlignToScaffold
    (
      const ConformationGraphConverter::AtomComparisonType &ATOM_TYPE,
      const ConfigurationalBondTypeData::Data &BOND_TYPE,
      const size_t &MIN_ISO_SIZE,
      const graph::CommonSubgraphIsomorphismBase::SolutionType &SOLUTION_TYPE
    ) :
        m_AtomType( ATOM_TYPE),
        m_BondType( BOND_TYPE),
        m_MinIsoSize( MIN_ISO_SIZE),
        m_SolutionType( SOLUTION_TYPE)
    {
    }

    //! @brief clone constructor
    FragmentAlignToScaffold *FragmentAlignToScaffold::Clone() const
    {
      return new FragmentAlignToScaffold( *this);
    }

    /////////////////
    // data access //
    /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FragmentAlignToScaffold::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the atom type used for comparison
    const ConformationGraphConverter::AtomComparisonType &FragmentAlignToScaffold::GetAtomType() const
    {
      return m_AtomType;
    }

    //! @brief return the bond type used for comparison
    const ConfigurationalBondTypeData::Data &FragmentAlignToScaffold::GetBondType() const
    {
      return m_BondType;
    }

    //! @brief return the size of the minimum allowed isomorphism
    const size_t FragmentAlignToScaffold::GetMinIsoSize() const
    {
      return m_MinIsoSize;
    }

    //! @brief return the solution type used for comparison
    const graph::CommonSubgraphIsomorphismBase::SolutionType &FragmentAlignToScaffold::GetSolutionType() const
    {
      return m_SolutionType;
    }

    //! @brief set the atom type used for comparison
    void FragmentAlignToScaffold::SetAtomType( const ConformationGraphConverter::AtomComparisonType &ATOM_TYPE)
    {
      m_AtomType = ATOM_TYPE;
    }

    //! @brief set the bond type used for comparison
    void FragmentAlignToScaffold::SetBondType( const ConfigurationalBondTypeData::Data &BOND_TYPE)
    {
      m_BondType = BOND_TYPE;
    }

    // @brief set the minimum isomorphism size
    void FragmentAlignToScaffold::SetMinIsoSize( const size_t MIN_ISO_SIZE)
    {
      m_MinIsoSize = MIN_ISO_SIZE;
    }

    //! @brief set the solution type used for comparison
    void FragmentAlignToScaffold::SetSolutionType( const graph::CommonSubgraphIsomorphismBase::SolutionType &SOLUTION_TYPE)
    {
      m_SolutionType = SOLUTION_TYPE;
    }

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //! @brief align small molecules by substructure
    //! @param TARGET_MOL the molecule to be aligned
    //! @param SCAFFOLD_MOL the molecule to which TARGET_MOL will be aligned
    //! @return return false if the isomorphism is less
    bool FragmentAlignToScaffold::AlignToScaffold
    (
      FragmentComplete &TARGET_MOL,
      const FragmentComplete &SCAFFOLD_MOL,
      const storage::Vector< size_t> &TARGET_MOL_INDICES,
      const storage::Vector< size_t> &SCAFFOLD_MOL_INDICES
    ) const
    {

      // extract the fragment for which the MCS search will be performed
      FragmentComplete target_mol( ExtractFragmentByIndices( TARGET_MOL, TARGET_MOL_INDICES));
      FragmentComplete scaffold_mol( ExtractFragmentByIndices( SCAFFOLD_MOL, SCAFFOLD_MOL_INDICES));

      // Store the coordinates of the fragments for alignment
      util::SiPtrVector< const linal::Vector3D> mol_b_coords( scaffold_mol.GetHeavyAtomCoordinates()), mol_a_coords( target_mol.GetHeavyAtomCoordinates());

      // get the isomorphism between the target and scaffold
      graph::CommonSubgraphIsomorphism< size_t, size_t> csi_substructure( FindMaximumCommonSubstructure( target_mol, scaffold_mol ) );
      csi_substructure.FindIsomorphism( csi_substructure.EstimateUpperBounds());
      storage::Map< size_t, size_t> isomorphism( csi_substructure.GetIsomorphism());

      // alignment fails if there is no acceptable isomorphism
      if( isomorphism.GetSize() < m_MinIsoSize)
      {
        return false;
      }

      // trim the fragments being compared to just their matched subgraph components
      mol_a_coords.Reorder( isomorphism.GetKeysAsVector());
      mol_b_coords.Reorder( isomorphism.GetMappedValues());

      // Generate transformation matrix based on common isomorphism between scaffold and current small molecule
      math::TransformationMatrix3D transform
      (
        quality::RMSD::SuperimposeCoordinates( mol_b_coords, mol_a_coords)
      );

      //Store the atom information of the small molecule
      storage::Vector< sdf::AtomInfo> atom_vector( TARGET_MOL.GetAtomInfo());

      //Transform the coordinates of the small molecule atoms based on the transformation matrix
      for
      (
        storage::Vector< sdf::AtomInfo>::iterator itr( atom_vector.Begin()), itr_end( atom_vector.End());
        itr != itr_end;
        ++itr
      )
      {
        linal::Vector3D temp( itr->GetCoordinates());
        itr->SetCoordinates( temp.Transform( transform));
      }

      // Store transformed coordinates and atom information in a new molecule and add it to output ensemble
      FragmentComplete new_molecule
      (
        AtomVector< AtomComplete>( atom_vector, TARGET_MOL.GetBondInfo()),
        TARGET_MOL.GetName()
      );

      // add back the original properties to the re-aligned structure
      new_molecule.StoreProperties( TARGET_MOL);
      TARGET_MOL = new_molecule;
      return true;
    }

    //! @brief maximum common substructure alignment of small molecules with pose-dependent scoring
    //! @param TARGET_MOL the molecule to be aligned
    //! @param SCAFFOLD_MOL the molecule against which the target is aligned
    //! @param BFACTORS vector of values indicating per-residue flexibility (higher values are more flexible)
    //! @param MDL the SDF file MDL property specifying the binding pocket filename
    //! @param BINDING_POCKET_FILENAME name of protein/protein binding pocket PDB file
    storage::Pair< bool, float> FragmentAlignToScaffold::PoseSensitiveAlignToScaffold
    (
      FragmentComplete &TARGET_MOL,
      const FragmentComplete &SCAFFOLD_MOL,
      const descriptor::CheminfoProperty &SCORE,
      const std::string &MDL,
      const std::string &BINDING_POCKET_FILENAME
    )
    {
      // set the necessary descriptors if needed
      if
      (
          !TARGET_MOL.GetStoredPropertiesNonConst().GetMDLProperty( MDL).size() ||
          TARGET_MOL.GetStoredPropertiesNonConst().GetMDLProperty( MDL) != BINDING_POCKET_FILENAME
      )
      {
        TARGET_MOL.GetStoredPropertiesNonConst().SetMDLProperty( MDL, BINDING_POCKET_FILENAME);
      }

      // align the molecule
      bool aligned( AlignToScaffold( TARGET_MOL, SCAFFOLD_MOL));

      // return the score of the new pose
      return std::make_pair( aligned, SCORE->SumOverObject( TARGET_MOL)( 0));
    }

    //! @brief maximum common substructure alignment of small molecule conformer ensembles
    //! @param TARGET_ENS the molecule ensemble to be aligned
    //! @param SCAFFOLD_MOL the molecule against which the targets are aligned
    //! @param COMPARER the metric to be used to compare alignments
    //! @return true if alignment occurs, false if the isomorphism size
    //! is below the minimum size allowed or the ensemble is empty
    storage::Vector< storage::Pair< bool, float> > FragmentAlignToScaffold::AlignEnsembleToScaffold
    (
      FragmentEnsemble &TARGET_ENS,
      const FragmentComplete &SCAFFOLD_MOL,
      const util::Implementation< ConformationComparisonInterface> &COMPARER
    ) const
    {
      // if empty ensemble then return false
      if( !TARGET_ENS.GetSize())
      {
        return false;
      }

      // output alignment success and MolAlign score
      storage::Vector< storage::Pair< bool, float> > scores( TARGET_ENS.GetSize());

      // iterate over conformational ensemble and score confs
      size_t score_index( 0);

      // align all conformers
      for
      (
        auto conf_itr( TARGET_ENS.Begin()), conf_itr_end( TARGET_ENS.End());
        conf_itr != conf_itr_end;
        ++conf_itr, ++score_index
      )
      {
        // align and score
        // TODO: if this is a conformer ensemble of the same molecule then graph comparison only needs to be done once
        scores( score_index).First() = ( AlignToScaffold( *conf_itr, SCAFFOLD_MOL));
        scores( score_index).Second() = COMPARER->operator ()( *conf_itr, SCAFFOLD_MOL);
        conf_itr->StoreProperty( COMPARER->GetAlias(), linal::Vector< double>( 1, scores( score_index).Second()));
        conf_itr->StoreProperty( "ConformerIndex", linal::Vector< double>( 1, score_index));
      }

      // end successfully
      return scores;
    }

    //! @brief build a new conformer of the target molecule starting from the conformation of the largest shared substructure with a scaffold molecule
    //! @param TARGET_MOL the molecule for which a new conformer will be generated
    //! @param SCAFFOLD_MOL the molecule whose MCS with TARGET_MOL will be the core of the new conformation for TARGET_MOL
    //! @param TARGET_MOL_INDICES atoms in TARGET_MOL that will be searched for a common substructure
    //! @param SCAFFOLD_MOL_INDICES atoms in SCAFFOLD_MOL that will be searched for a common substructure
    //! @return true if the conformer is generated successfully, false otherwise
    bool FragmentAlignToScaffold::GenerateConformerBasedOnScaffold
    (
      FragmentComplete &TARGET_MOL,
      const FragmentComplete &SCAFFOLD_MOL,
      const storage::Vector< size_t> &TARGET_MOL_INDICES,
      const storage::Vector< size_t> &SCAFFOLD_MOL_INDICES
    ) const
    {
      // extract the fragment for which the MCS search will be performed
      FragmentComplete target_mol( ExtractFragmentByIndices( TARGET_MOL, TARGET_MOL_INDICES));
      FragmentComplete scaffold_mol( ExtractFragmentByIndices( SCAFFOLD_MOL, SCAFFOLD_MOL_INDICES));

      // get the isomorphism between the target and scaffold
      graph::CommonSubgraphIsomorphism< size_t, size_t> csi_substructure( FindMaximumCommonSubstructure( target_mol, scaffold_mol ) );
      csi_substructure.FindIsomorphism( csi_substructure.EstimateUpperBounds());
      storage::Map< size_t, size_t> isomorphism( csi_substructure.GetIsomorphism( ) );

      // get the subgraph isomorphism for the target molecule
      graph::Subgraph< size_t, size_t> target_subgraph
      (
        util::OwnPtr< const graph::ConstGraph< size_t, size_t> >( &csi_substructure.GetGraphA(), false),
        csi_substructure.GetIsomorphism().GetKeysAsVector()
      );

      // alignment fails if there is no acceptable isomorphism
      if( isomorphism.GetSize() < m_MinIsoSize)
      {
        return false;
      }

      // for each atom in the isomorphism, set the target molecule atom coordinates to be those of the scaffold molecule
      storage::Vector< sdf::AtomInfo> target_atoms( TARGET_MOL.GetAtomInfo( ) );
      const storage::Vector< sdf::AtomInfo> &scaffold_atoms( SCAFFOLD_MOL.GetAtomInfo());
      for( auto iso_itr( isomorphism.Begin()), iso_itr_end( isomorphism.End()); iso_itr != iso_itr_end; ++iso_itr)
      {
        target_atoms( iso_itr->first).SetCoordinates( scaffold_atoms( iso_itr->second).GetCoordinates());
      }

      // get inverted subgraph of the new mol
      graph::Subgraph< size_t, size_t> target_subgraph_complement( target_subgraph.GetComplement( ) );
      if( target_subgraph_complement.GetSize())
      {
        // generate a new 3D conformer without sampling the coordinates of the subgraph isomorphism
        storage::Vector< size_t> moveable_atoms( target_subgraph_complement.GetVertexIndices());
        storage::Set< size_t> moveable_atoms_set( moveable_atoms.Begin(), moveable_atoms.End());
        FragmentMapConformer cleaner( "None", false, moveable_atoms);

        // add all the adjacent edges to our unique subgraph vertices
        storage::Set< size_t> all_adjacent_indices( cleaner.MapSubgraphAdjacentAtoms( target_subgraph_complement, 2) );
        moveable_atoms_set.InsertElements( all_adjacent_indices.Begin(), all_adjacent_indices.End());

        // add bad geometry atoms
        AtomVector< AtomComplete> temp_vec( target_atoms, TARGET_MOL.GetBondInfo());
        FragmentComplete temp_mol( temp_vec, TARGET_MOL.GetName());
        storage::Vector< size_t> bad_geo_atoms( temp_mol.GetAtomsWithBadGeometry());
        moveable_atoms_set.InsertElements( bad_geo_atoms.Begin(), bad_geo_atoms.End());

        ConformationGraphConverter::t_AtomGraph molecule_graph( ConformationGraphConverter::CreateGraphWithAtoms( temp_mol));
        FragmentSplitRings splitter( true, size_t( 3));
        storage::List< storage::Vector< size_t> > ring_components( splitter.GetComponentVertices( temp_mol, molecule_graph));
        storage::Vector< storage::Vector< size_t>> ring_components_vec( ring_components.Begin(), ring_components.End());
        for( size_t ring_i( 0), n_rings( ring_components_vec.GetSize()); ring_i < n_rings; ++ring_i) {
          storage::Set< size_t> ring_atoms( ring_components_vec( ring_i).Begin(), ring_components_vec( ring_i).End());
          for( auto bad_geo_atoms_itr( bad_geo_atoms.Begin()), bad_geo_atoms_itr_end( bad_geo_atoms.End()); bad_geo_atoms_itr != bad_geo_atoms_itr_end; ++bad_geo_atoms_itr)
          {
            if (ring_atoms.Find( *bad_geo_atoms_itr) != ring_atoms.End())
            {
              // if found, then this ring contains bad geometry atoms and we need to add the whole ring to the moveable atoms set
              moveable_atoms_set.InsertElements( ring_atoms.Begin(), ring_atoms.End());
              break;
            }
          }
        }

        // collect unique atoms to be mobile
        moveable_atoms = storage::Vector< size_t>( moveable_atoms_set.Begin(), moveable_atoms_set.End());
        cleaner.SetMoveableAtomIndices( moveable_atoms);

        // create cleaned molecule
        util::ShPtr< FragmentComplete> new_molecule
        (
          cleaner.Clean
          (
            AtomVector< AtomComplete>( target_atoms, TARGET_MOL.GetBondInfo() ),
            FragmentComplete(),
            "None",
            true
          )
        );

        if( new_molecule.IsDefined())
        {
          // re-align the subgraph isomorphism to the desired coordinates; if the subgraph is not the largest substructure then
          // it will not be re-positioned in real-space by default
          util::SiPtrVector< const linal::Vector3D> scaffold_subgraph_coords( scaffold_mol.GetHeavyAtomCoordinates()), target_subraph_coords( new_molecule->GetHeavyAtomCoordinates());
          target_subraph_coords.Reorder( isomorphism.GetKeysAsVector());
          scaffold_subgraph_coords.Reorder( isomorphism.GetMappedValues());

          math::TransformationMatrix3D transform
          (
            quality::RMSD::SuperimposeCoordinates( scaffold_subgraph_coords, target_subraph_coords)
          );

          storage::Vector< sdf::AtomInfo> new_molecule_atom_vector( new_molecule->GetAtomInfo());
          for
          (
            storage::Vector< sdf::AtomInfo>::iterator itr( new_molecule_atom_vector.Begin()),
            itr_end( new_molecule_atom_vector.End());
            itr != itr_end;
            ++itr
          )
          {
            linal::Vector3D coords( itr->GetCoordinates());
            itr->SetCoordinates( coords.Transform( transform));
          }

          // add back the original properties to the re-aligned structure
          AtomVector< AtomComplete> atoms_transformed( new_molecule_atom_vector, new_molecule->GetBondInfo());
          FragmentComplete new_molecule_transformed( atoms_transformed, TARGET_MOL.GetName());
          new_molecule_transformed.StoreProperties( TARGET_MOL);
          new_molecule_transformed.GetStoredPropertiesNonConst().SetMDLProperty( "SampleByParts", moveable_atoms);
          TARGET_MOL = new_molecule_transformed;
          return true;
        }
      }

      // it is possible for the entire substructure to be contained within the scaffold such that there are no complement subgraph atoms
      else
      {
        // clean atoms
        AtomVector< AtomComplete> clean_atoms(
          FragmentMapConformer::CleanAtoms
          (
            AtomVector< AtomComplete>( target_atoms, TARGET_MOL.GetBondInfo() ),
            "None",
            true
          )
        );
        if( clean_atoms.GetSize())
        {
          FragmentComplete new_molecule( clean_atoms, TARGET_MOL.GetName());
          new_molecule.StoreProperties( TARGET_MOL);
          new_molecule.GetStoredPropertiesNonConst().SetMDLProperty( "SampleByParts", storage::Vector< size_t>() );
          TARGET_MOL = new_molecule;
          return true;
        }
      }
      return false;
    }

    //////////////////////
    // helper functions //
    //////////////////////

    //! @brief extract a fragment from a molecule based on its indices
    FragmentComplete FragmentAlignToScaffold::ExtractFragmentByIndices( const FragmentComplete &MOLECULE, const storage::Vector< size_t> &INDICES)
    {
      AtomVector< AtomComplete> atoms( MOLECULE.GetAtomVector());
      if( INDICES.GetSize())
      {
        atoms.Reorder( INDICES);
      }
      FragmentComplete molecule( atoms, MOLECULE.GetName());
      return molecule;
    }

    //! @brief get the maximum common substructure between two molecules
    graph::CommonSubgraphIsomorphism< size_t, size_t> FragmentAlignToScaffold::FindMaximumCommonSubstructure( const FragmentComplete &MOL_A, const FragmentComplete &MOL_B) const
    {
      const ConformationGraphConverter arbitrary_graph_converter( m_AtomType, m_BondType, true );
      graph::ConstGraph< size_t, size_t>
          mol_a_graph( arbitrary_graph_converter( MOL_A)),
          mol_b_graph( arbitrary_graph_converter( MOL_B));

      graph::CommonSubgraphIsomorphism< size_t, size_t> csi_substructure( m_SolutionType);
      csi_substructure.SetGraphs( mol_a_graph, mol_b_graph);
      return csi_substructure;
    }

    //! @brief get the common substructures between two molecules
//    storage::Vector< storage::Map< size_t, size_t> > FragmentAlignToScaffold::FindAllCommonSubstructures( const FragmentComplete &MOL_A, const FragmentComplete &MOL_B) const
//    {
//      const ConformationGraphConverter arbitrary_graph_converter( m_AtomType, m_BondType, true);
//      graph::ConstGraph< size_t, size_t>
//          mol_a_graph( arbitrary_graph_converter( MOL_A)),
//          mol_b_graph( arbitrary_graph_converter( MOL_B));
//
//      graph::SubgraphIsomorphism< size_t, size_t> csi_substructure;
//      csi_substructure.SetGraphs( mol_a_graph, mol_b_graph);
//      csi_substructure.FindAllIsomorphisms();  // TODO add options for these
//      storage::Vector< storage::Map< size_t, size_t> >isomorphisms( csi_substructure.GetIsomorphisms( ) );
//      return isomorphisms;
//    }

    // initialize static
    sched::Mutex &FragmentAlignToScaffold::GetMutex()
    {
      static sched::Mutex s_Mutex;
      return s_Mutex;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FragmentAlignToScaffold::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &FragmentAlignToScaffold::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
