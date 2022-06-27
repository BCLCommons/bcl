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
#include "chemistry/bcl_chemistry_fragment_make_conformers.h"
#include "chemistry/bcl_chemistry_fragment_split_largest_component.h"
#include "chemistry/bcl_chemistry_fragment_split_rings_with_unsaturated_substituents.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "chemistry/bcl_chemistry_sample_conformations.h"
#include "chemistry/bcl_chemistry_stereocenters_handler.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "graph/bcl_graph_common_subgraph_isomorphism_base.h"
#include "graph/bcl_graph_const_graph.h"
#include "linal/bcl_linal_vector_3d.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "quality/bcl_quality_rmsd.h"
#include "random/bcl_random_distribution_interface.h"
#include "random/bcl_random_uniform_distribution.h"
#include "sdf/bcl_sdf_atom_info.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

#undef AddAtom
#undef ATOMS

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    // initialize static
    sched::Mutex &FragmentAlignToScaffold::GetMutex()
    {
      static sched::Mutex s_Mutex;
      return s_Mutex;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FragmentAlignToScaffold::FragmentAlignToScaffold() :
        m_AtomType( ConformationGraphConverter::AtomComparisonType::e_ElementType),
        m_BondType( ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness),
        m_MinIsoSize( size_t( 3))
    {
    }

    //! @brief constructor with binding pocket to be created from MDL property
    //! @param PROPERTY_SCORER enables pose-dependent optimization of score
    //! @param BFACTORS vector of values indicating per-residue flexibility (higher values are more flexible)
    //! @param MDL the SDF file MDL property specifying the binding pocket filename
    //! @param BINDING_POCKET_FILENAME name of protein/protein binding pocket PDB file
    FragmentAlignToScaffold::FragmentAlignToScaffold
    (
      const ConformationGraphConverter::AtomComparisonType &ATOM_TYPE,
      const ConfigurationalBondTypeData::Data &BOND_TYPE,
      const size_t &MIN_ISO_SIZE
    ) :
        m_AtomType( ATOM_TYPE),
        m_BondType( BOND_TYPE),
        m_MinIsoSize( MIN_ISO_SIZE)
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
      ConformationGraphConverter::AtomComparisonType atm_cmp( m_AtomType);
      ConfigurationalBondTypeData::Data bnd_cmp( m_BondType);
//      ConformationGraphConverter::AtomComparisonType atm_cmp( ConformationGraphConverter::e_ElementType);
//      ConfigurationalBondTypeData::Data bnd_cmp( ConfigurationalBondTypeData::e_Identity);
      ConformationGraphConverter arbitrary_graph_converter( atm_cmp, bnd_cmp, true);

      // get the atoms of each mol
      AtomVector< AtomComplete> atom_v_a( TARGET_MOL.GetAtomVector());
      AtomVector< AtomComplete> atom_v_b( SCAFFOLD_MOL.GetAtomVector());

      // Remove everything except specified indices
      if( TARGET_MOL_INDICES.GetSize())
      {
        atom_v_a.Reorder( TARGET_MOL_INDICES);
      }
      if( SCAFFOLD_MOL_INDICES.GetSize())
      {
        atom_v_b.Reorder( SCAFFOLD_MOL_INDICES);
      }

      // convert back to molecules
      FragmentComplete target_mol( atom_v_a, "");
      FragmentComplete scaffold_mol( atom_v_b, "");

      graph::ConstGraph< size_t, size_t>
          mol_b_graph( arbitrary_graph_converter( scaffold_mol)),
          mol_a_graph( arbitrary_graph_converter( target_mol));

      // Store the coordinates of the scaffold atoms for alignment
      auto mol_b_coords( scaffold_mol.GetHeavyAtomCoordinates()), mol_a_coords( target_mol.GetHeavyAtomCoordinates());

      graph::CommonSubgraphIsomorphism< size_t, size_t> csi_substructure
      (
        graph::CommonSubgraphIsomorphismBase::e_GreedyUnconnected
      );

      csi_substructure.SetGraphA
      (
        util::OwnPtr< graph::ConstGraph< size_t, size_t> >( &mol_a_graph, false)
      );
      csi_substructure.SetGraphB
      (
        util::OwnPtr< graph::ConstGraph< size_t, size_t> >( &mol_b_graph, false)
      );

      //Find isomorphism between the scaffold graph and the small molecule graph
      csi_substructure.FindIsomorphism( csi_substructure.EstimateUpperBounds());
      storage::Map< size_t, size_t> isomorphism( csi_substructure.GetIsomorphism());

      if( isomorphism.GetSize() < m_MinIsoSize)
      {
        return false;
      }

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

      //Store transformed coordinates and atom information in a new molecule and add it to output ensemble
      FragmentComplete new_molecule
      (
        AtomVector< AtomComplete>( atom_vector, TARGET_MOL.GetBondInfo()),
        TARGET_MOL.GetName()
      );

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
        scores( score_index).First() = ( AlignToScaffold( *conf_itr, SCAFFOLD_MOL));
        scores( score_index).Second() = COMPARER->operator ()( *conf_itr, SCAFFOLD_MOL);
        conf_itr->StoreProperty( COMPARER->GetAlias(), linal::Vector< double>( 1, scores( score_index).Second()));
        conf_itr->StoreProperty( "ConformerIndex", linal::Vector< double>( 1, score_index));
      }

      // end successfully
      return scores;
    }

  //////////////////////
  // helper functions //
  //////////////////////

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
