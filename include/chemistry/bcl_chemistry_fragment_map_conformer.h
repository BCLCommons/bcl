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

#ifndef BCL_CHEMISTRY_FRAGMENT_MAP_CONFORMER_H_
#define BCL_CHEMISTRY_FRAGMENT_MAP_CONFORMER_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically
#include "find/bcl_find.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_conformational_interface.h"
#include "bcl_chemistry_collector_valence.h"
#include "bcl_chemistry_fragment_align_to_scaffold.h"
#include "bcl_chemistry_fragment_complete.h"
#include "bcl_chemistry_fragment_constitution_shared.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "descriptor/bcl_descriptor_base.h"
#include "find/bcl_find_pick_interface.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_result.h"
#include "util/bcl_util_function_interface.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_si_ptr_list.h"
// external includes - sorted alphabetically

#undef AddAtom
#undef RemoveAtom
#undef ATOMS

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FragmentMapConformer
    //! @brief Used to map the 3D coordinates of a derived structure to its parent
    //!
    //! @see @link example_chemistry_fragment_map_conformer.cpp @endlink
    //! @author brownbp1, mendenjl
    //! @date Dec 13, 2019
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentMapConformer :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! type of drug likeness filter to apply during molecule cleaning
      std::string m_DrugLikenessType;

      //! MDL property label specifying path to protein binding pocket
      std::string m_MDL;

      //! filename for protein or protein binding pocket
      std::string m_BindingPocketFilename;

      //! scoring property to be used during the Clean phase
      descriptor::CheminfoProperty m_PropertyScorer;

      //! if true, sample conformations to resolve clashes with protein following mutate
      bool m_ResolveClashes;

      //! per-residue flexibility (lower numbers less rigid, higher numbers more rigid)
      storage::Vector< float> m_BFactors;

      //! alignment object
      FragmentAlignToScaffold m_Aligner;

      //! 3D VDW cutoff
      float m_VDWClashCutoff;

      //! enables corina conformer generation
      bool m_Corina;

      //! bypass the mapconformer step and directly set the moveable indices
      mutable storage::Vector< size_t> m_MoveableIndices;

      //! perform a quick substructure-based ensemble align and choose best conformer based on ChargeRMSD
      bool m_ChooseBestAlignedConf;

      //! Open conformational sampling to any bad geometry atoms/bonds
      bool m_FixGeometry;

      //! Open conformational sampling to neighbors (out to bond distance of size_t) of new atoms
      size_t m_AdjacentNeighbors;

      //! Open conformational sampling to new ring atoms
      bool m_MapSubgraphRingAtoms;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FragmentMapConformer();

      //! @brief druglikeness constructor
      //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
      FragmentMapConformer
      (
        const std::string &DRUG_LIKENESS_TYPE,
        const bool CORINA_CONFS,
        const storage::Vector< size_t> &MOVEABLE_INDICES = storage::Vector< size_t>()
      );

      //! @brief pose resolver constructor
      //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
      //! @param MDL property label containing path to protein binding pocket PDB file
      //! @param BINDING_POCKET_FILENAME name of protein/protein binding pocket PDB file
      //! @param PROPERTY_SCORER enables pose-dependent optimization of score
      //! @param RESOLVE_CLASHES true will seek to resolve clashes with pocket by changing the ligand conformer
      //! @param BFACTORS vector of values indicating per-residue flexibility (higher values are more flexible)
      FragmentMapConformer
      (
        const std::string &DRUG_LIKENESS_TYPE,
        const std::string &MDL,
        const std::string &BINDING_POCKET_FILENAME,
        const descriptor::CheminfoProperty &PROPERTY_SCORER,
        const bool RESOLVE_CLASHES,
        const storage::Vector< float> &BFACTORS,
        const bool CORINA_CONFS = false,
        const storage::Vector< size_t> &MOVEABLE_INDICES = storage::Vector< size_t>(),
        const bool CHOOSE_BEST_ALIGNED_CONF = false,
        const bool FIX_GEOMETRY = true,
        const size_t ADJACENT_NBRS = size_t( 1),
        const bool MAP_SUBGRAPH_RINGS = true
      );

      //! @brief clone constructor
      FragmentMapConformer *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the MDL property associated with this object
      //! @return the MDL string
      const std::string &GetMDL() const;

      //! @brief return the pocket filename associated with this object
      //! @return the PDB filename
      const std::string &GetPocketFilename() const;

      //! brief return the bfactors associated with this object
      //! return the bfactors for each atom in the pocket
      const storage::Vector< float> &GetBFactors() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief virtual operator taking an fragment and generating a new fragment by growing on a valence
      //! @param FRAGMENT small molecule of interest
      //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
      //! @return AtomVector< AtomComplete> of a new set of clean atoms following the mutation
      AtomVector< AtomComplete> CleanAtoms
      (
        const AtomVector< AtomComplete> &ATOM_VEC,
        const std::string &DRUG_LIKENESS_TYPE = "None",
        const bool &SKIP_NEUT = true,
        const bool &SKIP_SATURATE_H = false
      ) const;

      //! @brief virtual operator taking an fragment and generating a new fragment by growing on a valence
      //! @param FRAGMENT small molecule of interest
      //! @param REFERENCE_MOL the scaffold molecule for the substructure-based alignment of the new 3D conformer
      //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
      //! @return util::ShPtr< FragmentComplete> of a new 3D molecule following the mutation
      util::ShPtr< FragmentComplete> Clean
      (
        const AtomVector< AtomComplete> &ATOM_VEC,
        const FragmentComplete &REFERENCE_MOL,
        const std::string &DRUG_LIKENESS_TYPE = "None",
        const bool &SKIP_NEUT = true
      ) const;

      //! @brief preserve conformational information from starting molecule in new molecule
      //! @param STARTING_MOL the starting molecule
      //! @param NEW_MOL the resulting molecule post-design
      //! @param ATOM_COMPARISON the atom comaprison type to use for subgraph isomorphism
      //! @param BOND_COMPARISON the bond comaprison type to use for subgraph isomorphism
      //! @param COMPLEMENT if true, the mapped atom indices returned are the subgraph complement,
      //! if false then the returned indices are of the common subgraph
      //! @return the NEW_MOL indices mapped to STARTING_MOL
      storage::Set< size_t> MapAtoms
      (
        const FragmentComplete &STARTING_MOL,
        const FragmentComplete &NEW_MOL,
        const ConformationGraphConverter::AtomComparisonType &ATOM_COMPARISON = ConformationGraphConverter::AtomComparisonType::e_ElementType,
        const ConfigurationalBondTypeData::Data &BOND_COMPARISON = ConfigurationalBondTypeData::Data::e_BondOrderAmideOrAromaticWithRingness,
        const bool &COMPLEMENT = true
      ) const;

      //! @brief preserve conformational information from starting molecule in new molecule
      //! @param STARTING_MOL the starting molecule
      //! @param NEW_MOL the resulting molecule post-design
      //! @param ATOM_COMPARISON the atom comaprison type to use for subgraph isomorphism
      //! @param BOND_COMPARISON the bond comaprison type to use for subgraph isomorphism
      //! @param COMPLEMENT if true, the mapped atom indices returned are the subgraph complement,
      //! if false then the returned indices are of the common subgraph
      //! @return the NEW_MOL indices mapped to STARTING_MOL
      storage::Set< size_t> MapSubgraphRingAtoms
      (
        const FragmentComplete &STARTING_MOL,
        const FragmentComplete &NEW_MOL,
        const ConformationGraphConverter::AtomComparisonType &ATOM_COMPARISON = ConformationGraphConverter::AtomComparisonType::e_ElementType,
        const ConfigurationalBondTypeData::Data &BOND_COMPARISON = ConfigurationalBondTypeData::Data::e_BondOrderAmideOrAromaticWithRingness,
        const bool &COMPLEMENT = true
      ) const;

      //! @brief preserve conformational information from starting molecule in new molecule
      //! @param SUBGRAPH the subgraph containing rings
      //! @param RING_COMPONENTS the rings separated as component sets of indices
      //! @return the indices of ring component atoms in the subgraph
      storage::Set< size_t> MapSubgraphRingAtoms
      (
        const graph::Subgraph< size_t, size_t> &SUBGRAPH,
        const storage::List< storage::Vector< size_t> > &RING_COMPONENTS
      ) const;

      //! @brief preserve conformational information from starting molecule in new molecule
      //! @param SUBGRAPH the subgraph containing rings
      //! @param N_ADJACENT_NBRS the number of adjacent neighbors to include
      //! @return the indices of adjacent atoms in the subgraph
      storage::Set< size_t> MapSubgraphAdjacentAtoms
      (
        const graph::Subgraph< size_t, size_t> &SUBGRAPH,
        const size_t N_ADJACENT_NBRS
      ) const;

      //! @brief clean and generate a 3D structure of a molecule without perturbing atom indices
      //! @param MOL the molecule of interest
      //! @param ATOM_COMPARISON the atom comaprison type to use for subgraph isomorphism
      //! @param BOND_COMPARISON the bond comaprison type to use for subgraph isomorphism
      //! @return a 3D conformer of the original molecule with preserved atom indices
      FragmentComplete Clean3DCoords
      (
        const FragmentComplete &MOL,
        const ConformationGraphConverter::AtomComparisonType &ATOM_COMPARISON = ConformationGraphConverter::AtomComparisonType::e_ElementType,
        const ConfigurationalBondTypeData::Data &BOND_COMPARISON = ConfigurationalBondTypeData::Data::e_BondOrderAmideOrAromaticWithRingness
      ) const;

      //! @brief reset the moveable indices member data to an empty vector
      void ResetMoveableIndices() const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief get a mutex
      static sched::Mutex &GetMutex();

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class FragmentMapConformer

  } // namespace chemistry
} // namespace bcl

#endif //BCL_CHEMISTRY_FRAGMENT_MAP_CONFORMER_H_
