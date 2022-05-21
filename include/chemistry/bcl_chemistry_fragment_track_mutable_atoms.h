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

#ifndef BCL_CHEMISTRY_FRAGMENT_TRACK_MUTABLE_ATOMS_H_
#include "graph/bcl_graph_common_subgraph_isomorphism_base.h"

#define BCL_CHEMISTRY_FRAGMENT_TRACK_MUTABLE_ATOMS_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically
#include "find/bcl_find.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_conformational_interface.h"
#include "bcl_chemistry_collector_valence.h"
#include "bcl_chemistry_fragment_complete.h"
#include "bcl_chemistry_fragment_constitution_shared.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "find/bcl_find_pick_interface.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_result.h"
#include "util/bcl_util_function_interface.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_si_ptr_list.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FragmentTrackMutableAtoms
    //! @brief Used to map the 3D coordinates of a derived structure to its parent
    //!
    //! @see @link example_chemistry_fragment_track_mutable_atoms.cpp @endlink
    //! @author brownbp1
    //! @date Dec 13, 2019
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentTrackMutableAtoms :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! mutable fragment specified through indices or fragment
      static FragmentComplete s_MutableFragment;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FragmentTrackMutableAtoms();

      //! @brief clone constructor
      FragmentTrackMutableAtoms *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      // TODO: legacy; will be deprecated at some point in favor of GetAtomFromMutable and related functions
      //! @brief a function that returns the non-mutable base fragment provided a mutable fragment/atoms
      //! @param FRAGMENT the small molecule of interest
      //! @param MUTABLE_FRAGMENT sub-fragment of small molecule that can be mutated
      //! @param MUTABLE_ATOMS the atoms within the sub-fragment that can be mutated
      //! @param COMPLEMENT get the complement indices instead
      //! @return a fragment that can be modified by mutate interface classes
      static FragmentComplete GetBaseFragment
      (
        const FragmentComplete &FRAGMENT,
        const FragmentComplete &MUTATABLE_FRAGMENT,
        const storage::Vector< size_t> &MUTABLE_ATOMS,
        const bool COMPLEMENT = true, // legacy; needed for backwards compatibility
        const ConformationGraphConverter &GRAPH_MAKER = ConformationGraphConverter
        (
          ConformationGraphConverter::e_ElementType,
          ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness
        ),
        const graph::CommonSubgraphIsomorphismBase::SolutionType &GRAPH_SOLUTION_TYPE =
            graph::CommonSubgraphIsomorphismBase::SolutionType::e_GreedyUnconnected
      );

      //! @brief a function that randomly picks a fragment to be the mutable region and returns the non-mutable base
      //! @param FRAGMENT the small molecule of interest
      //! @param SPLITTER the splitter to use to make the fragmnts
      //! @return a fragment that can be modified by mutate interface classes
      static FragmentComplete GetRandomBaseFragment
      (
        const FragmentComplete &FRAGMENT,
        const util::Implementation< FragmentSplitInterface> &SPLITTER
      );

      //! @brief a function that randomly picks a fragment from an allowed fragment/atoms
      //! to be the mutable region and returns the non-mutable base
      //! @param FRAGMENT the small molecule of interest
      //! @param SPLITTER the splitter to use to make the fragmnts
      //! @param MUTABLE_FRAGMENT sub-fragment of small molecule that can be mutated
      //! @param MUTABLE_ATOMS the atoms within the sub-fragment that can be mutated
      //! @return a fragment that can be modified by mutate interface classes
      static FragmentComplete GetRandomRestrictedBaseFragment
      (
        const FragmentComplete &FRAGMENT,
        const FragmentSplitInterface &SPLITTER,
        const FragmentComplete &MUTATABLE_FRAGMENT,
        const storage::Vector< size_t> &MUTABLE_ATOMS
      );

      //! @brief a function that randomly picks an atom from a list of allowed element types for transformation
      //! @param FRAGMENT the small molecule of interest
      //! @param MUTABLE_ATOMS the atoms indicating the fragment that can be mutated
      //! @return a pointer to the selected atom
      static storage::Vector< size_t> SetMutableElements
      (
        const FragmentComplete &FRAGMENT,
        const storage::Vector< ElementType> &MUTABLE_ELEMENTS
      );

      //! @brief a function that randomly picks an atom from a molecular substructure for transformation
      //! @param FRAGMENT the small molecule of interest
      //! @param REF_FRAGMENT the starting molecule
      //! @return a pointer to the selected atom
      static storage::Vector< size_t> SetMutableFragments
      (
        const FragmentComplete &FRAGMENT,
        const FragmentEnsemble &MUTABLE_FRAGMENTS,
        const bool COMPLEMENT = false,
        const ConformationGraphConverter &GRAPH_MAKER = ConformationGraphConverter
        (
          ConformationGraphConverter::e_ElementType,
          ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness
        ),
        const graph::CommonSubgraphIsomorphismBase::SolutionType &GRAPH_SOLUTION_TYPE =
            graph::CommonSubgraphIsomorphismBase::SolutionType::e_GreedyUnconnected
      );

      //! @brief a function that randomly picks an atom for transformation
      //! @param FRAGMENT the small molecule of interest
      //! @param MUTABLE_ATOMS the atoms indicating the fragment that can be mutated
      //! @return a pointer to the selected atom
      static util::SiPtr< const AtomConformationalInterface> GetAtomFromMutable
      (
        const FragmentComplete &FRAGMENT,
        const bool ALL_ATOMS = true,
        const storage::Vector< size_t> &MUTABLE_ATOMS = storage::Vector< size_t>(),
        const storage::Vector< ElementType> &MUTABLE_ELEMENTS = storage::Vector< ElementType>(),
        const FragmentEnsemble &MUTABLE_FRAGMENTS = FragmentEnsemble(),
        const storage::Vector< size_t> &FIXED_ATOMS = storage::Vector< size_t>(),
        const storage::Vector< ElementType> &FIXED_ELEMENTS = storage::Vector< ElementType>(),
        const FragmentEnsemble &FIXED_FRAGMENTS = FragmentEnsemble(),
        const bool COMPLEMENT_MUTABLE_FRAGMENTS = false,
        const bool COMPLEMENT_FIXED_FRAGMENTS = false,
        const ConformationGraphConverter &MUTABLE_GRAPH_MAKER = ConformationGraphConverter
        (
          ConformationGraphConverter::e_ElementType,
          ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness
        ),
        const ConformationGraphConverter &FIXED_GRAPH_MAKER = ConformationGraphConverter
        (
          ConformationGraphConverter::e_ElementType,
          ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness
        )
      );

    //////////////////////
    // helper functions //
    //////////////////////

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

    }; // class FragmentTrackMutableAtoms

  } // namespace chemistry
} // namespace bcl

#endif //BCL_CHEMISTRY_FRAGMENT_TRACK_MUTABLE_ATOMS_H_
