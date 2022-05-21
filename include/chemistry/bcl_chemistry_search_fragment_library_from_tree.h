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

#ifndef BCL_CHEMISTRY_SEARCH_FRAGMENT_LIBRARY_FROM_TREE_H_
#define BCL_CHEMISTRY_SEARCH_FRAGMENT_LIBRARY_FROM_TREE_H_

// include the namespace header
#include "bcl_chemistry.h"
#include "bcl_chemistry_bond_angle_assignment.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_fragment_graph_marker.h"
#include "bcl_chemistry_rotamer_library_interface.h"
#include "bcl_chemistry_small_molecule_fragment_isomorphism.h"
#include "graph/bcl_graph_const_graph.h"
#include "sched/bcl_sched_mutex.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SearchFragmentLibraryFromTree
    //! @brief Searches Scaffolds( fragments) that have isomorphism with the given molecule.
    //! @details creates SmallMoleculeFragmentIsomorphism object for each of the scaffolds that are found for molecule
    //!
    //! @see @link example_chemistry_search_fragment_library_from_tree.cpp @endlink
    //! @author kothiwsk
    //! @date Jul 09, 2014
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SearchFragmentLibraryFromTree :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! share pointer to ensemble of fragments
      mutable util::Implementation< RotamerLibraryInterface> m_FragmentLibrary;

      //! graph that denotes the substructure tree of constitutions
      mutable graph::ConstGraph< size_t, size_t> m_SubstructureTree;

      //! Vector of muteces; one for each graph. Only the mutex for the graph currently desired is locked
      mutable storage::Vector< sched::Mutex> m_Muteces;

      //! string to keep track of which constitutions have already been seen so graphs are not re-created for them.
      mutable std::string m_ConstitutionsSeen;

      //! graphs of unique constitutions in the rotamer library
      mutable storage::Vector< graph::ConstGraph< size_t, size_t> > m_ConstitutionGraphs;

      //! mapping between constitutions and configurations
      mutable storage::Vector< storage::Set< size_t> > m_ConfigurationMapping;

      //! number of root nodes
      mutable storage::List< size_t> m_RootNodes;

      //! bond angle map
      mutable typename RotamerLibraryInterface::t_BondAngleMap m_BondAngleMap;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SearchFragmentLibraryFromTree();

      //! @brief constructor
      //! @param FRAGMENTLIBRARY pointer to ensemble of fragments
      explicit SearchFragmentLibraryFromTree( const RotamerLibraryInterface &FRAGMENT_LIBRARY);

      //! @brief Clone function
      //! @return pointer to new SearchFragmentLibraryFromTree
      SearchFragmentLibraryFromTree *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the rotamer library
      //! @return reference to the rotamer library
      const util::Implementation< RotamerLibraryInterface> &GetRotamerLibrary() const
      {
        return m_FragmentLibrary;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief find fragments that have isomorphism with molecule
      //! @param MOLECULE molecule of interest whose fragments need to be searched from the given library of fragments
      //! @return vector of SmallMoleculeFragmentIsomorphism
      util::ShPtrVector< SmallMoleculeFragmentIsomorphism> FindFragmentsOfMolecule
      (
        const ConformationInterface &MOLECULE
      ) const;

      //! @brief find bond angle information objects for a given molecule
      //! @param MOLECULE molecule of interest whose fragments need to be searched from the given library of fragments
      //! @return vector of SmallMoleculeFragmentIsomorphism
      util::ShPtrVector< BondAngleAssignment> GetBondAngleAssignments
      (
        const ConformationInterface &MOLECULE,
        const bool &ENFORCE_CHIRALITY
      ) const;

      //! @brief checks whether fragment graph is part of molecule graph
      //! @param MOLECULE_GRAPH molecule graph for which fragments need to be checked
      //! @param FRAGMENT_GRAPH fragment graph which needs to be checked if contained in molecule
      static bool CheckSubStructure
      (
        const graph::ConstGraph< size_t, size_t> &MOLECULE_GRAPH,
        const graph::ConstGraph< size_t, size_t> &FRAGMENT_GRAPH
      );

      //! @brief creates smallmolecule fragment isomorphism object
      //! @param MOLECULE_GRAPH molecule graph for which fragments need to be checked
      //! @param FRAGMENT_GRAPH fragment graph which needs to be checked if contained in molecule
      //! @param FRAGMENT fragment which is part of molecule
      //! @return return a small molecule fragment isomorphism object
      static util::ShPtr< SmallMoleculeFragmentIsomorphism> CreateFragmentIsomorphismObject
      (
        const FragmentComplete &MOLECULE,
        const graph::ConstGraph< size_t, size_t> &MOLECULE_GRAPH,
        const graph::ConstGraph< size_t, size_t> &FRAGMENT_GRAPH,
        const FragmentComplete &FRAGMENT
      );

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
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

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief load the rotamer library information. It is necessary to call this function before the first call to
      //!        any of the operations. It should not, however, be called in the constructor, because this class may be
      //!        constructed with a RotamerLibraryInterface during static initialization, and where that happens, the
      //!        rotamer library information does not exist
      void LoadRotamerLibraryInformation() const;

      //! @brief acquire mutex for acquiring an isomorphism object
      //! @return mutex used in AcquireIsomorphism and ReleaseIsomorphism
      static sched::Mutex &GetIsomorphismMutex();

      //! @brief acquire an isomorphism object
      //! @return iterator to the isomorphism object
      static storage::List< graph::SubgraphIsomorphism< size_t, size_t> >::iterator AcquireIsomorphism();

      //! @brief release a given isomorphism object
      //! @param ITR iterator to the isomorphism object
      static void ReleaseIsomorphism
      (
        const storage::List< graph::SubgraphIsomorphism< size_t, size_t> >::iterator &ITR
      );

      //! Pool accessors of isomorphism-computing objects, stored because they are very large and excessive memory reallocations
      //! can be avoided by maintaining such a pool, without making it a data member directly, which would make this
      //! class's functions not thread-safe

      //! @brief get the list of isomorphisms that have been allocated
      static storage::List< graph::SubgraphIsomorphism< size_t, size_t> > &GetAllocatedIsomorphisms();

      //! @brief get the list of isomorphisms that are available
      static storage::List< graph::SubgraphIsomorphism< size_t, size_t> > &GetAvailableIsomorphisms();

      //! @brief create molecule from empty molecule containing atom vector and bond vector information in properties
      //! @param molecule which contains atomvector and bondvector information in properties
      //! @return molecule with all the atom type, bond type and connectivity information
      static FragmentComplete CreateMoleculeFromRotlib
      (
        FragmentComplete MOLECULE
      );

      //! @brief search bond angle library for a given atom
      typename RotamerLibraryInterface::t_BondAngleMap::const_iterator SearchForBondAngles
      (
        const AtomConformationalInterface &ATOM,
        const ConformationGraphConverter::AtomComparisonType &ATOM_COMPARISON
      ) const;

    }; // class SearchFragmentLibraryFromTree

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_SEARCH_FRAGMENT_LIBRARY_FROM_TREE_H_
