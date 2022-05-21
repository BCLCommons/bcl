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
#include "chemistry/bcl_chemistry_fragment_split_gadd_fragments.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "graph/bcl_graph_connectivity.h"
#include "graph/bcl_graph_subgraph.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> FragmentSplitGADDFragments::s_Instance
    (
      util::Enumerated< FragmentSplitInterface>::AddInstance( new FragmentSplitGADDFragments)
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor
    //! @param GET_RINGS true if rings are desired, false if chains are desired
    //! @param MIN_SIZE get the minimum size of ring or chain that is desired
    FragmentSplitGADDFragments::FragmentSplitGADDFragments()
    {
    }

    //! virtual copy constructor
    FragmentSplitGADDFragments *FragmentSplitGADDFragments::Clone() const
    {
      return new FragmentSplitGADDFragments( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! the class name as const ref std::string
    const std::string &FragmentSplitGADDFragments::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &FragmentSplitGADDFragments::GetAlias() const
    {
      static const std::string s_name( "GADDFragments");
      return s_name;
    }

    //! @brief Get a description for what this class does (used when writing help)
    //! @return a description for what this class does (used when writing help)
    const std::string &FragmentSplitGADDFragments::GetClassDescription() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the minimum size of fragments
    //! @return the minimum size of fragments
    const size_t FragmentSplitGADDFragments::GetMinSize() const
    {
      return 0;
    }

  /////////////////
  //  operations //
  /////////////////

    //! @brief returns list of ring or chain fragments
    //! @param MOLECULE molecule of interest
    //! @param MOLECULE_GRAPH graph of molecule of interest
    //! @return list of ring or chain fragments
    storage::List< storage::Vector< size_t> > FragmentSplitGADDFragments::GetComponentVertices
    (
      const ConformationInterface &MOLECULE,
      ConformationGraphConverter::t_AtomGraph &MOLECULE_GRAPH
    ) const
    {

      // Determine which atoms are considered "heteroatoms"
      size_t atom_no( 0);
      iterate::Generic< const AtomConformationalInterface> itr_atom( MOLECULE.GetAtomsIterator());
      storage::Vector< size_t> heteroatoms;

      // If the atom is a heteroatom, or if it is a carbon multiple-bonded to a heteroatom, count it as a heteroatom
      const ElementType &carbon( GetElementTypes().e_Carbon);
      for( ; itr_atom.NotAtEnd(); ++itr_atom, ++atom_no)
      {
        if( itr_atom->GetElementType() != carbon)
        {
          heteroatoms.PushBack( atom_no);
        }
        else
        {

          bool is_heteroatom( false);
          const storage::Vector< BondConformational> &bonds( itr_atom->GetBonds());

          // If the carbon is multiple-bonded to any heteroatoms it should be considered a "heteroatom" itself
          for( size_t bond_no( 0); bond_no < bonds.GetSize(); ++bond_no)
          {
            if
            (
              bonds( bond_no).GetTargetAtom().GetElementType() != carbon
              && (
                   bonds( bond_no).GetBondType()->GetNumberOfElectrons() > 2
                   || bonds( bond_no).GetBondType()->GetConjugation() == ConstitutionalBondTypeData::e_Aromatic
                 )
            )
            {
              is_heteroatom = true;
              break;
            }
          }

          if( is_heteroatom)
          {
            heteroatoms.PushBack( atom_no);
          }
        }
      }

      util::OwnPtr< ConformationGraphConverter::t_AtomGraph> mol_graph_ptr( &MOLECULE_GRAPH, false);

      graph::Subgraph< util::SiPtr< const AtomConformationalInterface>, size_t> heteroatom_subgraph
      (
        mol_graph_ptr,
        heteroatoms
      );

      // Remove edges from the graph to isolate heteroatom and carbon fragments
      storage::List< storage::Pair< size_t, size_t> > adj_edges( heteroatom_subgraph.GetAdjacentEdgeIndices());
      for
      (
        storage::List< storage::Pair< size_t, size_t> >::const_iterator itr_edge( adj_edges.Begin()), itr_edge_end( adj_edges.End());
        itr_edge != itr_edge_end;
        ++itr_edge
      )
      {
        MOLECULE_GRAPH.RemoveEdge( itr_edge->First(), itr_edge->Second());
      }

      // Get the connected components, i.e. the individual carbon and heteroatom fragments
      storage::List< storage::Vector< size_t> > indices( graph::Connectivity::GetComponents( MOLECULE_GRAPH));

      return indices;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer FragmentSplitGADDFragments::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "splits molecules into their GA-based Drug Database fragments "
        "(see http://www.daylight.com/meetings/mug01/Yang/gadd/)"
      );
      return parameters;
    }

  } // namespace chemistry
} // namespace bcl
