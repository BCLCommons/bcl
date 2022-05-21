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
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_bond_isometry_handler.h"
#include "chemistry/bcl_chemistry_fragment_split_rings.h"
#include "chemistry/bcl_chemistry_priority_dihedral_angles.h"
#include "chemistry/bcl_chemistry_stereocenters_handler.h"
#include "graph/bcl_graph_subgraph_isomorphism.h"
#include "math/bcl_math_statistics.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  ///////////
  // Enums //
  ///////////

    //! @brief Data as string
    //! @param DATA the data whose name is desired
    //! @return the name as string
    const std::string &ConformationGraphConverter::GetAtomComparisonType( const AtomComparisonType &DATA)
    {
      static const std::string s_Names[ 1 + s_NumberAtomComparisonTypes] =
      {
        "Identity",
        "ElementType",
        "AtomType",
        "AtomTypeAndChirality",
        "AtomTypeAndComplexRingChirality",
        "AtomTypeAndSymmetry",
        "AtomTypeAndHasSymmetry",
        "AtomTypeAndNumberHydrogens",
        "AtomTypeAndNumberHydrogensOnRingAtoms",
        "AtomTypeAndNumberHydrogensOnRingsAndDistinguishHydrogens",
        "AtomTypeAndDistinguishHydrogens",
        "CIPPriorityHighToLow",
        "CouldHaveSubstituents",
        GetStaticClassName< AtomComparisonType>()
      };
      return s_Names[ DATA];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor, compares by atom type, chirality, and e_BondOrderAmideWithIsometryOrAromaticWithRingness
    ConformationGraphConverter::ConformationGraphConverter() :
      m_AtomRepresentation( e_AtomTypeAndChirality),
      m_BondRepresentation( ConfigurationalBondTypeData::e_BondOrderAmideWithIsometryOrAromaticWithRingness),
      m_RemoveH( false)
    {
    }

    //! @brief constructor from coloring schemes
    //! @param ATOM_COMPARISON_TYPE means by which to color vertices of a molecule
    //! @param BOND_TYPE_INFO the desired information use for the bonds
    //! @param REMOVE_H true to remove H from any graphs that are created
    ConformationGraphConverter::ConformationGraphConverter
    (
      const AtomComparisonType &ATOM_COMPARISON_TYPE,
      const ConfigurationalBondTypeData::Data &BOND_TYPE_INFO,
      const bool &REMOVE_H
    ) :
      m_AtomRepresentation( ATOM_COMPARISON_TYPE),
      m_BondRepresentation( BOND_TYPE_INFO),
      m_RemoveH( REMOVE_H)
    {
    }

    //! @brief virtual copy constructor
    ConformationGraphConverter *ConformationGraphConverter::Clone() const
    {
      return new ConformationGraphConverter( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ConformationGraphConverter::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief Create a graph of a conformation
    //! @param CONFORMATION a conformation
    //! @return The conformation converted into a graph with the given coloring schemes
    graph::ConstGraph< size_t, size_t>
      ConformationGraphConverter::operator()( const ConformationInterface &CONFORMATION) const
    {
      if( m_RemoveH && CONFORMATION.GetNumberHydrogens())
      {
        // if hydrogens are to be ignored, construct a new fragment with the info, and remove H, then create graph
        FragmentComplete temp_fragment
        (
          AtomVector< AtomComplete>( CONFORMATION.GetAtomInfo(), CONFORMATION.GetBondInfo()),
          ""
        );
        temp_fragment.RemoveH();
        return operator()( temp_fragment);
      }
      storage::Vector< size_t> vertices;
      vertices.AllocateMemory( CONFORMATION.GetNumberAtoms());

      // get data for each atom
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr( CONFORMATION.GetAtomsIterator());
        itr.NotAtEnd();
        ++itr
      )
      {
        vertices.PushBack( ConvertAtomData( *itr, m_AtomRepresentation.GetEnum()));
      }
      if( m_AtomRepresentation == e_AtomTypeAndSymmetry)
      {
        storage::Vector< size_t> symmetry_multiplicities( PriorityDihedralAngles::CalculateMultiplicity( CONFORMATION));
        size_t natypes( GetAtomTypes().GetEnumCount());
        for( size_t i( 0), sz( CONFORMATION.GetSize()); i < sz; ++i)
        {
          vertices( i) += natypes * ( symmetry_multiplicities( i) - 1);
        }
      }
      else if( m_AtomRepresentation == e_AtomTypeAndHasSymmetry)
      {
        storage::Vector< size_t> symmetry_multiplicities( PriorityDihedralAngles::CalculateMultiplicity( CONFORMATION));
        size_t natypes( GetAtomTypes().GetEnumCount());
        for( size_t i( 0), sz( CONFORMATION.GetSize()); i < sz; ++i)
        {
          vertices( i) += natypes * ( std::min( symmetry_multiplicities( i), size_t( 2)) - 1);
        }
      }
      else if( m_AtomRepresentation == e_CIPPriorityHighToLow)
      {
        vertices = PriorityDihedralAngles::GetPriority( CONFORMATION);
      }
      else if
      (
        m_AtomRepresentation == e_AtomTypeAndNumberHydrogensOnRingsAndDistinguishHydrogens
        || m_AtomRepresentation == e_AtomTypeAndDistinguishHydrogens
      )
      {
        auto atom_itr( CONFORMATION.GetAtomsIterator());
        for( size_t i( 0), sz( CONFORMATION.GetSize()); i < sz; ++i, ++atom_itr)
        {
          if( atom_itr->GetAtomType() != GetAtomTypes().H_S && atom_itr->GetElementType() != GetElementTypes().e_Fluorine)
          {
            continue;
          }
          const storage::Vector< BondConformational> &bonded_atom_bonds
          (
            atom_itr->GetBonds()( 0).GetTargetAtom().GetBonds()
          );
          size_t n_h( 0), n_f( 0);
          for( auto itr_bnd( bonded_atom_bonds.Begin()); 1; ++itr_bnd)
          {
            if( itr_bnd->GetTargetAtom().GetAtomType() == GetAtomTypes().H_S)
            {
              ++n_h;
            }
            if( itr_bnd->GetTargetAtom().GetElementType() == GetElementTypes().e_Fluorine)
            {
              ++n_f;
            }
            if( &itr_bnd->GetTargetAtom() == &*atom_itr)
            {
              vertices( i) += atom_itr->GetAtomType() == GetAtomTypes().H_S ? n_h : n_f;
              break;
            }
          }
        }
      }

      // construct the graph and return it
      return graph::ConstGraph< size_t, size_t>
             (
               vertices,
               CONFORMATION.GetAdjacencyList( m_BondRepresentation),
               GetConfigurationalBondTypes().e_Undefined
             );
    }

    //! @brief Create a tree of a conformation, vertex data indicates the number of atoms represented by the point
    //!        1 for atoms in a chain, > 2 for atoms in rings
    //! @param CONFORMATION a conformation
    //! @return The conformation converted into a graph with the given atom/bond representations,
    //!         along with a vector, each index corresponds to atom in the molecule, value in vector is the index in the
    //!         graph corresponding to that point
    storage::Pair< graph::ConstGraph< size_t, size_t>, storage::Vector< size_t> >
      ConformationGraphConverter::ToTree( const ConformationInterface &CONFORMATION) const
    {
      auto atom_graph( CreateGraphWithAtoms( CONFORMATION));
      storage::List< storage::Vector< size_t> > rings
      (
        FragmentSplitRings( true).GetComponentVertices( CONFORMATION, atom_graph)
      );
      const size_t n_rings( rings.GetSize());
      storage::Vector< size_t> new_id( CONFORMATION.GetSize(), util::GetUndefinedSize_t());
      size_t new_atom_id( 0);
      for( auto itr( CONFORMATION.GetAtomsIterator()); itr.NotAtEnd(); ++itr)
      {
        new_id( itr.GetPosition()) =
          itr->CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, size_t( 1)) < size_t( 2)
          ? new_atom_id++
          : util::GetUndefinedSize_t();
      }
      const size_t nr_chain_atoms( new_atom_id);
      const size_t nr_ring_atoms( CONFORMATION.GetSize() - nr_chain_atoms);
      const size_t graph_sz( n_rings + nr_chain_atoms);
      storage::Vector< size_t> vertex_size( graph_sz, size_t( 1));
      for( auto itr_ring( rings.Begin()), itr_ring_end( rings.End()); itr_ring != itr_ring_end; ++itr_ring)
      {
        const size_t ring_sz( itr_ring->GetSize());
        for( auto itr_itr_ring( itr_ring->Begin()), itr_itr_ring_end( itr_ring->End()); itr_itr_ring != itr_itr_ring_end; ++itr_itr_ring)
        {
          new_id( *itr_itr_ring) = new_atom_id;
        }
        vertex_size( new_atom_id) = ring_sz;
        ++new_atom_id;
      }
      const size_t tree_size( new_atom_id);

      storage::Vector< graph::UndirectedEdge< size_t> > initial_edges
      (
        CONFORMATION.GetAdjacencyList( m_BondRepresentation)
      );
      if( !util::IsDefined( new_id( 0)))
      {
        CONFORMATION.WriteMDL( util::GetLogger());
      }

      auto itr_edges_store( initial_edges.Begin());
      for
      (
        auto itr_edges( initial_edges.Begin()), itr_edges_end( initial_edges.End());
        itr_edges != itr_edges_end;
        ++itr_edges
      )
      {
        if( new_id( itr_edges->GetIndexLow()) != new_id( itr_edges->GetIndexHigh()))
        {
          *itr_edges_store =
            graph::UndirectedEdge< size_t>
            (
              new_id( itr_edges->GetIndexLow()),
              new_id( itr_edges->GetIndexHigh()),
              itr_edges->GetEdgeData()
            );
          ++itr_edges_store;
        }
      }
      initial_edges.Resize( std::distance( initial_edges.Begin(), itr_edges_store));
      return
        storage::Pair< graph::ConstGraph< size_t, size_t>, storage::Vector< size_t> >
        (
          graph::ConstGraph< size_t, size_t>
          (
            vertex_size,
            initial_edges,
            GetConfigurationalBondTypes().e_Undefined
          ),
          new_id
        );
    }

    //! @brief Create a tree of a conformation, vertex data indicates the number of atoms represented by the point
    //!        1 for atoms in a chain, > 2 for atoms in rings
    //! @param CONFORMATION a conformation
    //! @return The conformation converted into a graph with the given atom/bond representations,
    //!         along with a vector, each index corresponds to atom in the molecule, value in vector is the index in the
    //!         graph corresponding to that point
    size_t ConformationGraphConverter::CountNonRingVariantIsomorphisms( const ConformationInterface &CONFORMATION) const
    {
      FragmentComplete confcopy( CONFORMATION);
      confcopy.Canonicalize();
      graph::ConstGraph< size_t, size_t> base_graph
      (
        ConformationGraphConverter( e_AtomType, ConfigurationalBondTypeData::e_IsInRing)( CONFORMATION)
      );
      graph::ConstGraph< size_t, size_t> real_graph
      (
        ConformationGraphConverter( e_AtomType, ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness)( CONFORMATION)
      );
      graph::ConstGraph< size_t, size_t> base_grapho( base_graph);
      graph::ConstGraph< size_t, size_t> real_grapho( real_graph);

      for( size_t vertexn( 0), sz( confcopy.GetSize()); vertexn < sz; ++vertexn)
      {
        size_t nbonds( base_graph.GetNeighborData( vertexn).GetSize());
        if( nbonds < size_t( 3))
        {
          continue;
        }
        // count number of ring bonds
        size_t n_ring_bonds
        (
          math::Statistics::Sum
          (
            base_graph.GetNeighborData( vertexn).Begin(),
            base_graph.GetNeighborData( vertexn).End()
          )
        );
        if( n_ring_bonds < size_t( 2))
        {
          continue;
        }
        // break the highest priority ring bonds
        size_t n_to_remove( n_ring_bonds - size_t( 1));
        size_t ringpos( 0);
        for( size_t i( 0); i < n_to_remove && ringpos < base_graph.GetNeighborData( vertexn).GetSize(); ++i)
        {
          if( base_graph.GetNeighborData( vertexn)( ringpos) != size_t( 1))
          {
            ++ringpos;
            --i;
            continue;
          }
          size_t indx_rmv( base_graph.GetNeighborIndices( vertexn)( ringpos));
          size_t n_ring_bondsb
          (
            math::Statistics::Sum
            (
              base_graph.GetNeighborData( indx_rmv).Begin(),
              base_graph.GetNeighborData( indx_rmv).End()
            )
          );
          if( n_ring_bondsb < size_t( 2))
          {
            ++ringpos;
            --i;
            continue;
          }
          base_graph.RemoveEdge( vertexn, indx_rmv);
          if( !graph::Connectivity::IsConnected( base_graph))
          {
            base_graph.AddEdge( vertexn, indx_rmv, size_t( 1));
            ++ringpos;
            --i;
            continue;
          }
          real_graph.RemoveEdge( vertexn, indx_rmv);
        }
      }
      graph::SubgraphIsomorphism< size_t, size_t> isos;
      isos.SetGraphExternalOwnership( real_graph);
      isos.SetSubgraphExternalOwnership( real_graph);
      isos.FindAllIsomorphisms();
      if( isos.GetIsomorphisms().GetSize() == size_t( 0))
      {
        BCL_Debug( real_graph.GetBasicConnectivity());
        BCL_Debug( base_graph.GetBasicConnectivity());
        BCL_Debug( real_grapho.GetBasicConnectivity());
        BCL_Debug( base_grapho.GetBasicConnectivity());
      }
      return isos.GetIsomorphisms().GetSize();
    }

    //! @brief Create an atom vector from a graph with atoms
    //! @param GRAPH a graph created by this converter with CreateGraphWithAtoms
    //! @return the generated atom vector
    AtomVector< AtomComplete> ConformationGraphConverter::CreateAtomsFromGraph
    (
      const t_AtomGraph &GRAPH,
      const bool &RECALCULATE_CONFIGURATION
    )
    {
      // recreate atom info and bond info vectors
      storage::Vector< sdf::AtomInfo> atom_info( GRAPH.GetSize());
      storage::Vector< sdf::BondInfo> bond_info;
      bond_info.AllocateMemory( GRAPH.NumEdges() / 2);

      // iterate over atoms
      for( size_t atom_id( 0), number_atoms( atom_info.GetSize()); atom_id < number_atoms; ++atom_id)
      {
        // set the atom info up properly
        atom_info( atom_id) = GRAPH.GetVertexData( atom_id)->GetAtomInfo();

        // get the bonds of this atom and the bond data
        const storage::Vector< size_t> &neighbor_indices( GRAPH.GetNeighborIndices( atom_id));
        const storage::Vector< size_t> &neighbor_data( GRAPH.GetNeighborData( atom_id));

        for
        (
          storage::Vector< size_t>::const_iterator
            itr_index( neighbor_indices.Begin()), itr_data( neighbor_data.Begin()), itr_data_end( neighbor_data.End());
          itr_data != itr_data_end;
          ++itr_data, ++itr_index
        )
        {
          // only need one bond info per bond, so skip if the first atom index is >= than *itr_atom_index
          if( atom_id >= *itr_index)
          {
            continue;
          }

          bond_info.PushBack
          (
            sdf::BondInfo( atom_id, *itr_index, ConfigurationalBondType( *itr_data))
          );
        }
      }
      AtomVector< AtomComplete> atoms( atom_info, bond_info);

      if( RECALCULATE_CONFIGURATION)
      {
        // add isometry info, which may be different for a subgraph
        BondIsometryHandler::AddIsometryInformation( atoms, true);

        // add stereocenter information
        StereocentersHandler::AddChiralityFromConformation( atoms);
      }

      return atoms;
    }

    //! @brief Given a small molecule, instantiate its graphical representation
    //! @param CONFORMATION a conformation
    //! @return The resulting graph after the conversion.
    ConformationGraphConverter::t_AtomGraph
      ConformationGraphConverter::CreateGraphWithAtoms( const ConformationInterface &CONFORMATION)
    {
      // construct the graph and return it
      iterate::Generic< const AtomConformationalInterface> itr( CONFORMATION.GetAtomsIterator());
      return
       t_AtomGraph
       (
         storage::Vector< util::SiPtr< const AtomConformationalInterface> > // si-ptrs to the atoms involved
         (
           itr,
           itr.End()
         ),
         CONFORMATION.GetAdjacencyList( ConfigurationalBondTypeData::e_ConfigurationalBondType),
         GetConfigurationalBondTypes().e_Undefined                   // value for undefined bonds
       );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ConformationGraphConverter::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_AtomRepresentation, ISTREAM);
      io::Serialize::Read( m_BondRepresentation, ISTREAM);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &ConformationGraphConverter::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_AtomRepresentation, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_BondRepresentation, OSTREAM);
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief helper function to get AtomComparisonType from an atom
    //! @param ATOM atom to retrieve data from
    //! @param TYPE atom comparison type
    //! @return data for the atome atom
    size_t ConformationGraphConverter::ConvertAtomData( const AtomConformationalInterface &ATOM, const AtomComparisonType &TYPE)
    {
      switch( TYPE)
      {
        case e_Identity:
          return size_t( 1);
        case e_ElementType:
          return ATOM.GetElementType().GetIndex();
        case e_AtomType:
        case e_AtomTypeAndSymmetry:
          return ATOM.GetAtomType().GetIndex();
        case e_AtomTypeAndChirality:
          // pack the atom type index together with R/S, if it is known
          return s_NumberChiralities * ATOM.GetAtomType().GetIndex() + size_t( ATOM.GetChirality());
        case e_AtomTypeAndComplexRingChirality:
          // pack the atom type index together with R/S, if it is known
          return s_NumberChiralities * ATOM.GetAtomType().GetIndex()
                 + size_t
                   (
                     ATOM.CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, size_t( 1)) >= size_t( 4)
                     ? ATOM.GetChirality()
                     : 0
                   );
        case e_AtomTypeAndNumberHydrogens:
          return 9 * ATOM.GetAtomType().GetIndex()
                 + ATOM.GetNumberCovalentlyBoundHydrogens();
        case e_AtomTypeAndNumberHydrogensOnRings:
          return 9 * ATOM.GetAtomType().GetIndex()
                 +
                 (
                   ATOM.CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, size_t( 1)) >= size_t( 2)
                   && ATOM.CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsAromatic, size_t( 1)) < size_t( 2)
                   ? ATOM.GetNumberCovalentlyBoundHydrogens()
                   : 0
                 );
        case e_AtomTypeAndNumberHydrogensOnRingsAndDistinguishHydrogens:
          return 90 * ATOM.GetAtomType().GetIndex()
                 +
                 9 *
                 (
                   ATOM.CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, size_t( 1)) >= size_t( 2)
                   && ATOM.CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsAromatic, size_t( 1)) < size_t( 2)
                   ? ATOM.GetNumberCovalentlyBoundHydrogens()
                   : 0
                 );
        case e_AtomTypeAndDistinguishHydrogens:
          return 9 * ATOM.GetAtomType().GetIndex();
        case e_CIPPriorityHighToLow:
          return 0;
        case e_CouldHaveSubstituents:
          return ( ATOM.GetElementType()->GetMainGroup() != 6 && ATOM.GetElementType()->GetMainGroup() != 2);
        case s_NumberAtomComparisonTypes:
        default:
          return util::GetUndefined< size_t>();
      }
      return util::GetUndefined< size_t>();
    }

    //! @brief get the atom's representation as a size_t using m_AtomRepresentation
    //! @param ATOM atom type to retrieve data from
    //! @param TYPE the type of atom comparison that will be used
    //! @return the atom's representation as a size_t using m_AtomRepresentation
    size_t ConformationGraphConverter::ConvertAtomTypeData( const AtomType &ATOM, const AtomComparisonType &TYPE)
    {
      switch( TYPE)
      {
        case e_Identity:
          return size_t( 1);
        case e_ElementType:
          return ATOM->GetElementType().GetIndex();
        case e_AtomType:
          return ATOM.GetIndex();
        default:
          return util::GetUndefined< size_t>();
      }
      return util::GetUndefined< size_t>();
    }
  } // namespace chemistry
} // namespace bcl

