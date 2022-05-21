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
#include "chemistry/bcl_chemistry_fragment_split_rigid.h"

// includes from bcl - sorted alphabetically
#include "graph/bcl_graph_connectivity.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> FragmentSplitRigid::s_Instance
    (
      util::Enumerated< FragmentSplitInterface>::AddInstance( new FragmentSplitRigid())
    );

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> FragmentSplitRigid::s_AmideInstance
    (
      util::Enumerated< FragmentSplitInterface>::AddInstance( new FragmentSplitRigid( 1, false))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor
    //! @param MIN_SIZE the minimum size of rigid fragment that is desired
    FragmentSplitRigid::FragmentSplitRigid( const size_t MIN_SIZE, const bool &CONSIDER_AMIDE_RIGID) :
      m_MinSize( MIN_SIZE),
      m_ConsiderAmideRigid( CONSIDER_AMIDE_RIGID)
    {
    }

    //! virtual copy constructor
    FragmentSplitRigid *FragmentSplitRigid::Clone() const
    {
      return new FragmentSplitRigid( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! the class name as const ref std::string
    const std::string &FragmentSplitRigid::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &FragmentSplitRigid::GetAlias() const
    {
      static const std::string s_name_ring( "Rigid"), s_name( "RigidSansAmide");
      return m_ConsiderAmideRigid ? s_name_ring : s_name;
    }

    //! get the minimum size of a component of interest
    const size_t FragmentSplitRigid::GetMinSize() const
    {
      return m_MinSize;
    }

  /////////////////
  //  operations //
  /////////////////

    namespace
    {
      enum BondRigidity
      {
        e_Rigid,           // Ring, 2x, 3x, amide, or links to atom with triple bond
        e_AttatchedToRing, // Adjacent to ring, but not itself rigid
        e_Terminal,        // Bond is to a terminal atom or atom that is otherwise only connected to H
        e_Flexible,        // All others
        s_NumberBondRigidity
      };

      //! @brief helper function that always returns "Nameless" no matter what was passed to it
      const std::string &Nameless( const BondRigidity &ENUM)
      {
        static const std::string s_nameless( "Nameless");
        return s_nameless;
      }

      typedef util::WrapperEnum< BondRigidity, &Nameless, s_NumberBondRigidity> BondRigidityEnum;

    }

    //! @brief returns list of ring or chain fragments
    //! @param MOLECULE molecule of interest
    //! @param MOLECULE_GRAPH graph of molecule of interest
    //! @return list of ring or chain fragments
    storage::List< storage::Vector< size_t> > FragmentSplitRigid::GetComponentVertices
    (
      const ConformationInterface &MOLECULE,
      ConformationGraphConverter::t_AtomGraph &MOLECULE_GRAPH
    ) const
    {
      // determine whether each atom is terminal
      const size_t n_atoms( MOLECULE_GRAPH.GetSize());
      linal::Vector< size_t> is_terminal( n_atoms, size_t( 0));
      linal::Vector< size_t> is_in_ring( n_atoms, size_t( 0));
      linal::Vector< size_t> has_triple_bond( n_atoms, size_t( 0));
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atom( MOLECULE.GetIterator());
        itr_atom.NotAtEnd();
        ++itr_atom
      )
      {
        size_t n_nonterminal_bonds( 0), n_uniq_term_bond_types( 0);
        if( itr_atom->GetBonds().GetSize() > size_t( 1))
        {
          ElementType last_terminal_bond_type( GetElementTypes().e_Undefined);
          for
          (
            storage::Vector< BondConformational>::const_iterator
              itr_bonds( itr_atom->GetBonds().Begin()), itr_bonds_end( itr_atom->GetBonds().End());
            itr_bonds != itr_bonds_end;
            ++itr_bonds
          )
          {
            const AtomConformationalInterface &tatom( itr_bonds->GetTargetAtom());
            if( tatom.GetAtomType()->GetNumberBonds() == size_t( 1))
            {
              if( !last_terminal_bond_type.IsDefined())
              {
                last_terminal_bond_type = tatom.GetElementType();
              }
              else if( last_terminal_bond_type != tatom.GetElementType())
              {
                n_uniq_term_bond_types = 2;
                break;
              }
            }
            else if( ++n_nonterminal_bonds > size_t( 1))
            {
              break;
            }
          }
        }
        if
        (
          itr_atom->GetAtomType()->GetNumberBonds() <= size_t( 1)
          || ( n_nonterminal_bonds <= size_t( 1) && n_uniq_term_bond_types <= size_t( 1))
        )
        {
          is_terminal( itr_atom.GetPosition()) = 1;
        }
        else
        {
          if( itr_atom->CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, 1))
          {
            is_in_ring( itr_atom.GetPosition()) = 1;
          }
          if( itr_atom->CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_BondOrder, 3))
          {
            for
            (
              storage::Vector< BondConformational>::const_iterator
                itr_bonds( itr_atom->GetBonds().Begin()), itr_bonds_end( itr_atom->GetBonds().End());
              itr_bonds != itr_bonds_end;
              ++itr_bonds
            )
            {
              const AtomConformationalInterface &tatom( itr_bonds->GetTargetAtom());
              if
              (
                itr_bonds->GetBondType()->GetNumberOfElectrons() == size_t( 6)
                &&
                (
                  tatom.GetAtomType()->GetNumberBonds() == size_t( 1)
                  || tatom.GetNumberCovalentlyBoundHydrogens() + tatom.GetNumberofValenceBondsWithOrder( 1)
                     == tatom.GetAtomType()->GetNumberBonds() - 1
                )
              )
              {
                has_triple_bond( itr_atom.GetPosition()) = 1;
                break;
              }
            }
          }
        }
      }
      storage::Vector< sdf::BondInfo> bond_info( MOLECULE.GetBondInfo());

      storage::Vector< BondRigidityEnum> bond_rigidity( bond_info.GetSize(), BondRigidityEnum( e_Flexible));
      storage::Vector< BondRigidityEnum>::iterator itr_bond_rigidity( bond_rigidity.Begin());
      storage::Vector< storage::Vector< sdf::BondInfo> > removed_bonds( n_atoms);
      storage::Vector< size_t> n_terminal_bonds( n_atoms, size_t( 0));

      for
      (
        storage::Vector< sdf::BondInfo>::const_iterator
          itr_bond( bond_info.Begin()), itr_bond_end( bond_info.End());
        itr_bond != itr_bond_end;
        ++itr_bond, ++itr_bond_rigidity
      )
      {
        if
        (
          itr_bond->GetConfigurationalBondType()->IsBondInRing()
          || itr_bond->GetConfigurationalBondType()->GetNumberOfElectrons() > size_t( 2)
        )
        {
          *itr_bond_rigidity = e_Rigid;
          continue;
        }
        const size_t atom_id_a( itr_bond->GetAtomIndexLow());
        const size_t atom_id_b( itr_bond->GetAtomIndexHigh());
        if( has_triple_bond( atom_id_a) || has_triple_bond( atom_id_b))
        {
          *itr_bond_rigidity = e_Rigid;
          continue;
        }

        // check for amide
        //      const AtomConformationalInterface &atom_a( *MOLECULE_GRAPH.GetVertexData( atom_id_a));
        //     const AtomConformationalInterface &atom_b( *MOLECULE_GRAPH.GetVertexData( atom_id_b));
        const AtomType &at_a( MOLECULE_GRAPH.GetVertexData( atom_id_a)->GetAtomType());
        const AtomType &at_b( MOLECULE_GRAPH.GetVertexData( atom_id_b)->GetAtomType());
        // C_TrTrTrPi - C_TrTrTrPi is almost always rather rigid; usually stabilized by resonance
        if
        (
          at_a->GetNumberBonds() == size_t( 3) && at_b->GetNumberBonds() == size_t( 3)
          && at_a->GetNumberElectronsInBonds() >= 4
          && at_b->GetNumberElectronsInBonds() >= 4
          && !is_in_ring( atom_id_a) && !is_in_ring( atom_id_b)
        )
        {
          //BCL_MessageStd( "Haha Conj");
          *itr_bond_rigidity = e_Rigid;
          continue;
        }
        else if
        (
          !is_in_ring( atom_id_a) && !is_in_ring( atom_id_b)
          && at_a->GetHybridOrbitalType() == GetHybridOrbitalTypes().e_SP2
          && at_b->GetHybridOrbitalType() == GetHybridOrbitalTypes().e_SP2
        )
        {
          //BCL_MessageStd( "Haha OTrTr");
          *itr_bond_rigidity = e_Rigid;
          continue;
        }
        if( itr_bond->GetConfigurationalBondType() == GetConfigurationalBondTypes().e_AmideSingleBond && m_ConsiderAmideRigid)
        {
          *itr_bond_rigidity = e_Rigid;
          continue;
        }

        if( is_terminal( atom_id_a) || is_terminal( atom_id_b))
        {
          ++n_terminal_bonds( atom_id_a);
          ++n_terminal_bonds( atom_id_b);
          *itr_bond_rigidity = e_Terminal;
          continue;
        }
        // remove the bond from the molecule graph for simplicity
        removed_bonds( atom_id_a).PushBack( *itr_bond);
        removed_bonds( atom_id_b).PushBack( *itr_bond);
        MOLECULE_GRAPH.RemoveEdge( atom_id_a, atom_id_b);
      }
      storage::List< storage::Vector< size_t> > components( graph::Connectivity::GetComponents( MOLECULE_GRAPH));

      // walk through components, add adjacent bonds that were previously removed
      storage::Vector< size_t> is_in_a_component( n_atoms, size_t( 0));
      for
      (
        storage::List< storage::Vector< size_t> >::iterator
          itr( components.Begin()), itr_end( components.End());
        itr != itr_end;
        ++itr
      )
      {
        storage::Set< size_t> vertices_to_have;
        for
        (
          storage::Vector< size_t>::const_iterator itr_c( itr->Begin()), itr_c_end( itr->End());
          itr_c != itr_c_end;
          ++itr_c
        )
        {
          vertices_to_have.Insert( *itr_c);
          is_in_a_component( *itr_c) = 1;
          for
          (
            storage::Vector< sdf::BondInfo>::const_iterator
              itr_removed( removed_bonds( *itr_c).Begin()), itr_removed_end( removed_bonds( *itr_c).End());
            itr_removed != itr_removed_end;
            ++itr_removed
          )
          {
            const size_t atom_to_add
            (
              *itr_c != itr_removed->GetAtomIndexLow()
              ? itr_removed->GetAtomIndexLow()
              : itr_removed->GetAtomIndexHigh()
            );
            vertices_to_have.Insert( atom_to_add);
          }
        }
        *itr = storage::Vector< size_t>( vertices_to_have.Begin(), vertices_to_have.End());
      }

      size_t atom_index_a( 0);
      for
      (
        storage::Vector< storage::Vector< sdf::BondInfo> >::const_iterator
          itr( removed_bonds.Begin()), itr_end( removed_bonds.End());
        itr != itr_end;
        ++itr, ++atom_index_a
      )
      {
        for
        (
          storage::Vector< sdf::BondInfo>::const_iterator
            itr_removed( itr->Begin()), itr_removed_end( itr->End());
          itr_removed != itr_removed_end;
          ++itr_removed
        )
        {
          if( atom_index_a == itr_removed->GetAtomIndexLow())
          {
            const size_t atom_index_b( itr_removed->GetAtomIndexHigh());
            MOLECULE_GRAPH.AddEdge( atom_index_a, atom_index_b, itr_removed->GetConfigurationalBondType().GetIndex());
            if( is_in_a_component( atom_index_b) || is_in_a_component( atom_index_a))
            {
              continue;
            }
            components.PushBack( storage::Vector< size_t>::Create( atom_index_a, atom_index_b));
          }
        }
      }
      return components;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer FragmentSplitRigid::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Rigid components; defined by breaking all single bonds that are not in a ring, amide, or which connect to a "
        "terminal atom (disregarding H)"
      );
      return parameters;
    }

  } // namespace chemistry
} // namespace bcl
