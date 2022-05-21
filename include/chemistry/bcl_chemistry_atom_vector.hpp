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

#ifndef BCL_CHEMISTRY_ATOM_VECTOR_HPP_
#define BCL_CHEMISTRY_ATOM_VECTOR_HPP_

// include the header of this class
#include "bcl_chemistry_atom_vector.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! Default constructor
    template< typename t_Atom>
    AtomVector< t_Atom>::AtomVector() :
      m_Size( 0),
      m_NumberBonds( 0),
      m_Atoms( NULL)
    {
    }

    //! @brief Constructor from a vector of objects that can be used to construct a t_Atom
    //! @param ATOM_INFO atom information about the molecule
    //! @param BONDS bond connectivities of molecule
    template< typename t_Atom>
    AtomVector< t_Atom>::AtomVector
    (
      const storage::Vector< sdf::AtomInfo> &ATOM_INFO,
      const storage::Vector< sdf::BondInfo> &BONDS
    ) :
      m_Size( ATOM_INFO.GetSize()),
      m_NumberBonds( 0),
      m_Atoms( new t_Atom[ m_Size])
    {
      // iterate over all atoms and set them in m_Data
      for( size_t i( 0); i < m_Size; ++i)
      {
        m_Atoms[ i] = t_Atom( ATOM_INFO( i));
      }
      // initialize a vector of maps for each atom, where each map describes the target atom and type of bond
      storage::Vector< storage::Map< size_t, ConfigurationalBondType> > atom_to_bonds( m_Size);
      for
      (
        storage::Vector< sdf::BondInfo>::const_iterator
          itr( BONDS.Begin()), itr_end( BONDS.End());
        itr != itr_end;
        ++itr
      )
      {
        // set bonds in both directions
        atom_to_bonds( itr->GetAtomIndexLow())[ itr->GetAtomIndexHigh()] = itr->GetConfigurationalBondType();
        atom_to_bonds( itr->GetAtomIndexHigh())[ itr->GetAtomIndexLow()] = itr->GetConfigurationalBondType();
      }

      // set the bonds on the atoms
      // compute the # of bonds here too, since it is possible that the bonds vector had non-unique bonds, we cannot
      // trust that # bonds is given by the vector
      for( size_t i( 0); i < m_Size; ++i)
      {
        storage::Vector< typename t_Atom::t_Bond> bonds;

        const size_t number_bonds( atom_to_bonds( i).GetSize());
        bonds.AllocateMemory( number_bonds);

        for
        (
          typename storage::Map< size_t, ConfigurationalBondType>::const_iterator
            itr_bond_map( atom_to_bonds( i).Begin()), itr_bond_map_end( atom_to_bonds( i).End());
          itr_bond_map != itr_bond_map_end;
          ++itr_bond_map
        )
        {
          bonds.PushBack( typename t_Atom::t_Bond( m_Atoms[ itr_bond_map->first], itr_bond_map->second));
        }
        m_Atoms[ i].SetBonds( bonds);
        m_NumberBonds += number_bonds;
      }
      // divide # bonds by 2 since each bond was counted twice
      m_NumberBonds /= 2;
    }

    //! Copy constructor
    template< typename t_Atom>
    AtomVector< t_Atom>::AtomVector( const AtomVector< t_Atom> &BASE) :
      m_Size( BASE.m_Size),
      m_NumberBonds( BASE.m_NumberBonds),
      m_Atoms( new t_Atom[ m_Size])
    {
      // set all elements
      for( size_t i( 0); i < m_Size; ++i)
      {
        m_Atoms[ i] = BASE.m_Atoms[ i];
      }
    }

    //! destructor
    template< typename t_Atom>
    AtomVector< t_Atom>::~AtomVector()
    {
      delete [] m_Atoms;
    }

    //! virtual copy constructor
    template< typename t_Atom>
    AtomVector< t_Atom> *AtomVector< t_Atom>::Clone() const
    {
      return new AtomVector< t_Atom>( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get the current class name
    template< typename t_Atom>
    const std::string &AtomVector< t_Atom>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! return number of elements
    template< typename t_Atom>
    const size_t &AtomVector< t_Atom>::GetSize() const
    {
      return m_Size;
    }

    //! return number of bonds
    template< typename t_Atom>
    const size_t &AtomVector< t_Atom>::GetNumberBonds() const
    {
      return m_NumberBonds;
    }

    //! return reference to changeable element ( POS)
    template< typename t_Atom>
          t_Atom &AtomVector< t_Atom>::operator()( const size_t POS)
    {
      AssertIsValidPosition( POS);
      return m_Atoms[ POS];
    }

    //! return copy of element ( POS)
    template< typename t_Atom>
    const t_Atom &AtomVector< t_Atom>::operator()( const size_t POS) const
    {
      AssertIsValidPosition( POS);
      return m_Atoms[ POS];
    }

    //! return reference to changeable atom from the interface
    template< typename t_Atom>
    t_Atom &AtomVector< t_Atom>::operator()( const typename t_Atom::t_Interface &ATOM)
    {
      return m_Atoms[ GetAtomIndex( ATOM)];
    }

    //! return reference to the atom given the interface
    template< typename t_Atom>
    const t_Atom &AtomVector< t_Atom>::operator()( const typename t_Atom::t_Interface &ATOM) const
    {
      return m_Atoms[ GetAtomIndex( ATOM)];
    }

    //! C-style data access with [] gives a pointer on the element
    template< typename t_Atom>
          t_Atom *AtomVector< t_Atom>::operator[]( const size_t POS)
    {
      AssertIsValidPosition( POS);
      return m_Atoms + POS;
    }

    //! C-style data access with [] gives a pointer on the element
    template< typename t_Atom>
    const t_Atom *AtomVector< t_Atom>::operator[]( const size_t POS) const
    {
      AssertIsValidPosition( POS);
      return m_Atoms + POS;
    }

    //! return pointer on begin
    template< typename t_Atom>
          t_Atom *AtomVector< t_Atom>::Begin()
    {
      return m_Atoms;
    }

    //! return const pointer on begin
    template< typename t_Atom>
    const t_Atom *AtomVector< t_Atom>::Begin() const
    {
      return m_Atoms;
    }

    //! return pointer on end
    template< typename t_Atom>
          t_Atom *AtomVector< t_Atom>::End()
    {
      return m_Atoms + m_Size;
    }

    //! return const pointer on end
    template< typename t_Atom>
    const t_Atom *AtomVector< t_Atom>::End() const
    {
      return m_Atoms + m_Size;
    }

    //! return const reference to first element
    template< typename t_Atom>
    const t_Atom &AtomVector< t_Atom>::First() const
    {
      BCL_Assert( !IsEmpty(), "no first element in empty object");
      return m_Atoms[ 0];
    }

    //! return reference to first element
    template< typename t_Atom>
          t_Atom &AtomVector< t_Atom>::First()
    {
      BCL_Assert( !IsEmpty(), "no first element in empty object");
      return m_Atoms[ 0];
    }

    //! return const reference to last element
    template< typename t_Atom>
    const t_Atom &AtomVector< t_Atom>::Last() const
    {
      BCL_Assert( !IsEmpty(), "no last element in empty object");
      return m_Atoms[ m_Size - 1];
    }

    //! return reference to last element
    template< typename t_Atom>
          t_Atom &AtomVector< t_Atom>::Last()
    {
      BCL_Assert( !IsEmpty(), "no last element in empty object");
      return m_Atoms[ m_Size - 1];
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief set all elements of vector to one given VALUE
    //! @param VALUE value that needs to be set
    template< typename t_Atom>
    AtomVector< t_Atom> &AtomVector< t_Atom>::operator =( const AtomVector< t_Atom> &VALUE)
    {
      // set all elements
      if( this != &VALUE)
      {
        if( m_Size != VALUE.m_Size)
        {
          m_Size = VALUE.m_Size;
          delete [] m_Atoms;
          if( m_Size != 0)
          {
            m_Atoms = new t_Atom[ m_Size];
          }
          else
          {
            m_Atoms = NULL;
          }
        }
        for( size_t i( 0); i < m_Size; ++i)
        {
          m_Atoms[ i] = VALUE.m_Atoms[ i];
        }
        m_NumberBonds = VALUE.m_NumberBonds;
      }

      //end
      return *this;
    }

    //! @brief get index of an atom
    //! @param ATOM atom whose index is required
    //! @return index of the required atom
    template< typename t_Atom>
    size_t AtomVector< t_Atom>::GetAtomIndex( const typename t_Atom::t_Interface &ATOM) const
    {
      const size_t index( ( const t_Atom *)( &ATOM) - m_Atoms);
      BCL_Assert( index < m_Size, "tried to get index of atom not in vector");
      return index;
    }

    //! @brief get info on all atoms in the molecule
    //! @return all the atom info about the molecule
    template< typename t_Atom>
    storage::Vector< sdf::AtomInfo> AtomVector< t_Atom>::GetAtomInfo() const
    {
      storage::Vector< sdf::AtomInfo> atom_info( m_Size);

      for( size_t i( 0); i < m_Size; ++i)
      {
        atom_info( i) = m_Atoms[ i].GetAtomInfo();
      }
      return atom_info;
    }

    //! @brief get bonds of molecule
    //! @return all the connectivities of a molecule
    template< typename t_Atom>
    storage::Vector< sdf::BondInfo> AtomVector< t_Atom>::GetBondInfo() const
    {
      storage::Vector< sdf::BondInfo> bond_info;
      bond_info.AllocateMemory( m_NumberBonds);

      for( size_t i( 0); i < m_Size; ++i)
      {
        for
        (
          typename storage::Vector< typename t_Atom::t_Bond>::const_iterator
            itr_bond( m_Atoms[ i].GetBonds().Begin()), itr_bond_end( m_Atoms[ i].GetBonds().End());
          itr_bond != itr_bond_end;
          ++itr_bond
        )
        {
          const size_t target_atom_id( GetAtomIndex( itr_bond->GetTargetAtom()));

          // record the bond info; to prevent duplicates, only pushback if the target atom id is larger
          if( target_atom_id > i)
          {
            bond_info.PushBack( sdf::BondInfo( i, target_atom_id, itr_bond->GetBondType()));
          }
        }
      }
      return bond_info;
    }

    //! @brief Reorder the atoms in this atom vector
    //! @param NEW_ORDER Indices to map the indices of this vector onto
    template< typename t_Atom>
    void AtomVector< t_Atom>::Reorder( const storage::Vector< size_t> &NEW_ORDER)
    {
      // compute the reverse reordering so that bonds can be updated appropriately
      storage::Vector< size_t> new_order_inverse( m_Size, util::GetUndefined< size_t>());
      const size_t new_size( NEW_ORDER.GetSize());
      for( size_t i( 0); i < new_size; ++i)
      {
        BCL_Assert( NEW_ORDER( i) < m_Size, "Reorder was called with bad order vector (index out of range)");
        BCL_Assert
        (
          !util::IsDefined( new_order_inverse( NEW_ORDER( i))),
          "Reorder was called with bad order vector (indices non-unique)"
        );
        new_order_inverse( NEW_ORDER( i)) = i;
      }
      t_Atom *new_data( new t_Atom[ NEW_ORDER.GetSize()]);

      // reset # of bonds, which may change here
      m_NumberBonds = 0;

      // compute the difference between the current vector and the new vector
      ptrdiff_t difference( ( char *)( new_data) - ( char *)( m_Atoms));
      for( size_t i( 0), new_size( NEW_ORDER.GetSize()); i < new_size; ++i)
      {
        // copy the atom vector, accounting for the actual difference in atom vector
        new_data[ i].Copy( m_Atoms[ NEW_ORDER( i)], difference);

        storage::Map< size_t, typename t_Atom::t_Bond> bond_indices_to_bonds_map;

        for
        (
          typename storage::Vector< typename t_Atom::t_Bond>::const_iterator
            itr_bond( new_data[ i].GetBonds().Begin()), itr_bond_end( new_data[ i].GetBonds().End());
          itr_bond != itr_bond_end;
          ++itr_bond
        )
        {
          const typename t_Atom::t_Bond &bond_to_copy( *itr_bond);

          // get the id of the atom prior to the remapping
          const size_t target_atom_id( ( const t_Atom *)( &bond_to_copy.GetTargetAtom()) - new_data);

          // get the new target atom id
          const size_t new_target_atom_id( new_order_inverse( target_atom_id));

          // skip bonds to atoms that will no longer be in the molecule
          if( !util::IsDefined( new_target_atom_id))
          {
            continue;
          }

          // set the bond up to point to the newly-mapped atom.
          // the position of remapped atoms is given by the new_order_inverse function
          bond_indices_to_bonds_map[ new_target_atom_id] =
            typename t_Atom::t_Bond( new_data[ new_target_atom_id], bond_to_copy.GetBondType());
        }

        new_data[ i].SetBonds( bond_indices_to_bonds_map.GetMappedValues());
        m_NumberBonds += new_data[ i].GetBonds().GetSize();
      }

      // delete m_Data and set it equal to new data, copy size too
      delete[] m_Atoms;
      m_Atoms = new_data;
      m_Size = new_size;
      // divide # bonds by 2 since each bond was counted twice
      m_NumberBonds /= 2;
    }

    //! @brief return the adjacency list
    //! @param BOND_SCHEME how to represent bond data as a size_t
    //! @return the adjacency list
    template< typename t_Atom>
    storage::Vector< graph::UndirectedEdge< size_t> >
      AtomVector< t_Atom>::GetAdjacencyList( const typename t_Atom::t_BondData &BOND_SCHEME) const
    {
      // iterate through the bonds and use the data to populate the connectivity list
      storage::Vector< graph::UndirectedEdge< size_t> > adjacency_list;
      adjacency_list.AllocateMemory( m_NumberBonds);

      size_t atom_index_a( 0);
      for( const t_Atom *itr( m_Atoms), *itr_end( m_Atoms + m_Size); itr != itr_end; ++itr, ++atom_index_a)
      {
        for
        (
          typename storage::Vector< typename t_Atom::t_Bond>::const_iterator
            itr_bond( itr->GetBonds().Begin()), itr_bond_end( itr->GetBonds().End());
          itr_bond != itr_bond_end;
          ++itr_bond
        )
        {
          // get the index of the atom attached to *itr
          const size_t atom_index_b( GetAtomIndex( itr_bond->GetTargetAtom()));

          // check atom indices so that undirected edges are not added twice
          if( atom_index_a <= atom_index_b)
          {
            // create the new edge
            graph::UndirectedEdge< size_t> new_edge
            (
              atom_index_a,
              atom_index_b,
              itr_bond->GetBondType()->GetBondData( BOND_SCHEME)
            );
            adjacency_list.PushBack( new_edge);
          }
        }
      }
      return adjacency_list;
    }

    //! @brief add atoms to the atom vector, connecting them where specified
    //! @param NEW_ATOMS atoms to add to this molecule
    //! @param CONNECTING_BONDS bonds connecting atoms in this object to atoms in the incoming vector. Note:
    //! atom indices in the BondInfos must reflect what they would be AFTER appending the new vector
    //! e.g. to join atom 0 of an object which has 5 atoms to atom 0 of an incoming vector, the bondinfo
    //! would be made with indices 0 (atom 0 of current obj) and 5 (for atom 0 of the incoming vector)
    template< typename t_Atom>
    void AtomVector< t_Atom>::AddAtomsWithConnectivity
    (
      const AtomVector< t_Atom> &NEW_ATOMS,
      const storage::Vector< sdf::BondInfo> &CONNECTING_BONDS
    )
    {
      if( NEW_ATOMS.m_Size == size_t( 0))
      {
        // no atoms to add, so just return
        return;
      }

      // make a pointer to the atoms already held by this vector
      t_Atom *current_atoms( m_Atoms);

      // allocate a new vector large enough to hold all atoms together
      m_Atoms = new t_Atom[ m_Size + NEW_ATOMS.m_Size];

      // copy the atoms previously held
      for( size_t i( 0); i < m_Size; ++i)
      {
        m_Atoms[ i] = current_atoms[ i];
      }

      // copy the atoms held in the new atoms vector
      for( size_t i( 0); i < NEW_ATOMS.m_Size; ++i)
      {
        m_Atoms[ m_Size + i] = NEW_ATOMS.m_Atoms[ i];
      }

      delete [] current_atoms;

      m_Size += NEW_ATOMS.m_Size;

      // if there are no new bonds to create, just return
      if( CONNECTING_BONDS.IsEmpty())
      {
        return;
      }

      // if there are bonds between the atoms, then we need to recreate the bond vectors

      // initialize a vector of maps for each atom, where each map describes the target atom and type of bond
      storage::Vector< storage::Map< size_t, ConfigurationalBondType> > atom_to_bonds( m_Size);
      for
      (
        storage::Vector< sdf::BondInfo>::const_iterator
          itr( CONNECTING_BONDS.Begin()), itr_end( CONNECTING_BONDS.End());
        itr != itr_end;
        ++itr
      )
      {
        // set bonds in both directions
        atom_to_bonds( itr->GetAtomIndexLow())[ itr->GetAtomIndexHigh()] = itr->GetConfigurationalBondType();
        atom_to_bonds( itr->GetAtomIndexHigh())[ itr->GetAtomIndexLow()] = itr->GetConfigurationalBondType();
      }

      // recompute # of bonds
      m_NumberBonds = 0;

      // set the bonds on the atoms
      for( size_t i( 0); i < m_Size; ++i)
      {
        if( atom_to_bonds( i).IsEmpty()) // if there are no bonds to add, continue
        {
          m_NumberBonds += m_Atoms[ i].GetBonds().GetSize();
          continue;
        }

        storage::Vector< typename t_Atom::t_Bond> bonds;
        bonds.AllocateMemory( m_Atoms[ i].GetBonds().GetSize() + atom_to_bonds( i).GetSize());
        bonds.Append( m_Atoms[ i].GetBonds());

        for
        (
          storage::Map< size_t, ConfigurationalBondType>::const_iterator
            itr_bond_map( atom_to_bonds( i).Begin()), itr_bond_map_end( atom_to_bonds( i).End());
          itr_bond_map != itr_bond_map_end;
          ++itr_bond_map
        )
        {
          bonds.PushBack
          (
            typename t_Atom::t_Bond( m_Atoms[ itr_bond_map->first], itr_bond_map->second)
          );
        }
        m_Atoms[ i].SetBonds( bonds);
        m_NumberBonds += bonds.GetSize();
      }
      // divide # bonds by 2 since each bond was counted twice
      m_NumberBonds /= 2;
    }

  //////////////////////
  // input and output //
  //////////////////////

    template< typename t_Atom>
    std::istream &AtomVector< t_Atom>::Read( std::istream &ISTREAM)
    {
      // read atom info
      storage::Vector< sdf::AtomInfo> atom_info;
      io::Serialize::Read( atom_info, ISTREAM);

      // read bond info
      storage::Vector< sdf::BondInfo> bond_info;
      io::Serialize::Read( bond_info, ISTREAM);

      // create atom vector
      *this = AtomVector< t_Atom>( atom_info, bond_info);

      // end
      return ISTREAM;
    }

    //! write AtomVector to std::ostream
    template< typename t_Atom>
    std::ostream &AtomVector< t_Atom>::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( GetAtomInfo(), OSTREAM, INDENT) << '\n';
      io::Serialize::Write( GetBondInfo(), OSTREAM, INDENT);
      return OSTREAM;
    }

  ////////////////
  // operations //
  ////////////////

    //! check whether position is valid
    template< typename t_Atom>
    void AtomVector< t_Atom>::AssertIsValidPosition( const size_t &POS) const
    {
      BCL_Assert
      (
        POS < m_Size,
        "cannot access atom outside range! " + util::Format()( POS) + " >= " + util::Format()( m_Size)
      );
    }

    //! check whether object is empty
    //! @return true if size == 0
    template< typename t_Atom>
    bool AtomVector< t_Atom>::IsEmpty() const
    {
      return m_Size == size_t( 0);
    }

  } // namespace chemistry
} // namespace bcl

#endif //BCL_CHEMISTRY_ATOM_VECTOR_HPP_
