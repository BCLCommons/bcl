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

#ifndef BCL_CHEMISTRY_ATOM_VECTOR_H_
#define BCL_CHEMISTRY_ATOM_VECTOR_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_complete.h"
#include "bcl_chemistry_atom_configurational_shared.h"
#include "bcl_chemistry_atom_conformational_shared.h"
#include "bcl_chemistry_atom_constitutional_shared.h"
#include "graph/bcl_graph_undirected_edge.h"
#include "iterate/bcl_iterate_generic.h"
#include "sdf/bcl_sdf_atom_info.h"
#include "sdf/bcl_sdf_bond_info.h"
#include "storage/bcl_storage_vector.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AtomVector
    //! @brief base class handling memory for atom classes
    //! @details This class is a template base class handling different atom classes
    //!
    //! @tparam t_Atom the type of the elements in the vector
    //!
    //! @see @link example_chemistry_atom_vector.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date Dec 21, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Atom>
    class AtomVector :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      size_t  m_Size;        //!< size of vector == # atoms
      size_t  m_NumberBonds; //!< # bonds
      t_Atom *m_Atoms;       //!< allocated vector of atoms

    public:

      // iteration typedefs
      typedef t_Atom*       iterator;
      typedef const t_Atom* const_iterator;
      typedef t_Atom&       reference;
      typedef const t_Atom& const_reference;

      //! Default constructor
      AtomVector();

      //! @brief Constructor from a vector of objects that can be used to construct a t_Atom
      //! @param ATOM_INFO atom information about the molecule
      //! @param BONDS bond connectivities of molecule
      AtomVector
      (
        const storage::Vector< sdf::AtomInfo> &ATOM_INFO,
        const storage::Vector< sdf::BondInfo> &BONDS
      );

      //! Copy constructor
      AtomVector( const AtomVector< t_Atom> &BASE);

      //! abstract virtual copy constructor
      virtual AtomVector< t_Atom> *Clone() const;

      //! destructor
      virtual ~AtomVector();

    /////////////////
    // data access //
    /////////////////

      //! @brief get the current class name
      virtual const std::string &GetClassIdentifier() const;

      //! return number of atoms
      const size_t &GetSize() const;

      //! return number of bonds
      const size_t &GetNumberBonds() const;

      //! return reference to changeable element ( POS)
      virtual       t_Atom &operator()( const size_t POS);

      //! return copy of element ( POS)
      virtual const t_Atom &operator()( const size_t POS) const;

      //! return reference to changeable atom from the interface
      virtual       t_Atom &operator()( const typename t_Atom::t_Interface &ATOM);

      //! return reference to the atom given the interface
      virtual const t_Atom &operator()( const typename t_Atom::t_Interface &ATOM) const;

      //! C-style data access with [] gives a pointer on the element
      virtual       t_Atom *operator[]( const size_t POS);

      //! C-style data access with [] gives a pointer on the element
      virtual const t_Atom *operator[]( const size_t POS) const;

      //! return pointer on begin
      virtual       t_Atom *Begin();

      //! return const pointer on begin
      virtual const t_Atom *Begin() const;

      //! return pointer on end
      virtual       t_Atom *End();

      //! return const pointer on end
      virtual const t_Atom *End() const;

      //! return const reference to first element
      virtual const t_Atom &First() const;

      //! return reference to first element
      virtual       t_Atom &First();

      //! return const reference to last element
      virtual const t_Atom &Last() const;

      //! return reference to last element
      virtual       t_Atom &Last();

      //! check whether object is empty
      //! @return true if size == 0
      bool IsEmpty() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief link to a lower level layer
      //! @param LAYER atoms of the next lower layer
      template< typename t_OtherAtomType>
      void LinkToLayer( iterate::Generic< const t_OtherAtomType> LAYER)
      {
        BCL_Assert
        (
          LAYER.GetSize() == m_Size, "Incorrect size for layer, this layer size :"
          + util::Format()( m_Size) + " given layer size : " + util::Format()( LAYER.GetSize())
        );
        for( size_t i( 0); LAYER.NotAtEnd(); ++LAYER, ++i)
        {
          m_Atoms[ i].SetAtom( *LAYER);
        }
      }

      //! @brief set all elements of vector to one given VALUE
      //! @param VALUE value that needs to be set
      AtomVector< t_Atom> &operator =( const AtomVector< t_Atom> &VALUE);

      //! @brief get info on all atoms in the molecule
      //! @return all the atom info about the molecule
      storage::Vector< sdf::AtomInfo> GetAtomInfo() const;

      //! @brief get bonds of molecule
      //! @return all the connectivities of a molecule
      storage::Vector< sdf::BondInfo> GetBondInfo() const;

      //! @brief get index of an atom
      //! @param ATOM atom whose index is required
      //! @return index of the required atom
      size_t GetAtomIndex( const typename t_Atom::t_Interface &ATOM) const;

      //! @brief Reorder the atoms in this atom vector
      //! @param NEW_ORDER Indices to map the indices of this vector onto
      void Reorder( const storage::Vector< size_t> &NEW_ORDER);

      //! @brief return the adjacency list
      //! @param BOND_SCHEME how to represent bond data as a size_t
      //! @return the adjacency list
      storage::Vector< graph::UndirectedEdge< size_t> >
        GetAdjacencyList( const typename t_Atom::t_BondData &BOND_SCHEME) const;

      //! @brief add atoms to the atom vector, connecting them where specified
      //! @param NEW_ATOMS atoms to add to this molecule
      //! @param CONNECTING_BONDS bonds connecting atoms in this object to atoms in the incoming vector. Note:
      //! atom indices in the BondInfos must reflect what they would be AFTER appending the new vector
      //! e.g. to join atom 0 of this object (which has 5 atoms) to atom 0 of the incoming vector, the bondinfo
      //! would be made with indices 0 (atom 0 of this obj) and 5 (for atom 0 of the incoming vector)
      void AddAtomsWithConnectivity
      (
        const AtomVector< t_Atom> &NEW_ATOMS,
        const storage::Vector< sdf::BondInfo> &CONNECTING_BONDS
      );

    protected:

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! check whether position is valid
      void AssertIsValidPosition( const size_t &POS) const;

    }; //end template class AtomVector

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API AtomVector< AtomConstitutionalShared>;
    BCL_EXPIMP_TEMPLATE template class BCL_API AtomVector< AtomConfigurationalShared>;
    BCL_EXPIMP_TEMPLATE template class BCL_API AtomVector< AtomConformationalShared>;
    BCL_EXPIMP_TEMPLATE template class BCL_API AtomVector< AtomComplete>;

  } // namespace chemistry
} // namespace bcl

#endif //BCL_CHEMISTRY_ATOM_VECTOR_H_
