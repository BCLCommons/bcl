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
#include "chemistry/bcl_chemistry_ring_fragment_map.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> RingFragmentMap::s_Instance
    (
      GetObjectInstances().AddInstance( new RingFragmentMap())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    RingFragmentMap::RingFragmentMap()
    {
    }

    //! @brief constructor given atoms molecule, fragment atom that need to be connected
    //! @params MOLECULE_ATOM connecting atom in molecule adjacent to the fragment atom
    //! @params FRAGMENT_ATOM fragment atom that is connected to another atom of molecule
    //! @params MOLECULE_ATOM_IN_RING whether MOLECULE_ATOM is in a ring
    RingFragmentMap::RingFragmentMap
    (
      const size_t MOLECULE_ATOM,
      bool MOLECULE_ATOM_IN_RING,
      const storage::Vector< size_t> &CONNECTED_VERTICES,
      const storage::Vector< size_t> &CONNECTED_VERTICES_FOR_ROTATION
    ) :
      m_MoleculeAtom( MOLECULE_ATOM),
      m_MoleculeAtomInRing( MOLECULE_ATOM_IN_RING),
      m_ConnectedVertices( CONNECTED_VERTICES),
      m_ConnectedVerticesForRotation( CONNECTED_VERTICES_FOR_ROTATION)
    {
    }

    //! @brief Clone function
    //! @return pointer to new RingFragmentMap
    RingFragmentMap *RingFragmentMap::Clone() const
    {
      return new RingFragmentMap( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &RingFragmentMap::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns map between indices of Node to Parent
    //! @return the map between indices of Node to Parent
    const size_t &RingFragmentMap::GetMoleculeAtom() const
    {
      return m_MoleculeAtom;
    }

    //! @brief returns map between indices of Parent to This
    //! @return the map between indices of Parent to This
    const storage::Vector< size_t> &RingFragmentMap::GetConnectedVertices() const
    {
      return m_ConnectedVertices;
    }

    //! @brief returns map between indices of Parent to This
    //! @return the map between indices of Parent to This
    const storage::Vector< size_t> &RingFragmentMap::GetConnectedVerticesForRotation() const
    {
      return m_ConnectedVerticesForRotation;
    }

    //! @brief returns map between indices of Node to Parent
    //! @return the map between indices of Node to Parent
    const bool &RingFragmentMap::GetMoleculeAtomInRing() const
    {
      return m_MoleculeAtomInRing;
    }

  ////////////////
  // operations //
  ////////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RingFragmentMap::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_MoleculeAtom, ISTREAM);
      io::Serialize::Read( m_MoleculeAtomInRing, ISTREAM);
      io::Serialize::Read( m_ConnectedVertices, ISTREAM);

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &RingFragmentMap::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_MoleculeAtom, OSTREAM, INDENT) << '\t' << m_MoleculeAtomInRing << '\n';
      io::Serialize::Write( m_ConnectedVertices, OSTREAM, INDENT) << '\n';
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace chemistry
} // namespace bcl
