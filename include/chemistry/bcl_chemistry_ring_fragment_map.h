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

#ifndef BCL_CHEMISTRY_RING_FRAGMENT_MAP_H_
#define BCL_CHEMISTRY_RING_FRAGMENT_MAP_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_conformational_interface.h"
#include "bcl_chemistry_fragment_complete.h"
#include "graph/bcl_graph_const_graph.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RingFragmentMap
    //! @brief a helper class for the FragmentMolecule class
    //! @details stores a
    //!
    //! @see @link example_chemistry_ring_fragment_map.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date May 28, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RingFragmentMap :
      public util::ObjectInterface
    {

    private:

        size_t                                m_MoleculeAtom;

        bool                                  m_MoleculeAtomInRing;

        storage::Vector< size_t>              m_ConnectedVertices;
        storage::Vector< size_t>              m_ConnectedVerticesForRotation;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      RingFragmentMap();

      //! @brief constructor given atoms molecule, fragment atom that need to be connected
      //! @params MOLECULE_ATOM connecting atom in molecule adjacent to the fragment atom
      //! @params FRAGMENT_ATOM fragment atom that is connected to another atom of molecule
      //! @params MOLECULE_ATOM_IN_RING whether MOLECULE_ATOM is in a ring
      RingFragmentMap
      (
        const size_t MOLECULE_ATOM,
        bool MOLECULE_ATOM_IN_RING,
        const storage::Vector< size_t> &CONNECTED_VERTICES,
        const storage::Vector< size_t> &CONNECTED_VERTICES_FOR_ROTATION
      );

      //! @brief Clone function
      //! @return pointer to new RingFragmentMap
      RingFragmentMap *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns vector of indices from this RingFragmentMap to the original molecule
      //! @return the vector of indices from this RingFragmentMap to the original molecule
      const size_t &GetMoleculeAtom() const;

      //! @brief returns vector of indices from this RingFragmentMap to the parent molecule
      //! @return the vector of indices from this RingFragmentMap to the parent molecule
      const storage::Vector< size_t> &GetConnectedVertices() const;

      //! @brief returns vector of indices from this RingFragmentMap to the parent molecule
      //! @return the vector of indices from this RingFragmentMap to the parent molecule
      const storage::Vector< size_t> &GetConnectedVerticesForRotation() const;

      //! @brief returns set of indices from this RingFragmentMap to the original molecule
      //! @return the set of indices from this RingFragmentMap to the original molecule
      const bool &GetMoleculeAtomInRing() const;

    ////////////////
    // operations //
    ////////////////

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
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class RingFragmentMap

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_RING_FRAGMENT_MAP_H_
