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

#ifndef BCL_CHEMISTRY_VALENCE_HANDLER_H_
#define BCL_CHEMISTRY_VALENCE_HANDLER_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_complete.h"
#include "bcl_chemistry_atom_vector.h"
#include "sdf/bcl_sdf_molecule_reading_pref.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ValenceHandler
    //! @brief Standardizes FragmentComplete object
    //!
    //! @see @link example_chemistry_valence_handler.cpp @endlink
    //! @author kothiwsk
    //! @date Oct 19, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API ValenceHandler :
      public util::ObjectInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new ValenceHandler
      ValenceHandler *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief determine the coordinates of all missing hydrogens for a given atom
      //! @param ATOM an AtomComplete that should have both type and position defined
      //! @return a storage vector containing the positions of all missing hydrogens
      static storage::Vector< linal::Vector3D> DetermineCoordinates( const AtomConformationalInterface &ATOM);

      //! @brief get the idealized geometry for a particular hybrid orbital type and bond length
      //! @param ORBITAL_TYPE the type of orbital associated with the geometry
      //! @param BOND_LENGTH the length of the bonds in the geometry
      static storage::Vector< linal::Vector3D> GetIdealizedGeometry
      (
        const HybridOrbitalType &ORBITAL_TYPE,
        const double &BOND_LENGTH
      );

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

    }; // class ValenceHandler

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_VALENCE_HANDLER_H_

