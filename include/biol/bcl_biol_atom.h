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

#ifndef BCL_BIOL_ATOM_H_
#define BCL_BIOL_ATOM_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_biol_atom_types.h"
#include "coord/bcl_coord_movable_interface.h"
#include "linal/bcl_linal_vector_3d.h"
#include "util/bcl_util_si_ptr.h"

#undef FindAtom

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Atom
    //! @brief this class represents a biological atom
    //! @details This representation excludes chemical information that are not required for biological classes.
    //! The only information it stores is the AtomType, pdb id and a b-factor and the coordinates
    //!
    //! @see @link example_biol_atom.cpp @endlink
    //! @author karakam, woetzen
    //! @date 11/26/2007
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Atom :
      public coord::MovableInterface
    {

    private:

    //////////
    // data //
    //////////

      linal::Vector3D m_Coordinates; //!< atom coordinates
      AtomType       m_Type;        //!< Atom type related data
      size_t         m_PdbID;       //!< ID provided in the PDB file for this atom
      double         m_BFactor;     //!< BFactor value defined in the PDB

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct undefined atom
      Atom();

      //! @brief construct atom from a single given ATOM_TYPE
      //! @param ATOM_TYPE AtomType of interest
      Atom( const AtomType &ATOM_TYPE);

      //! @brief construct atom from atom type ( string name can also be used instead), pdb id and b factor
      //! @param ATOM_TYPE AtomType of interest
      //! @param PDB_ID pdb id of atom
      //! @param B_FACTOR B factor of atom
      Atom
      (
        const AtomType &ATOM_TYPE,
        const size_t PDB_ID,
        const double B_FACTOR = util::GetUndefinedDouble()
      );

      //! @brief construct atom from coordinates, atom type (string name can also be used instead), pdb id and b factor
      //! @param COORDINATES 3D coordinates of the atom
      //! @param ATOM_TYPE AtomType of interest
      //! @param PDB_ID pdb id of atom
      //! @param B_FACTOR B factor of atom
      Atom
      (
        const linal::Vector3D &COORDINATES,
        const AtomType &ATOM_TYPE,
        const size_t PDB_ID = util::GetUndefinedSize_t(),
        const double B_FACTOR = util::GetUndefinedDouble()
      );

      //! @brief virtual copy constructor
      Atom *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the atom type
      //! @return the atom type
      const AtomType &GetType() const
      {
        return m_Type;
      }

      //! @brief return the atom pdb id
      //! @return the atom pdb id
      size_t GetPdbID() const
      {
        return m_PdbID;
      }

      //! @brief return the atom B factor
      //! @return atom B factor
      double GetBFactor() const
      {
        return m_BFactor;
      }

      //! return coordinates as Vector3D
      const linal::Vector3D &GetCoordinates() const;

      //! change coordinates and increase state
      void SetCoordinates( const linal::Vector3D &COORDINATES);

    ////////////////
    // operations //
    ////////////////

      //! @brief static function to find an atom of the specified ATOM_TYPE from given ATOM_LIST
      //! @param ATOM_LIST SiPtrVector of Atoms
      //! @param ATOM_TYPE AtomType of interest
      //! @return SiPtr to requested ATOM_TYPE from given ATOM_LIST  return undefined SiPtr if not found
      static const util::SiPtr< const Atom> &FindAtom
      (
        const util::SiPtrVector< const Atom> &ATOM_LIST,
        const AtomType &ATOM_TYPE
      );

      //! @brief AllCoordinatesDefined determines if the x,y,and z coordinates for the atom are defined
      //! @return returns a bool - true if all coordinates are defined, false otherwise
      bool AllCoordinatesDefined() const;

      //! transforms the coordinates with a transformationmatrix3D
      void Transform( const math::TransformationMatrix3D &TRANSFORMATIONMATRIX3D);

      //! translate coordinates
      void Translate( const linal::Vector3D &TRANSLATION);

      //! rotates coordinates
      void Rotate( const math::RotationMatrix3D &ROTATIONMATRIX3D);

      //! returns the Center
      linal::Vector3D GetCenter() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream
      //! @return ostream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; //class Atom

  ////////////////
  // operations //
  ////////////////

    //! returns square distance between two atoms
    BCL_API double SquareDistance( const Atom &ATOM_A, const Atom &ATOM_B);

    //! returns distance between two atoms
    BCL_API double Distance( const Atom &ATOM_A, const Atom &ATOM_B);

    //! returns angle between three atoms
    BCL_API double ProjAngle( const Atom &ATOM_A, const Atom &ATOM_B, const Atom &ATOM_C);

    //! returns dihedral between four atoms
    BCL_API double Dihedral( const Atom &ATOM_A, const Atom &ATOM_B, const Atom &ATOM_C, const Atom &ATOM_D);

  } // namespace biol
} // namespace bcl

#endif // BCL_BIOL_ATOM_H_

