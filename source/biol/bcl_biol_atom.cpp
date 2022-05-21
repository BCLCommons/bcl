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
#include "biol/bcl_biol_atom.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d_operations.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Atom::s_Instance
    (
      GetObjectInstances().AddInstance( new Atom())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct undefined atom
    Atom::Atom() :
      m_Coordinates( util::GetUndefined< double>()),
      m_Type( util::GetUndefined< AtomType>()),
      m_PdbID( util::GetUndefined< size_t>()),
      m_BFactor( util::GetUndefined< double>())
    {
    }

    //! @brief construct atom from a single given ATOM_TYPE
    //! @param ATOM_TYPE AtomType of interest
    Atom::Atom( const AtomType &ATOM_TYPE) :
      m_Coordinates( util::GetUndefined< double>()),
      m_Type( ATOM_TYPE),
      m_PdbID( util::GetUndefined< size_t>()),
      m_BFactor( util::GetUndefined< double>())
    {
    }

    //! @brief construct atom from atom type ( string name can also be used instead), pdb id and b factor
    //! @param ATOM_TYPE AtomType of interest
    //! @param PDB_ID pdb id of atom
    //! @param B_FACTOR B factor of atom
    Atom::Atom
    (
      const AtomType &ATOM_TYPE,
      const size_t PDB_ID,
      const double B_FACTOR
    ) :
      m_Coordinates( util::GetUndefined< double>()),
      m_Type( ATOM_TYPE),
      m_PdbID( PDB_ID),
      m_BFactor( B_FACTOR)
    {
    }

    //! @brief construct atom from coordinates, atom type (string name can also be used instead), pdb id and b factor
    //! @param COORDINATES 3D coordinates of the atom
    //! @param ATOM_TYPE AtomType of interest
    //! @param PDB_ID pdb id of atom
    //! @param B_FACTOR B factor of atom
    Atom::Atom
    (
      const linal::Vector3D &COORDINATES,
      const AtomType &ATOM_TYPE,
      const size_t PDB_ID,
      const double B_FACTOR
    ) :
      m_Coordinates( COORDINATES),
      m_Type( ATOM_TYPE),
      m_PdbID( PDB_ID),
      m_BFactor( B_FACTOR)
    {
    }

    //! @brief virtual copy constructor
    Atom *Atom::Clone() const
    {
      return new Atom( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Atom::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! return coordinates as Vector3D
    const linal::Vector3D &Atom::GetCoordinates() const
    {
      return m_Coordinates;
    }

    //! change coordinates and increase state
    void Atom::SetCoordinates( const linal::Vector3D &COORDINATES)
    {
      m_Coordinates = COORDINATES;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief static function to find an atom of the specified ATOM_TYPE from given ATOM_LIST
    //! @param ATOM_LIST SiPtrVector of Atoms
    //! @param ATOM_TYPE AtomType of interest
    //! @return SiPtr to requested ATOM_TYPE from given ATOM_LIST  return undefined SiPtr if not found
    const util::SiPtr< const Atom> &Atom::FindAtom
    (
      const util::SiPtrVector< const Atom> &ATOM_LIST,
      const AtomType &ATOM_TYPE
    )
    {
      // initialize undefined atom
      static const Atom s_undefined_atom;

      // initialize pointer to undefined atom
      static const util::SiPtr< const Atom> s_undefined_atom_pointer( s_undefined_atom);

      // iterate over atoms in the list
      for
      (
        util::SiPtrVector< const Atom>::const_iterator atom_itr( ATOM_LIST.Begin()), atom_itr_end( ATOM_LIST.End());
        atom_itr != atom_itr_end;
        ++atom_itr
      )
      {
        // if the type of this atom matches the specified type, return it
        if( ( *atom_itr)->GetType() == ATOM_TYPE)
        {
          return *atom_itr;
        }
      }

      // if the atom was not found return the simple pointer to undefined atom
      return s_undefined_atom_pointer;
    }

    //! transforms the coordinates with a transformationmatrix3D
    void Atom::Transform( const math::TransformationMatrix3D &TRANSFORMATIONMATRIX3D)
    {
      m_Coordinates.Transform( TRANSFORMATIONMATRIX3D);
    }

    //! translate coordinates
    void Atom::Translate( const linal::Vector3D &TRANSLATION)
    {
      m_Coordinates.Translate( TRANSLATION);
    }

    //! rotates coordinates
    void Atom::Rotate( const math::RotationMatrix3D &ROTATIONMATRIX3D)
    {
      m_Coordinates.Rotate( ROTATIONMATRIX3D);
    }

    //! returns the Center
    linal::Vector3D Atom::GetCenter() const
    {
      return m_Coordinates;
    }

    //! @brief AllCoordinatesDefined determines if the x,y,and z coordinates for the atom are defined
    //! @return returns a bool - true if all coordinates are defined, false otherwise
    bool Atom::AllCoordinatesDefined() const
    {
      return m_Coordinates.IsDefined();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! read Atom from io::IFStream
    std::istream &Atom::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Coordinates, ISTREAM);
      io::Serialize::Read( m_Type       , ISTREAM);
      io::Serialize::Read( m_PdbID      , ISTREAM);
      io::Serialize::Read( m_BFactor    , ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &Atom::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Coordinates, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Type       , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PdbID      , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_BFactor    , OSTREAM, INDENT);

      // end and return the ostream
      return OSTREAM;
    }

  ////////////////
  // operations //
  ////////////////

    //! returns square distance between two atoms
    double SquareDistance( const Atom &ATOM_A, const Atom &ATOM_B)
    {
      return linal::SquareDistance
      (
        ATOM_A.GetCoordinates(),
        ATOM_B.GetCoordinates()
      );
    }

    //! returns distance between two atoms
    double Distance( const Atom &ATOM_A, const Atom &ATOM_B)
    {
      return linal::Distance
      (
        ATOM_A.GetCoordinates(),
        ATOM_B.GetCoordinates()
      );
    }

    //! returns angle between three atoms
    double ProjAngle( const Atom &ATOM_A, const Atom &ATOM_B, const Atom &ATOM_C)
    {
      return linal::ProjAngle
      (
        ATOM_A.GetCoordinates(),
        ATOM_B.GetCoordinates(),
        ATOM_C.GetCoordinates()
      );
    }

    //! returns dihedral between four atoms
    double Dihedral( const Atom &ATOM_A, const Atom &ATOM_B, const Atom &ATOM_C, const Atom &ATOM_D)
    {
      return linal::Dihedral
      (
        ATOM_A.GetCoordinates(),
        ATOM_B.GetCoordinates(),
        ATOM_C.GetCoordinates(),
        ATOM_D.GetCoordinates()
      );
    }

  } // namespace biol
} // namespace bcl
