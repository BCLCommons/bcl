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
#include "biol/bcl_biol_aa_ca_cb.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_classes.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AACaCb::s_Instance
    (
      GetObjectInstances().AddInstance( new AACaCb())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AACaCb::AACaCb() :
      AABase(),
      m_CA( GetAtomTypes().CA),
      m_FirstSidechainAtom(),
      m_Atoms( util::SiPtrVector< const Atom>::Create( m_CA, m_FirstSidechainAtom))
    {
    }

    //! @brief construct AACaCb from util::ShPtr to AAData
    //! @brief SP_AA_DATA ShPTr to AAData to be copied
    AACaCb::AACaCb( const util::ShPtr< AAData> &SP_AA_DATA) :
      AABase( SP_AA_DATA),
      m_CA( GetAtomTypes().CA),
      m_FirstSidechainAtom( SP_AA_DATA->GetType()->GetFirstSidechainAtomType()),
      m_Atoms( util::SiPtrVector< const Atom>::Create( m_CA, m_FirstSidechainAtom))
    {
    }

    //! @brief constructor from AABase
    //! @param AA_BASE AABase
    AACaCb::AACaCb( const AABase &AA_BASE) :
      AABase( AA_BASE),
      m_CA( GetAtomTypes().CA),
      m_FirstSidechainAtom( GetType()->GetFirstSidechainAtomType()),
      m_Atoms( util::SiPtrVector< const Atom>::Create( m_CA, m_FirstSidechainAtom))
    {
    }

    //! @brief construct AACaCb from AABase and atoms CA, CB
    //! @param AA_BASE AABase
    //! @param ATOM_CA CA atom
    //! @param ATOM_CB CB atom
    AACaCb::AACaCb
    (
      const AABase &AA_BASE,
      const Atom &ATOM_CA,
      const Atom &ATOM_CB
    ) :
      AABase( AA_BASE),
      m_CA( ATOM_CA),
      m_FirstSidechainAtom( ATOM_CB),
      m_Atoms( util::SiPtrVector< const Atom>::Create( m_CA, m_FirstSidechainAtom))
    {
    }

    //! @brief construct AACaCb from AABASE and util::SiPtrVector< const Atom> ATOMS
    //! @param AA_BASE AABase
    //! @param ATOMS SiPtrVector of Atoms
    AACaCb::AACaCb
    (
      const AABase &AA_BASE,
      const util::SiPtrVector< const Atom> &ATOMS
    ) :
      AABase( AA_BASE),
      m_CA( Atom::FindAtom( ATOMS, GetAtomTypes().CA)->GetType().IsDefined() ? *Atom::FindAtom( ATOMS, GetAtomTypes().CA) : GetAtomTypes().CA),
      m_FirstSidechainAtom( Atom::FindAtom( ATOMS, GetType()->GetFirstSidechainAtomType())->GetType().IsDefined() ? *Atom::FindAtom( ATOMS, GetType()->GetFirstSidechainAtomType()) : GetType()->GetFirstSidechainAtomType()),
      m_Atoms( util::SiPtrVector< const Atom>::Create( m_CA, m_FirstSidechainAtom))
    {
    }

    //! @brief copy constructor
    //! @param AACACB AACaCb to be copied
    AACaCb::AACaCb( const AACaCb &AACACB) :
      AABase( AACACB),
      m_CA( AACACB.m_CA),
      m_FirstSidechainAtom( AACACB.m_FirstSidechainAtom),
      m_Atoms( util::SiPtrVector< const Atom>::Create( m_CA, m_FirstSidechainAtom))
    {
    }

    //! @brief virtual copy constructor
    AACaCb *AACaCb::Clone() const
    {
      return new AACaCb( *this);
    }

    //! @brief virtual empty constructor with AAData
    //! @param SP_AA_DATA AAData object with information for the new AA
    //! this function is designed to be used in cases where AAClass is used for production multiple AA's
    //! which do not share one common AAData
    AACaCb *AACaCb::Empty( const util::ShPtr< AAData> &SP_AA_DATA) const
    {
      return new AACaCb( SP_AA_DATA);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AACaCb::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get all atoms of specified types
    //! @param ATOM_TYPES AtomTypes of interest
    //! @return SiPtrVector of Atoms with specified type
    util::SiPtrVector< const Atom> AACaCb::GetAtoms( const storage::Set< AtomType> &ATOM_TYPES) const
    {
      // create and initialize a siptr vector
      util::SiPtrVector< const Atom> atoms;

      // allocate memory
      atoms.AllocateMemory( ATOM_TYPES.GetSize());

      // iterate over atom types
      for
      (
        storage::Set< AtomType>::const_iterator atom_itr( ATOM_TYPES.Begin()), atom_itr_end( ATOM_TYPES.End());
        atom_itr != atom_itr_end; ++atom_itr
      )
      {
        // get the atom
        const Atom &this_atom( GetAtom( *atom_itr));

        // only if the atom is defined
        if( this_atom.GetType().IsDefined())
        {
          // pushback the atom
          atoms.PushBack( util::ToSiPtr( this_atom));
        }
      }
      // return coordinates
      return atoms;
    }

    //! @brief set the specified atom
    //! @param ATOM Atom to be set
    void AACaCb::SetAtom( const Atom &ATOM)
    {
      // if ATOM is CA set CA
      if( ATOM.GetType() == GetAtomTypes().CA)
      {
        m_CA = ATOM;
      }
      // else if ATOM is first side chain atom
      else if( ATOM.GetType() == GetType()->GetFirstSidechainAtomType())
      {
        m_FirstSidechainAtom = ATOM;
      }

      // if else issue warning
      BCL_MessageStd
      (
        "This amino acid does not have the specified atom type " + ATOM.GetType().GetName()
      );
    }

    //! @brief set all atoms
    //! @param ATOMS SiPtrVector of Atoms to be set
    void AACaCb::SetAtoms( const util::SiPtrVector< const Atom> &ATOMS)
    {
      // if CA is found assign CA
      const util::SiPtr< const Atom> sp_ca( Atom::FindAtom( ATOMS, GetAtomTypes().CA));
      if( sp_ca->GetType().IsDefined())
      {
        m_CA = *sp_ca;
      }
      // if first side chain atom is found, assign first side chain atom
      const util::SiPtr< const Atom> sp_cb( Atom::FindAtom( ATOMS, GetType()->GetFirstSidechainAtomType()));
      if( sp_cb->GetType().IsDefined())
      {
        m_FirstSidechainAtom = *sp_cb;
      }
    }

    //! @brief get AAClass
    //! @return AAClass
    const AAClass &AACaCb::GetAAClass() const
    {
      return GetAAClasses().e_AACaCb;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get all atom coordinates
    //! @return all atom coordinates
    util::SiPtrVector< const linal::Vector3D> AACaCb::GetAtomCoordinates() const
    {

      // create and initialize coordinates vector
      util::SiPtrVector< const linal::Vector3D> coordinates;
      coordinates.PushBack( util::ToSiPtr( m_CA.GetCoordinates()));
      coordinates.PushBack( util::ToSiPtr( m_FirstSidechainAtom.GetCoordinates()));

      // return coordinates
      return coordinates;
    }

    //! @brief get all atom coordinates for the specified atom types
    //! @param ATOM_TYPES AtomTypes of interest
    //! @return coordinates of atoms of types specified in ATOM_TYPES
    util::SiPtrVector< const linal::Vector3D> AACaCb::GetAtomCoordinates
    (
      const storage::Set< AtomType> &ATOM_TYPES
    ) const
    {
      // create and initialize a siptr vector
      util::SiPtrVector< const linal::Vector3D> coordinates;
      coordinates.AllocateMemory( ATOM_TYPES.GetSize());

      for
      (
        storage::Set< AtomType>::const_iterator atom_itr( ATOM_TYPES.Begin()), atom_itr_end( ATOM_TYPES.End());
        atom_itr != atom_itr_end; ++atom_itr
      )
      {
        // get the atom
        const Atom &this_atom( GetAtom( *atom_itr));

        // only if the atom is defined
        if( this_atom.GetType().IsDefined())
        {
          // pushback the coordinates
          coordinates.PushBack( util::ToSiPtr( this_atom.GetCoordinates()));
        }
      }
      // return coordinates
      return coordinates;
    }

    //! @brief set the aminoacid to the ideal conformation according to the SS_TYPE with given TRANSFORMATION_MATRIX_3D
    //! @param SS_TYPE SSType this AABase derived class is in
    //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
    void AACaCb::SetToIdealConformation
    (
      const SSType &SS_TYPE,
      const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D
    )
    {
      // if ss type is not helix or strand skip
      if( !SS_TYPE->IsStructured())
      {
        return;
      }

      // set CA coordinates to corresponding to SSType
      m_CA.SetCoordinates
      (
        linal::Vector3D
        (
          GetAtomTypes().CA->GetCoordinates( SS_TYPE).GetCartesianCoordinates()
        ).Transform( TRANSFORMATION_MATRIX_3D)
      );

      // set first side chain atom coordinates to corresponding to SSType
      m_FirstSidechainAtom.SetCoordinates
      (
        linal::Vector3D
        (
          GetType()->GetFirstSidechainAtomType()->GetCoordinates( SS_TYPE).GetCartesianCoordinates()
        ).Transform( TRANSFORMATION_MATRIX_3D)
      );
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief assignment operator for AACaCb class
    //! @param AACACB_RHS AACaCb object to be copied
    //! @return this after assignment to AACACB_RHS is done
    AACaCb &AACaCb::operator =( const AACaCb &AACACB_RHS)
    {
      // assign base class and data members
      AABase::operator =( AACACB_RHS);
      m_CA = AACACB_RHS.m_CA;
      m_FirstSidechainAtom = AACACB_RHS.m_FirstSidechainAtom;

      // return
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AACaCb::Read( std::istream &ISTREAM)
    {
      // read base class
      AABase::Read( ISTREAM);

      // read member
      io::Serialize::Read( m_CA, ISTREAM);
      io::Serialize::Read( m_FirstSidechainAtom, ISTREAM);

      //return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &AACaCb::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write base class
      AABase::Write( OSTREAM, INDENT) << '\n';

      // write members
      io::Serialize::Write( m_CA, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_FirstSidechainAtom, OSTREAM, INDENT);

      //return
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
