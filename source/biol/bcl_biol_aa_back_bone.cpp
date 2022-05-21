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
#include "biol/bcl_biol_aa_back_bone.h"

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
    const util::SiPtr< const util::ObjectInterface> AABackBone::s_Instance
    (
      GetObjectInstances().AddInstance( new AABackBone())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AABackBone::AABackBone() :
      AABase(),
      m_N( GetAtomTypes().N),
      m_CA( GetAtomTypes().CA),
      m_C( GetAtomTypes().C),
      m_O( GetAtomTypes().O),
      m_FirstSidechainAtom(),
      m_Atoms( util::SiPtrVector< const Atom>::Create( m_N, m_CA, m_C, m_O, m_FirstSidechainAtom))
    {
    }

    //! @brief construct AABackBone from util::ShPtr to AAData
    //! @brief SP_AA_DATA ShPTr to AAData to be copied
    AABackBone::AABackBone( const util::ShPtr< AAData> &SP_AA_DATA) :
      AABase( SP_AA_DATA),
      m_N( GetAtomTypes().N),
      m_CA( GetAtomTypes().CA),
      m_C( GetAtomTypes().C),
      m_O( GetAtomTypes().O),
      m_FirstSidechainAtom( GetType()->GetFirstSidechainAtomType()),
      m_Atoms( util::SiPtrVector< const Atom>::Create( m_N, m_CA, m_C, m_O, m_FirstSidechainAtom))
    {
    }

    //! @brief constructor from AABase
    //! @param AA_BASE AABase
    AABackBone::AABackBone( const AABase &AA_BASE) :
      AABase( AA_BASE),
      m_N( GetAtomTypes().N),
      m_CA( GetAtomTypes().CA),
      m_C( GetAtomTypes().C),
      m_O( GetAtomTypes().O),
      m_FirstSidechainAtom( GetType()->GetFirstSidechainAtomType()),
      m_Atoms( util::SiPtrVector< const Atom>::Create( m_N, m_CA, m_C, m_O, m_FirstSidechainAtom))
    {
    }

    //! @brief construct AACaCb from AABASE and util::SiPtrVector< const Atom> ATOMS
    //! @param AA_BASE AABase
    //! @param ATOMS SiPtrVector of Atoms
    AABackBone::AABackBone( const AABase &AA_BASE, const util::SiPtrVector< const Atom> &ATOMS) :
      AABase( AA_BASE),
      m_N(  Atom::FindAtom( ATOMS, GetAtomTypes().N)->GetType().IsDefined()  ? *Atom::FindAtom( ATOMS, GetAtomTypes().N)  : GetAtomTypes().N),
      m_CA( Atom::FindAtom( ATOMS, GetAtomTypes().CA)->GetType().IsDefined() ? *Atom::FindAtom( ATOMS, GetAtomTypes().CA) : GetAtomTypes().CA),
      m_C(  Atom::FindAtom( ATOMS, GetAtomTypes().C)->GetType().IsDefined()  ? *Atom::FindAtom( ATOMS, GetAtomTypes().C)  : GetAtomTypes().C),
      m_O(  Atom::FindAtom( ATOMS, GetAtomTypes().O)->GetType().IsDefined()  ? *Atom::FindAtom( ATOMS, GetAtomTypes().O)  : GetAtomTypes().O),
      m_FirstSidechainAtom( Atom::FindAtom( ATOMS, GetType()->GetFirstSidechainAtomType())->GetType().IsDefined() ? *Atom::FindAtom( ATOMS, GetType()->GetFirstSidechainAtomType()) : GetType()->GetFirstSidechainAtomType()),
      m_Atoms( util::SiPtrVector< const Atom>::Create( m_N, m_CA, m_C, m_O, m_FirstSidechainAtom))
    {
    }

    //! @brief copy constructor
    //! @param AA_BACKBONE AABackBone to be copied
    AABackBone::AABackBone( const AABackBone &AA_BACKBONE) :
      AABase( AA_BACKBONE),
      m_N( AA_BACKBONE.m_N),
      m_CA( AA_BACKBONE.m_CA),
      m_C( AA_BACKBONE.m_C),
      m_O( AA_BACKBONE.m_O),
      m_FirstSidechainAtom( AA_BACKBONE.m_FirstSidechainAtom),
      m_Atoms( util::SiPtrVector< const Atom>::Create( m_N, m_CA, m_C, m_O, m_FirstSidechainAtom))
    {
    }

    //! @brief constructor from CaCb amino acids
    //! @param AA_CACB AACaCb to be copied
    AABackBone::AABackBone( const AACaCb &AA_CACB) :
      AABase( AA_CACB),
      m_N( AA_CACB.GetAtom( GetAtomTypes().N)),
      m_CA( AA_CACB.GetCA()),
      m_C( AA_CACB.GetAtom( GetAtomTypes().C)),
      m_O( AA_CACB.GetAtom( GetAtomTypes().O)),
      m_FirstSidechainAtom( AA_CACB.GetFirstSidechainAtom()),
      m_Atoms( util::SiPtrVector< const Atom>::Create( m_N, m_CA, m_C, m_O, m_FirstSidechainAtom))
    {
    }

    //! @brief virtual copy constructor
    AABackBone *AABackBone::Clone() const
    {
      return new AABackBone( *this);
    }

    //! @brief virtual empty constructor with AAType
    //! @param SP_AA_DATA AAData object with information for the new AA
    //! this function is designed to be used in cases where AAClass is used for production multiple AA's
    //! which do not share one common AAData
    AABackBone *AABackBone::Empty( const util::ShPtr< AAData> &SP_AA_DATA) const
    {
      return new AABackBone( SP_AA_DATA);
    }

    //! @brief destructor
    AABackBone::~AABackBone()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AABackBone::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get types of atoms
    //! @return set of AtomTypes
    const storage::Set< AtomType> &AABackBone::GetTypesOfAtoms() const
    {
      // initialize a static set of AtomType from all backbone atoms
      static storage::Set< AtomType> s_atom_type_set
      (
        storage::Set< AtomType>::Create
        (
          GetAtomTypes().N,
          GetAtomTypes().CA,
          GetAtomTypes().C,
          GetAtomTypes().O,
          GetAtomTypes().CB,
          GetAtomTypes().HA2
        )
      );

      // return
      return s_atom_type_set;
    }

    //! @brief get all atoms of specified types
    //! @param ATOM_TYPES AtomTypes of interest
    //! @return SiPtrVector of Atoms with specified type
    util::SiPtrVector< const Atom> AABackBone::GetAtoms( const storage::Set< AtomType> &ATOM_TYPES) const
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
    void AABackBone::SetAtom( const Atom &ATOM)
    {
      // switch over type of ATOM and assign it to corresponding backbone atom with same type
           if( ATOM.GetType() == GetAtomTypes().N)  { m_N  = ATOM;}
      else if( ATOM.GetType() == GetAtomTypes().CA) { m_CA = ATOM;}
      else if( ATOM.GetType() == GetAtomTypes().C)  { m_C  = ATOM;}
      else if( ATOM.GetType() == GetAtomTypes().O)  { m_O  = ATOM;}
      else if( ATOM.GetType() == GetType()->GetFirstSidechainAtomType()) { m_FirstSidechainAtom = ATOM;}
      else
      {
        // if no match was found issue warning
        BCL_MessageStd
        (
          "This amino acid does not have the specified atom type " + ATOM.GetType().GetName()
        );
      }
    }

    //! @brief set all atoms
    //! @param ATOMS SiPtrVector of Atoms to be set
    void AABackBone::SetAtoms( const util::SiPtrVector< const Atom> &ATOMS)
    {
      // search for each individual backbone atom, and if it is found assign
      const util::SiPtr< const Atom> sp_n( Atom::FindAtom( ATOMS, GetAtomTypes().N));
      if( sp_n->GetType().IsDefined())
      {
        m_N = *sp_n;
      }
      const util::SiPtr< const Atom> sp_ca( Atom::FindAtom( ATOMS, GetAtomTypes().CA));
      if( sp_ca->GetType().IsDefined())
      {
        m_CA = *sp_ca;
      }
      const util::SiPtr< const Atom> sp_c( Atom::FindAtom( ATOMS, GetAtomTypes().C));
      if( sp_c->GetType().IsDefined())
      {
        m_C = *sp_c;
      }
      const util::SiPtr< const Atom> sp_o( Atom::FindAtom( ATOMS, GetAtomTypes().O));
      if( sp_o->GetType().IsDefined())
      {
        m_O = *sp_o;
      }
      const util::SiPtr< const Atom> sp_first_sc_atom( Atom::FindAtom( ATOMS, GetType()->GetFirstSidechainAtomType()));
      if( sp_first_sc_atom->GetType().IsDefined())
      {
        m_FirstSidechainAtom = *sp_first_sc_atom;
      }
    }

    //! @brief get AAClass
    //! @return AAClass
    const AAClass &AABackBone::GetAAClass() const
    {
      return GetAAClasses().e_AABackBone;
    }

    //! @brief calculate Omega backbone angle
    //! @param PREVIOUS_CA previous CA atom
    //! @param PREVIOUS_C previous carbon atom
    //! @return omega angle
    double AABackBone::CalculateOmega( const Atom &PREVIOUS_CA, const Atom &PREVIOUS_C) const
    {
      // if the type of PREVIOUS_CA is not CA
      if( PREVIOUS_CA.GetType() != GetAtomTypes().CA)
      {
        // issue warning
        BCL_MessageCrt
        (
          "given PREVIOUS_CA is not of type CA. but: " + util::Format()( PREVIOUS_CA.GetType())
        );

        // return undefined value
        return util::GetUndefinedDouble();
      }

      // if PREVIOUS_C is not a carbon
      if( PREVIOUS_C.GetType() != GetAtomTypes().C)
      {
        // issue warning
        BCL_MessageCrt
        (
          "given PREVIOUS_C is not of type C. but: " + util::Format()( PREVIOUS_C.GetType())
        );

        // return undefined
        return util::GetUndefinedDouble();
      }

      // calculate the dihedral angle and return it
      return Dihedral( PREVIOUS_CA, PREVIOUS_C, m_N, m_CA);
    }

    //! @brief calculate phi backbone angle if hydrogen is part of this aa
    //! @return phi angle
    double AABackBone::Phi() const
    {
      return util::GetUndefined< double>();
    }

    //! @brief calculate phi backbone angle
    //! @param PREVIOUS_C previous carbon atom
    //! @return phi angle
    double AABackBone::CalculatePhi( const Atom &PREVIOUS_C) const
    {
      // if PREVIOUS_C is not of type Carbon
      if( PREVIOUS_C.GetType() != GetAtomTypes().C)
      {
        // issue warning
        BCL_MessageCrt
        (
          "given PREVIOUS_C is not of type C. but: " + util::Format()( PREVIOUS_C.GetType())
        );

        // return undefined
        return util::GetUndefinedDouble();
      }

      // calculate dihedral and return it
      return Dihedral( PREVIOUS_C, m_N, m_CA, m_C);
    }

    //! @brief calculate psi backbone angle
    //! @param FOLLOWING_N following nitrogen atom
    //! @return psi angle
    double AABackBone::CalculatePsi( const Atom &FOLLOWING_N) const
    {
      // if FOLLOWING_N is not of type Nitrogen
      if( FOLLOWING_N.GetType() != GetAtomTypes().N)
      {
        // issue warning
        BCL_MessageCrt
        (
          "given FOLLOWING_N is not of type N. but: " + util::Format()( FOLLOWING_N.GetType())
        );

        // return undefined
        return util::GetUndefinedDouble();
      }

      // calculate dihedral and return it
      return Dihedral( m_N, m_CA, m_C, FOLLOWING_N);
    }

    //! @brief calculate psi backbone angle if oxygen is part of this aa
    //! @return psi angle
    double AABackBone::Psi() const
    {
      double psi( Dihedral( m_N, m_CA, m_C, m_O));

      // rotate by 180, since the oxygen is opposite to the nitrogen usually used for psi
      if( psi > 0.0)
      {
        psi -= math::g_Pi;
      }
      else
      {
        psi += math::g_Pi;
      }

      return psi;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get all atom coordinates
    //! @return all atom coordinates
    util::SiPtrVector< const linal::Vector3D> AABackBone::GetAtomCoordinates() const
    {
      // create and initialize coordinates vector
      util::SiPtrVector< const linal::Vector3D> coordinates;
      coordinates.PushBack( util::ToSiPtr( m_N.GetCoordinates()));
      coordinates.PushBack( util::ToSiPtr( m_CA.GetCoordinates()));
      coordinates.PushBack( util::ToSiPtr( m_C.GetCoordinates()));
      coordinates.PushBack( util::ToSiPtr( m_O.GetCoordinates()));
      coordinates.PushBack( util::ToSiPtr( m_FirstSidechainAtom.GetCoordinates()));

      // return coordinates
      return coordinates;
    }

    //! @brief get all atom coordinates for the specified atom types
    //! @param ATOM_TYPES AtomTypes of interest
    //! @return coordinates of atoms of types specified in ATOM_TYPES
    util::SiPtrVector< const linal::Vector3D> AABackBone::GetAtomCoordinates
    (
      const storage::Set< AtomType> &ATOM_TYPES
    ) const
    {
      // create and initialize a siptr vector
      util::SiPtrVector< const linal::Vector3D> coordinates;

      // allocate memory
      coordinates.AllocateMemory( ATOM_TYPES.GetSize());

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
          // pushback the coordinates
          coordinates.PushBack( util::ToSiPtr( this_atom.GetCoordinates()));
        }
      }
      // return coordinates
      return coordinates;
    }

    //! @brief translate the object along a given TRANSLATION vector
    //! @param TRANSLATION translation vector to be applied
    void AABackBone::Translate( const linal::Vector3D &TRANSLATION)
    {
      // transform all atoms
      m_N.Translate( TRANSLATION);
      m_CA.Translate( TRANSLATION);
      m_C.Translate( TRANSLATION);
      m_O.Translate( TRANSLATION);
      m_FirstSidechainAtom.Translate( TRANSLATION);
    }

    //! @brief transform the object by a given TRANSFORMATION_MATRIX_3D
    //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
    void AABackBone::Transform( const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D)
    {
      // transform all atoms
      m_N.Transform( TRANSFORMATION_MATRIX_3D);
      m_CA.Transform( TRANSFORMATION_MATRIX_3D);
      m_C.Transform( TRANSFORMATION_MATRIX_3D);
      m_O.Transform( TRANSFORMATION_MATRIX_3D);
      m_FirstSidechainAtom.Transform( TRANSFORMATION_MATRIX_3D);
    }

    //! @brief rotate the object by a given ROTATION_MATRIX_3D
    //! @param ROTATION_MATRIX_3D RotationMatrix3D to be applied
    void AABackBone::Rotate( const math::RotationMatrix3D &ROTATION_MATRIX_3D)
    {
      // transform all atoms
      m_N.Rotate( ROTATION_MATRIX_3D);
      m_CA.Rotate( ROTATION_MATRIX_3D);
      m_C.Rotate( ROTATION_MATRIX_3D);
      m_O.Rotate( ROTATION_MATRIX_3D);
      m_FirstSidechainAtom.Rotate( ROTATION_MATRIX_3D);
    }

    //! @brief returns the geometric center of the object
    //! @return geometric center of the object
    linal::Vector3D AABackBone::GetCenter() const
    {
      return coord::CenterOfMass( GetAtomCoordinates());
    }

    //! @brief set the aminoacid to the ideal conformation according to the SS_TYPE with given TRANSFORMATION_MATRIX_3D
    //! @param SS_TYPE SSType this AABase derived class is in
    //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
    void AABackBone::SetToIdealConformation
    (
      const SSType &SS_TYPE,
      const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D
    )
    {
      // if SS_TYPE is not helix nor strand, then skip
      if( !SS_TYPE->IsStructured())
      {
        return;
      }

      // set each atom to corresponding SS_TYPE specific coordinates
      m_N.SetCoordinates
      (
        linal::Vector3D
        (
          GetAtomTypes().N->GetCoordinates( SS_TYPE).GetCartesianCoordinates()
        ).Transform( TRANSFORMATION_MATRIX_3D)
      );
      m_CA.SetCoordinates
      (
        linal::Vector3D
        (
          GetAtomTypes().CA->GetCoordinates( SS_TYPE).GetCartesianCoordinates()
        ).Transform( TRANSFORMATION_MATRIX_3D)
      );
      m_C.SetCoordinates
      (
        linal::Vector3D
        (
          GetAtomTypes().C->GetCoordinates( SS_TYPE).GetCartesianCoordinates()
        ).Transform( TRANSFORMATION_MATRIX_3D)
      );
      m_O.SetCoordinates
      (
        linal::Vector3D
        (
          GetAtomTypes().O->GetCoordinates( SS_TYPE).GetCartesianCoordinates()
        ).Transform( TRANSFORMATION_MATRIX_3D)
      );
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

    //! @brief assignment operator for AABackBone class
    //! @param AA_BACKBONE_RHS AABackBone object to be copied
    //! @return this after assignment to AA_BACKBONE_RHS is done
    AABackBone &AABackBone::operator =( const AABackBone &AA_BACKBONE_RHS)
    {
      // assign base class AABase and data members
      AABase::operator =( AA_BACKBONE_RHS);
      m_N = AA_BACKBONE_RHS.m_N;
      m_CA = AA_BACKBONE_RHS.m_CA;
      m_C = AA_BACKBONE_RHS.m_C;
      m_O = AA_BACKBONE_RHS.m_O;
      m_FirstSidechainAtom = AA_BACKBONE_RHS.m_FirstSidechainAtom;

      // return
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AABackBone::Read( std::istream &ISTREAM)
    {
      // read base class
      AABase::Read( ISTREAM);

      // read members
      io::Serialize::Read( m_N, ISTREAM);
      io::Serialize::Read( m_CA, ISTREAM);
      io::Serialize::Read( m_C, ISTREAM);
      io::Serialize::Read( m_O, ISTREAM);
      io::Serialize::Read( m_FirstSidechainAtom, ISTREAM);

      //return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &AABackBone::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write base class
      AABase::Write( OSTREAM, INDENT) << '\n';

      // write members
      io::Serialize::Write( m_N, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_CA, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_C, OSTREAM, INDENT)  << '\n';
      io::Serialize::Write( m_O, OSTREAM, INDENT)  << '\n';
      io::Serialize::Write( m_FirstSidechainAtom, OSTREAM, INDENT);

      //return
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
