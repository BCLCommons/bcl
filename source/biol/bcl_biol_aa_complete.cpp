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
#include "biol/bcl_biol_aa_complete.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_back_bone.h"
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
    const util::SiPtr< const util::ObjectInterface> AAComplete::s_Instance
    (
      GetObjectInstances().AddInstance( new AAComplete())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AAComplete::AAComplete() :
      AABase(),
      m_Atoms( DefaultBackboneAtoms()),
      m_AtomList( m_Atoms.Begin(), m_Atoms.End())
    {
    }

    //! @brief construct AABase from util::ShPtr to AAData
    //! @brief SP_AA_DATA ShPTr to AAData to be copied
    AAComplete::AAComplete( const util::ShPtr< AAData> &SP_AA_DATA) :
      AABase( SP_AA_DATA),
      m_Atoms( DefaultBackboneAtomsWithFirstSideChainAtom( SP_AA_DATA->GetType())),
      m_AtomList( m_Atoms.Begin(), m_Atoms.End())
    {
    }

    //! @brief constructor from AABase
    //! @param AA_BASE AABase
    AAComplete::AAComplete( const AABase &AA_BASE) :
      AABase( AA_BASE),
      m_Atoms( DefaultBackboneAtomsWithFirstSideChainAtom( AA_BASE.GetType())),
      m_AtomList( m_Atoms.Begin(), m_Atoms.End())
    {
      SetAtoms( AA_BASE.GetAtoms());
    }

    //! @brief construct AACaCb from AABASE and util::SiPtrVector< const Atom> ATOMS
    //! @param AA_BASE AABase
    //! @param ATOMS SiPtrVector of Atoms
    AAComplete::AAComplete( const AABase &AA_BASE, const util::SiPtrVector< const Atom> &ATOMS) :
      AABase( AA_BASE),
      m_Atoms( DefaultBackboneAtomsWithFirstSideChainAtom( AA_BASE.GetType())),
      m_AtomList( m_Atoms.Begin(), m_Atoms.End())
    {
      SetAtoms( ATOMS);
    }

    //! @brief copy constructor
    //! @param AA_COMPLETE AAComplete to be copied
    AAComplete::AAComplete( const AAComplete &AA_COMPLETE) :
      AABase( AA_COMPLETE),
      m_Atoms( DefaultBackboneAtomsWithFirstSideChainAtom( AA_COMPLETE.GetType())),
      m_AtomList( m_Atoms.Begin(), m_Atoms.End())
    {
      SetAtoms( AA_COMPLETE.GetAtoms());
    }

    //! @brief copy constructor
    //! @param AA_BACKBONE AABackBone to be copied
    AAComplete::AAComplete( const AABackBone &AA_BACKBONE) :
      AABase( AA_BACKBONE),
      m_Atoms( DefaultBackboneAtomsWithFirstSideChainAtom( AA_BACKBONE.GetType())),
      m_AtomList( m_Atoms.Begin(), m_Atoms.End())
    {
      SetAtoms( AA_BACKBONE.GetAtoms());
    }

    //! @brief constructor from CaCb amino acids
    //! @param AA_CACB AACaCb to be copied
    AAComplete::AAComplete( const AACaCb &AA_CACB) :
      AABase( AA_CACB),
      m_Atoms( DefaultBackboneAtomsWithFirstSideChainAtom( AA_CACB.GetType())),
      m_AtomList( m_Atoms.Begin(), m_Atoms.End())
    {
      SetAtoms( AA_CACB.GetAtoms());
    }

    //! @brief virtual copy constructor
    AAComplete *AAComplete::Clone() const
    {
      return new AAComplete( *this);
    }

    //! @brief virtual empty constructor with AAData
    //! @param SP_AA_DATA AAData object with information for the new AA
    //! this function is designed to be used in cases where AAClass is used for production multiple AA's
    //! which do not share one common AAData
    AAComplete *AAComplete::Empty( const util::ShPtr< AAData> &SP_AA_DATA) const
    {
      return new AAComplete( SP_AA_DATA);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AAComplete::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get types of atoms
    //! @return set of AtomTypes
    const storage::Set< AtomType> &AAComplete::GetTypesOfAtoms() const
    {
      // return the atom types that belong to that complete amino acid
      return AABase::GetType()->GetAllowedAtomTypes();
    }

    //! @brief get all atoms of specified types
    //! @param ATOM_TYPES AtomTypes of interest
    //! @return SiPtrVector of Atoms with specified type
    util::SiPtrVector< const Atom> AAComplete::GetAtoms( const storage::Set< AtomType> &ATOM_TYPES) const
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

    //! @brief get CB atom
    //! @return CB atom
    const Atom &AAComplete::GetFirstSidechainAtom() const
    {
      // return cb, it amino acid type has cb, HA2 otherwise (for GLY type amino acids)
      return GetAtom( GetType()->GetFirstSidechainAtomType());
    }

    //! @brief set the specified atom
    //! @param ATOM Atom to be set
    void AAComplete::SetAtom( const Atom &ATOM)
    {
      const AtomType &given_atom_type( ATOM.GetType());

      // iterate over all atoms to find the one with the correct type
      for
      (
        storage::List< Atom>::iterator itr( m_Atoms.Begin()), itr_end( m_Atoms.End());
        itr != itr_end;
        ++itr
      )
      {
        // check if atom agrees
        if( itr->GetType() == given_atom_type)
        {
          *itr = ATOM;
          return;
        }
      }

      // if there was non found, check if the type does match the aa types atom types or is a terminal atom type and insert
      if( GetType()->DoesContainAtomType( given_atom_type) || GetAtomTypes().GetTerminalExtraAtomTypes().Has( given_atom_type))
      {
        m_Atoms.PushBack( ATOM);
      }
      // the atom cannot go into the amino acid
      else
      {
        // if no match was found issue warning
        BCL_MessageVrb
        (
          "This amino acid does not have the atom of the specified atom type " + ATOM.GetType().GetName() +
          "\nthe type of the amino acid type " + GetType().GetName() + " does not allow it"
        );

        // end
        return;
      }

      // update the siptr vector
      m_AtomList = util::SiPtrVector< const Atom>( m_Atoms.Begin(), m_Atoms.End());
    }

    //! @brief set all atoms
    //! @param ATOMS SiPtrVector of Atoms to be set
    void AAComplete::SetAtoms( const util::SiPtrVector< const Atom> &ATOMS)
    {
      // call set atom for each given atom
      for
      (
        util::SiPtrVector< const Atom>::const_iterator itr( ATOMS.Begin()), itr_end( ATOMS.End());
        itr != itr_end;
        ++itr
      )
      {
        SetAtom( **itr);
      }
    }

    //! @brief get AAClass
    //! @return AAClass
    const AAClass &AAComplete::GetAAClass() const
    {
      return GetAAClasses().e_AAComplete;
    }

    //! @brief calculate Omega backbone angle
    //! @param PREVIOUS_CA previous CA atom
    //! @param PREVIOUS_C previous carbon atom
    //! @return omega angle
    double AAComplete::CalculateOmega( const Atom &PREVIOUS_CA, const Atom &PREVIOUS_C) const
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
      return Dihedral( PREVIOUS_CA, PREVIOUS_C, GetAtom( GetAtomTypes().N), GetAtom( GetAtomTypes().CA));
    }

    //! @brief calculate phi backbone angle if hydrogen is part of this aa
    //! @return phi angle
    double AAComplete::Phi() const
    {
      // calculate dihedral
      double phi( Dihedral( GetAtom( GetAtomTypes().H), GetAtom( GetAtomTypes().N), GetAtom( GetAtomTypes().CA), GetAtom( GetAtomTypes().C)));
      if( phi > 0.0)
      {
        phi -= math::g_Pi;
      }
      else
      {
        phi += math::g_Pi;
      }

      // return
      return phi;
    }

    //! @brief calculate phi backbone angle
    //! @param PREVIOUS_C previous carbon atom
    //! @return phi angle
    double AAComplete::CalculatePhi( const Atom &PREVIOUS_C) const
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
      return Dihedral( PREVIOUS_C, GetAtom( GetAtomTypes().N), GetAtom( GetAtomTypes().CA), GetAtom( GetAtomTypes().C));
    }

    //! @brief calculate psi backbone angle if oxygen is part of this aa
    //! @return psi angle
    double AAComplete::Psi() const
    {
      // calculate dihedral
      double psi( Dihedral( GetAtom( GetAtomTypes().N), GetAtom( GetAtomTypes().CA), GetAtom( GetAtomTypes().C), GetAtom( GetAtomTypes().O)));
      if( psi > 0.0)
      {
        psi -= math::g_Pi;
      }
      else
      {
        psi += math::g_Pi;
      }

      // return
      return psi;
    }

    //! @brief calculate psi backbone angle
    //! @param FOLLOWING_N following nitrogen atom
    //! @return psi angle
    double AAComplete::CalculatePsi( const Atom &FOLLOWING_N) const
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
      return Dihedral( GetAtom( GetAtomTypes().N), GetAtom( GetAtomTypes().CA), GetAtom( GetAtomTypes().C), FOLLOWING_N);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get all atom coordinates
    //! @return all atom coordinates
    util::SiPtrVector< const linal::Vector3D> AAComplete::GetAtomCoordinates() const
    {
      // create and initialize coordinates vector
      util::SiPtrVector< const linal::Vector3D> coordinates;

      // iterate over all atoms and collect their coordinates
      for
      (
        storage::List< Atom>::const_iterator itr( m_Atoms.Begin()), itr_end( m_Atoms.End());
        itr != itr_end;
        ++itr
      )
      {
        coordinates.PushBack( util::ToSiPtr( itr->GetCoordinates()));
      }

      // return coordinates
      return coordinates;
    }

    //! @brief get all atom coordinates for the specified atom types
    //! @param ATOM_TYPES AtomTypes of interest
    //! @return coordinates of atoms of types specified in ATOM_TYPES
    util::SiPtrVector< const linal::Vector3D> AAComplete::GetAtomCoordinates
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
        atom_itr != atom_itr_end;
        ++atom_itr
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
    void AAComplete::Translate( const linal::Vector3D &TRANSLATION)
    {
      // transform all atoms
      for
      (
        storage::List< Atom>::iterator itr( m_Atoms.Begin()), itr_end( m_Atoms.End());
        itr != itr_end;
        ++itr
      )
      {
        itr->Translate( TRANSLATION);
      }
    }

    //! @brief transform the object by a given TRANSFORMATION_MATRIX_3D
    //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
    void AAComplete::Transform( const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D)
    {
      // transform all atoms
      for
      (
        storage::List< Atom>::iterator itr( m_Atoms.Begin()), itr_end( m_Atoms.End());
        itr != itr_end;
        ++itr
      )
      {
        itr->Transform( TRANSFORMATION_MATRIX_3D);
      }
    }

    //! @brief rotate the object by a given ROTATION_MATRIX_3D
    //! @param ROTATION_MATRIX_3D RotationMatrix3D to be applied
    void AAComplete::Rotate( const math::RotationMatrix3D &ROTATION_MATRIX_3D)
    {
      // transform all atoms
      for
      (
        storage::List< Atom>::iterator itr( m_Atoms.Begin()), itr_end( m_Atoms.End());
        itr != itr_end;
        ++itr
      )
      {
        itr->Rotate( ROTATION_MATRIX_3D);
      }
    }

    //! @brief returns the geometric center of the object
    //! @return geometric center of the object
    linal::Vector3D AAComplete::GetCenter() const
    {
      return coord::CenterOfMass( GetAtomCoordinates());
    }

    //! @brief set the aminoacid to the ideal conformation according to the SS_TYPE with given TRANSFORMATION_MATRIX_3D
    //! @param SS_TYPE SSType this AABase derived class is in
    //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
    void AAComplete::SetToIdealConformation
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

      // this function will be called whenever an SSE is constructed from an aa sequence when the main geometry is
      // calculated, so this is probably not a terribly useful message
      BCL_MessageDbg
      (
        "SetToIdealConformation will strip off side chain atoms except CB for " + GetStaticClassName( *this)
      );

      storage::List< Atom>::iterator atom_itr( m_Atoms.Begin()), atom_itr_end( m_Atoms.End());

      // delete non backbone atoms
      while( atom_itr != atom_itr_end)
      {
        if( atom_itr->GetType()->IsBackBone() || atom_itr->GetType() == GetAtomTypes().CB)
        {
          // set the backbone to ideal conformation
          atom_itr->SetCoordinates
          (
            linal::Vector3D
            (
              atom_itr->GetType()->GetCoordinates( SS_TYPE).GetCartesianCoordinates()
            ).Transform( TRANSFORMATION_MATRIX_3D)
          );

          // got to next atom
          ++atom_itr;
        }
        else
        {
          atom_itr = m_Atoms.Remove( atom_itr);
        }
      }

      // reinitialize the atom pointer vector
      m_AtomList = util::SiPtrVector< const Atom>( m_Atoms.Begin(), m_Atoms.End());
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief assignment operator for AAComplete class
    //! @param AA_COMPLETE AAComplete object to be copied
    //! @return this after assignment to m_Atoms
    AAComplete &AAComplete::operator =( const AAComplete &AA_COMPLETE)
    {
      // assign base class AABase and data members
      AABase::operator =( AA_COMPLETE);

      // assign the atoms
      m_Atoms = AA_COMPLETE.m_Atoms;

      // return
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AAComplete::Read( std::istream &ISTREAM)
    {
      // read base class
      AABase::Read( ISTREAM);

      // read member
      io::Serialize::Read( m_Atoms, ISTREAM);

      // reinitialize the siptrvector
      m_AtomList = util::SiPtrVector< const Atom>( m_Atoms.Begin(), m_Atoms.End());

      //return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @return ostream which was written to
    std::ostream &AAComplete::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write base class
      AABase::Write( OSTREAM, INDENT) << '\n';

      // write members
      io::Serialize::Write( m_Atoms, OSTREAM, INDENT);

      //return
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief create a list of default backbone atoms, except CB (firstsidechainatom)
    //! @return List with  atoms, CA, C, N, O
    storage::List< Atom> AAComplete::DefaultBackboneAtoms()
    {
      // list of backbone atoms
      storage::List< Atom> backbone_atoms;

      // insert an atom for each backbone atom type
      for
      (
        storage::Set< AtomType>::const_iterator
          itr( GetAtomTypes().GetBackBoneAtomTypes().Begin()), itr_end( GetAtomTypes().GetBackBoneAtomTypes().End());
        itr != itr_end;
        ++itr
      )
      {
        backbone_atoms.PushBack( Atom( *itr));
      }

      // end
      return backbone_atoms;
    }

    //! @brief create a list of default backbone atoms
    //! @return List with  atoms, CA, C, N, O
    storage::List< Atom> AAComplete::DefaultBackboneAtomsWithFirstSideChainAtom( const AAType &AA_TYPE)
    {
      // get the backbone atoms
      storage::List< Atom> backbone_atoms( DefaultBackboneAtoms());

      // if type is Glycine
      if( AA_TYPE == GetAATypes().GLY)
      {
        backbone_atoms.PushBack( Atom( GetAtomTypes().HA2));
      }
      else
      {
        backbone_atoms.PushBack( Atom( GetAtomTypes().CB));
      }

      // end
      return backbone_atoms;
    }

  } // namespace biol
} // namespace bcl
