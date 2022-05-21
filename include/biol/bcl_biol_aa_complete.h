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

#ifndef BCL_BIOL_AA_COMPLETE_H_
#define BCL_BIOL_AA_COMPLETE_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_biol_aa_base.h"
#include "bcl_biol_atom.h"
#include "storage/bcl_storage_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AAComplete
    //! @brief class for holding all backbone and side chain atoms for an amino acid
    //! @details Unlike AABackbone or AACaCb, this class represents each amino acid with all the atoms that AAType has.
    //! Since these types differ from AAType to AAType, the internal storage is managed with a list of instead of
    //! individual atom data members.
    //!
    //! @see @link example_biol_aa_complete.cpp @endlink
    //! @author woetzen
    //! @date Sep 4, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AAComplete :
      public AABase
    {

    private:

    //////////
    // data //
    //////////

      storage::List< Atom>           m_Atoms;    //!< atoms for that amino acid
      util::SiPtrVector< const Atom> m_AtomList; //!< siptr vector for easier access

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AAComplete();

      //! @brief construct AABase from util::ShPtr to AAData
      //! @brief SP_AA_DATA ShPTr to AAData to be copied
      AAComplete( const util::ShPtr< AAData> &SP_AA_DATA);

      //! @brief constructor from AABase
      //! @param AA_BASE AABase
      AAComplete( const AABase &AA_BASE);

      //! @brief construct AACaCb from AABASE and util::SiPtrVector< const Atom> ATOMS
      //! @param AA_BASE AABase
      //! @param ATOMS SiPtrVector of Atoms
      AAComplete( const AABase &AA_BASE, const util::SiPtrVector< const Atom> &ATOMS);

      //! @brief copy constructor
      //! @param AA_COMPLETE AAComplete to be copied
      AAComplete( const AAComplete &AA_COMPLETE);

      //! @brief copy constructor
      //! @param AA_BACKBONE AABackBone to be copied
      AAComplete( const AABackBone &AA_BACKBONE);

      //! @brief constructor from CaCb amino acids
      //! @param AA_CACB AACaCb to be copied
      AAComplete( const AACaCb &AA_CACB);

      //! @brief Clone function
      //! @return pointer to new AAComplete
      AAComplete *Clone() const;

      //! @brief virtual empty constructor with AAData
      //! @param SP_AA_DATA AAData object with information for the new AA
      //! this function is designed to be used in cases where AAClass is used for production multiple AA's
      //! which do not share one common AAData
      AAComplete *Empty( const util::ShPtr< AAData> &SP_AA_DATA) const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get number of atoms
      //! @return number of atoms
      size_t GetNumberOfAtoms() const
      {
        return m_AtomList.GetSize();
      }

      //! @brief get types of atoms
      //! @return set of AtomTypes
      const storage::Set< AtomType> &GetTypesOfAtoms() const;

      //! @brief get all atoms
      //! @return SiPtrVector of Atoms
      const util::SiPtrVector< const Atom> &GetAtoms() const
      {
        // return
        return m_AtomList;
      }

      //! @brief get all atoms of specified types
      //! @param ATOM_TYPES AtomTypes of interest
      //! @return SiPtrVector of Atoms with specified type
      util::SiPtrVector< const Atom> GetAtoms( const storage::Set< AtomType> &ATOM_TYPES) const;

      //! @brief get the specified atom
      //! @brief ATOM_TYPE AtomType of interest
      //! @return atom with type ATOM_TYPE
      const Atom &GetAtom( const AtomType &ATOM_TYPE) const
      {
        return *Atom::FindAtom( m_AtomList, ATOM_TYPE);
      }

      //! @brief get CA atom
      //! @return CA atom
      const Atom &GetCA() const
      {
        return GetAtom( GetAtomTypes().CA);
      }

      //! @brief get CB atom
      //! @return CB atom
      const Atom &GetFirstSidechainAtom() const;

      //! @brief set the specified atom
      //! @param ATOM Atom to be set
      void SetAtom( const Atom &ATOM);

      //! @brief set all atoms
      //! @param ATOMS SiPtrVector of Atoms to be set
      void SetAtoms( const util::SiPtrVector< const Atom> &ATOMS);

      //! @brief get AAClass
      //! @return AAClass
      const AAClass &GetAAClass() const;

      //! @brief calculate Omega backbone angle
      //! @param PREVIOUS_CA previous CA atom
      //! @param PREVIOUS_C previous carbon atom
      //! @return omega angle
      double CalculateOmega( const Atom &PREVIOUS_CA, const Atom &PREVIOUS_C) const;

      //! @brief calculate phi backbone angle if hydrogen is part of this aa
      //! @return phi angle
      double Phi() const;

      //! @brief calculate phi backbone angle
      //! @param PREVIOUS_C previous carbon atom
      //! @return phi angle
      double CalculatePhi( const Atom &PREVIOUS_C) const;

      //! @brief calculate psi backbone angle if oxygen is part of this aa
      //! @return psi angle
      double Psi() const;

      //! @brief calculate psi backbone angle
      //! @param FOLLOWING_N following nitrogen atom
      //! @return psi angle
      double CalculatePsi( const Atom &FOLLOWING_N) const;

    ////////////////
    // operations //
    ////////////////

      //! @brief get all atom coordinates
      //! @return all atom coordinates
      util::SiPtrVector< const linal::Vector3D> GetAtomCoordinates() const;

      //! @brief get all atom coordinates for the specified atom types
      //! @param ATOM_TYPES AtomTypes of interest
      //! @return coordinates of atoms of types specified in ATOM_TYPES
      util::SiPtrVector< const linal::Vector3D> GetAtomCoordinates( const storage::Set< AtomType> &ATOM_TYPES) const;

      //! @brief translate the object along a given TRANSLATION vector
      //! @param TRANSLATION translation vector to be applied
      void Translate( const linal::Vector3D &TRANSLATION);

      //! @brief transform the object by a given TRANSFORMATION_MATRIX_3D
      //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
      void Transform( const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D);

      //! @brief rotate the object by a given ROTATION_MATRIX_3D
      //! @param ROTATION_MATRIX_3D RotationMatrix3D to be applied
      void Rotate( const math::RotationMatrix3D &ROTATION_MATRIX_3D);

      //! @brief returns the geometric center of the object
      //! @return geometric center of the object
      linal::Vector3D GetCenter() const;

      //! @brief set the aminoacid to the ideal conformation according to the SS_TYPE with given TRANSFORMATION_MATRIX_3D
      //! @param SS_TYPE SSType this AABase derived class is in
      //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
      void SetToIdealConformation
      (
        const SSType &SS_TYPE,
        const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D
      );

    ///////////////
    // operators //
    ///////////////

      //! @brief assignment operator for AAComplete class
      //! @param AA_COMPLETE AAComplete object to be copied
      //! @return this after assignment to AA_COMPLETE_RHS is done
      AAComplete &operator =( const AAComplete &AA_COMPLETE);

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

    private:

      //! @brief create a list of default backbone atoms
      //! @return List with  atoms, CA, C, N, O
      static storage::List< Atom> DefaultBackboneAtoms();

      //! @brief create a list of default backbone atoms
      //! @return List with  atoms, CA, C, N, O
      static storage::List< Atom> DefaultBackboneAtomsWithFirstSideChainAtom( const AAType &AA_TYPE);

    }; // class AAComplete

  } // namespace biol
} // namespace bcl

#endif // BCL_BIOL_AA_COMPLETE_H_
