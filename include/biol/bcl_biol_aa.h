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

#ifndef BCL_BIOL_AA_H_
#define BCL_BIOL_AA_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_biol_aa_base.h"
#include "bcl_biol_atom.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AA
    //! @brief This is a class for the simplest amino acid without atoms derived from the AABase interface
    //! @details This AA class only has the information that comes from AAData and is suitable for use in applications
    //! where atom coordinates are not required such as alignment.
    //!
    //! @see @link example_biol_aa.cpp @endlink
    //! @author meilerj, staritrd, woetzen, karakam
    //! @date 21.11.2005
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AA :
      public AABase
    {

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
      AA();

      //! @brief construct AAData from util::ShPtr to AA
      //! @param SP_AA_DATA ShPtr to AAData to be copied
      AA( const util::ShPtr< AAData> &SP_AA_DATA);

      //! @brief constructor from AABase
      //! @param AA_BASE AABase
      AA( const AABase &AA_BASE);

      //! @brief copy constructor - makes just a soft copy of m_Data
      //! @param AMINO_ACID AA to be copied
      AA( const AA &AMINO_ACID);

      //! @brief virtual copy constructor
      AA *Clone() const;

      //! @brief virtual empty constructor with AAData
      //! @param SP_AA_DATA AAData object with information for the new AA
      //! this function is designed to be used in cases where AAClass is used for production multiple AA's
      //! which do not share one common AAData
      AA *Empty( const util::ShPtr< AAData> &SP_AA_DATA) const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get number of atoms
      //! @return number of atoms
      size_t GetNumberOfAtoms() const
      {
        return size_t( 0);
      }

      //! @brief get types of atoms
      //! @return set of AtomTypes
      const storage::Set< AtomType> &GetTypesOfAtoms() const;

      //! @brief get all atoms
      //! @return SiPtrVector of Atoms
      const util::SiPtrVector< const Atom> &GetAtoms() const;

      //! @brief get all atoms of specified types
      //! @param ATOM_TYPES AtomTypes of interest
      //! @return SiPtrVector of Atoms with specified type
      util::SiPtrVector< const Atom> GetAtoms( const storage::Set< AtomType> &ATOM_TYPES) const
      {
        return GetAtoms();
      }

      //! @brief get the specified atom
      //! @brief ATOM_TYPE AtomType of interest
      //! @return atom with type ATOM_TYPE
      const Atom &GetAtom( const AtomType &ATOM_TYPE) const;

      //! @brief get CA atom
      //! @return CA atom
      const Atom &GetCA() const
      {
        // initialize undefined atom
        static const Atom s_undefined_atom;

        // return
        return s_undefined_atom;
      }

      //! @brief get CB atom
      //! @return CB atom
      const Atom &GetFirstSidechainAtom() const
      {
        // initialize undefined atom
        static const Atom s_undefined_atom;

        // return
        return s_undefined_atom;
      }

      //! @brief set the specified atom
      //! @param ATOM Atom to be set
      void SetAtom( const Atom &ATOM)
      {
      }

      //! @brief set all atoms
      //! @param ATOMS SiPtrVector of Atoms to be set
      void SetAtoms( const util::SiPtrVector< const Atom> &ATOMS)
      {
      }

      //! @brief get AAClass
      //! @return AAClass
      const AAClass &GetAAClass() const;

      //! @brief calculate Omega backbone angle
      //! @param PREVIOUS_CA previous CA atom
      //! @param PREVIOUS_C previous carbon atom
      //! @return omega angle
      double CalculateOmega( const Atom &PREVIOUS_CA, const Atom &PREVIOUS_C) const
      {
        // undefined omega since does not have atoms
        return util::GetUndefinedDouble();
      }

      //! @brief calculate phi backbone angle if hydrogen is part of this aa
      //! @return phi angle
      double Phi() const
      {
        return util::GetUndefined< double>();
      }

      //! @brief calculate phi backbone angle
      //! @param PREVIOUS_C previous carbon atom
      //! @return phi angle
      double CalculatePhi( const Atom &PREVIOUS_C) const
      {
        // undefined phi since does not have atoms
        return util::GetUndefinedDouble();
      }

      //! @brief calculate psi backbone angle if oxygen is part of this aa
      //! @return psi angle
      double Psi() const
      {
        return util::GetUndefined< double>();
      }

      //! @brief calculate psi backbone angle
      //! @param FOLLOWING_N following nitrogen atom
      //! @return psi angle
      double CalculatePsi( const Atom &FOLLOWING_N) const
      {
        // undefined psi since does not have atoms
        return util::GetUndefinedDouble();
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief get all atom coordinates
      //! @return all atom coordinates
      util::SiPtrVector< const linal::Vector3D> GetAtomCoordinates() const
      {
        // initialize vector of undefined coordinates
        static const util::SiPtrVector< const linal::Vector3D> s_undefined_coordinates;

        // retrun
        return s_undefined_coordinates;
      }

      //! @brief get all atom coordinates for the specified atom types
      //! @param ATOM_TYPES AtomTypes of interest
      //! @return coordinates of atoms of types specified in ATOM_TYPES
      util::SiPtrVector< const linal::Vector3D> GetAtomCoordinates( const storage::Set< AtomType> &ATOM_TYPES) const
      {
        // initalize vector of undefined coordinates
        static const util::SiPtrVector< const linal::Vector3D> s_undefined_coordinates;

        // return
        return s_undefined_coordinates;
      }

      //! @brief translate the object along a given TRANSLATION vector
      //! @param TRANSLATION translation vector to be applied
      void Translate( const linal::Vector3D &TRANSLATION)
      {
      }

      //! @brief transform the object by a given TRANSFORMATION_MATRIX_3D
      //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
      void Transform( const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D)
      {
      }

      //! @brief rotate the object by a given ROTATION_MATRIX_3D
      //! @param ROTATION_MATRIX_3D RotationMatrix3D to be applied
      void Rotate( const math::RotationMatrix3D &ROTATION_MATRIX_3D)
      {
      }

      //! @brief returns the geometric center of the object
      //! @return geometric center of the object
      linal::Vector3D GetCenter() const
      {
        // initialize center of undefined coordinates
        static const linal::Vector3D s_undefined_center( util::GetUndefined< double>());

        // return
        return s_undefined_center;
      }

      //! @brief set the aminoacid to the ideal conformation according to the SSTYPE with given TRANSFORMATION_MATRIX_3D
      //! @param SS_TYPE SSType this AABase derived class is in
      //! @param TRANSFORMATION_MATRIX3D TransformationMatrix3D to be applied
      void SetToIdealConformation
      (
        const SSType &SS_TYPE,
        const math::TransformationMatrix3D &TRANSFORMATION_MATRIX3D
      )
      {
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief assignment operator for AA class
      //! @param AA_RHS AA object to be copied
      //! @return this after assignment to AA_RHS is done
      AA &operator =( const AA &AA_RHS);

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

    }; // class AA

  } // namespace biol
} // namespace bcl

#endif //BCL_BIOL_AA_H_

