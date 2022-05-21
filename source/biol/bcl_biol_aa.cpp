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
#include "biol/bcl_biol_aa.h"

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
    const util::SiPtr< const util::ObjectInterface> AA::s_Instance
    (
      GetObjectInstances().AddInstance( new AA())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AA::AA() :
      AABase()
    {
    }

    //! @brief construct AAData from util::ShPtr to AA
    //! @param SP_AA_DATA ShPtr to AAData to be copied
    AA::AA( const util::ShPtr< AAData> &SP_AA_DATA) :
      AABase( SP_AA_DATA)
    {
    }

    //! @brief constructor from AABase
    //! @param AA_BASE AABase
    AA::AA( const AABase &AA_BASE) :
      AABase( AA_BASE)
    {
    }

    //! @brief copy constructor - makes just a soft copy of m_Data
    //! @param AMINO_ACID AA to be copied
    AA::AA( const AA &AMINO_ACID) :
      AABase( AMINO_ACID)
    {
    }

    //! @brief virtual copy constructor
    AA *AA::Clone() const
    {
      return new AA( *this);
    }

    //! @brief virtual empty constructor with AAData
    //! @param SP_AA_DATA AAData object with information for the new AA
    //! this function is designed to be used in cases where AAClass is used for production multiple AA's
    //! which do not share one common AAData
    AA *AA::Empty( const util::ShPtr< AAData> &SP_AA_DATA) const
    {
      return new AA( SP_AA_DATA);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AA::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get types of atoms
    //! @return set of AtomTypes
    const storage::Set< AtomType> &AA::GetTypesOfAtoms() const
    {
      // initialize set of undefined atom types
      static const storage::Set< AtomType> s_atom_type_set;

      // return
      return s_atom_type_set;
    }

    //! @brief get all atoms
    //! @return SiPtrVector of Atoms
    const util::SiPtrVector< const Atom> &AA::GetAtoms() const
    {
      static const util::SiPtrVector< const Atom> s_undefined_atoms;
      return s_undefined_atoms;
    }

    //! @brief get the specified atom
    //! @brief ATOM_TYPE AtomType of interest
    //! @return atom with type ATOM_TYPE
    const Atom &AA::GetAtom( const AtomType &ATOM_TYPE) const
    {
      // initialize undefined atom
      static const Atom s_undefined_atom;

      // return
      return s_undefined_atom;
    }

    //! @brief get AAClass
    //! @return AAClass
    const AAClass &AA::GetAAClass() const
    {
      return GetAAClasses().e_AA;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief assignment operator for AA class
    //! @param AA_RHS AA object to be copied
    //! @return this after assignment to AA_RHS is done
    AA &AA::operator =( const AA &AA_RHS)
    {
      // assign base class
      AABase::operator =( AA_RHS);

      // return
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AA::Read( std::istream &ISTREAM)
    {
      // read AABase
      AABase::Read( ISTREAM);

      //end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &AA::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write bases
      AABase::Write( OSTREAM, INDENT);

      //end
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
