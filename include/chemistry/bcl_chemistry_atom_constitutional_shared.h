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

#ifndef BCL_CHEMISTRY_ATOM_CONSTITUTIONAL_SHARED_H_
#define BCL_CHEMISTRY_ATOM_CONSTITUTIONAL_SHARED_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_bond_constitutional.h"

// external includes - sorted alphabetically
#include <cstddef>

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AtomConstitutionalShared
    //! @brief This class contains atom constitution information
    //! @details  stores the atom type, bond type of bonds that atom has and valence bonds of the atom
    //!
    //! @see @link example_chemistry_atom_constitutional_shared.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date Dec 05, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AtomConstitutionalShared :
      public AtomConstitutionalInterface
    {
    private:

    /////////////
    // friends //
    /////////////

      friend class AtomVector< AtomConstitutionalShared>;

    //////////////
    // typedefs //
    //////////////

      // these typedefs are needed so that AtomVector can be templated on a single type (AtomConstitutionShared), rather
      //than having to pass all these types in as well
      typedef BondConstitutional                    t_Bond;      //!< Type of bonds this class contains
      typedef ConstitutionalBondTypeData::DataEnum  t_BondData;  //!< enum of data that can be obtained from a bond
      typedef AtomConstitutionalInterface           t_Interface; //!< atom interface class used by the bond

    //////////
    // data //
    //////////

      AtomType                                 m_AtomType;     //!< Atom type of atom
      storage::Vector< BondConstitutional>     m_Bonds;        //!< Bonds of this atom

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      AtomConstitutionalShared();

      //! @brief constructor from initializer, needed by atom vector
      //! @param ATOM_INFO atom information
      AtomConstitutionalShared( const sdf::AtomInfo &ATOM_INFO);

      //! virtual copy constructor
      AtomConstitutionalShared *Clone() const;

      //! @brief copy constructor
      //! @param ATOM atom constitution object of atom of interest
      AtomConstitutionalShared( const AtomConstitutionalShared &ATOM);

      //! @brief assignment constructor
      //! @param ATOM atom constitution object of atom of interest
      AtomConstitutionalShared &operator =( const AtomConstitutionalShared &ATOM);

    /////////////////
    // data access //
    /////////////////

    public:

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns atom attribute
      //! @return a reference to atom attribute
      const AtomType &GetAtomType() const;

      //! @brief returns bond configuration
      //! @return a reference to bond configuration objects
      const storage::Vector< BondConstitutional> &GetBonds() const;

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

    private:

      //! @brief set bonds constitution of atom
      //! @param BONDS the bonds that the atom connects
      void SetBonds( const storage::Vector< BondConstitutional> &BONDS);

      //! @brief copy the atom and bonds using the specified difference in pointer address
      //! @param ATOM the atom to copy
      //! @param DIFF difference in begin of vector that atoms are located in
      void Copy( const AtomConstitutionalShared &ATOM, const ptrdiff_t &DIFF);

    }; // class AtomConstitutionalShared

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_ATOM_CONSTITUTIONAL_SHARED_H_
