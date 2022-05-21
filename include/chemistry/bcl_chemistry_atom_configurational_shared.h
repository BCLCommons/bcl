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

#ifndef BCL_CHEMISTRY_ATOM_CONFIGURATIONAL_SHARED_H_
#define BCL_CHEMISTRY_ATOM_CONFIGURATIONAL_SHARED_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_configurational_interface.h"
#include "bcl_chemistry_atom_constitutional_shared.h"

// external includes - sorted alphabetically
#include <cstddef>

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AtomConfigurationalShared
    //! @brief This class contains store atom configurational information
    //! @details stores atom chirality, atom constitution and configurations of bonds attached to the atom .
    //!
    //! @see @link example_chemistry_atom_configurational_shared.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date Dec 05, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AtomConfigurationalShared :
      public AtomConfigurationalInterface
    {
    private:

    /////////////
    // friends //
    /////////////

      friend class AtomVector< AtomConfigurationalShared>;
      friend class FragmentConfigurationShared;

    //////////////
    // typedefs //
    //////////////

      // these typedefs are needed so that AtomVector can be templated on a single type (AtomConfigurationShared),rather
      //than having to pass all these types in as well
      typedef BondConfigurational               t_Bond;      //!< Type of bonds this class contains
      typedef ConfigurationalBondTypeData::Data t_BondData;  //!< enum of data that can be obtained from a bond
      typedef AtomConfigurationalInterface      t_Interface; //!< atom interface class used by the bond

      //! data needed to initialize this class
      typedef Chirality                         t_Initializer;

    //////////
    // data //
    //////////

      ChiralityEnum                                m_Chirality;    //!< chirality of atom
      util::SiPtr< const AtomConstitutionalShared> m_Constitution; //!< constitution of atom
      storage::Vector< BondConfigurational>        m_Bonds;        //!< Bonds of this atom

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      AtomConfigurationalShared();

      //! @brief constructor from initializer, needed by atom vector
      //! @param ATOM_INFO all info about the atom; only chirality will be used
      AtomConfigurationalShared( const sdf::AtomInfo &ATOM_INFO);

      //! virtual copy constructor
      AtomConfigurationalShared *Clone() const;

      //! @brief copy constructor
      //! @param ATOM atom configuration object of atom of interest
      AtomConfigurationalShared( const AtomConfigurationalShared &ATOM);

      //! @brief assignment constructor
      //! @param ATOM atom configuration object of atom of interest
      AtomConfigurationalShared &operator =( const AtomConfigurationalShared &ATOM);

    /////////////////
    // data access //
    /////////////////

    public:

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns chirality
      //! @return a reference to chirality attribute
      const Chirality &GetChirality() const;

      //! @brief returns atom type
      //! @return a reference of stom type data attribute
      const AtomType &GetAtomType() const;

      //! @brief returns bonds
      //! @return the bond configuration of bond connected the atom
      const storage::Vector< BondConfigurational> &GetBonds() const;

      //! @brief return a simple pointer to atom constitutional object
      //! @return pointer to atom constitutional object
      const util::SiPtr< const AtomConstitutionalShared> GetConstitution() const;

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

    private:

      //! @brief set bonds
      //! @param BONDS the bonds that the atom connects
      void SetBonds( const storage::Vector< BondConfigurational> &BONDS);

      //! @param set atom
      //! @param ATOM_CONSTITUTIONAL atom constitutional
      void SetAtom( const AtomConstitutionalInterface &ATOM_CONSTITUTIONAL)
      {
        m_Constitution = util::ToSiPtr( ATOM_CONSTITUTIONAL);
      }

      //! @brief copy the atom and bonds using the specified difference in pointer address
      //! @param ATOM the atom to copy
      //! @param DIFF difference in begin of vector that atoms are located in
      void Copy( const AtomConfigurationalShared &ATOM, const ptrdiff_t &DIFF);

    }; // class AtomConfigurationalShared

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_ATOM_CONFIGURATIONAL_SHARED_H_
