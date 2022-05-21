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

#ifndef BCL_CHEMISTRY_ATOM_CONFORMATIONAL_SHARED_H_
#define BCL_CHEMISTRY_ATOM_CONFORMATIONAL_SHARED_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_configurational_shared.h"
#include "bcl_chemistry_atom_conformational_interface.h"
#include "bcl_chemistry_bond_conformational.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AtomConformationalShared
    //! @brief This class contains atom conformation information
    //! @details stores atom configuration, 3D coordinates of atoms and conformation of bonds attached to atom.
    //!
    //! @see @link example_chemistry_atom_conformational_shared.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date Dec 05, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AtomConformationalShared :
      public AtomConformationalInterface
    {
    private:

    /////////////
    // friends //
    /////////////

      friend class ConformationalFactory;
      friend class AtomVector< AtomConformationalShared>;

    //////////////
    // typedefs //
    //////////////

      // these typedefs are needed so that AtomVector can be templated on a single type (AtomConformationShared), rather
      //than having to pass all these types in as well
      typedef BondConformational                t_Bond;      //!< Type of bonds this class contains
      typedef ConfigurationalBondTypeData::Data t_BondData;  //!< enum of data that can be obtained from a bond
      typedef AtomConformationalInterface       t_Interface; //!< atom interface class used by the bond

    //////////
    // data //
    //////////

      util::SiPtr< const AtomConfigurationalShared> m_Configuration; //!< configuration of atom
      storage::Vector< BondConformational>          m_Bonds;         //!< Bonds of this atom
      linal::Vector3D                               m_Coordinates;   //!< 3D coordinates of atom

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      AtomConformationalShared();

      //! @brief constructor from initializer, needed by atom vector
      //! @param ATOM_INFO all info about the atom, only coordinates will be used
      AtomConformationalShared( const sdf::AtomInfo &ATOM_INFO);

      //! virtual copy constructor
      AtomConformationalShared *Clone() const;

      //! @brief copy constructor
      //! @param ATOM atom conformation object of atom of interest
      AtomConformationalShared( const AtomConformationalShared &ATOM);

      //! @brief assignment constructor
      //! @param ATOM atom conformation object of atom of interest
      AtomConformationalShared &operator =( const AtomConformationalShared &ATOM);

    /////////////////
    // data access //
    /////////////////

    public:

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns bonds
      //! @return the bond conformations of the atom
      const storage::Vector< BondConformational> &GetBonds() const;

      //! @brief returns chirality
      //! @return reference of chirality attribute
      const Chirality &GetChirality() const;

      //! @brief returns atom type
      //! @return a reference of stom type data attribute
      const AtomType &GetAtomType() const;

      //! @brief returns coordinates
      //! @return a reference to the position vector
      const linal::Vector3D &GetPosition() const;

      //! @brief return a simple pointer to atom configurational object
      //! @return pointer to atom configurational object
      const util::SiPtr< const AtomConfigurationalShared> GetConfiguration() const;

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

      //! @param set atom
      //! @param ATOM_CONFIGURATIONAL atom configurational
      void SetAtom( const AtomConfigurationalInterface &ATOM_CONFIGURATIONAL)
      {
        m_Configuration = util::ToSiPtr( ATOM_CONFIGURATIONAL);
      }

      //! @brief set the position
      //! @param POSITION the new coordinates
      void SetPosition( const linal::Vector3D &POSITION)
      {
        m_Coordinates = POSITION;
      }

      //! @brief set bonds
      //! @param BONDS the bonds that the atom connects
      void SetBonds( const storage::Vector< BondConformational> &BONDS);

      //! @brief copy the atom and bonds using the specified difference in pointer address
      //! @param ATOM the atom to copy
      //! @param DIFF difference in begin of vector that atoms are located in
      void Copy( const AtomConformationalShared &ATOM, const ptrdiff_t &DIFF);

    }; // class AtomConformationalShared

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_ATOM_CONFORMATIONAL_SHARED_H_
