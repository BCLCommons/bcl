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

#ifndef BCL_CHEMISTRY_BOND_CONFORMATIONAL_H_
#define BCL_CHEMISTRY_BOND_CONFORMATIONAL_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_constitutional_bond_types.h"
#include "util/bcl_util_si_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BondConformational
    //! @brief This class contains bond conformational information, i.e. atom and bond type involved.
    //! @details information about bond configuration and atom conformation to which bond belongs
    //!
    //! @see @link example_chemistry_bond_conformational.cpp @endlink
    //! @author mendenjl, kothiwsk
    //! @date Dec 02, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API BondConformational :
      public util::ObjectInterface
    {

    private:

    /////////////
    // friends //
    /////////////

      friend class AtomComplete; //!< AtomComplete needs to be able to set bond types
      friend class FragmentCompleteStandardizer; //!< needs to be able to set the bond types

    //////////
    // data //
    //////////

      util::SiPtr< const AtomConformationalInterface> m_AtomConformational; //!< The atom to which this bond connects
      ConfigurationalBondType m_BondType;                                   //!< Type of the bond

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      BondConformational();

      //! @brief constructor from AtomConformational and BondType charge
      //! @param ATOM_CONFORMATIONAL atom involved in the bond
      //! @param BASE_BOND bondtype of the bond
      BondConformational
      (
        const AtomConformationalInterface &ATOM_CONFORMATIONAL,
        const ConfigurationalBondType &BASE_BOND
      );

      //! virtual copy constructor
      BondConformational *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get a reference to bond configurational attribute
      //! @return a reference to bond configurational attribute
      const ConfigurationalBondType &GetBondType() const;

      //! @brief get a reference to atom conformational to which bond belongs
      //! @return a reference to atom conformational to which bond belongs
      const AtomConformationalInterface &GetTargetAtom() const;

      //! @brief test whether the pointer is to a defined atom
      //! @return true if this bond is defined
      bool IsDefined() const;

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

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief set bond type of the bond
      //! @param BOND_TYPE the desired bond type
      void SetBondType( const ConfigurationalBondType &BOND_TYPE);

    }; // class BondConformational

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_BOND_CONFORMATIONAL_H_
