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

#ifndef BCL_CHEMISTRY_BOND_CONFIGURATIONAL_H_
#define BCL_CHEMISTRY_BOND_CONFIGURATIONAL_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_configurational_interface.h"
#include "bcl_chemistry_configurational_bond_types.h"
#include "util/bcl_util_si_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BondConfigurational
    //! @brief This class contains bond configurational information.
    //! @details information about bond configuration and atom configuration to which bond belongs
    //!
    //! @see @link example_chemistry_bond_configurational.cpp @endlink
    //! @author kothiwsk
    //! @date Dec 02, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API BondConfigurational :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      util::SiPtr< const AtomConfigurationalInterface> m_AtomConfigurational; //!< The atom to which this bond connects
      ConfigurationalBondType m_BondType;                                     //!< Type of the bond

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      BondConfigurational();

      //! @brief constructor from AtomConfigurational and BondType charge
      //! @param ATOM_CONFIGURATIONAL atom involved in the bond
      //! @param BASE_BOND bondtype of the bond
      BondConfigurational
      (
        const AtomConfigurationalInterface &ATOM_CONFIGURATIONAL,
        const ConfigurationalBondType &BASE_BOND
      );

      //! virtual copy constructor
      BondConfigurational *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get a reference to bond configurational attribute
      //! @return a reference to bond configurational attribute
      const ConfigurationalBondType &GetBondType() const;

      //! @brief get a reference to atom configuration to which bond belongs
      //! @return a reference to atom configuration to which bond belongs
      const AtomConfigurationalInterface &GetTargetAtom() const;

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

    }; // class BondConfigurational

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_BOND_CONFIGURATIONAL_H_
