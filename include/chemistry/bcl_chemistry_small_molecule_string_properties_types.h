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

#ifndef BCL_CHEMISTRY_SMALL_MOLECULE_STRING_PROPERTIES_TYPES_H_
#define BCL_CHEMISTRY_SMALL_MOLECULE_STRING_PROPERTIES_TYPES_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_configuration_interface.h"
#include "bcl_chemistry_conformation_interface.h"
#include "bcl_chemistry_string_property_interface.h"
#include "function/bcl_function_member_const.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SmallMoleculeStringPropertiesTypes
    //! @brief retrieves a string from a molecule by calling a molecule member function
    //!
    //! @see @link example_chemistry_small_molecule_string_properties_types.cpp @endlink
    //! @author mendenjl
    //! @date Nov 29, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SmallMoleculeStringPropertiesTypes :
      public StringPropertyInterface
    {

    private:

    //////////
    // data //
    //////////

      //! Instance that gets atom types of all atoms in the molecule
      static const util::SiPtr< const util::ObjectInterface> s_GetAtomTypesInstance;

      //! Instance that gets bond types from the molecule
      static const util::SiPtr< const util::ObjectInterface> s_GetBondTypeInstance;

      //! Instance that gets chirality types of all atoms in the molecule
      static const util::SiPtr< const util::ObjectInterface> s_GetChiralityInstance;

      //! Instance that gets sum formula of all atoms in the molecule
      static const util::SiPtr< const util::ObjectInterface> s_GetSumFormulaInstance;

      //! alias
      std::string m_Alias;

      //! class description
      std::string m_ClassDescription;

      //! function that actually retrieves a string from a small molecule
      function::MemberConst< ConformationInterface, std::string> m_StringRetrievalFunction;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from members
      //! @param ALIAS alias used for this property over the command line
      //! @param DESCRIPTION description shown for this class in the help
      //! @param FUNCTION function used to retrieve a string from a small molecule
      SmallMoleculeStringPropertiesTypes
      (
        const std::string &ALIAS,
        const std::string &DESCRIPTION,
        const function::MemberConst< ConformationInterface, std::string> &FUNCTION
      );

      //! virtual copy constructor
      SmallMoleculeStringPropertiesTypes *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief operator the implements the assignment operation on the two arguments returning a result
      //! @param MOLECULE the molecule to calculate the property for
      //! @return the property as a string
      std::string operator()( const ConformationInterface &MOLECULE) const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class SmallMoleculeStringPropertiesTypes

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_SMALL_MOLECULE_STRING_PROPERTIES_TYPES_H_
