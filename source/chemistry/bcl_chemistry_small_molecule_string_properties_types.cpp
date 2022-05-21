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
#include "chemistry/bcl_chemistry_small_molecule_string_properties_types.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // Instance that gets atom types of all atoms in the molecule
    const util::SiPtr< const util::ObjectInterface> SmallMoleculeStringPropertiesTypes::s_GetAtomTypesInstance
    (
      util::Enumerated< StringPropertyInterface>::AddInstance
      (
        new SmallMoleculeStringPropertiesTypes
        (
          "AtomTypes",
          "Retrieves the atom types of all atoms in a molecule",
          function::MemberFunction( &ConformationInterface::GetAtomTypesString)
        )
      )
    );

    // Instance that gets bond types from the molecule
    const util::SiPtr< const util::ObjectInterface> SmallMoleculeStringPropertiesTypes::s_GetBondTypeInstance
    (
      util::Enumerated< StringPropertyInterface>::AddInstance
      (
        new SmallMoleculeStringPropertiesTypes
        (
          "BondTypes",
          "Retrieves the bond types of all bonds in a molecule",
          function::MemberFunction( &ConformationInterface::GetBondTypesString)
        )
      )
    );

    // Instance that gets chirality types of all atoms in the molecule
    const util::SiPtr< const util::ObjectInterface> SmallMoleculeStringPropertiesTypes::s_GetChiralityInstance
    (
      util::Enumerated< StringPropertyInterface>::AddInstance
      (
        new SmallMoleculeStringPropertiesTypes
        (
          "Chirality",
          "Retrieves the chirality of all atoms in a molecule",
          function::MemberFunction( &ConformationInterface::GetChiralityString)
        )
      )
    );

    // Instance that gets chirality types of all atoms in the molecule
    const util::SiPtr< const util::ObjectInterface> SmallMoleculeStringPropertiesTypes::s_GetSumFormulaInstance
    (
      util::Enumerated< StringPropertyInterface>::AddInstance
      (
        new SmallMoleculeStringPropertiesTypes
        (
          "SumFormula",
          "Retrieves the sum formula a molecule",
          function::MemberFunction( &ConformationInterface::GetSumFormula)
        )
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from members
    //! @param ALIAS alias used for this property over the command line
    //! @param DESCRIPTION description shown for this class in the help
    //! @param FUNCTION function used to retrieve a string from a small molecule
    SmallMoleculeStringPropertiesTypes::SmallMoleculeStringPropertiesTypes
    (
      const std::string &ALIAS,
      const std::string &DESCRIPTION,
      const function::MemberConst< ConformationInterface, std::string> &FUNCTION
    ) :
      m_Alias( ALIAS),
      m_ClassDescription( DESCRIPTION),
      m_StringRetrievalFunction( FUNCTION)
    {
    }

    //! virtual copy constructor
    SmallMoleculeStringPropertiesTypes *SmallMoleculeStringPropertiesTypes::Clone() const
    {
      return new SmallMoleculeStringPropertiesTypes( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SmallMoleculeStringPropertiesTypes::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &SmallMoleculeStringPropertiesTypes::GetAlias() const
    {
      return m_Alias;
    }

    //! @brief operator the implements the assignment operation on the two arguments returning a result
    //! @param MOLECULE the molecule to calculate the property for
    //! @return the property as a string
    std::string SmallMoleculeStringPropertiesTypes::operator()( const ConformationInterface &MOLECULE) const
    {
      // return the output of the function of interest
      return util::TrimString( m_StringRetrievalFunction( MOLECULE));
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SmallMoleculeStringPropertiesTypes::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( m_ClassDescription);
      return parameters;
    }
  } // namespace chemistry
} // namespace bcl
