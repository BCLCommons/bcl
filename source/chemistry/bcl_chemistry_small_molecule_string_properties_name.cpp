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
#include "chemistry/bcl_chemistry_small_molecule_string_properties_name.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    //! Instance that gets the (trimmed) name of the molecule
    const util::SiPtr< const util::ObjectInterface> SmallMoleculeStringPropertiesName::s_Instance
    (
      util::Enumerated< StringPropertyInterface>::AddInstance
      (
        new SmallMoleculeStringPropertiesName()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! virtual copy constructor
    SmallMoleculeStringPropertiesName *SmallMoleculeStringPropertiesName::Clone() const
    {
      return new SmallMoleculeStringPropertiesName( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SmallMoleculeStringPropertiesName::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &SmallMoleculeStringPropertiesName::GetAlias() const
    {
      static const std::string s_name( "Name");
      return s_name;
    }

    //! @brief operator the implements the assignment operation on the two arguments returning a result
    //! @param MOLECULE the molecule to calculate the property for
    //! @return the property as a string
    std::string SmallMoleculeStringPropertiesName::operator()( const ConformationInterface &MOLECULE) const
    {
      std::string original_name( util::TrimString( MOLECULE.GetName()));
      std::replace( original_name.begin(), original_name.end(), '\n', '\t');
      // return the name, with newlines replaced by tabs
      return original_name;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SmallMoleculeStringPropertiesName::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Retrieves the name of the molecule"
      );
      return parameters;
    }
  } // namespace chemistry
} // namespace bcl
