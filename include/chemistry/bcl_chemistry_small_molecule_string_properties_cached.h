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

#ifndef BCL_CHEMISTRY_SMALL_MOLECULE_STRING_PROPERTIES_CACHED_H_
#define BCL_CHEMISTRY_SMALL_MOLECULE_STRING_PROPERTIES_CACHED_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_interface.h"
#include "bcl_chemistry_string_property_interface.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SmallMoleculeStringPropertiesCached
    //! @brief retrieves a string from a molecule by calling a molecule member function
    //!
    //! @see @link example_chemistry_small_molecule_string_properties_cached.cpp @endlink
    //! @author mendenjl
    //! @date Nov 29, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SmallMoleculeStringPropertiesCached :
      public StringPropertyInterface
    {

    private:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! name of a cached property to retrieve
      util::ObjectDataLabel m_CachedPropertyName;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SmallMoleculeStringPropertiesCached();

      //! virtual copy constructor
      SmallMoleculeStringPropertiesCached *Clone() const;

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

      //! @brief Set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read of this class
      //! @param ERROR_STREAM stream with which to write errors
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    }; // class SmallMoleculeStringPropertiesCached

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_SMALL_MOLECULE_STRING_PROPERTIES_CACHED_H_
