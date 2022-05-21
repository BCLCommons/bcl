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

#ifndef BCL_CHEMISTRY_SMALL_MOLECULE_STRING_PROPERTIES_MAPPED_H_
#define BCL_CHEMISTRY_SMALL_MOLECULE_STRING_PROPERTIES_MAPPED_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_interface.h"
#include "bcl_chemistry_string_property_interface.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_function_interface.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_si_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SmallMoleculeStringPropertiesMapped
    //! @brief uses one property value as a key to return a value for that key that is given in a file
    //! Input files must have at least one space/tab per row, with the key before the first space/tab
    //! Everything after the first space/tab on each line will be the value (expect additional spaces)
    //! in the first row, e.g.
    //! ID  Value
    //! 2   5.0 6.0 12.0
    //! 1   3.0 14.0 18.0
    //!
    //! All properties can be strings too, e.g.:
    //! Letter  NextLetter
    //! Alpha   Beta
    //! Beta    Gamma
    //! Gamma   Phi
    //!
    //! It is an error if the value property is omitted for a particular key
    //!
    //! @see @link example_chemistry_small_molecule_string_properties_mapped.cpp @endlink
    //! @author mendenjl
    //! @date Mar 16, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SmallMoleculeStringPropertiesMapped :
      public StringPropertyInterface
    {

    private:

    //////////
    // data //
    //////////

      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      std::string m_Filename;       //!< name of the file containing the table

      std::string m_Delimiters;     //!< Set of allowed delimiters

      // implementation of the property that retrieves the key property for the molecule
      util::Implementation< StringPropertyInterface> m_KeyProperty;

      storage::Map< std::string, std::string> m_Map; //!< Map from key to value

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SmallMoleculeStringPropertiesMapped() :
        m_Delimiters( " \t,;=")
      {
      }

      //! virtual copy constructor
      SmallMoleculeStringPropertiesMapped *Clone() const;

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

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief given a line, read the next key and value into passed-in strings
      //! @param LINE the next line
      //! @param KEY storage for the key
      //! @param VALUE storage for the value
      //! @return true iff a value and key were each read
      bool ReadKeyAndValue( const std::string &LINE, std::string &KEY, std::string &VALUE) const;

    }; // class SmallMoleculeStringPropertiesMapped

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_SMALL_MOLECULE_STRING_PROPERTIES_MAPPED_H_
