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

#ifndef BCL_OPENCL_EXTENSION_DATA_H_
#define BCL_OPENCL_EXTENSION_DATA_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ExtensionData
    //! @brief this class contains the data for an opencl extension, which is always a combined string of the extension
    //!        vendor and the name of the extensions. Enumerated, it provides easy handling of opencl extensions.
    //!
    //! @see @link example_opencl_extension_data.cpp @endlink
    //! @author woetzen
    //! @date Aug 24, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ExtensionData :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      std::string m_Vendor; //!< vendor name for extension e.g. "khr", "amd", "nv" ...
      std::string m_Name;   //!< name of the extension e.g. "fp64"

    public:

    //////////
    // data //
    //////////

      //! @brief prefix to all extensions "cl"
      //! @return prefix string for opencl extensions
      static const std::string &GetPrefix();

      //! @brief separator for prefix, vendor and name
      //! @return separator '_'
      static char GetSeparator();

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ExtensionData();

      //! @brief construct from extension string
      //! @param EXTENSION_STRING string for extension with prefix, vendor and name, e.g. "cl_amd_fp64"
      ExtensionData( const std::string &EXTENSION_STRING);

      //! @brief Clone function
      //! @return pointer to new ExtensionData
      ExtensionData *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the vendor, eg. "amd"
      //! @return name of vendor
      const std::string &GetVendor() const;

      //! @brief return the name of extension, eg. "nv"
      //! @return name of extension
      const std::string &GetName() const;

      //! @brief get the extension string made up of prefix, vendor and name
      //! @return full extension string, e.g. "cl_amd_fp64"
      std::string GetString() const;

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
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief set vendor name and extension name form extension string
      //! @param EXTENSION_STRING string for extension with prefix, vendor and name, e.g. "cl_amd_fp64"
      //! @return true, if successful
      bool SetMembersFromExtensionString( const std::string &EXTENSION_STRING);

    }; // class ExtensionData

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_EXTENSION_DATA_H_ 
