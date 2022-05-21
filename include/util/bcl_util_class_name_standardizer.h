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

#ifndef BCL_UTIL_CLASS_NAME_STANDARDIZER_H_
#define BCL_UTIL_CLASS_NAME_STANDARDIZER_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically
#include <string>

namespace bcl
{
  namespace util
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ClassNameStandardizer
    //! @brief standardizes the output of __PRETTY_FUNCTION__ names across MSVS and GCC
    //!
    //! @see @link example_util_class_descriptor.cpp @endlink
    //! @author mendenjl
    //! @date Apr 26, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ClassNameStandardizer
    {
    public:

      //! @brief Standardize a class name that has already been extracted from any template
      //! @param NAME the name from GetStaticClassName<>
      //! @return the standardized class name
      static std::string Standardize( const std::string &NAME);

      //! @brief Standardize a class name that came from a function template (e.g. GetStaticClassName<>)
      //! @param NAME the name from the GetStaticClassName<>
      //! @return the standardized class name
      static std::string StandardizeTemplateFunctionParameter( const std::string &NAME);

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief a helper function for standardizing a class name
      //! @param NAME should be the name passed only from Standardize class name
      //! @param WORD should be the word to migrate forward in the type declaration to the last sequence character (*>,)
      static void MigrateWordForward( std::string &NAME, const std::string &WORD);

      //! @brief a helper function for standardizing a class name
      //! @param NAME should be the name passed only from any of the standardize class name functions
      static void AbbreviateEnums( std::string &NAME);

      //! @brief Remove struct/class/union/enum artifacts put into __PRETTY_NAME__ by the MSVS compiler
      //! @param MSVS_NAME the output of the __PRETTY_NAME__
      //! @return MSVS_NAME without the MSVS artifacts.  Should not be used if the string came from GCC
      static std::string HandleMSVSDataTypeString( const std::string MSVS_NAME);

      //! @brief extracts the template list
      //! The '__PRETTY_FUNCTION__' macro provided by the compiler should be used.
      //! If the Line looks like "virtual const std::string& bcl::Test< T>::GetName() const [with T = double]" the
      //! string "double" will be extracted
      //! @param PRETTY_FUNCTION_NAME should be passed from the  '__PRETTY_FUNCTION__' macro
      //! @return the template parameter in PRETTY_FUNCTION_NAME
      static std::string ExtractTemplateParameter( const std::string &PRETTY_FUNCTION_NAME);

      //! @brief remove any ul ui, l, i, or other suffix from numerical template parameters
      static void StandardizeNumericTemplateParameters( std::string &NAME);

    }; // class ClassNameStandardizer

  } // namespace util
} // namespace bcl

#endif // BCL_UTIL_CLASS_NAME_STANDARDIZER_H_ 

