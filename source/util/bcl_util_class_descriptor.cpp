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
#include "util/bcl_util_class_descriptor.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_name_standardizer.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    //! @brief Standardizes the output of __PRETTY_FUNCTION__ between compilers / architectures
    //! This function should only be called from a (non-member) template function (like GetStaticClassName)
    //! @param PRETTY_FUNCTION_NAME should be passed from the  '__PRETTY_FUNCTION__' macro
    //! @return the standardized class, struct, union, or enum name as std::string
    std::string ClassNameViaTemplatedFunction( const std::string &PRETTY_FUNCTION_NAME)
    {
      return ClassNameStandardizer::StandardizeTemplateFunctionParameter( PRETTY_FUNCTION_NAME);
    }

    //! @brief extracts the class name
    //! The '__PRETTY_FUNCTION__' macro provided by the compiler should be used.
    //! If the Line looks like "virtual const std::string& bcl::Test< T>::GetName() const [with T = double]" the
    //! string "double" will be extracted
    //! @param PRETTY_FUNCTION_NAME should be passed from the  '__PRETTY_FUNCTION__' macro
    //! @return the Class name as std::string
    std::string StandardizeClassName( const std::string &PRETTY_FUNCTION_NAME)
    {
      return ClassNameStandardizer::Standardize( PRETTY_FUNCTION_NAME);
    }

    //! @brief extracts the namespace name before the "::"
    //! The '__PRETTY_FUNCTION__' macro provided by the compiler should be used.
    //! Example: ExtractNamespaceIdentifier( "virtual const std::string& bcl::Test::GetName() const) " the string "bcl::Test" will be extracted
    //! Algorithm: Finds the first '(' after the first
    //! @param PRETTY_FUNCTION_NAME should be passed from the  '__PRETTY_FUNCTION__' macro
    //! @return the Class name as std::string
    const std::string ExtractNamespaceIdentifier( const std::string &PRETTY_FUNCTION_NAME)
    {
      // standardize name first
      const std::string standardized_name( StandardizeClassName( PRETTY_FUNCTION_NAME));

      // find the first :: in the string
      const size_t first_scope( standardized_name.find( "::"));
      if( first_scope == std::string::npos)
      {
        return standardized_name; // no scope; the user probably passed a string that was already a namespace
      }

      // find the first of (, and < after the first :: in the string.  These cannot be in a namespace name
      const size_t function_parens( standardized_name.find_first_of( "(<", first_scope));
      if( function_parens == std::string::npos)
      {
        return standardized_name; // no function start, the string was likely already a namespace
      }

      const size_t last_pos( standardized_name.rfind( "::", function_parens));
      if( function_parens == std::string::npos)
      {
        return std::string(); // no scope, just a function name, so no namespace
      }

      // find the first non-alphanumeric character other than ':' and '_', starting with the character before (
      size_t first_pos( function_parens - 1);
      while
      (
        first_pos > 0
        &&
        (
          isalnum( standardized_name[ first_pos - 1])
          || standardized_name[ first_pos - 1] == '_'
          || standardized_name[ first_pos - 1] == ':'
        )
      )
      {
        --first_pos;
      }

      //return the substring between begin and end of class name
      return standardized_name.substr( first_pos, last_pos - first_pos);
    }

  } // namespace util
} // namespace bcl

