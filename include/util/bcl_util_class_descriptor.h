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

#ifndef BCL_UTIL_CLASS_DESCRIPTOR_H_
#define BCL_UTIL_CLASS_DESCRIPTOR_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically
#include <string>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//!
//! @file bcl_util_class_descriptor.h
//! @brief functions for converting the __PRETTY_FUNCTION__ emitted by the compiler into a standardized, readable string
//!
//! @author mendenjl
//! @date April 01, 2010
//! @see @link example_util_class_descriptor.cpp @endlink
//!
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace bcl
{
  namespace util
  {
    //! @brief Standardizes the output of __PRETTY_FUNCTION__ between compilers / architectures
    //! This function should only be called from a (non-member) template function (like GetStaticClassName)
    //! @param PRETTY_FUNCTION_NAME should be passed from the  '__PRETTY_FUNCTION__' macro
    //! @return the standardized class, struct, union, or enum name as std::string
    BCL_API std::string ClassNameViaTemplatedFunction( const std::string &PRETTY_FUNCTION_NAME);

    //! @brief Standardize a class name
    //! @param PRETTY_FUNCTION_NAME should be passed from ClassNameViaTemplatedFunction or from a file
    //! @return the standardized class, struct, union, or enum name as std::string
    BCL_API std::string StandardizeClassName( const std::string &PRETTY_FUNCTION_NAME);

    //! @brief extracts the class name before the "::"
    //! The '__PRETTY_FUNCTION__' macro provided by the compiler should be used.
    //! If the Line looks like "virtual const std::string& bcl::Test::GetName() const" the string "bcl::Test" will be extracted
    //! @param PRETTY_FUNCTION_NAME should be passed from the  '__PRETTY_FUNCTION__' macro
    //! @return the Class name as std::string
    BCL_API const std::string ExtractNamespaceIdentifier( const std::string &IDENTIFIER);

  } // namespace util

  //! @brief Get the standardized name of a class/struct/union/or enum without a member of the class
  //! example usage: GetStaticClassName< MyClass< double> >();
  //! @see @link example_util_class_descriptor.cpp @endlink
  //! @return the standardized class, struct, union, or enum name as std::string
  //! @note does not follow inheritance hierarchy, so every (non-abstract) class should still have:
  //! virtual const std::string &GetClassIdentifier() const
  //! {
  //!   return GetStaticClassName( *this);
  //! }
  template< typename t_ClassType>
  inline
  const std::string &GetStaticClassName()
  {
    static const std::string s_name( util::ClassNameViaTemplatedFunction( __PRETTY_FUNCTION__));
    return s_name;
  }

  //! @brief Get the standardized name of a class/struct/union/or enum without a member of the class
  //! example usage: MyClass< double> a; GetStaticClassName(a);
  //! @see @link example_util_class_descriptor.cpp @endlink
  //! @return the standardized class, struct, union, or enum name as std::string
  //! @note non-abstract classes need:
  //! const std::string &GetClassIdentifier() const
  //! {
  //!   return GetStaticClassName( *this);
  //! }
  template< typename t_ClassType>
  inline
  const std::string &GetStaticClassName( t_ClassType const &)
  {
    return GetStaticClassName< t_ClassType>();
  }

} // namespace bcl

#endif // BCL_UTIL_CLASS_DESCRIPTOR_H_ 

