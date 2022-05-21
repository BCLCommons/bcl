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

#ifndef BCL_UTIL_ENUM_HPP_
#define BCL_UTIL_ENUM_HPP_

// include the namespace header
#include "bcl_util_enum.hpp"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_assert.h"
#include "bcl_util_class_descriptor.h"
#include "bcl_util_enum_data.h"
#include "bcl_util_object_interface.h"
#include "bcl_util_undefined.h"
#include "io/bcl_io_serialization_interface.h"
#include "io/bcl_io_serialize.h"

// external includes - sorted alphabetically
#include <vector>

#ifndef BCL_ENUM_DECLARATION_ONLY

#ifndef BCL_UTIL_ENUM_IMPLEMENTATION_H_
#define BCL_UTIL_ENUM_IMPLEMENTATION_H_

namespace bcl
{
  namespace util
  {

    //! read from std::istream
    //! @param ISTREAM stream to read from
    template< typename t_DataType, typename t_Derived>
    std::istream &Enum< t_DataType, t_Derived>::Read( std::istream &ISTREAM)
      {
        // read name from stream
        std::string name;
        io::Serialize::Read( name, ISTREAM);

        // find the object belonging to the name and reassign the reference
        *this = t_Derived::GetEnums().GetEnumFromName( name);

        // end
        return ISTREAM;
      }

    //! @brief set the value of the corresponding member based on the label
      //! @param LABEL label that is used to set the string
      //! @param ERROR_STREAM stream to write errors to
      //! @return bool - if the parameter was set (parameters are not set if they are not allowed parameters)
      template< typename t_DataType, typename t_Derived>
      bool Enum< t_DataType, t_Derived>::TryRead( const ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
      {
        return t_Derived::GetEnums().SetEnumFromName( *this, LABEL.GetValue(), ERROR_STREAM);
      }

      //! @brief writes the help for the label
      //! @param OSTREAM the stream to which the help is written to
      //! @param INDENT the amount of indent
      //! @return the given stream to which the help was written to
      template< typename t_DataType, typename t_Derived>
      std::ostream &Enum<t_DataType, t_Derived>::WriteHelp( std::ostream &OSTREAM, const size_t INDENT) const
      {
        OSTREAM << "Choose from the following: ";
        return t_Derived::GetEnums().WriteList( OSTREAM);
      }

  } // namespace util

} // namespace bcl

#endif // BCL_UTIL_ENUM_IMPLEMENTATION_H_

#endif // BCL_ENUM_DECLARATION_ONLY

#endif // BCL_UTIL_ENUM_HPP_

