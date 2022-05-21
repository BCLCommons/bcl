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

#ifndef BCL_UTIL_OBJECT_INTERFACE_H_
#define BCL_UTIL_OBJECT_INTERFACE_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically
#include <iosfwd>

namespace bcl
{
  namespace util
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ObjectInterface
    //! @brief bcl standard interface. every class has to be derived from
    //! @details This class provides functions necessary for each bcl object:
    //! a virtual Clone and Empty function that guarantees to get the correct object when you are pointing to a class and
    //! want to retrieve a copy without slicing. An Empty that helps to build up a repository of those objects.
    //! A GetClassIdentifier() function, that gives Users the possibility to know which object is actually behind a pointer.
    //!
    //! @see @link example_util_object_interface.cpp @endlink
    //! @author woetzen
    //! @date 09.08.2007
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ObjectInterface
    {
    /////////////
    // friends //
    /////////////

      template< typename t_DataType> friend class ShPtr;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! A pointer to a base class interface will return a new copy of the actual object that is pointed to.
      //! In order to make this work derived objects have to overwrite this with:
      //! virtual {ClassName} *Clone() const
      //! {
      //!   return new {ClassName}( *this);
      //! }
      //! @return new Pointer to a copy of the actual object behind the pointer
      virtual ObjectInterface *Clone() const = 0;

      //! @brief virtual destructor
      virtual ~ObjectInterface();

      //! @brief returns class name of the object behind a pointer or the current object
      //! To guarantee the hierarchy it is necessary to overwrite this function in every class with this code:
      //! virtual const std::string &GetClassIdentifier() const
      //! {
      //!   return GetStaticClassName( *this);
      //! }
      //! @return the class name if this function is overwritten
      virtual const std::string &GetClassIdentifier() const = 0;

    //////////////////////
    // input and output //
    //////////////////////

    private:

      //! @brief read from std::istream
      //! this can only be called, if the identifier is already extracted from the stream, as the ShPtr will do it
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &ReadWithOutIdentifier( std::istream &ISTREAM)
      {
        return Read( ISTREAM);
      }

    public:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      virtual std::istream &Read( std::istream &ISTREAM) = 0;

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return output stream which was written to
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const = 0;

      //! @brief read the identifier
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &ReadIdentifier( std::istream &ISTREAM) const;

      //! @brief write the identifier
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return output stream which was written to
      std::ostream &WriteIdentifier( std::ostream &OSTREAM, const size_t INDENT) const;

      //! @brief read from std::istream
      //! This is the main function to be used for reading an object. It will read the first line assuming there is the
      //! class name as identifier - will be asserted, and after that it will call the overwritten Read function for the
      //! object.
      //! @param ISTREAM input stream
      //! @return istream which was read from
      virtual std::istream &ReadObject( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! This is the main function to be used for writing an object. It will write the classname as identifier,
      //! and after that it will call the overwritten Write function for the object
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return output stream which was written to
      virtual std::ostream &WriteObject( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief extracts the identifier string for the given ISTREAM (removes it)
      //! @param ISTREAM input stream where identifier is extracted from
      //! @return identifier string
      static std::string ExtractIdentifier( std::istream &ISTREAM);

      //! @brief extract namespace name from the output of an arbitrary __PRETTY_FUNCTION_NAME__ macro
      //! @param IDENTIFIER object identifier such as "bcl::util::ClassName"
      //! @return namespace name without the "bcl::" part such as "util" or "assemble"
      static const std::string ExtractNamespaceName( const std::string &PRETTY_FUNCTION_NAME);

    }; // class Object Interface

    //! @brief operator >> to read ReadWrite derived Object from std::istream
    //! @param ISTREAM input stream the object is read from
    //! @param OBJECT the object derived from ReadWrite
    //! @return the input stream where the object was read from
    BCL_API std::istream &operator >>( std::istream &ISTREAM, ObjectInterface &OBJECT);

    //! @brief operator >> to read const ReadWrite derived Object from std::istream which will assert, since this is not possible
    //! @param ISTREAM input stream the object is read from
    //! @param OBJECT the object derived from ReadWrite
    //! @return the input stream where the object was read from
    BCL_API std::istream &operator >>( std::istream &ISTREAM, const ObjectInterface &OBJECT);

    //! @brief write ReadWrite derived Object to std::ostream
    //! @param OSTREAM output stream the object is written to
    //! @param OBJECT object derived from ReadWrite
    //! @return output stream the object was written to
    BCL_API std::ostream &operator <<( std::ostream &OSTREAM, const ObjectInterface &OBJECT);

  } // namespace util
} // namespace bcl

#endif //BCL_UTIL_OBJECT_INTERFACE_H_
