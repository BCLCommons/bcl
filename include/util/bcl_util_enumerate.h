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

#ifndef BCL_UTIL_ENUMERATE_H_
#define BCL_UTIL_ENUMERATE_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_enum.h"

// external includes - sorted alphabetically
#include <list>

namespace bcl
{
  namespace util
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Enumerate
    //! @brief provides a framework for enumerating different instances of a class
    //! @details provides a framework for enumerating different instances of a class where only one particular
    //! instance is necessary. t_DataType is the class type that is enumerated. E.g. this can be ElementTypeData as for
    //! each elementtype the properties are constant and known a priori. The Enumerate now provides the framework, to
    //! to hold the instances of this class with an additional name to uniquely identifies that instance and an index,
    //! that starts with 0 and is automatically incremented for each additional instances, that is added using the
    //! AddEnum function.
    //!
    //! class EnumEx:
    //!   public Enumerate< biol::ElementTypeData, EnumEx>
    //! {
    //! private:
    //!   const EnumType e_Hydrogen;
    //!   const EnumType e_Helium;
    //!   //...code
    //! }
    //!
    //! note that one template type is the derived class itself. When deriving an enumerator every enum will be one
    //! datamember. The will be one private default constructor that is required by the GetEnums function implemented in
    //! the Enumerate and is the singleton that provides the only access to the one and only instance of the derived
    //! Enumerator. That default constructor initializes all EnumType members properly, by calling
    //! AddEnum( NAME, t_DataType( ...)). The unique name for that t_DataType instance has to be determined, the index
    //! is auto incremented and automatically assigned.
    //!
    //! @see @link example_util_enumerate.cpp @endlink
    //! @author heinzes1, woetzen, karakam
    //! @date May 08, 2007
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_Derived>
    class Enumerate :
      public ObjectInterface
    {
    public:

      //! @brief EnumDataType underlying data type
      typedef EnumData< t_DataType> EnumDataType;

      //! @brief EnumType (pointer to the data type)
      typedef Enum< t_DataType, t_Derived> EnumType;

      //! @brief const iterator on the Enum vector
      typedef typename std::vector< EnumType>::const_iterator const_iterator;

      //! @brief iterator on the Enum vector
      typedef typename std::vector< EnumType>::iterator iterator;

    private:

    //////////
    // data //
    //////////

      //! this list contains all EnumData
      std::list< EnumDataType> m_DataList;

      //! this list contains a copy of all enums, pointing to their according EnumData in the m_DataList
      std::vector< EnumType> m_EnumVector;

    public:

      //! this is an instance of the UndefinedEnum
      EnumType e_Undefined;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    protected:

      //! @brief construct Enums - only used by derived classes
      //! constructs a default constructed Enum as very first and undefined Enum
      //! @param IS_READ_WRITABLE declare that Enum derived instance as read and writable for the user over the
      //!        EnumsInstances and ObjectInstances mechanisms
      Enumerate( const bool IS_READ_WRITABLE = true);

    public:

      //! virtual copy constructor
      Enumerate< t_DataType, t_Derived> *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief access to the only instance of the derived enumerators
      //! @return reference to singleton instance of t_Derived
      static t_Derived &GetEnums();

      //! @brief get strings of all enums in a vector
      //! @return std::vector filled with strings of all enums, sorted by enum index
      std::vector< std::string> GetEnumStrings() const;

    ////////////////
    // operations //
    ////////////////

      //! protected constructor, no creation of objects outside of the class
      //! @param NAME   name of the current enum
      //! @param OBJECT object to be enumerated
      virtual EnumType &AddEnum
      (
        const std::string &NAME,
        const t_DataType &OBJECT
      );

      //! begin iterator for iterating over enums
      //! @return const_iterator to begin of all enums
      const_iterator Begin() const;

      //! begin iterator for iterating over enums
      //! @return iterator to begin of all enums
      iterator Begin();

      //! end iterator for iterating over enums
      //! @return const_iterator to end of all enums
      const_iterator End() const;

      //! end iterator for iterating over enums
      //! @return iterator to end of all enums
      iterator End();

      //! get number of enum objects
      //! @return number of enums
      size_t GetEnumCount() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief get iterator to enum object from index of enum
      //! @param INDEX index of the enum to be searched for
      //! @return iterator to the enum with the INDEX - will point to end if there is no such enum with the INDEX
      iterator GetEnumIteratorFromIndex( const size_t INDEX);

      //! @brief get const_iterator to enum object from index of enum
      //! @param INDEX index of the enum to be searched for
      //! @return iterator to the enum with the INDEX - will point to end if there is no such enum with the INDEX
      const_iterator GetEnumIteratorFromIndex( const size_t INDEX) const;

      //! @brief get iterator to enum from name
      //! @param NAME name of the enum to be searched for
      //! @return iterator to the enum with the NAME - will point to end if there is no such enum with the NAME
      iterator GetEnumIteratorFromName( const std::string &NAME);

      //! @brief get const_iterator to enum from name
      //! @param NAME name of the enum to be searched for
      //! @return iterator to the enum with the NAME - will point to end if there is no such enum with the NAME
      const_iterator GetEnumIteratorFromName( const std::string &NAME) const;

      //! @brief set an enum to the enum with a given name
      //! @param ENUM enum to set
      //! @param NAME name of the enum to be searched for
      //! @param ERR_STREAM error output stream, for if name is not an enum
      //! @return true on success
      bool SetEnumFromName( EnumType &ENUM, const std::string &NAME, std::ostream &ERR_STREAM) const;

      //! @brief Get an undefined EnumData of this type
      //! @return const reference to undefined Data
      static const EnumDataType &GetUndefinedData();

      //! @brief Get an undefined Enum of this type
      //! @brief return reference to undefined enum, which is e_Undefined
      const EnumType &GetUndefined() const;

      //! @brief get iterator to enum from enum instance
      //! @param ENUM enum to be searched for
      //! @return iterator to the provided enum
      iterator GetEnumIteratorFromEnum( const EnumType &ENUM);

      //! @brief get iterator to enum from enum instance
      //! @param ENUM enum to be searched for
      //! @return const_iterator to the provided enum
      const_iterator GetEnumIteratorFromEnum( const EnumType &ENUM) const;

      //! @brief get enum from index
      //! @param INDEX index of the enum
      //! @return const enum reference with the index, undefined if non such enum exists
      const EnumType &GetEnumFromIndex( const size_t INDEX) const;

      //! @brief get enum from index
      //! @param INDEX index of the enum
      //! @return enum reference with the index, undefined if non such enum exists
      EnumType &GetEnumFromIndex( const size_t INDEX);

      //! @brief get enum from name
      //! @param NAME name of the enum to be searched for
      //! @return the const reference to enum with the NAME - will give undefined enum if no enum with NAME was found
      const EnumType &GetEnumFromName( const std::string &NAME) const;

      //! @brief get enum from name
      //! @param NAME name of the enum to be searched for
      //! @return the reference to enum with the NAME - will give undefined enum if no enum with NAME was found
      EnumType &GetEnumFromName( const std::string &NAME);

      //! @brief check if an enum with that INDEX exists
      //! @param INDEX index of suspected enum
      //! @return true if an enum with INDEX exists, false if there is no such enum or if INDEX is undefined
      bool HaveEnumWithIndex( const size_t INDEX) const;

      //! @brief check if an enum with the NAME exists
      //! @param NAME name of suspected enum
      //! @return true if an enum with NAME exists alos, if it is the name of the undefined enum, false otherwise
      bool HaveEnumWithName( const std::string &NAME) const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief writes the list of enums
      //! @param OSTREAM the stream to which the help is written to
      //! @return the given stream to which the list was written to
      //! Virtual to allow derived classes alter how the help is displayed without overriding Enum
      virtual std::ostream &WriteList( std::ostream &OSTREAM) const;

    protected:

      //! @brief write Enums to std::ostream
      //! @param OSTREAM ostream to write to
      //! @param INDENT the indent to be used for each line
      //! @return ref to the stream written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! @brief read Enums from io::IFStream
      //! @param ISTREAM istream to read from
      //! @return ref to the stream written to
      std::istream &Read( std::istream &ISTREAM);

    private:

      // To prevent copying Enumerate classes, declare a copy constructor and an assignment operator
      // by declaring them private, there is no way to copy construct Enumerate-derived classes

      //! @brief Undefined copy constructor, used to prevent copying this class
      Enumerate( const Enumerate &);

      //! Undefined assignment operator, used to prevent copying this class
      Enumerate &operator =( const Enumerate &);

    }; // template class Enumerate

  } // namespace util
} // namespace bcl

#endif // BCL_UTIL_ENUMERATE_H_
