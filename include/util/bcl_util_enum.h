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

#ifndef BCL_UTIL_ENUM_H_
#define BCL_UTIL_ENUM_H_

// include the namespace header
#include "bcl_util.h"

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

namespace bcl
{
  namespace util
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Enum
    //! @brief Enum points is a pointer class to EnumData that provides convenience functions to access the data.
    //! @details Besides providing access to the index and name to the enum data, it has operators to convert to string
    //!          or size_t and the iterator on the vector of Enums can be accessed.
    //!
    //! @see @link example_util_enum.cpp @endlink
    //! @author heinzes1, woetzen, karakam
    //! @date Nov 4, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_Derived>
    class Enum :
      public io::SerializationInterface
    {
    private:

      //! reference to associated object
      EnumData< t_DataType> *m_ObjectRef;

    public:

      //! @brief for iterator on the t_Derived Enum vector
      typedef typename std::vector< Enum< t_DataType, t_Derived> >::const_iterator const_iterator;

      //! @brief for iterator on the t_Derived Enum vector
      typedef typename std::vector< Enum< t_DataType, t_Derived> >::iterator iterator;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! @details points to the undefined enum
      Enum() :
        m_ObjectRef( &const_cast< EnumData< t_DataType> &>( EnumData< t_DataType>::GetUndefinedData()))
      {
      }

      //! @brief copy constructor
      //! @param ENUM the enum to be copied using its m_ObejectRef
      Enum( const Enum &ENUM) :
        m_ObjectRef( const_cast< EnumData< t_DataType> *>( ENUM.m_ObjectRef))
      {
      }

      //! @brief construct an undefined enum
      //! @param UNDEFINED an undefined object
      explicit Enum( const UndefinedObject UNDEFINED) :
        m_ObjectRef( &const_cast< EnumData< t_DataType> &>( EnumData< t_DataType>::GetUndefinedData()))
      {
      }

      //! @brief constructor from EnumData
      //! @param ENUM_DATA enum data object that this enum will point to
      Enum( EnumData< t_DataType> &ENUM_DATA) :
        m_ObjectRef( &ENUM_DATA)
      {
      }

      //! @brief construct enum from index
      //! @param INDEX index of enum to be constructed
      explicit Enum( const size_t &INDEX) :
        m_ObjectRef( t_Derived::GetEnums().GetEnumFromIndex( INDEX).m_ObjectRef)
      {
      }

      //! @brief construct enum from name
      //! @param NAME name of enum to be constructed
      explicit Enum( const std::string &NAME) :
        m_ObjectRef( t_Derived::GetEnums().GetEnumFromName( NAME).m_ObjectRef)
      {
      }

      //! @brief construct enum from name
      //! @param NAME name of enum to be constructed
      explicit Enum( const char *NAME) :
        m_ObjectRef( t_Derived::GetEnums().GetEnumFromName( std::string( NAME)).m_ObjectRef)
      {
      }

      //! copy function
      Enum *Clone() const
      {
        return new Enum( *this);
      }

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get index for that enum object
      //! @return index of current enum
      //! @return index of enum object, undefined size_t for undefined enum
      size_t GetIndex() const
      {
        return m_ObjectRef->m_Index;
      }

      //! @brief get name of enum, "Undefined name" for undefined enum
      //! @return name of current enum
      const std::string &GetName() const
      {
        return m_ObjectRef->m_Name;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief conversion of the enum to its index
      //! @return the index of this enum
      operator size_t() const
      {
        return m_ObjectRef->m_Index;
      }

      //! @brief conversion of the enum to its string nmae
      //! @return the string of this enum
      operator std::string() const
      {
        return m_ObjectRef->m_Name;
      }

      //! @brief convert the enum to the const enumerated data object
      //! @return the const data object for this enum
      operator const EnumData< t_DataType>() const
      {
        return *m_ObjectRef;
      }

      //! @brief convert the enum to the enumerated data object
      //! @return the data object for this enum
      operator EnumData< t_DataType>()
      {
        return *m_ObjectRef;
      }

      //! @brief operator -> that returns the const data object
      //! @return pointer to the const data object for convinient access to the data
      const EnumData< t_DataType> *operator ->() const
      {
        return m_ObjectRef;
      }

      //! @brief operator -> that returns the data object
      //! @return pointer to the data object for convinient access to the data
      t_DataType *operator ->()
      {
        return m_ObjectRef;
      }

      //! @brief operator * that returns the const data object
      //! @return const reference to the data object for access
      const t_DataType &operator *() const
      {
        return *m_ObjectRef;
      }

      //! @brief operator * that returns the data object
      //! @return reference to the data object for access
      t_DataType &operator *()
      {
        return *m_ObjectRef;
      }

      //! @brief access to the const_iterator for that enum
      //! @return const_iterator to the enum
      const_iterator GetIterator() const
      {
        return t_Derived::GetEnums().GetEnumIteratorFromIndex( m_ObjectRef->m_Index);
      }

      //! @brief access to the iterator for that enum
      //! @return iterator to the enum
      iterator GetIterator()
      {
        return t_Derived::GetEnums().GetEnumIteratorFromIndex( m_ObjectRef->m_Index);
      }

      //! @brief Get the label for an object of the given type
      //! @param WITH_DATA whether to include any data members, else, only include initialization members
      ObjectDataLabel GetLabel( const bool &WITH_DATA = false) const
      {
        return ObjectDataLabel( "", GetName());
      }

      //! @brief check if this enum is defined
      //! @return true if this enum is defined as indicated by a defined index
      bool IsDefined() const
      {
        return m_ObjectRef->m_Index != GetUndefined< size_t>();
      }

      //! @brief comparison operator equal
      //! @param RHS right hand side enum to compare to
      //! @return true if both enums point to the same object, false otherwise
      bool operator ==( const Enum &RHS) const
      {
        return m_ObjectRef == RHS.m_ObjectRef;
      }

      //! @brief comparison operator not equal
      //! @param RHS right hand side enum to compare to
      //! @return false if both enums point to the same object, true otherwise
      bool operator !=( const Enum &RHS) const
      {
        return m_ObjectRef != RHS.m_ObjectRef;
      }

      //! @brief comparison operator less than (<) for enums
      //! @param RHS right hand side enum to compare to
      //! @return true, of this index is smaller than RHS index
      bool operator <( const Enum &RHS) const
      {
        return m_ObjectRef->m_Index < RHS.m_ObjectRef->m_Index;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! write to std::ostream
      //! @param OSTREAM stream to write to
      //! @param INDENT indentation of output
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write m_Name
        io::Serialize::Write( m_ObjectRef->m_Name, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

      //! read from std::istream
      //! @param ISTREAM stream to read from
      std::istream &Read( std::istream &ISTREAM);
      // {
      //   // read name from stream
      //   std::string name;
      //   io::Serialize::Read( name, ISTREAM);

      //   // find the object belonging to the name and reassign the reference
      //   *this = t_Derived::GetEnums().GetEnumFromName( name);

      //   // end
      //   return ISTREAM;
      // }

    public:

      //! @brief set the value of the corresponding member based on the label
      //! @param LABEL label that is used to set the string
      //! @param ERROR_STREAM stream to write errors to
      //! @return bool - if the parameter was set (parameters are not set if they are not allowed parameters)
      bool TryRead( const ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);
      // {
      //   return t_Derived::GetEnums().SetEnumFromName( *this, LABEL.GetValue(), ERROR_STREAM);
      // }

      //! @brief writes the help for the label
      //! @param OSTREAM the stream to which the help is written to
      //! @param INDENT the amount of indent
      //! @return the given stream to which the help was written to
      std::ostream &WriteHelp( std::ostream &OSTREAM, const size_t INDENT = 0) const;
      // {
      //   OSTREAM << "Choose from the following: ";
      //   return t_Derived::GetEnums().WriteList( OSTREAM);
      // }
    }; // class Enum

    //! @brief overload IsDefiend for Enum derived classes
    //! @param ENUMERATOR an enum to be checked if is defined
    //! @return true if current ENUM is defined
    template< typename t_DataType, typename t_Derived>
    inline bool IsDefined( const Enum< t_DataType, t_Derived> &ENUMERATOR)
    {
      // return IsDefined for the argument ENUMERATOR
      return ENUMERATOR.IsDefined();
    }

  } // namespace util

  //! @brief overload < operator for Enum derived classes
  //! @param ENUMERATOR_L left handside enumerator
  //! @param ENUMERATOR_R right handside enumerator
  //! @return true if ENUMERATOR_L is smaller than ENUMERATOR_R by index
  template< typename t_DataType, typename t_Derived>
  inline bool operator <
  (
    const util::Enum< t_DataType, t_Derived> &ENUMERATOR_L,
    const util::Enum< t_DataType, t_Derived> &ENUMERATOR_R
  )
  {
    // return
    return ENUMERATOR_L.GetIndex() < ENUMERATOR_R.GetIndex();
  }

  //! @brief overload IsDefiend for Enum derived classes
  //! @param ENUMERATOR an enum to be checked if is defined
  //! @return true if current ENUM is defined
  template< typename t_DataType, typename t_Derived>
  inline bool IsDefined( const util::Enum< t_DataType, t_Derived> &ENUMERATOR)
  {
    // return IsDefined for the argument ENUMERATOR
    return ENUMERATOR.IsDefined();
  }

} // namespace bcl

// Include implementation header
#include "bcl_util_enum.hpp"

#endif //BCL_UTIL_ENUM_H_

