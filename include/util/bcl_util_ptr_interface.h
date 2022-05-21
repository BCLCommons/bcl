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

#ifndef BCL_UTIL_PTR_INTERFACE_H_
#define BCL_UTIL_PTR_INTERFACE_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_assert.h"
#include "bcl_util_class_descriptor.h"
#include "bcl_util_object_interface.h"
#include "type/bcl_type_enable_if.h"
#include "type/bcl_type_is_a.h"

// external includes - sorted alphabetically
#include <memory>

namespace bcl
{
  namespace util
  {
    //! @brief Writes out a warning message that a pointer cast failed between the two types
    //! @param PTR_TYPE the pointer type that was received
    //! @param CAST_TYPE the type the pointer was cast to
    //! @note this function provides a convenient hook for debuggers, which need only set a single breakpoint in the cpp
    BCL_API void NotifyUserBadPointerCast( const std::string &PTR_TYPE, const std::string &CAST_TYPE);

    //! @brief get the string for an empty pointer, e.g. "NULL"
    //! @return the string for an empty pointer, e.g. "NULL"
    BCL_API const std::string &GetNullDescriptor();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PtrInterface
    //! @brief derives all bcl pointer classes
    //! @details This interface provides basic functions, that are necessary for pointer classes, like operator-> and ptr comparisons
    //!
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date May 20, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class PtrInterface :
      public ObjectInterface
    {

    public:

      typedef t_DataType &reference_type;
      typedef t_DataType value_type;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new PtrInterface
      virtual PtrInterface< t_DataType> *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief return the pointer
      //! every object in the bcl should be derived from ObjectInterface, so this function returns the address to that
      //! ObjectInterface
      //! @return pointer to ObjectInterface
      virtual const t_DataType *GetPointer() const = 0;

      //! @brief return the const pointer
      //! every object in the bcl should be derived from ObjectInterface, so this function returns the address to that
      //! ObjectInterface
      //! @return pointer to const ObjectInterface
      virtual t_DataType *GetPointer() = 0;

    ////////////////
    // operations //
    ////////////////

      //! check whether Ptr is different from NULL
      virtual bool IsDefined() const = 0;

      //! @brief assert that the pointer does not point to NULL
      void AssertIsDefined() const
      {
        BCL_Assert( IsDefined(), "Invalid operation for NULL pointer of type: " + GetStaticClassName< t_DataType>());
      }

      //! @brief checks that a dynamic cast of a pointer to t_DataType succeeds
      //! @param PTR the pointer to try to convert into t_DataType
      template< typename t_OtherDataType>
      static typename type::EnableIf< type::IsA< t_OtherDataType, t_DataType>::value, t_DataType *>::Type
        CheckValidPointerConversion( t_OtherDataType *const PTR)
      {
        return static_cast< t_DataType *>( PTR);
      }

      //! @brief checks that a dynamic cast of a pointer to t_DataType succeeds
      //! @param PTR the pointer to try to convert into t_DataType
      template< typename t_OtherDataType>
      static typename type::EnableIf< !type::IsA< t_OtherDataType, t_DataType>::value, t_DataType *>::Type
        CheckValidPointerConversion( t_OtherDataType *const PTR)
      {
        t_DataType *cast_ptr( dynamic_cast< t_DataType *>( PTR));

        // check whether the conversion worked.
        // Invalid conversions cause cast_ptr to be set to NULL, but if it was NULL to begin with, it succeeds
        if( PTR != NULL && cast_ptr == NULL)
        {
          NotifyUserBadPointerCast( GetStaticClassName< t_OtherDataType>(), GetStaticClassName< t_DataType>());
        }
        return cast_ptr;
      }

      //! @brief checks that a dynamic cast of a pointer to t_DataType succeeds
      //! @param PTR the pointer to try to convert into t_DataType
      template< template< typename> class t_PtrType, typename t_OtherDataType>
      static typename type::EnableIf< type::IsA< t_OtherDataType, t_DataType>::value, t_PtrType< t_DataType> >::Type
        CheckValidSmartPointerConversion( const t_PtrType< t_OtherDataType> &PTR)
      {
        return std::static_pointer_cast< t_DataType>( PTR);
      }

      //! @brief checks that a dynamic cast of a pointer to t_DataType succeeds
      //! @param PTR the pointer to try to convert into t_DataType
      template< template< typename> class t_PtrType, typename t_OtherDataType>
      static typename type::EnableIf< !type::IsA< t_OtherDataType, t_DataType>::value, t_PtrType< t_DataType> >::Type
        CheckValidSmartPointerConversion( const t_PtrType< t_OtherDataType> &PTR)
      {
        t_PtrType< t_DataType> cast_ptr( std::dynamic_pointer_cast< t_DataType>( PTR));

        // check whether the conversion worked.
        // Invalid conversions cause cast_ptr to be set to NULL, but if it was NULL to begin with, it succeeds
        if( PTR != nullptr && cast_ptr == nullptr)
        {
          NotifyUserBadPointerCast( GetStaticClassName< t_OtherDataType>(), GetStaticClassName< t_DataType>());
        }
        return cast_ptr;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator ->
      //! @return pointer to ObjectInterface
      virtual const t_DataType *operator->() const = 0;

      //! overload operator -> to return pointer to changeable object t_DataType
      virtual t_DataType *operator->() = 0;

      //! overload operator * to return const reference of object t_DataType
      virtual const t_DataType &operator *() const = 0;

      //! overload operator * to reference to changeable object t_DataType
      virtual t_DataType &operator *() = 0;

    }; // template class PtrInterface

    //! @brief perform equal comparison of two Ptr
    //! @param POINTER_LEFT  left hand argument
    //! @param POINTER_RIGHT right hand argument
    //! @return bool - true if addresses are equal, false if addresses do not match
    template< typename t_DataTypeRight, typename t_DataTypeLeft>
    inline bool operator ==
    (
      const PtrInterface< t_DataTypeRight> &POINTER_LEFT,
      const PtrInterface< t_DataTypeLeft> &POINTER_RIGHT
    )
    {
      return POINTER_LEFT.GetPointer() == POINTER_RIGHT.GetPointer();
    }

    //! @brief perform not equal comparison of two Ptr
    //! @param POINTER_LEFT  left hand argument
    //! @param POINTER_RIGHT right hand argument
    //! @return true iff addresses are not equal or either pointer is undefined
    //! This emulates the result of comparing double NaNs
    template< typename t_DataTypeRight, typename t_DataTypeLeft>
    inline bool operator !=
    (
      const PtrInterface< t_DataTypeRight> &POINTER_LEFT,
      const PtrInterface< t_DataTypeLeft> &POINTER_RIGHT
    )
    {
      return POINTER_LEFT.GetPointer() != POINTER_RIGHT.GetPointer();
    }

    //! @brief perform smaller than  comparison of two Ptr
    //! @param POINTER_LEFT  left hand argument
    //! @param POINTER_RIGHT right hand argument
    //! @return bool - true if addresses of left pointer is smaller than addresses of right pointer
    template< typename t_DataTypeRight, typename t_DataTypeLeft>
    inline bool operator <
    (
      const PtrInterface< t_DataTypeRight> &POINTER_LEFT,
      const PtrInterface< t_DataTypeLeft> &POINTER_RIGHT
    )
    {
      return POINTER_LEFT.GetPointer() < POINTER_RIGHT.GetPointer();
    }

    //! @brief perform larger than  comparison of two Ptr
    //! @param POINTER_LEFT  left hand argument
    //! @param POINTER_RIGHT right hand argument
    //! @return bool - false if addresses of left pointer is smaller or equal than addresses of right pointer
    template< typename t_DataTypeRight, typename t_DataTypeLeft>
    inline bool operator >
    (
      const PtrInterface< t_DataTypeRight> &POINTER_LEFT,
      const PtrInterface< t_DataTypeLeft> &POINTER_RIGHT
    )
    {
      return POINTER_LEFT.GetPointer() > POINTER_RIGHT.GetPointer();
    }

    //! @brief perform smaller equal than  comparison of two Ptr
    //! @param POINTER_LEFT  left hand argument
    //! @param POINTER_RIGHT right hand argument
    //! @return bool - true if addresses of left pointer is smaller or equal than addresses of right pointer
    template< typename t_DataTypeRight, typename t_DataTypeLeft>
    inline bool operator <=
    (
      const PtrInterface< t_DataTypeRight> &POINTER_LEFT,
      const PtrInterface< t_DataTypeLeft> &POINTER_RIGHT
    )
    {
      return POINTER_LEFT.GetPointer() <= POINTER_RIGHT.GetPointer();
    }

    //! @brief perform larger equal than  comparison of two Ptr
    //! @param POINTER_LEFT  left hand argument
    //! @param POINTER_RIGHT right hand argument
    //! @return bool - false if addresses of left pointer is smaller than addresses of right pointer
    template< typename t_DataTypeRight, typename t_DataTypeLeft>
    inline bool operator >=
    (
      const PtrInterface< t_DataTypeRight> &POINTER_LEFT,
      const PtrInterface< t_DataTypeLeft> &POINTER_RIGHT
    )
    {
      return POINTER_LEFT.GetPointer() >= POINTER_RIGHT.GetPointer();
    }

    //! @brief check if prt is defined
    //! @param PTR to be checked if it is defined
    //! @return true if pointer is defined, false otherwise
    template< typename t_DataType>
    inline bool IsDefined( const PtrInterface< t_DataType> &PTR)
    {
      return PTR.IsDefined();
    }

  } // namespace util
} // namespace bcl

#endif // BCL_UTIL_PTR_INTERFACE_H_
