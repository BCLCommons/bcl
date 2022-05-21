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

#ifndef BCL_UTIL_SI_PTR_H_
#define BCL_UTIL_SI_PTR_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_ptr_interface.h"
#include "io/bcl_io_serialize.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SiPtr
    //! @brief This class is a template for SiPtrs
    //! @details SiPtr< t_DataType> is a non-ownership pointer. It should be used to pass arguments to functions. It should
    //! completely replace *t_DataType to avoid confusion about ownership of t_DataType.
    //! It can be created from a ShPtr< t_DataType>, from a SiPtr< t_DataType> or from t_DataType &. Note that
    //! SiPtr< t_DataType> is NOT for storing t_DataType (use ShPtr< t_DataType>).
    //! SiPtr< t_DataType> are ONLY for temporary use! They will be undefined if OBJECT gets deleted.
    //! Since SiPtr< t_DataType> is not supposed to own an t_DataType, no hard_copy operations are implemented.
    //!
    //! @see @link example_util_si_ptr.cpp @endlink
    //! @author meilerj, woetzen
    //! @date 02.09.2005
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class SiPtr :
      public PtrInterface< t_DataType>
    {
    /////////////
    // friends //
    /////////////

      template< typename t_OtherDataType> friend class SiPtr;

    private:

    //////////
    // data //
    //////////

      t_DataType *m_Pointer; //!< pointer to object

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      SiPtr() :
        m_Pointer( NULL)
      {
      }

      //! construct from reference to t_DataType
      SiPtr( t_DataType &OBJECT) :
        m_Pointer( &OBJECT)
      {
      }

      //! construct from t_DataType pointer
      explicit SiPtr( t_DataType *const POINTER) :
        m_Pointer( POINTER)
      {
      }

      //! construct from t_OtherDataType pointer
      template< typename t_OtherDataType>
      explicit SiPtr( t_OtherDataType *const POINTER) :
        m_Pointer( PtrInterface< t_DataType>::CheckValidPointerConversion( POINTER))
      {
      }

      //! construct from SiPtr< t_OtherDataType>
      template< typename t_OtherDataType>
      explicit SiPtr( const SiPtr< t_OtherDataType> &POINTER) :
        m_Pointer( PtrInterface< t_DataType>::CheckValidPointerConversion( POINTER.m_Pointer))
      {
      }

      //! construct from PtrInterface< t_DataType>
      SiPtr( PtrInterface< t_DataType> &POINTER) :
        m_Pointer( POINTER.GetPointer())
      {
      }

      //! construct from PtrInterface< t_DataType>
      SiPtr( const PtrInterface< t_DataType> &POINTER) :
        m_Pointer( CheckValidPointerConversion( POINTER.GetPointer()))
      {
      }

      //! construct from PtrInterface< t_OtherDataType>
      template< typename t_OtherDataType>
      SiPtr( PtrInterface< t_OtherDataType> &POINTER) :
        m_Pointer( PtrInterface< t_DataType>::CheckValidPointerConversion( POINTER.GetPointer()))
      {
      }

      //! construct from PtrInterface< t_OtherDataType>
      template< typename t_OtherDataType>
      SiPtr( const PtrInterface< t_OtherDataType> &POINTER) :
        m_Pointer( PtrInterface< t_DataType>::CheckValidPointerConversion( POINTER.GetPointer()))
      {
      }

      //! copy constructor from SiPtr< t_DataType>
      SiPtr( const SiPtr< t_DataType> &POINTER) :
        m_Pointer( POINTER.m_Pointer)
      {
      }

      //! clone
      SiPtr< t_DataType> *Clone() const
      {
        return new SiPtr< t_DataType>( *this);
      }

    ///////////////
    // operators //
    ///////////////

      //! allows to set SiPtr from t_DataType pointer
      SiPtr< t_DataType> &operator  =( t_DataType *const POINTER)
      {
        m_Pointer = POINTER;
        return *this;
      }

      //! allows to set SiPtr from t_OtherDataType pointer
      template< typename t_OtherDataType>
      SiPtr< t_DataType> &operator  =( t_OtherDataType *const POINTER)
      {
        // assign m_Pointer
        m_Pointer = PtrInterface< t_DataType>::CheckValidPointerConversion( POINTER);
        return *this;
      }

      //! allows to set SiPtr< t_DataType> from SiPtr< t_DataType>
      SiPtr< t_DataType> &operator =( const SiPtr< t_DataType> &POINTER)
      {
        m_Pointer = POINTER.m_Pointer;
        return *this;
      }

      //! allows to set SiPtr< t_DataType> from SiPtr< t_OtherDataType>
      template< typename t_OtherDataType>
      SiPtr< t_DataType> &operator =( const SiPtr< t_OtherDataType> &POINTER)
      {
        // assign m_Pointer
        m_Pointer = PtrInterface< t_DataType>::CheckValidPointerConversion( POINTER.m_Pointer);
        return *this;
      }

      //! allows to set SiPtr< t_DataType> from PtrInterface< t_DataType>
      SiPtr< t_DataType> &operator =( const PtrInterface< t_DataType> &POINTER)
      {
        m_Pointer = POINTER.GetPointer();
        return *this;
      }

      //! allows to set SiPtr< t_DataType> from PtrInterface< t_OtherDataType>
      template< typename t_OtherDataType>
      SiPtr< t_DataType> &operator =( const PtrInterface< t_OtherDataType> &POINTER)
      {
        // assign m_Pointer
        m_Pointer = PtrInterface< t_DataType>::CheckValidPointerConversion( POINTER.GetPointer());
        return *this;
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! overload GetPointer to return const pointer on object t_DataType
      const t_DataType *GetPointer() const
      {
        return m_Pointer;
      }

      //! overload GetPointer to return pointer on object t_DataType
      t_DataType *GetPointer()
      {
        return m_Pointer;
      }

      //! overload operator -> to return const pointer on object t_DataType
      const t_DataType *operator->() const
      {
        this->AssertIsDefined();
        return m_Pointer;
      }

      //! overload operator -> to return pointer to changeable object t_DataType
      t_DataType       *operator->()
      {
        this->AssertIsDefined();
        return m_Pointer;
      }

      //! overload operator * to return const reference of object t_DataType
      const t_DataType &operator *() const
      {
        this->AssertIsDefined();
        return *m_Pointer;
      }

      //! overload operator * to reference to changeable object t_DataType
      t_DataType       &operator *()
      {
        this->AssertIsDefined();
        return *m_Pointer;
      }

      //! @brief check whether SiPtr is different from NULL
      //! @return true, if pointer != NULL
      bool IsDefined() const
      {
        return ( m_Pointer != NULL);
      }

    ////////////////
    // operations //
    ////////////////

      //! sets the pointer to NULL
      void Reset()
      {
        m_Pointer = NULL;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! reads Object in SiPtr from istream
      std::istream &Read( std::istream &ISTREAM)
      {
        BCL_Assert( IsDefined(), "Cannot read into a null si-ptr pointer!");
        return io::Serialize::Read( *m_Pointer, ISTREAM);
      }

      //! writes SiPtr and Object to ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        if( !IsDefined())
        {
          // write a Null descriptor
          io::Serialize::InsertIndent( OSTREAM, INDENT);
          OSTREAM << GetNullDescriptor();
          return OSTREAM;
        }
        return io::Serialize::Write( *m_Pointer, OSTREAM, INDENT);
      }

    }; // template class SiPtr

    //! redirect reading to reading function of object
    template< typename t_DataType>
    inline
    std::istream &
    operator >>( std::istream &ISTREAM, SiPtr< const t_DataType> &)
    {
      // it was necessary to disable this function, because if you would like to have a SiPtr< const t_DataType> then
      // there is no function that reads to a const argument, so this is generally not possible
      // there is also danger in reading to a NULL pointer, since the object is there anywhere, you should use the >> operator
      // for the object directly and dereference the SiPtr
      BCL_Exit( "not possible to read to a SiPtr< const t_DataType> because you are not allowed to alter t_DataType!", -1);
      return ISTREAM;
    }

    //! redirect reading to reading function of object
    template< typename t_DataType>
    inline
    std::istream &
    operator >>( std::istream &ISTREAM, const SiPtr< const t_DataType> &)
    {
      // it was necessary to disable this function, because if you would like to have a SiPtr< const t_DataType> then
      // there is no function that reads to a const argument, so this is generally not possible
      // there is also danger in reading to a NULL pointer, since the object is there anywhere, you should use the >> operator
      // for the object directly and dereference the SiPtr
      BCL_Exit( "not possible to read to a const SiPtr< const t_DataType> because you are not allowed to alter t_DataType!", -1);
      return ISTREAM;
    }

    //! @brief function to create a (const) si ptr to reference an object; saves writing long typenames
    //! @param OBJECT the object to create a si ptr from
    //! @return a const si ptr to the object
    template< typename t_DataType>
    inline SiPtr< const t_DataType> ToSiPtr( const t_DataType &OBJECT)
    {
      return SiPtr< const t_DataType>( OBJECT);
    }

    //! @brief function to create a (const) si ptr to reference an object; saves writing long typenames
    //! @param OBJECT the object to create a si ptr from
    //! @return a const si ptr to the object
    template< typename t_DataType>
    SiPtr< const t_DataType> ToSiPtr( const ShPtr< const t_DataType> &OBJECT);

    //! @brief function to create a (non-const) si ptr to reference an object; saves writing long typenames
    //! @param OBJECT the object to create a si ptr from
    //! @return a si ptr to the object
    template< typename t_DataType>
    inline SiPtr< t_DataType> ToSiPtrNonConst( t_DataType &OBJECT)
    {
      return SiPtr< t_DataType>( OBJECT);
    }
  } // namespace util
} // namespace bcl

#endif //BCL_UTIL_SI_PTR_H_
