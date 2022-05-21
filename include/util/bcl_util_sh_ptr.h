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

#ifndef BCL_UTIL_SH_PTR_H_
#define BCL_UTIL_SH_PTR_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_message.h"
#include "bcl_util_object_instances.h"
#include "bcl_util_own_ptr.h"
#include "type/bcl_type_enable_if.h"
#include "type/bcl_type_is_a.h"

// external includes - sorted alphabetically
#include <memory>

namespace bcl
{
  namespace util
  {
    //! returns hard copy of object - if it is a null pointer it just returns a null pointer
    //! This implementation handles BCL objects
    template< typename t_DataType>
    typename type::EnableIf< type::IsA< t_DataType, ObjectInterface>::value, ShPtr< t_DataType> >::Type
      PtrHardCopy( const t_DataType *const &PTR);

    //! returns hard copy of object - if it is a null pointer it just returns a null pointer
    //! This implementation handles non-bcl objects (e.g. double)
    template< typename t_DataType>
    typename type::EnableIf< !type::IsA< t_DataType, ObjectInterface>::value, ShPtr< t_DataType> >::Type
      PtrHardCopy( const t_DataType *const &PTR);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ShPtr
    //! @brief This class is a template for ShPtrs
    //! @details This pointer is only to be used as an ownership pointer. In turn it should always obtain ownership
    //! of the object at the point of initialization. Hence, if ShPtr is initialized from a normal pointer
    //! ALWAYS use ShPtr< t_DataType> OBJECT( new t_DataType).
    //! Note that the ShPtr will delete the OBJECT for you at the end of its scope. so NEVER
    //! call delete( OBJECT). Also, do not attempt to create a ShPtr from an existing OBJECT, *OBJECT,
    //! or an existing SiPtr< OBJECT>.
    //!
    //! @see @link example_util_sh_ptr.cpp @endlink
    //! @author meilerj, woetzen, mendenjl
    //! @date 02.09.2005
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class ShPtr :
      public PtrInterface< t_DataType>
    {

    /////////////
    // friends //
    /////////////

      template< typename t_OtherDataType> friend class ShPtr;

    private:

    //////////
    // data //
    //////////

      std::shared_ptr< t_DataType> m_Pointer; //!< pointer to object

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      ShPtr() :
        m_Pointer( nullptr)
      {
      }

      //! construct from pointer, always pass new t_DataType( POINTER) to this function (explicit assures usage of new)
      explicit ShPtr( t_DataType *const POINTER) :
        m_Pointer( POINTER)
      {
      }

      //! construct from pointer, always pass new t_OtherDataType( POINTER) to this function
      template< typename t_OtherDataType>
      explicit ShPtr( t_OtherDataType *const POINTER) :
        m_Pointer( PtrInterface< t_DataType>::CheckValidPointerConversion( POINTER))
      {
      }

      //! copy constructor from ShPtr
      ShPtr( const ShPtr< t_DataType> &SH_PTR) :
        m_Pointer( SH_PTR.m_Pointer)
      {
      }

      //! copy constructor from ShPtr< t_OtherDataType
      template< typename t_OtherDataType>
      ShPtr( const ShPtr< t_OtherDataType> &SH_PTR) :
        m_Pointer( PtrInterface< t_DataType>::CheckValidSmartPointerConversion( SH_PTR.m_Pointer))
      {
      }

      //! clone
      ShPtr< t_DataType> *Clone() const
      {
        return new ShPtr< t_DataType>( *this);
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

    ///////////////
    // operators //
    ///////////////

      //! soft copy, returns new ShPtr
      ShPtr< t_DataType> &operator =( const ShPtr< t_DataType> &SH_PTR)
      {
        // assign
        m_Pointer = SH_PTR.m_Pointer;

        // return *this
        return *this;
      }

      //! soft copy, returns *this if conversion succeeded
      template< typename t_OtherDataType>
      ShPtr< t_DataType> &operator =( const ShPtr< t_OtherDataType> &SH_PTR)
      {
        return operator =( ShPtr< t_DataType>( SH_PTR));
      }

      //! returns a constant reference to the undefined state for an ShPtr of t_DataType
      void Reset()
      {
        m_Pointer.reset();
      }

      //! returns hard copy of object - if it is a null pointer it just returns a null pointer
      //! This implementation handles BCL objects
      ShPtr< t_DataType> HardCopy() const
      {
        return PtrHardCopy( m_Pointer.get());
      }

    /////////////////
    // data access //
    /////////////////

      //! overload operator () to return const pointer on object t_DataType
      const t_DataType *GetPointer() const
      {
        return m_Pointer.get();
      }

      //! overload operator () to return pointer on object t_DataType
      t_DataType *GetPointer()
      {
        return m_Pointer.get();
      }

      //! use GetConstPointer to return const pointer to an object in a function that would otherwise
      //! allow it to be modified
      const t_DataType *GetConstPointer() const
      {
        this->AssertIsDefined();
        return m_Pointer.get();
      }

      //! overload operator -> to return const pointer on object t_DataType
      const t_DataType *operator ->() const
      {
        this->AssertIsDefined();
        return m_Pointer.get();
      }

      //! overload operator -> to return changeable pointer on object t_DataType
      t_DataType *operator ->()
      {
        this->AssertIsDefined();
        return m_Pointer.get();
      }

      //! overload operator * to return const reference of object t_DataType
      const t_DataType &operator *() const
      {
        this->AssertIsDefined();
        return *m_Pointer;
      }

      //! overload operator * to return reference to changeable object t_DataType
      t_DataType &operator *()
      {
        this->AssertIsDefined();
        return *m_Pointer;
      }

      //! check whether ShPtr is different from NULL
      bool IsDefined() const
      {
        return m_Pointer != nullptr;
      }

      //! Count the number of ShPtrs that also link to this particular ShPtr
      size_t GetSharedCommunitySize() const
      {
        return m_Pointer.use_count();
      }

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief returns a ShPtr to an object from a stream, that only contains an object
      //! @param ISTREAM the stream, with the object identifier and the content of interest
      static ShPtr< t_DataType> GetShPtrToNewObjectFromStream( std::istream &ISTREAM)
      {
        ShPtr< t_DataType> new_ptr;
        GetPtrToNewObjectFromStream( ISTREAM, new_ptr);
        return new_ptr;
      }

    protected:

      //! reads Object in ShPtr from istream
      std::istream &Read( std::istream &ISTREAM)
      {
        // reset this ShPtr
        Reset();

        // read address
        size_t address( 0);
        io::Serialize::Read( address, ISTREAM);

        if( address)
        {
          // read the object in from the stream
          GetPtrToNewObjectFromStream( ISTREAM, *this);
        }

        // end
        return ISTREAM;
      }

      //! writes ShPtr and Object to ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        WriteShPtr( m_Pointer.get(), OSTREAM, INDENT);
        //end
        return OSTREAM;
      }

    private:

      template< typename t_OtherDataType>
      typename type::EnableIf< type::IsDerivedFrom< t_OtherDataType, ObjectInterface>::value, std::ostream &>::Type
      WriteShPtr( t_OtherDataType *PTR, std::ostream &OSTREAM, const size_t INDENT) const
      {
        // call the object instances version of write
        GetObjectInstances().WriteShPtr( PTR, OSTREAM, INDENT);
        return OSTREAM;
      }

      template< typename t_OtherDataType>
      typename type::EnableIf< !type::IsDerivedFrom< t_OtherDataType, ObjectInterface>::value, std::ostream &>::Type
      WriteShPtr( t_OtherDataType *PTR, std::ostream &OSTREAM, const size_t INDENT) const
      {
        GetObjectInstances().WritePtrAddress( size_t( PTR), OSTREAM, INDENT);
        if( PTR)
        {
          io::Serialize::Write( *PTR, OSTREAM, INDENT);
        }
        return OSTREAM;
      }

      //! @brief initializes a passed in pointer given a stream
      //! @param ISTREAM the stream, with the object identifier and the content of interest
      //! @param PTR pointer to initialize
      //! @note overload for non-object interface derived types
      template< typename t_OtherDataType>
      static typename type::EnableIf< type::IsA< t_OtherDataType, ObjectInterface>::value, void>::Type
      GetPtrToNewObjectFromStream( std::istream &ISTREAM, ShPtr< t_OtherDataType> &PTR)
      {
        // read identifier for object
        const std::string identifier( ObjectInterface::ExtractIdentifier( ISTREAM));

        if( identifier == GetNullDescriptor())
        {
          PTR = ShPtr< t_OtherDataType>();
          return;
        }

        // try to case the object in the object instances map to t_DataType
        PTR =
          ShPtr< t_OtherDataType>
          (
            PtrInterface< t_DataType>::CheckValidPointerConversion
            (
              GetObjectInstances().GetNewPtrToObjectFromIdentifier( identifier)
            )
          );

        //read object
        BCL_Assert
        (
          PTR.IsDefined(),
          "requested object of name \"" + identifier +
          "\" is not a derived class of " + GetStaticClassName< t_DataType>()
        );

        PTR->ReadWithOutIdentifier( ISTREAM);
      }

      //! @brief initializes a passed in pointer given a stream
      //! @param ISTREAM the stream, with the object identifier and the content of interest
      //! @param PTR pointer to initialize
      //! @note overload for non-object interface derived types
      template< typename t_OtherDataType>
      static typename type::EnableIf< !type::IsA< t_OtherDataType, ObjectInterface>::value, void>::Type
      GetPtrToNewObjectFromStream( std::istream &ISTREAM, ShPtr< t_OtherDataType> &PTR)
      {
        // read identifier for object
        const std::string identifier( ObjectInterface::ExtractIdentifier( ISTREAM));

        if( identifier == GetNullDescriptor())
        {
          PTR = ShPtr< t_OtherDataType>();
        }
        else
        {
          PTR = ShPtr< t_OtherDataType>( new t_OtherDataType);
          std::stringstream stream( identifier);
          stream >> *PTR;
        }
      }

    }; // template class ShPtr

    //! returns hard copy of object - if it is a null pointer it just returns a null pointer
    //! This implementation handles BCL objects
    template< typename t_DataType>
    typename type::EnableIf< type::IsA< t_DataType, ObjectInterface>::value, ShPtr< t_DataType> >::Type
      PtrHardCopy( const t_DataType *const &PTR)
    {
      return ShPtr< t_DataType>( PTR ? PTR->Clone() : nullptr);
    }

    //! returns hard copy of object - if it is a null pointer it just returns a null pointer
    //! This implementation handles non-bcl objects (e.g. double)
    template< typename t_DataType>
    typename type::EnableIf< !type::IsA< t_DataType, ObjectInterface>::value, ShPtr< t_DataType> >::Type
      PtrHardCopy( const t_DataType *const &PTR)
    {
      return ShPtr< t_DataType>( PTR ? new t_DataType( *PTR) : nullptr);
    }

    //! @brief clone the given object into a ShPtr
    //! this can be used to create a clone copy of an object hold by a shared pointer, if a function requires a ShPtr
    //! to an object but one only has the reference to one - it is as good as writing ShPtr< t_DataType>( OBJECT.Clone())
    //! but it does not require to know t_DataType
    //! @param OBJECT the object to clone
    //! @return ShPtr to a new object, which is a of OBJECT
    template< typename t_DataType>
    ShPtr< t_DataType> CloneToShPtr( const t_DataType &OBJECT)
    {
      return ShPtr< t_DataType>( OBJECT.Clone());
    }

    BCL_EXPIMP_TEMPLATE template class BCL_API ShPtr< ObjectInterface>;

  } // namespace util
} // namespace bcl

#endif //BCL_UTIL_SH_PTR_H_
