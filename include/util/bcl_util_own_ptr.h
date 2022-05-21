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

#ifndef BCL_UTIL_OWN_PTR_H_
#define BCL_UTIL_OWN_PTR_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_ptr_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class OwnPtr
    //! @brief generalization of a simple pointer.  OwnPtr can optionally "own" (be responsible for deleting) its pointer
    //!
    //! @see @link example_util_own_ptr.cpp @endlink
    //! @author woetzen
    //! @date 01.05.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class OwnPtr :
      public PtrInterface< t_DataType>
    {
    private:

    //////////
    // data //
    //////////

      t_DataType * m_Pointer; //!< pointer to an t_DataType
      bool m_IsOwner;         //!< if IsOwner is true, call to Delete will delete the pointer

      template< typename t_OtherDataType>
      friend class OwnPtr;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      OwnPtr() :
        m_Pointer( NULL),
        m_IsOwner( false)
      {
      }

      //! @brief construct from pointer to t_DataType and optionally indicate, if OwnPtr gets ownership
      //! @param OBJECT_PTR pointer to t_DataType
      //! @param OWNERSHIP true - has ownership and needs to delete pointer; false - no ownership, nothing to be done
      explicit OwnPtr( t_DataType *const OBJECT_PTR, const bool OWNERSHIP = true) :
        m_Pointer( OBJECT_PTR),
        m_IsOwner( OWNERSHIP)
      {
      }

      //! @brief construct from pointer to t_DataType and optionally indicate, if OwnPtr gets ownership
      //! @param OBJECT_PTR pointer to t_DataType
      //! @param OWNERSHIP true - has ownership and needs to delete pointer; false - no ownership, nothing to be done
      template< typename t_OtherDataType>
      explicit OwnPtr( t_OtherDataType *const OBJECT_PTR, const bool OWNERSHIP = true) :
        m_Pointer( PtrInterface< t_DataType>::CheckValidPointerConversion( OBJECT_PTR)),
        m_IsOwner( OWNERSHIP)
      {
      }

      //! @brief copy constructor
      OwnPtr( const OwnPtr< t_DataType> &OWN_PTR) :
        m_Pointer
        (
          OWN_PTR.m_IsOwner && OWN_PTR.IsDefined()
          ? OWN_PTR.m_Pointer->Clone()
          : OWN_PTR.m_Pointer
        ),
        m_IsOwner( OWN_PTR.m_IsOwner)
      {
      }

      //! @brief construct from OwnPtr to derived class
      template< typename t_OtherDataType>
      OwnPtr( const OwnPtr< t_OtherDataType> &OWN_PTR) :
        m_Pointer
        (
          PtrInterface< t_DataType>::CheckValidPointerConversion
          (
            OWN_PTR.m_IsOwner && OWN_PTR.IsDefined()
            ? OWN_PTR.m_Pointer->Clone()
            : OWN_PTR.m_Pointer
          )
        ),
        m_IsOwner( OWN_PTR.m_IsOwner)
      {
      }

      //! @brief virtual copy constructor
      OwnPtr< t_DataType> *Clone() const
      {
        return new OwnPtr< t_DataType>( *this);
      }

      //! @brief destructor
      virtual ~OwnPtr()
      {
        // only delete if OwnPtr is actually owner
        // deleting a NULL pointer is fine
        if( m_IsOwner)
        {
          delete m_Pointer;
        }
      }

      //! Assignment operator
      OwnPtr< t_DataType> &operator=( const OwnPtr< t_DataType> &ORIGINAL)
      {
        // return immediately if self-assignment
        if( m_Pointer == ORIGINAL.m_Pointer)
        {
          return *this;
        }

        // delete the existing object if it was owned
        if( m_IsOwner)
        {
          delete m_Pointer;
        }

        // copy the ownership permissions
        m_IsOwner = ORIGINAL.m_IsOwner;

        // copy the pointer itself; clone it if this pointer will also own the object.
        m_Pointer = m_IsOwner && ORIGINAL.m_Pointer
                    ? ORIGINAL.m_Pointer->Clone()
                    : ORIGINAL.m_Pointer;

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

      //! @brief return the const pointer
      //! @return pointer to const t_DataType
      const t_DataType *GetPointer() const
      {
        return m_Pointer;
      }

      //! @brief return the pointer
      //! @return pointer to t_DataType
      t_DataType *GetPointer()
      {
        return m_Pointer;
      }

      //! @brief check whether this own pointer owns the pointed-to object
      //! @return true if this pointer owns its object, and will delete it
      bool IsOwner() const
      {
        return m_IsOwner;
      }

    ////////////////
    // operations //
    ////////////////

      //! check whether OwnPtr is different from NULL
      bool IsDefined() const
      {
        return ( m_Pointer != NULL);
      }

    ///////////////
    // operators //
    ///////////////

      //! overload operator -> to return pointer to object t_DataType
      const t_DataType *operator ->() const
      {
        this->AssertIsDefined();
        return m_Pointer;
      }

      //! overload operator -> to return pointer to object t_DataType
      t_DataType *operator ->()
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

      //! overload operator * to return reference to changeable object t_DataType
      t_DataType &operator *()
      {
        this->AssertIsDefined();
        return *m_Pointer;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // return
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return ostream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // return
        return OSTREAM;
      }

    }; // template class OwnPtr

  } // namespace util
} // namespace bcl

#endif // BCL_UTIL_OWN_PTR_H_
