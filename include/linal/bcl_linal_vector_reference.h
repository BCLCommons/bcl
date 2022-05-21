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

#ifndef BCL_LINAL_VECTOR_REFERENCE_H_
#define BCL_LINAL_VECTOR_REFERENCE_H_

// include the namespace header
#include "bcl_linal.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_linal_vector_interface.h"
#include "io/bcl_io_serialize.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class VectorReference
    //! @brief reference a particular vector range
    //!
    //! @see @link example_linal_vector_reference.cpp @endlink
    //! @author mendenjl
    //! @date Dec 18, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class VectorReference :
      public VectorInterface< t_DataType>
    {

    private:

    //////////
    // data //
    //////////

      //! length of VectorReference
      size_t m_Size;

      //! range of memory of size m_Size
      t_DataType *m_Data;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      VectorReference() :
        m_Size( 0),
        m_Data( NULL)
      {
      }

      //! @brief construct from dimension and pointer to data
      //! @param SIZE number of elements in DATA to reference
      //! @param DATA pointer to field of data
      VectorReference( const size_t SIZE, t_DataType *DATA) :
        m_Size( SIZE),
        m_Data( DATA)
      {
        // check if allocation was successful
        BCL_Assert
        (
          m_Data != NULL || m_Size == size_t( 0),
          "cannot reference NULL data elements of type: " + GetStaticClassName< t_DataType>()
        );
      }

      //! @brief copy constructor
      VectorReference( const VectorReference< t_DataType> &REFERENCE) :
        m_Size( REFERENCE.m_Size),
        m_Data( REFERENCE.m_Data)
      {
      }

      //! @brief initialize from VectorInterface
      //! @param VECTOR vector interface to be copied from
      explicit VectorReference( VectorInterface< t_DataType> &VECTOR) :
        m_Size( VECTOR.GetSize()),
        m_Data( VECTOR.Begin())
      {
      }

      //! @brief Clone function
      //! @return pointer to new VectorReference< t_DataType>
      VectorReference< t_DataType> *Clone() const
      {
        return new VectorReference< t_DataType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief size of VectorReference
      //! @return size of VectorReference
      size_t GetSize() const
      {
        return m_Size;
      }

      //! @brief pointer to First Element
      //! @return const pointer to first element in range containing all elements of VectorReference
      t_DataType *Begin()
      {
        return m_Data;
      }

      //! @brief pointer to end of range
      //! @return const pointer to address one after last element in VectorReference
      t_DataType *End()
      {
        return m_Data + m_Size;
      }

      //! @brief pointer to First Element
      //! @return const pointer to first element in range containing all elements of VectorReference
      const t_DataType *Begin() const
      {
        return m_Data;
      }

      //! @brief pointer to end of range
      //! @return const pointer to address one after last element in VectorReference
      const t_DataType *End() const
      {
        return m_Data + m_Size;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief create subVectorReference from this VectorReference
      //! @param SIZE number of elements in subVectorReference
      //! @param POS starting position in VectorReference - default 0
      //! @return VectorReference, that is a subVectorReference of this
      VectorReference< t_DataType> CreateSubVectorReference
      (
        const size_t SIZE,
        const size_t POS = 0
      ) const
      {
        // check that there are enough cols and rows
        BCL_Assert
        (
          SIZE + POS <= m_Size,
          "this VectorReference ( size = " + util::Format()( m_Size)
          + ") is too small for requested reference [" + util::Format()( POS)
          + "," + util::Format()( POS + SIZE) + ")"
        );

        // if size is equal (pos have to be 0 as asserted)
        if( SIZE == m_Size)
        {
          return *this;
        }

        // create VectorReference
        return VectorReference< t_DataType>( SIZE, m_Data + POS);
      }

      //! @brief fill the vector reference with a particular value
      //! @param VALUE value to set each element of the reference to
      VectorReference< t_DataType> &CopyValues( const VectorConstInterface< t_DataType> &VALUES)
      {
        BCL_Assert
        (
          VALUES.GetSize() == m_Size,
          "Tried to CopyValues from vector of size " + util::Format()( VALUES.GetSize())
          + " to a vector reference of size: " + util::Format()( m_Size)
        );
        std::copy( VALUES.Begin(), VALUES.End(), m_Data);
        return *this;
      }

      //! @brief make this reference point to the given vector
      //! @param VECTOR vector to reference
      VectorReference< t_DataType> &Reference( VectorInterface< t_DataType> &VALUES)
      {
        m_Size = VALUES.GetSize();
        m_Data = VALUES.Begin();
        return *this;
      }

    ///////////////
    // operators //
    ///////////////

      //! return reference of element ( POS)
      const t_DataType &operator()( const size_t POS) const
      {
        AssertValidPosition( POS);
        return m_Data[ POS];
      }

      //! return pointer to element
      const t_DataType *operator[]( const size_t POS) const
      {
        AssertValidPosition( POS);
        return m_Data + POS;
      }

      //! return reference of element ( POS)
      t_DataType &operator()( const size_t POS)
      {
        AssertValidPosition( POS);
        return m_Data[ POS];
      }

      //! return pointer to element
      t_DataType *operator[]( const size_t POS)
      {
        AssertValidPosition( POS);
        return m_Data + POS;
      }

      //! @brief fill the vector reference with a particular value
      //! @param VALUE value to set each element of the reference to
      VectorReference< t_DataType> &operator =( const t_DataType &VALUE)
      {
        std::fill( m_Data, m_Data + m_Size, VALUE);
        return *this;
      }

    //////////////////////
    // input and output //
    //////////////////////

      //! write VectorReference to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const util::Format &FORMAT) const
      {
        // write size
        OSTREAM << m_Size << '\n';

        //write data
        for( const t_DataType *ptr( Begin()), *ptr_end( End()); ptr != ptr_end; ++ptr)
        {
          OSTREAM << ' ' << FORMAT( *ptr);
        }

        // end
        return OSTREAM;
      }

    protected:

      //! write VectorReference to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write size
        io::Serialize::Write( m_Size, OSTREAM, INDENT) << '\n';
        io::Serialize::InsertIndent( OSTREAM, INDENT);

        //write data
        for( const t_DataType *ptr( Begin()), *ptr_end( End()); ptr != ptr_end; ++ptr)
        {
          io::Serialize::Write( *ptr, OSTREAM) << '\t';
        }

        // end
        return OSTREAM;
      }

      //! read VectorReference from std::istream
      std::istream &Read( std::istream &ISTREAM)
      {
        BCL_Exit( "You can't read in a vector reference because it does not have ownership!", -1);

        // end
        return ISTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

      //! check whether position is valid
      void AssertValidPosition( const size_t POS) const
      {
        BCL_Assert
        (
          POS < m_Size,
          "cannot access element outside range! " + util::Format()( POS) + " >= " + util::Format()( m_Size)
        );
      }

    private:

      //! @brief equal operator
      //! @param VECTOR source VectorReference
      //! @return reference to this assigned VectorReference
      //! undefined because it is unclear whether operator = would be shallow (copy pointers) or deep (copy values)
      //! Normally, copying references would be effectively a deep copy, but this class cannot always act like a true
      //! reference, since it has no way to dynamically resize the memory block
      VectorReference< t_DataType> &operator =( const VectorReference< t_DataType> &VECTOR);

    }; // template class VectorReference

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_VECTOR_REFERENCE_H_
