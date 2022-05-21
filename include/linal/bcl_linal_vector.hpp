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

#ifndef BCL_LINAL_VECTOR_HPP_
#define BCL_LINAL_VECTOR_HPP_

// include header of this class
#include "bcl_linal_vector.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_linal.h"
#include "io/bcl_io_serialize.h"
#include "math/bcl_math_statistics.h"
#include "util/bcl_util_assert.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
  //////////
  // data //
  //////////

    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> Vector< t_DataType>::s_Instance
    (
      GetObjectInstances().AddInstance( new Vector< t_DataType>())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    template< typename t_DataType>
    Vector< t_DataType>::Vector() :
      m_Size( 0),
      m_Data( NULL)
    {
    }

    //! @brief construct from size and possible filler
    //! @param SIZE number fo elements in Vector
    //! @param FILL_VALUE assign every element to that value
    template< typename t_DataType>
    Vector< t_DataType>::Vector( const unsigned long &SIZE, const t_DataType &FILL_VALUE) :
      m_Size( SIZE),
      m_Data( new t_DataType[ SIZE])
    {
      // check if allocation was successful
      BCL_Assert
      (
        m_Data,
        "unable to allocate memory for " + util::Format()( m_Size) + " elements of type: " +
        GetStaticClassName< t_DataType>()
      );

      // set all values to FILL_VALUE
      std::fill( Begin(), End(), FILL_VALUE);
    }

    //! @brief construct from size and possible filler
    //! @param SIZE number fo elements in Vector
    //! @param FILL_VALUE assign every element to that value
    template< typename t_DataType>
    Vector< t_DataType>::Vector( const unsigned int &SIZE, const t_DataType &FILL_VALUE) :
      m_Size( SIZE),
      m_Data( new t_DataType[ SIZE])
    {
      // check if allocation was successful
      BCL_Assert
      (
        m_Data,
        "unable to allocate memory for " + util::Format()( m_Size) + " elements of type: " +
        GetStaticClassName< t_DataType>()
      );

      // set all values to FILL_VALUE
      std::fill( Begin(), End(), FILL_VALUE);
    }

    //! @brief construct from size and possible filler
    //! @param SIZE number fo elements in Vector
    //! @param FILL_VALUE assign every element to that value
    template< typename t_DataType>
    Vector< t_DataType>::Vector( const int &SIZE, const t_DataType &FILL_VALUE) :
      m_Size( SIZE),
      m_Data( new t_DataType[ SIZE])
    {
      // check if allocation was successful
      BCL_Assert
      (
        m_Data,
        "unable to allocate memory for " + util::Format()( m_Size) + " elements of type: " +
        GetStaticClassName< t_DataType>()
      );

      // set all values to FILL_VALUE
      std::fill( Begin(), End(), FILL_VALUE);
    }

    //! @brief construct from length and pointer to data
    template< typename t_DataType>
    Vector< t_DataType>::Vector( const size_t SIZE, const t_DataType *DATA) :
      m_Size( SIZE),
      m_Data( new t_DataType[ SIZE])
    {
      // check if allocation was successful
      BCL_Assert
      (
        m_Data,
        "unable to allocate memory for " + util::Format()( m_Size) + " elements of type: " +
        GetStaticClassName< t_DataType>()
      );

      std::copy( DATA, DATA + SIZE, m_Data);
    }

    //! @brief copy constructor
    //! @param VECTOR copy the given Vector
    template< typename t_DataType>
    Vector< t_DataType>::Vector( const Vector< t_DataType> &VECTOR) :
      m_Size( VECTOR.m_Size),
      m_Data( new t_DataType[ m_Size])
    {
      // check if allocation was successful
      BCL_Assert
      (
        m_Data,
        "unable to allocate memory for " + util::Format()( m_Size) + " elements of type: " +
        GetStaticClassName< t_DataType>()
      );

      std::copy( VECTOR.m_Data, VECTOR.m_Data + m_Size, m_Data);
    }

    //! @brief copy constructor
    //! @param VECTOR_INTERFACE copy the given Vector
    template< typename t_DataType>
    Vector< t_DataType>::Vector( const VectorConstInterface< t_DataType> &VECTOR_INTERFACE) :
      m_Size( VECTOR_INTERFACE.GetSize()),
      m_Data( new t_DataType[ m_Size])
    {
      // check if allocation was successful
      BCL_Assert
      (
        m_Data,
        "unable to allocate memory for " + util::Format()( m_Size) + " elements of type: " +
        GetStaticClassName< t_DataType>()
      );

      std::copy( VECTOR_INTERFACE.Begin(), VECTOR_INTERFACE.End(), m_Data);
    }

    //! @brief constructor from VectorInterface and padding
    //! @param VECTOR_INTERFACE copy the given Vector
    //! @param PADDING number of additional hidden elements to pad the vector
    //! @param FILL_VALUE assign every pad element to that value
    template< typename t_DataType>
    Vector< t_DataType>::Vector
    (
      const VectorConstInterface< t_DataType> &VECTOR_INTERFACE,
      const size_t PADDING,
      const t_DataType &FILL_VALUE
    ) :
      m_Size( VECTOR_INTERFACE.GetSize() + PADDING),
      m_Data( new t_DataType[ m_Size])
    {
      // check if allocation was successful
      BCL_Assert
      (
        m_Size == 0 || m_Data,
        "unable to allocate memory for " + util::Format()( m_Size) + " elements of type: " +
        GetStaticClassName< t_DataType>()
      );

      std::copy( VECTOR_INTERFACE.Begin(), VECTOR_INTERFACE.End(), m_Data);
      std::fill( m_Data + VECTOR_INTERFACE.GetSize(), m_Data + m_Size, FILL_VALUE);
    }

    //! construct Vector from storage::Vector< t_DataType>
    template< typename t_DataType>
    Vector< t_DataType>::Vector( const storage::Vector< t_DataType> &STORAGEVECTOR) :
      m_Size( STORAGEVECTOR.GetSize()),
      m_Data( new t_DataType[ m_Size])
    {
      std::copy( STORAGEVECTOR.Begin(), STORAGEVECTOR.End(), m_Data);
    }

    //! @brief move constructor
    template< typename t_DataType>
    Vector< t_DataType>::Vector( Vector< t_DataType> && VECTOR) :
      m_Size( VECTOR.m_Size),
      m_Data( VECTOR.m_Data)
    {
      VECTOR.m_Data = nullptr;
      VECTOR.m_Size = 0;
    }

    //! @brief Clone function
    //! @return pointer to new Vector< t_DataType>
    template< typename t_DataType>
    Vector< t_DataType> *Vector< t_DataType>::Clone() const
    {
      return new Vector< t_DataType>( *this);
    }

    //! @ brief destructor
    template< typename t_DataType>
    Vector< t_DataType>::~Vector()
    {
      delete[] m_Data;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType>
    const std::string &Vector< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief create subvector from this vector
    //! @param SIZE number of elements in subvector
    //! @param POS starting position in vector - default 0
    //! @return vector, that is a subvector of this
    template< typename t_DataType>
    Vector< t_DataType> Vector< t_DataType>::CreateSubVector
    (
      const size_t SIZE,
      const size_t POS
    ) const
    {
      // check that there are enough cols and rows
      BCL_Assert( SIZE + POS <= m_Size, "this vector is too small for subvector");

      // if size is equal (pos have to be 0 as asserted)
      if( SIZE == m_Size)
      {
        return *this;
      }

      // create vector
      return Vector< t_DataType>( SIZE, m_Data + POS);
    }

    //! @brief discard unnecessary points from the end of the vector (without reallocating it)
    //! @param SIZE new size of the vector
    template< typename t_DataType>
    void Vector< t_DataType>::Shrink( const size_t SIZE)
    {
      BCL_Assert( SIZE <= m_Size, "Cannot increase vector size with the Shrink function!");
      m_Size = SIZE;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief move operator
    //! @param VECTOR source vector
    //! @return reference to this assigned Vector
    template< typename t_DataType>
    Vector< t_DataType> &Vector< t_DataType>::operator =( Vector< t_DataType> && VECTOR)
    {
      // check that data is different
      if( m_Data != VECTOR.m_Data)
      {
        // compare sizes
        if( m_Size)
        {
          // delete data
          delete[] m_Data;
        }
        m_Size = VECTOR.m_Size;
        m_Data = VECTOR.m_Data;
        VECTOR.m_Data = nullptr;
        VECTOR.m_Size = 0;
      }
      return *this;
    }

    //! @brief equal operator
    //! @param VECTOR source vector
    //! @return reference to this assigned Vector
    template< typename t_DataType>
    Vector< t_DataType> &Vector< t_DataType>::operator =( const Vector< t_DataType> &VECTOR)
    {
      // check that data is different
      if( m_Data != VECTOR.m_Data)
      {
        // compare sizes
        if( m_Size != VECTOR.m_Size)
        {
          // delete data
          delete[] m_Data;

          // assign size and allocate
          m_Size = VECTOR.m_Size;
          m_Data = new t_DataType[ m_Size];

          // check for allocation
          BCL_Assert
          (
            m_Data,
            "unable to allocate memory for " + util::Format()( m_Size) + " elements of type: " +
            GetStaticClassName< t_DataType>()
          );
        }

        // copy all elements
        std::copy( VECTOR.m_Data, VECTOR.m_Data + m_Size, m_Data);
      }

      // return reference to this Vector
      return *this;
    }

    //! @brief equal operator
    //! @param VECTOR source vector
    //! @return reference to this assigned Vector
    template< typename t_DataType>
    Vector< t_DataType> &Vector< t_DataType>::operator =( const VectorConstInterface< t_DataType> &VECTOR)
    {
      // compare sizes
      if( m_Size != VECTOR.GetSize())
      {
        // delete data
        delete[] m_Data;

        // assign size and allocate
        m_Size = VECTOR.GetSize();
        m_Data = new t_DataType[ m_Size];
      }

      // copy all elements
      if( m_Data != VECTOR.Begin())
      {
        std::copy( VECTOR.Begin(), VECTOR.End(), m_Data);
      }

      // return reference to this Vector
      return *this;
    }

    //! @brief equal operator
    //! @param VALUE all elements are set to that value
    //! @return reference to this assigned Vector
    template< typename t_DataType>
    Vector< t_DataType> &Vector< t_DataType>::operator =( const t_DataType &VALUE)
    {
      // set all element to given VALUE
      std::fill( m_Data, m_Data + m_Size, VALUE);

      // return reference to this Vector
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! write Vector to std::ostream
    template< typename t_DataType>
    std::ostream &Vector< t_DataType>::Write( std::ostream &OSTREAM, const util::Format &FORMAT) const
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

    //! write Vector to std::ostream
    template< typename t_DataType>
    std::ostream &Vector< t_DataType>::Write( std::ostream &OSTREAM, const size_t INDENT) const
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

    //! read Vector from std::istream
    template< typename t_DataType>
    std::istream &Vector< t_DataType>::Read( std::istream &ISTREAM)
    {
      // read size
      BCL_Assert( ISTREAM >> m_Size, "unable to read size");

      delete [] m_Data;
      m_Data = new t_DataType[ m_Size];

      // check if allocation was successful
      BCL_Assert
      (
        m_Data,
        "unable to allocate memory for " + util::Format()( m_Size) + " elements of type: " +
        GetStaticClassName< t_DataType>()
      );

      ISTREAM >> std::ws;

      if( !m_Size)
      {
        // nothing to read
        return ISTREAM;
      }

      t_DataType *ptr( Begin());
      // read 1 line into a string stream, which is then read into a matrix.
      // This is > 3x faster if the stream is coming from a compressed file than reading directly due to the peeks and
      // functions employed in reading numbers directly from a data stream, operations which are not as fast on
      // compressed streams. It is 10-20% faster on uncompressed streams as well. Only if the data is another string
      // stream will this be marginally slower
      std::string templine;
      std::stringstream vec_stream;
      std::getline( ISTREAM, templine);
      vec_stream.str( templine);
      for( t_DataType *ptr_end( ptr + m_Size); ptr != ptr_end; ++ptr)
      {
        BCL_Assert
        (
          io::Serialize::Read( *ptr, vec_stream),
          "unable to read element " + util::Format()( size_t( ptr - m_Data)) + " of " + util::Format()( m_Size) +
          " elements, so far: " + util::Format()( *this)
        );
      }

      // end
      return ISTREAM;
    }

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_VECTOR_HPP_
