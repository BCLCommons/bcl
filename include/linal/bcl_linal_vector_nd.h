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

#ifndef BCL_LINAL_VECTOR_ND_H_
#define BCL_LINAL_VECTOR_ND_H_

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
    //! @class VectorND
    //! @brief template class VectorND TODO: add a brief comment
    //! @details TODO: add an general comment to this class
    //!
    //! @see @link example_linal_vector_nd.cpp @endlink
    //! @author woetzen
    //! @date Aug 11, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, size_t t_N>
    class VectorND :
      public VectorInterface< t_DataType>
    {

    private:

    //////////
    // data //
    //////////

      //! range of dynamically allocated memory of size t_N
      t_DataType m_Data[ t_N];

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
      VectorND()
      {
        std::fill_n( m_Data, t_N, t_DataType( 0));
      }

      //! @brief constrcut from element
      //! @param ELEMENT all elements will be set to that value
      VectorND( const t_DataType &ELEMENT)
      {
        std::fill_n( m_Data, t_N, ELEMENT);
      }

      //! @brief constrcut from element
      //! @param ELEMENT all elements will be set to that value
      VectorND( const t_DataType &A, const t_DataType &B)
      {
        BCL_Assert( t_N == size_t( 2), "Wrong # elements!");
        m_Data[ 0] = A;
        m_Data[ 1] = B;
      }

      //! @brief constrcut from element
      //! @param ELEMENT all elements will be set to that value
      VectorND( const t_DataType &A, const t_DataType &B, const t_DataType &C)
      {
        BCL_Assert( t_N == size_t( 3), "Wrong # elements!");
        m_Data[ 0] = A;
        m_Data[ 1] = B;
        m_Data[ 2] = C;
      }

      //! @brief copy constructor
      //! @param VECTOR_ND copy the given Vector_ND
      VectorND( const VectorND< t_DataType, t_N> &VECTOR_ND)
      {
        std::copy( VECTOR_ND.m_Data, VECTOR_ND.m_Data + t_N, m_Data);
      }

      //! @brief copy constructor
      //! @param VECTOR other vector to copy
      VectorND( const VectorConstInterface< t_DataType> &VECTOR)
      {
        // check that sizes match
        BCL_Assert
        (
          t_N == VECTOR.GetSize(),
          "trying to set Vectors of different sizes equal: " + util::Format()( t_N) +
          " != " + util::Format()( VECTOR.GetSize())
        );

        std::copy( VECTOR.Begin(), VECTOR.End(), m_Data);
      }

      //! @brief construct a Vector from iterator [FIRST, LAST) range
      //! @param FIRST iterator to the first element to be copied
      //! @param LAST iterator to the first element after the last copied element
      template< typename t_InputIterator>
      VectorND( const t_InputIterator &FIRST, const t_InputIterator &LAST)
      {
        // check that sizes match
        BCL_Assert
        (
          t_N == std::distance( FIRST, LAST),
          "trying to set Vectors of different sizes equal: " + util::Format()( t_N) +
          " != " + util::Format()( std::distance( FIRST, LAST))
        );

        std::copy( FIRST, LAST, m_Data);
      }

      //! @brief Clone function
      //! @return pointer to new VectorND< t_DataType>
      VectorND< t_DataType, t_N> *Clone() const
      {
        return new VectorND< t_DataType, t_N>( *this);
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

      //! @brief size of vector
      //! @return size of Vector
      size_t GetSize() const
      {
        return t_N;
      }

      //! @brief pointer to First Element
      //! @return const pointer to first element in range containing all elements of Vector
      const t_DataType *Begin() const
      {
        return m_Data;
      }

      //! @brief pointer to First Element
      //! @return pointer to first element in range containing all elements of Vector
      t_DataType *Begin()
      {
        return m_Data;
      }

      //! @brief pointer to end of range
      //! @return const pointer to address one after last element in Vector
      const t_DataType *End() const
      {
        return m_Data + t_N;
      }

      //! @brief pointer to end of range
      //! @return pointer to address one after last element in Vector
      t_DataType *End()
      {
        return m_Data + t_N;
      }

    ////////////////
    // operations //
    ////////////////

      //! return reference to changeable element ( POS)
      t_DataType &operator()( const size_t POS)
      {
        return m_Data[ POS];
      }

      //! return copy of element ( POS)
      const t_DataType &operator()( const size_t POS) const
      {
        return m_Data[ POS];
      }

      //! @brief return pointer to element ( POS)
      //! @param POS position of the element requested
      //! @return pointer to element ( POS)
      t_DataType *operator[]( const size_t POS)
      {
        return m_Data + POS;
      }

      //! @brief return pointer to const element ( POS)
      //! @param POS position of the element requested
      //! @return pointer to const element ( POS)
      const t_DataType *operator[]( const size_t POS) const
      {
        return m_Data + POS;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief assignment operator from scalar
      //! @param SCALAR the scalar to fill this vector with
      //! @return reference to this VectorND
      VectorND &operator =( const t_DataType &SCALAR)
      {
        // fill the vector with SCALAR
        std::fill( m_Data, m_Data + t_N, SCALAR);

        // return reference to this Vector
        return *this;
      }

      //! @brief assignment operator from other vector nd
      //! @return reference to this VectorND
      VectorND &operator =( const VectorND< t_DataType, t_N> &VECTOR_ND)
      {
        // copy all element
        if( m_Data != VECTOR_ND.m_Data)
        {
          std::copy( VECTOR_ND.m_Data, VECTOR_ND.m_Data + t_N, m_Data);
        }

        // return reference to this Vector
        return *this;
      }

      //! @brief assignment operator from other vector
      //! @return reference to this VectorND
      VectorND &operator =( const VectorConstInterface< t_DataType> &VECTOR)
      {
        // check that sizes match
        BCL_Assert
        (
          t_N == VECTOR.GetSize(),
          "trying to set Vectors of different sizes equal: " + util::Format()( t_N) +
          " != " + util::Format()( VECTOR.GetSize())
        );

        // copy all elements
        if( m_Data != VECTOR.Begin())
        {
          std::copy( VECTOR.Begin(), VECTOR.End(), m_Data);
        }

        // return reference to this Vector
        return *this;
      }

    //////////////////////
    // input and output //
    //////////////////////

      //! write VectorND to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const util::Format &FORMAT) const
      {
        //write data
        for( const t_DataType *ptr( Begin()), *ptr_end( End()); ptr != ptr_end; ++ptr)
        {
          OSTREAM << ' ' << FORMAT( *ptr);
        }

        return OSTREAM;
      }

    protected:

      //! write VectorND to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        io::Serialize::InsertIndent( OSTREAM, INDENT);
        //write data
        for( const t_DataType *ptr( Begin()), *ptr_end( End()); ptr != ptr_end; ++ptr)
        {
          io::Serialize::Write( *ptr, OSTREAM) << '\t';
        }

        // end
        return OSTREAM;
      }

      //! read VectorND from std::istream
      std::istream &Read( std::istream &ISTREAM)
      {
        // read each element from ISTREAM
        for( t_DataType *ptr( Begin()), *ptr_end( End()); ptr != ptr_end; ++ptr)
        {
          BCL_Assert( io::Serialize::Read( *ptr, ISTREAM), "unable to read element");
        }

        // end
        return ISTREAM;
      }

    }; // template class VectorND

    //! @brief helper function to construct a vector ND 2
    //! @param A, B the values
    //! @return a new vector 2D with the specified values
    template< typename t_DataType>
    VectorND< t_DataType, 2> MakeVectorND( const t_DataType &A, const t_DataType &B)
    {
      VectorND< t_DataType, 2> vector;
      vector( 0) = A;
      vector( 1) = B;
      return vector;
    }

    //! @brief helper function to construct a vector ND 3
    //! @param A, B, C the values
    //! @return a new vector 3D with the specified values
    template< typename t_DataType>
    VectorND< t_DataType, 3> MakeVectorND( const t_DataType &A, const t_DataType &B, const t_DataType &C)
    {
      VectorND< t_DataType, 3> vector;
      vector( 0) = A;
      vector( 1) = B;
      vector( 2) = C;
      return vector;
    }

    //! @brief helper function to construct a vector ND 3
    //! @param A, B, C, D the values
    //! @return a new vector 4D with the specified values
    template< typename t_DataType>
    VectorND< t_DataType, 4> MakeVectorND( const t_DataType &A, const t_DataType &B, const t_DataType &C, const t_DataType &D)
    {
      VectorND< t_DataType, 4> vector;
      vector( 0) = A;
      vector( 1) = B;
      vector( 2) = C;
      vector( 3) = D;
      return vector;
    }

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_VECTOR_ND_H_
