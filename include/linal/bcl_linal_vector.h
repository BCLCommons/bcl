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

#ifndef BCL_LINAL_VECTOR_H_
#define BCL_LINAL_VECTOR_H_

// include the namespace header
#include "bcl_linal.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_linal_vector_interface.h"
#include "util/bcl_util_assert.h"
#include "util/bcl_util_format.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Vector
    //! @brief template class Vector TODO: add a brief comment
    //! @details TODO: add an general comment to this class
    //!
    //! @see @link example_linal_vector.cpp @endlink
    //! @author woetzen
    //! @date Aug 11, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class Vector :
      public VectorInterface< t_DataType>
    {

    private:

    //////////
    // data //
    //////////

      //! length of vector
      size_t m_Size;

      //! range of dynamically allocated memory of size m_Size
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
      Vector();

      //! @brief construct from size and possible filler
      //! @param SIZE number fo elements in Vector
      //! @param FILL_VALUE assign every element to that value
      explicit Vector( const unsigned long &SIZE, const t_DataType &FILL_VALUE = t_DataType( 0));

      //! @brief construct from size and possible filler
      //! @param SIZE number fo elements in Vector
      //! @param FILL_VALUE assign every element to that value
      explicit Vector( const unsigned int &SIZE, const t_DataType &FILL_VALUE = t_DataType( 0));

      //! @brief construct from size and possible filler
      //! @param SIZE number fo elements in Vector
      //! @param FILL_VALUE assign every element to that value
      explicit Vector( const int &SIZE, const t_DataType &FILL_VALUE = t_DataType( 0));

      //! @brief construct from length and pointer to data
      Vector( const size_t SIZE, const t_DataType *DATA);

      //! construct Vector from storage::Vector< t_DataType>
      Vector( const storage::Vector< t_DataType> &STORAGEVECTOR);

      //! @brief construct from iterators
      template< typename t_InputIteratorType>
      Vector( const t_InputIteratorType &BEGIN, const t_InputIteratorType &END) :
        m_Size( std::distance( BEGIN, END)),
        m_Data( new t_DataType[ m_Size])
      {
        // copy the elements into memory
        std::copy( BEGIN, END, m_Data);
      }

      //! @brief copy constructor
      //! @param VECTOR copy the given Vector
      Vector( const Vector< t_DataType> &VECTOR);

      //! @brief copy constructor
      //! @param VECTOR copy the given Vector
      template< typename t_OtherDataType>
      Vector( const Vector< t_OtherDataType> &BASE) :
        m_Size( BASE.GetSize()),
        m_Data( new t_DataType[ m_Size])
      {
        // check if allocation was successful
        BCL_Assert
        (
          m_Data,
          "unable to allocate memory for " + util::Format()( m_Size) + " elements of type: " +
          GetStaticClassName< t_DataType>()
        );

        // copy elements, one by one
        typename Vector< t_OtherDataType>::const_iterator itr_other( BASE.Begin());
        for( t_DataType *itr( m_Data), *itr_end( m_Data + m_Size); itr < itr_end; ++itr, ++itr_other)
        {
          *itr = t_DataType( *itr_other);
        }
      }

      //! @brief copy constructor
      //! @param VECTOR_INTERFACE copy the given Vector
      Vector( const VectorConstInterface< t_DataType> &VECTOR_INTERFACE);

      //! @brief constructor from VectorInterface and padding
      //! @param VECTOR_INTERFACE copy the given Vector
      //! @param PADDING number of additional hidden elements to pad the vector
      //! @param FILL_VALUE assign every pad element to that value
      Vector
      (
        const VectorConstInterface< t_DataType> &VECTOR_INTERFACE,
        const size_t PADDING,
        const t_DataType &FILL_VALUE = t_DataType( 0)
      );

      //! @brief move constructor
      Vector( Vector< t_DataType> && VECTOR);

      //! @brief Clone function
      //! @return pointer to new Vector< t_DataType>
      Vector< t_DataType> *Clone() const;

      //! @ brief destructor
      ~Vector();

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief size of vector
      //! @return size of Vector
      size_t GetSize() const
      {
        return m_Size;
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
        return m_Data + m_Size;
      }

      //! @brief pointer to end of range
      //! @return pointer to address one after last element in Vector
      t_DataType *End()
      {
        return m_Data + m_Size;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief create subvector from this vector
      //! @param SIZE number of elements in subvector
      //! @param POS starting position in vector - default 0
      //! @return vector, that is a subvector of this
      Vector< t_DataType> CreateSubVector
      (
        const size_t SIZE,
        const size_t POS = 0
      ) const;

      //! @brief discard unnecessary points from the end of the vector (without reallocating it)
      //! @param SIZE new size of the vector
      void Shrink( const size_t SIZE);

    ///////////////
    // operators //
    ///////////////

      //! return reference to changeable element ( POS)
      t_DataType &operator()( const size_t POS)
      {
        AssertValidPosition( POS);
        return m_Data[ POS];
      }

      //! return reference of element ( POS)
      const t_DataType &operator()( const size_t POS) const
      {
        AssertValidPosition( POS);
        return m_Data[ POS];
      }

      //! @brief return pointer to element ( POS)
      //! @param POS position of the element requested
      //! @return pointer to element ( POS)
      t_DataType *operator[]( const size_t POS)
      {
        AssertValidPosition( POS);
        return m_Data + POS;
      }

      //! @brief return pointer to const element ( POS)
      //! @param POS position of the element requested
      //! @return pointer to const element ( POS)
      const t_DataType *operator[]( const size_t POS) const
      {
        AssertValidPosition( POS);
        return m_Data + POS;
      }

      //! @brief move operator
      //! @param VECTOR source vector
      //! @return reference to this assigned Vector
      Vector< t_DataType> &operator =( Vector< t_DataType> && VECTOR);

      //! @brief equal operator
      //! @param VECTOR source vector
      //! @return reference to this assigned Vector
      Vector< t_DataType> &operator =( const Vector< t_DataType> &VECTOR);

      //! @brief equal operator
      //! @param VECTOR source vector
      //! @return reference to this assigned Vector
      Vector< t_DataType> &operator =( const VectorConstInterface< t_DataType> &VECTOR);

      //! @brief equal operator
      //! @param VALUE all elements are set to that value
      //! @return reference to this assigned Vector
      Vector< t_DataType> &operator =( const t_DataType &VALUE);

    //////////////////////
    // input and output //
    //////////////////////

      //! write Vector to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const util::Format &FORMAT) const;

    protected:

      //! write Vector to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! read Vector from std::istream
      std::istream &Read( std::istream &ISTREAM);

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

    }; // template class Vector

    //! construct vector from one element
    template< typename t_DataType>
    inline Vector< t_DataType> MakeVector( const t_DataType &X)
    {
      return Vector< t_DataType>( 1, X);
    }

    //! construct vector from two elements
    template< typename t_DataType>
    inline Vector< t_DataType> MakeVector( const t_DataType &X, const t_DataType &Y)
    {
      Vector< t_DataType> newvector( 2);
      newvector( 0) = X;
      newvector( 1) = Y;

      return newvector;
    }

    //! construct vector from three elements
    template< typename t_DataType>
    inline Vector< t_DataType> MakeVector( const t_DataType &X, const t_DataType &Y, const t_DataType &Z)
    {
      Vector< t_DataType> newvector( 3);
      newvector( 0) = X;
      newvector( 1) = Y;
      newvector( 2) = Z;

      return newvector;
    }

    //! construct vector from four elements
    template< typename t_DataType>
    inline Vector< t_DataType> MakeVector
    (
      const t_DataType &X,
      const t_DataType &Y,
      const t_DataType &Z,
      const t_DataType &T
    )
    {
      Vector< t_DataType> newvector( 4);
      newvector( 0) = X;
      newvector( 1) = Y;
      newvector( 2) = Z;
      newvector( 3) = T;

      return newvector;
    }

    //! construct vector from four elements
    template< typename t_DataType>
    inline Vector< t_DataType> MakeVector
    (
      const t_DataType &A,
      const t_DataType &B,
      const t_DataType &C,
      const t_DataType &D,
      const t_DataType &E
    )
    {
      Vector< t_DataType> newvector( 5);
      newvector( 0) = A;
      newvector( 1) = B;
      newvector( 2) = C;
      newvector( 3) = D;
      newvector( 4) = E;

      return newvector;
    }

    //! construct vector from six elements
    template< typename t_DataType>
    inline Vector< t_DataType> MakeVector
    (
      const t_DataType &A,
      const t_DataType &B,
      const t_DataType &C,
      const t_DataType &D,
      const t_DataType &E,
      const t_DataType &F
    )
    {
      Vector< t_DataType> newvector( 6);
      newvector( 0) = A;
      newvector( 1) = B;
      newvector( 2) = C;
      newvector( 3) = D;
      newvector( 4) = E;
      newvector( 5) = F;

      return newvector;
    }

    //! construct vector from ten elements
    template< typename t_DataType>
    inline Vector< t_DataType> MakeVector
    (
      const t_DataType &A,
      const t_DataType &B,
      const t_DataType &C,
      const t_DataType &D,
      const t_DataType &E,
      const t_DataType &F,
      const t_DataType &G,
      const t_DataType &H,
      const t_DataType &I,
      const t_DataType &J
    )
    {
      Vector< t_DataType> newvector( 10);
      newvector( 0) = A;
      newvector( 1) = B;
      newvector( 2) = C;
      newvector( 3) = D;
      newvector( 4) = E;
      newvector( 5) = F;
      newvector( 6) = G;
      newvector( 7) = H;
      newvector( 8) = I;
      newvector( 9) = J;

      return newvector;
    }

    //! construct vector from start value, length and delta
    //! @param LENGTH number of elements in the vector
    //! @param START value of element 0
    //! @param DELTA element[ i] = START + i * DELTA
    //! @return vector that is filled with incrementing values
    template< typename t_DataType>
    inline Vector< t_DataType> FillVector( const size_t LENGTH, const t_DataType &START, const t_DataType &DELTA)
    {
      Vector< t_DataType> newvector( LENGTH, START);
      size_t i( 0);
      for( t_DataType *ptr( newvector.Begin()), *ptr_end( newvector.End()); ptr != ptr_end; ++ptr, ++i)
      {
        *ptr += i * DELTA;
      }

      return newvector;
    }

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Vector< double>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Vector< float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Vector< int>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Vector< unsigned int>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Vector< unsigned long>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Vector< unsigned long long>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Vector< char>;

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_VECTOR_H_
