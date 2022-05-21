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

#ifndef BCL_LINAL_VECTOR_2D_H_
#define BCL_LINAL_VECTOR_2D_H_

// include the namespace header
#include "bcl_linal.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_linal_vector_interface.h"
#include "util/bcl_util_assert.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Vector2D
    //! @brief a vector with 2 doubles
    //! @details this is a basic implementation for representing coordinates in a plane
    //!
    //! @see @link example_linal_vector_2d.cpp @endlink
    //! @author mendenjl
    //! @date Dec 05, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Vector2D :
      public VectorInterface< double>
    {
    private:

    //////////
    // data //
    //////////

      double m_X; //!< X
      double m_Y; //!< Y

      static const size_t s_Size = 2; //!< static size of the Vector2D

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
      Vector2D();

      //! @brief constructor from a single value
      //! @param VALUE common value for all the elements
      explicit Vector2D( const double &VALUE);

      //! @brief constructor from two values
      //! @param X, Y two elements
      Vector2D( const double &X, const double &Y);

      //! @brief constructor from a pointer to two elements
      //! @param PTR_DATA a pointer to two elements
      explicit Vector2D( const double *PTR_DATA);

      //! @brief constructor from VectorInterface
      //! @param VECTOR vector
      explicit Vector2D( const VectorConstInterface< double> &VECTOR);

      //! copy constructor
      Vector2D *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns true if all the coordinates are defined
      //! @return boolean true if all coordinates are defined - false otherwise
      bool IsDefined() const;

      //! @brief return changeable X
      //! @return changeable X
      double &X()
      {
        return m_X;
      }

      //! @brief return copy of X
      //! @return const X
      const double &X() const
      {
        return m_X;
      }

      //! @brief return changeable Y
      //! @return changeable Y
      double &Y()
      {
        return m_Y;
      }

      //! @brief return copy of Y
      //! @return const Y
      const double &Y() const
      {
        return m_Y;
      }

      //! @brief pointer to First Element
      //! @return const pointer to m_X
      const double *Begin() const
      {
        return &m_X;
      }

      //! @brief pointer to First Element
      //! @return pointer to m_X
      double *Begin()
      {
        return &m_X;
      }

      //! @brief pointer to end of range
      //! @return const pointer to address one after m_Y
      const double *End() const
      {
        return &m_X + s_Size;
      }

      //! @brief pointer to end of range
      //! @return pointer to address one after m_Y
      double *End()
      {
        return &m_X + s_Size;
      }

      //! @brief size of vector
      //! @return size of Vector
      size_t GetSize() const
      {
        return s_Size;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief return non-const reference to element ( POS)
      //! @param POS position of the element requested
      //! @return non-const reference to element ( POS)
      double &operator()( const size_t POS)
      {
        BCL_Assert( POS < s_Size, "Out of range access!");
        return ( &m_X)[ POS];
      }

      //! @brief return const reference to element ( POS)
      //! @param POS position of the element requested
      //! @return const reference to element ( POS)
      const double &operator()( const size_t POS) const
      {
        BCL_Assert( POS < s_Size, "Out of range access!");
        return ( &m_X)[ POS];
      }

      //! @brief return pointer to element ( POS)
      //! @param POS position of the element requested
      //! @return pointer to element ( POS)
      double *operator[]( const size_t POS)
      {
        BCL_Assert( POS < s_Size, "Out of range access!");
        return &m_X + POS;
      }

      //! @brief return pointer to const element ( POS)
      //! @param POS position of the element requested
      //! @return pointer to const element ( POS)
      const double *operator[]( const size_t POS) const
      {
        BCL_Assert( POS < s_Size, "Out of range access!");
        return &m_X + POS;
      }

      //! @brief set all elements to a common scalar value
      //! @param SCALAR scalar value to set this vector equal to
      //! @return reference to this, after setting all elements to SCALAR
      Vector2D &operator =( const double &SCALAR)
      {
        m_X = m_Y = SCALAR;
        return *this;
      }

      //! @brief set all elements to elements from another vector
      //! @param VECTOR vector, which must contain 2 elements
      //! @return reference to this, after copying vector
      Vector2D &operator =( const VectorConstInterface< double> &VECTOR);

      //! @brief set all elements to a vector 2d
      //! @param VECTOR the vector to copy into this vector
      //! @return reference to this, after copying vector
      Vector2D &operator =( const Vector2D &VECTOR);

    ////////////////
    // operations //
    ////////////////

      //! change all values
      Vector2D &Set( const double &X, const double &Y);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class Vector2D

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_VECTOR_2D_H_
