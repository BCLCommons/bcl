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

#ifndef BCL_LINAL_VECTOR_3D_H_
#define BCL_LINAL_VECTOR_3D_H_

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
    //! @class Vector3D
    //! @brief a vector with 3 doubles
    //! @details this is a basic implementation for representing coordinates and vectors
    //!
    //! @see @link example_linal_vector_3d.cpp @endlink
    //! @author karakam, woetzen
    //! @date Nov 6, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Vector3D :
      public VectorInterface< double>
    {
    private:

    //////////
    // data //
    //////////

      double m_X; //!< X
      double m_Y; //!< Y
      double m_Z; //!< Z

      static const size_t s_Size = 3; //!< static size of the Vector3D

    public:

    //////////
    // data //
    //////////

      // //! typedef for const_iterator
      /* typedef const double * const_iterator; */

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Vector3D();

      //! @brief constructor from a single value
      //! @param VALUE common value for all the elements
      explicit Vector3D( const double &VALUE);

      //! @brief constructor from three values
      //! @param X, Y, Z three elements
      Vector3D( const double &X, const double &Y, const double &Z);

      //! @brief constructor from a pointer to three elements
      //! @param PTR_DATA a pointer to three elements
      explicit Vector3D( const double *PTR_DATA);

      //! @brief constructor from VectorInterface
      //! @param VECTOR vector
      explicit Vector3D( const VectorConstInterface< double> &VECTOR);

      //! @brief construct a Vector from iterator [FIRST, LAST) range
      //! @param FIRST iterator to the first element to be copied
      //! @param LAST iterator to the first element after the last copied element
      template< typename t_InputIterator>
      Vector3D( const t_InputIterator &FIRST, const t_InputIterator &LAST)
      {
        BCL_Assert( std::distance( FIRST, LAST) == 3, "Vector3D requires exactly three coordinates");
        m_X = *FIRST;
        m_Y = *( FIRST + 1);
        m_Z = *( FIRST + 2);
      }

      //! copy constructor
      Vector3D *Clone() const;

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

      //! @brief return changeable Z
      //! @return changeable Z
      double &Z()
      {
        return m_Z;
      }

      //! @brief return copy of Z
      //! @return const Z
      const double &Z() const
      {
        return m_Z;
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
      //! @return const pointer to address one after m_Z
      const double *End() const
      {
        return &m_X + s_Size;
      }

      //! @brief pointer to end of range
      //! @return pointer to address one after m_Z
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
      Vector3D &operator =( const double &SCALAR)
      {
        m_X = m_Y = m_Z = SCALAR;
        return *this;
      }

      //! @brief set all elements to elements from another vector
      //! @param VECTOR vector, which must contain 3 elements
      //! @return reference to this, after copying vector
      Vector3D &operator =( const VectorConstInterface< double> &VECTOR);

      //! @brief set all elements to a vector 3d
      //! @param VECTOR the vector to copy into this vector
      //! @return reference to this, after copying vector
      Vector3D &operator =( const Vector3D &VECTOR);

      //! @brief add a vector to this one
      //! @param VECTOR the vector to add
      //! @return reference to this
      Vector3D &operator +=( const Vector3D &VECTOR);

      //! @brief add a scalar to this vector
      //! @param SCALAR the scalar to add
      //! @return reference to this
      Vector3D &operator +=( const double &SCALAR);

      //! @brief subtract a vector from this one
      //! @param VECTOR the vector to subtract
      //! @return reference to this
      Vector3D &operator -=( const Vector3D &VECTOR);

      //! @brief subtract a scalar from this vector
      //! @param SCALAR the value to subtract
      //! @return reference to this
      Vector3D &operator -=( const double &SCALAR);

    ////////////////
    // operations //
    ////////////////

      //! change all values
      Vector3D &Set( const double &X, const double &Y, const double &Z);

      //! translate vector with translation
      Vector3D &Translate( const Vector3D &TRANSLATE);

      //! translate vector with translation
      Vector3D &Translate( const double X, const double Y, const double Z);

      //! rotate vector with rotation matrix
      Vector3D &Rotate( const math::RotationMatrix3D &ROTATE);

      //! transform vector with transformation matrix
      Vector3D &Transform( const math::TransformationMatrix3D &TRANSFORM);

      //! return random translation vector equally distributed in a sphere of RADIUS
      Vector3D &SetRandomTranslation( const double &RADIUS);

      //! return random translation vector equally distributed in a ellipse of RADII
      //! @param RADII different radii for each directions in space
      Vector3D &SetRandomTranslation( const Vector3D &RADII);

      //! @brief normalize the vector (such that inner product == 1)
      //! Overrides function of the same name in VectorInterface for performance reasons
      Vector3D &Normalize();

      //! @brief normalize the vector (such that inner product == 1)
      //! Overrides function of the same name in VectorInterface for performance reasons
      Vector3D Normalized() const;

      //! @brief norm = length of vector
      //! @return length of vector
      double Norm() const;

      //! @brief square norm = square length of vector
      //! @return square length of vector
      double SquareNorm() const;

      //! @brief sum up all elements
      //! @return sum of all elements
      double Sum() const;

    //////////////////////
    // input and output //
    //////////////////////

      //! print in a user readable format
      //! @return the string containing the information
      std::string ToString() const;

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

    }; // class Vector3D

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_VECTOR_3D_H_
