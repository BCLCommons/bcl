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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "linal/bcl_linal_vector_2d.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_rotation_matrix_3d.h"
#include "math/bcl_math_transformation_matrix_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Vector2D::s_Instance
    (
      GetObjectInstances().AddInstance( new Vector2D())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Vector2D::Vector2D() :
      m_X( 0.0),
      m_Y( 0.0)
    {
    }

    //! @brief constructor from a single value
    //! @param VALUE common value for all the elements
    Vector2D::Vector2D( const double &VALUE) :
      m_X( VALUE),
      m_Y( VALUE)
    {
    }

    //! @brief constructor from two values
    //! @param X, Y, Z two elements
    Vector2D::Vector2D( const double &X, const double &Y) :
      m_X( X),
      m_Y( Y)
    {
    }

    //! @brief constructor from a pointer to two elements
    //! @param PTR_DATA a pointer to two elements
    Vector2D::Vector2D( const double *PTR_DATA) :
      m_X( PTR_DATA[ 0]),
      m_Y( PTR_DATA[ 1])
    {
    }

    //! @brief constructor from VectorInterface
    //! @param VECTOR vector
    Vector2D::Vector2D( const VectorConstInterface< double> &VECTOR)
    {
      BCL_Assert( VECTOR.GetSize() == s_Size, "given vector has wrong length");
      std::copy( VECTOR.Begin(), VECTOR.End(), &m_X);
    }

    //! copy constructor
    Vector2D *Vector2D::Clone() const
    {
      return new Vector2D( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Vector2D::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns true if all the coordinates are defined
    //! @return boolean true if all coordinates are defined - false otherwise
    bool Vector2D::IsDefined() const
    {
      return util::IsDefined( m_X) && util::IsDefined( m_Y);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief set all elements to elements from another vector
    //! @param VECTOR vector, which must contain 2 elements
    //! @return reference to this, after copying vector
    Vector2D &Vector2D::operator =( const VectorConstInterface< double> &VECTOR)
    {
      BCL_Assert( VECTOR.GetSize() == s_Size, "Tried to copy non-2d vector onto vector2d");
      m_X = VECTOR( 0);
      m_Y = VECTOR( 1);
      return *this;
    }

    //! @brief set all elements to a vector 2d
    //! @param VECTOR the vector to copy into this vector
    //! @return reference to this, after copying vector
    Vector2D &Vector2D::operator =( const Vector2D &VECTOR)
    {
      m_X = VECTOR.m_X;
      m_Y = VECTOR.m_Y;
      return *this;
    }

  ////////////////
  // operations //
  ////////////////

    //! change all values
    Vector2D &Vector2D::Set( const double &X, const double &Y)
    {
      m_X = X;
      m_Y = Y;
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Vector2D::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_X, ISTREAM);
      io::Serialize::Read( m_Y, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Vector2D::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_X, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_Y, OSTREAM, 0);

      // return the stream
      return OSTREAM;
    }

  } // namespace linal
} // namespace bcl
