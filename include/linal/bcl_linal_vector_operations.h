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

#ifndef BCL_LINAL_VECTOR_OPERATIONS_H_
#define BCL_LINAL_VECTOR_OPERATIONS_H_

// include the namespace header
#include "bcl_linal.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_linal_vector.h"
#include "math/bcl_math.h"

// external includes - sorted alphabetically
#include <algorithm>
#include <cmath>
#include <numeric>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @file bcl_linal_vector_operations.h
  //! @brief mathematical operators/operations that work on a vector
  //!
  //! @see @link example_linal_vector_operations.cpp @endlink
  //! @author woetzen
  //! @date Jul 22, 2010
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace linal
  {

  //////////////////////
  // binary operators //
  //////////////////////

    //! @brief add one vector to another
    //! @param VECTOR_LHS vector to add to
    //! @param VECTOR_RHS vector to add
    //! @return the changed lhs vector
    template< typename t_DataType>
    inline
    VectorInterface< t_DataType> &
    operator +=
    (
      VectorInterface< t_DataType> &VECTOR_LHS,
      const VectorConstInterface< t_DataType> &VECTOR_RHS
    )
    {
      BCL_Assert
      (
        VECTOR_LHS.GetSize() == VECTOR_RHS.GetSize(),
        "provided Vectors have different sizes; rhs != lhs : " + util::Format()( VECTOR_LHS.GetSize()) + " != " +
        util::Format()( VECTOR_RHS.GetSize())
      );

      std::transform
      (
        VECTOR_LHS.Begin(), VECTOR_LHS.End(),
        VECTOR_RHS.Begin(),
        VECTOR_LHS.Begin(),
        std::plus< t_DataType>()
      );

      return VECTOR_LHS;
    }

    //! @brief subtract one vector from another
    //! @param VECTOR_LHS vector to subtract from
    //! @param VECTOR_RHS vector to subtract
    //! @return the changed lhs vector
    template< typename t_DataType>
    inline
    VectorInterface< t_DataType> &
    operator -=
    (
      VectorInterface< t_DataType> &VECTOR_LHS,
      const VectorConstInterface< t_DataType> &VECTOR_RHS
    )
    {
      BCL_Assert
      (
        VECTOR_LHS.GetSize() == VECTOR_RHS.GetSize(),
        "provided Vectors have different sizes; rhs != lhs : " + util::Format()( VECTOR_LHS.GetSize()) + " != " +
        util::Format()( VECTOR_RHS.GetSize())
      );

      std::transform
      (
        VECTOR_LHS.Begin(), VECTOR_LHS.End(),
        VECTOR_RHS.Begin(),
        VECTOR_LHS.Begin(),
        std::minus< t_DataType>()
      );

      return VECTOR_LHS;
    }

    //! @brief divide one vector by another
    //! @param VECTOR_LHS vector to divided
    //! @param VECTOR_RHS vector to divide by
    //! @return the changed lhs vector
    template< typename t_DataType>
    inline
    VectorInterface< t_DataType> &
    operator /=
    (
      VectorInterface< t_DataType> &VECTOR_LHS,
      const VectorConstInterface< t_DataType> &VECTOR_RHS
    )
    {
      BCL_Assert
      (
        VECTOR_LHS.GetSize() == VECTOR_RHS.GetSize(),
        "provided Vectors have different sizes; rhs != lhs : " + util::Format()( VECTOR_LHS.GetSize()) + " != " +
        util::Format()( VECTOR_RHS.GetSize())
      );

      std::transform
      (
        VECTOR_LHS.Begin(), VECTOR_LHS.End(),
        VECTOR_RHS.Begin(),
        VECTOR_LHS.Begin(),
        std::divides< t_DataType>()
      );

      return VECTOR_LHS;
    }

    //! @brief add scalar to vector
    //! @param VECTOR_LHS vector to add to
    //! @param VALUE scalar to be added
    //! @return the changed lhs vector
    template< typename t_DataType>
    inline
    VectorInterface< t_DataType> &
    operator +=
    (
      VectorInterface< t_DataType> &VECTOR_LHS,
      const t_DataType &VALUE
    )
    {
      std::transform
      (
        VECTOR_LHS.Begin(), VECTOR_LHS.End(),
        VECTOR_LHS.Begin(),
        std::bind2nd( std::plus< t_DataType>(), VALUE)
      );

      return VECTOR_LHS;
    }

    //! @brief subtract scalar from vector
    //! @param VECTOR_LHS vector to subtract from
    //! @param VALUE scalar to be added
    //! @return the changed lhs vector
    template< typename t_DataType>
    inline
    VectorInterface< t_DataType> &
    operator -=
    (
      VectorInterface< t_DataType> &VECTOR_LHS,
      const t_DataType &VALUE
    )
    {
      std::transform
      (
        VECTOR_LHS.Begin(), VECTOR_LHS.End(),
        VECTOR_LHS.Begin(),
        std::bind2nd( std::minus< t_DataType>(), VALUE)
      );

      return VECTOR_LHS;
    }

    //! @brief multiply vector with scalar
    //! @param VECTOR_LHS vector to multiply to
    //! @param SCALAR scalar to be multiplied
    //! @return the changed lhs vector
    template< typename t_DataType, typename t_ScalarDataType>
    inline
    VectorInterface< t_DataType> &
    operator *=
    (
      VectorInterface< t_DataType> &VECTOR_LHS,
      const t_ScalarDataType &SCALAR
    )
    {
      std::transform
      (
        VECTOR_LHS.Begin(), VECTOR_LHS.End(),
        VECTOR_LHS.Begin(),
        std::bind2nd( std::multiplies< t_DataType>(), t_DataType( SCALAR))
      );

      return VECTOR_LHS;
    }

    //! @brief divide vector by scalar
    //! @param VECTOR_LHS vector to divide
    //! @param SCALAR scalar to divide by
    //! @return the changed lhs vector
    template< typename t_DataType>
    inline
    VectorInterface< t_DataType> &
    operator /=
    (
      VectorInterface< t_DataType> &VECTOR_LHS,
      const t_DataType &SCALAR
    )
    {
      std::transform
      (
        VECTOR_LHS.Begin(), VECTOR_LHS.End(),
        VECTOR_LHS.Begin(),
        std::bind2nd( std::divides< t_DataType>(), SCALAR)
      );

      return VECTOR_LHS;
    }

  //////////////////////////////
  // binary logical operators //
  //////////////////////////////

    //! @brief compare to vectors for equality
    //! @param VECTOR_LHS lhs vector
    //! @param VECTOR_RHS rhs vector
    //! @return true is they are equal in size and all pairs of items are equal
    template< typename t_DataType>
    inline
    bool
    operator ==
    (
      const VectorConstInterface< t_DataType>& VECTOR_LHS,
      const VectorConstInterface< t_DataType>& VECTOR_RHS
    )
    {
      // not equal if different size
      if( VECTOR_LHS.GetSize() != VECTOR_RHS.GetSize())
      {
        return false;
      }

      // check that all elements in bot ranges comparing each corresponding pair, is equal
      return std::equal( VECTOR_LHS.Begin(), VECTOR_LHS.End(), VECTOR_RHS.Begin());
    }

    //! @brief compare to vectors for inequality
    //! @param VECTOR_LHS lhs vector
    //! @param VECTOR_RHS rhs vector
    //! @return !( VECTOR_LHS == VECTOR_RHS)
    template< typename t_DataType>
    inline
    bool
    operator !=
    (
      const VectorConstInterface< t_DataType>& VECTOR_LHS,
      const VectorConstInterface< t_DataType>& VECTOR_RHS
    )
    {
      return !( VECTOR_LHS == VECTOR_RHS);
    }

    //! @brief compare if all items in vector are equal to a given VALUE
    //! @param VECTOR_LHS vector with values
    //! @param VALUE_RHS value that is compared against
    //! @return true if vector is empty are all elements in vector are equal to given VALUE
    template< typename t_DataType>
    inline
    bool
    operator ==
    (
      const VectorConstInterface< t_DataType> &VECTOR_LHS,
      const t_DataType &VALUE_RHS
    )
    {
      return std::find_if
             (
               VECTOR_LHS.Begin(), VECTOR_LHS.End(),
               std::bind2nd( std::not_equal_to< t_DataType>(), VALUE_RHS)
             ) == VECTOR_LHS.End();
    }

    //! @brief compare if all items in vector are equal to a given VALUE
    //! @param VALUE_LHS value that is compared against
    //! @param VECTOR_RHS vector with values
    //! @return true if vector is empty are all elements in vector are equal to given VALUE
    template< typename t_DataType>
    inline
    bool
    operator ==
    (
      const t_DataType &VALUE_LHS,
      const VectorConstInterface< t_DataType> &VECTOR_RHS
    )
    {
      return ( VECTOR_RHS == VALUE_LHS);
    }

    //! @brief compare if all items in vector are not equal to a given VALUE
    //! @param VECTOR_LHS vector with values
    //! @param VALUE_RHS value that is compared against
    //! @return false if vector is empty are all elements in vector are equal to given VALUE
    template< typename t_DataType>
    inline
    bool
    operator !=
    (
      const VectorConstInterface< t_DataType> &VECTOR_LHS,
      const t_DataType &VALUE_RHS
    )
    {
      return std::find_if
             (
               VECTOR_LHS.Begin(), VECTOR_LHS.End(),
               std::bind2nd( std::not_equal_to< t_DataType>(), VALUE_RHS)
             ) != VECTOR_LHS.End();
    }

    //! @brief compare if all items in vector are not equal to a given VALUE
    //! @param VALUE_LHS value that is compared against
    //! @param VECTOR_RHS vector with values
    //! @return false if vector is empty are all elements in vector are equal to given VALUE
    template< typename t_DataType>
    inline
    bool
    operator !=
    (
      const t_DataType &VALUE_LHS,
      const VectorConstInterface< t_DataType> &VECTOR_RHS
    )
    {
      return ( VECTOR_RHS != VALUE_LHS);
    }

    //! @brief operator < to compare to vectors
    //! @param VECTOR_LHS left hand side operand
    //! @param VECTOR_RHS right hand side operand
    //! @return true if the first position in VECTOR_LHS LHS that differs from VECTOR_RHS is less
    template< typename t_DataType>
    inline bool operator <
    (
      const VectorConstInterface< t_DataType> &VECTOR_LHS,
      const VectorConstInterface< t_DataType> &VECTOR_RHS
    )
    {
      // compare the first non-equal elements according to <
      for( size_t i( 0), vec_size( std::min( VECTOR_LHS.GetSize(), VECTOR_RHS.GetSize())); i < vec_size; ++i)
      {
        if( VECTOR_LHS( i) != VECTOR_RHS( i))
        {
          return VECTOR_LHS( i) < VECTOR_RHS( i);
        }
      }
      return VECTOR_LHS.GetSize() < VECTOR_RHS.GetSize();
    }

  //////////////////////
  // binary operators //
  //////////////////////

    //! @brief sum two vectors of equal size
    //! @param VECTOR_LHS lhs vector
    //! @param VECTOR_RHS rhs vector
    //! @return vector with all individual summed elements of lhs and rhs vector
    template< typename t_DataType>
    inline
    Vector< t_DataType>
    operator +
    (
      const VectorConstInterface< t_DataType> &VECTOR_LHS,
      const VectorConstInterface< t_DataType> &VECTOR_RHS
    )
    {
      BCL_Assert
      (
        VECTOR_LHS.GetSize() == VECTOR_RHS.GetSize(),
        "vectors of different sizes supplied: " +
        util::Format()( VECTOR_LHS.GetSize()) + " != " + util::Format()( VECTOR_RHS.GetSize())
      );

      Vector< t_DataType> new_vector( VECTOR_LHS);
      return new_vector += VECTOR_RHS;
    }

    //! @brief subtract two vectors of equal size
    //! @param VECTOR_LHS lhs vector
    //! @param VECTOR_RHS rhs vector
    //! @return vector with all individual subtracted elements of rhs from lhs vector
    template< typename t_DataType>
    inline
    Vector< t_DataType>
    operator -
    (
      const VectorConstInterface< t_DataType> &VECTOR_LHS,
      const VectorConstInterface< t_DataType> &VECTOR_RHS
    )
    {
      BCL_Assert
      (
        VECTOR_LHS.GetSize() == VECTOR_RHS.GetSize(),
        "vectors of different sizes supplied: " +
        util::Format()( VECTOR_LHS.GetSize()) + " != " + util::Format()( VECTOR_RHS.GetSize())
      );

      Vector< t_DataType> new_vector( VECTOR_LHS);
      return new_vector -= VECTOR_RHS;
    }

    //! @brief divide two vectors of equal size by dividing each corresponding pair of values
    //! @param VECTOR_LHS lhs vector
    //! @param VECTOR_RHS rhs vector
    //! @return vector with all individual divided elements of rhs from lhs vector
    template< typename t_DataType>
    inline
    Vector< t_DataType>
    operator /
    (
      const VectorConstInterface< t_DataType> &VECTOR_LHS,
      const VectorConstInterface< t_DataType> &VECTOR_RHS
    )
    {
      BCL_Assert
      (
        VECTOR_LHS.GetSize() == VECTOR_RHS.GetSize(),
        "vectors of different sizes supplied: " +
        util::Format()( VECTOR_LHS.GetSize()) + " != " + util::Format()( VECTOR_RHS.GetSize())
      );

      Vector< double> new_vector( VECTOR_LHS.GetSize());

      std::transform
      (
        VECTOR_LHS.Begin(), VECTOR_LHS.End(),
        VECTOR_RHS.Begin(),
        new_vector.Begin(),
        std::divides< t_DataType>()
      );

      return new_vector;
    }

    //! @brief multiply two vectors of equal size by building the inner product yielding the scalar product
    //! @param VECTOR_LHS lhs vector
    //! @param VECTOR_RHS rhs vector
    //! @return scalar representing root of inner product of the two ranges
    template< typename t_DataType>
    inline
    t_DataType
    operator *
    (
      const VectorConstInterface< t_DataType> &VECTOR_LHS,
      const VectorConstInterface< t_DataType> &VECTOR_RHS
    )
    {
      BCL_Assert
      (
        VECTOR_LHS.GetSize() == VECTOR_RHS.GetSize(),
        "VectorInterfaces have different sizes of " + util::Format()( VECTOR_LHS.GetSize()) + " and " +
        util::Format()( VECTOR_RHS.GetSize()) + " !"
      );

      // square root of inner product
      return std::inner_product( VECTOR_LHS.Begin(), VECTOR_LHS.End(), VECTOR_RHS.Begin(), t_DataType( 0));
    }

    //! @brief add value to vector
    //! @param VECTOR_LHS lhs vector
    //! @param VALUE_RHS rhs value to be added
    //! @return vector that has the value added to each value of the lhs given vector
    template< typename t_DataType>
    inline
    Vector< t_DataType>
    operator +
    (
      const VectorConstInterface< t_DataType> &VECTOR_LHS,
      const t_DataType &VALUE_RHS
    )
    {
      Vector< t_DataType> new_vector( VECTOR_LHS);
      return new_vector += VALUE_RHS;
    }

    //! @brief add vector to value
    //! @param VALUE_LHS lhs value to be added
    //! @param VECTOR_RHS rhs vector
    //! @return vector that has the value added to each value of the lhs given vector
    template< typename t_DataType>
    inline
    Vector< t_DataType>
    operator +
    (
      const t_DataType &VALUE_LHS,
      const VectorConstInterface< t_DataType> &VECTOR_RHS
    )
    {
      Vector< t_DataType> new_vector( VECTOR_RHS);
      return new_vector += VALUE_LHS;
    }

    //! @brief subtract value from vector
    //! @param VECTOR_LHS lhs vector
    //! @param VALUE_RHS rhs value to be subtracted
    //! @return vector that has the value subtracted from each value of the lhs given vector
    template< typename t_DataType>
    inline
    Vector< t_DataType>
    operator -
    (
      const VectorConstInterface< t_DataType> &VECTOR_LHS,
      const t_DataType &VALUE_RHS
    )
    {
      Vector< t_DataType> new_vector( VECTOR_LHS);
      return new_vector -= VALUE_RHS;
    }

    //! @brief subtract vector from value
    //! @param VALUE_LHS rhs value to be subtracted
    //! @param VECTOR_RHS lhs vector
    //! @return vector that has the values in the vector subtracted from the value
    template< typename t_DataType>
    inline
    Vector< t_DataType>
    operator -
    (
      const t_DataType &VALUE_LHS,
      const VectorConstInterface< t_DataType> &VECTOR_RHS
    )
    {
      Vector< t_DataType> new_vector( VECTOR_RHS.GetSize(), VALUE_LHS);
      return new_vector -= VECTOR_RHS;
    }

    //! @brief multiply scalar with vector
    //! @param SCALAR_LHS lhs value to be multiplied
    //! @param VECTOR_RHS rhs vector
    //! @return vector that has the values multiplied with the scalar
    template< typename t_DataType>
    inline
    Vector< t_DataType>
    operator *
    (
      const t_DataType &SCALAR_LHS,
      const VectorConstInterface< t_DataType> &VECTOR_RHS
    )
    {
      Vector< t_DataType> new_vector( VECTOR_RHS);
      return new_vector *= SCALAR_LHS;
    }

    //! @brief multiply vector with scalar
    //! @param VECTOR_LHS lhs vector
    //! @param SCALAR_RHS rhs value to be multiplied
    //! @return vector that has the values multiplied with the scalar
    template< typename t_DataType>
    inline
    Vector< t_DataType>
    operator *
    (
      const VectorConstInterface< t_DataType> &VECTOR_LHS,
      const t_DataType &SCALAR_RHS
    )
    {
      Vector< t_DataType> new_vector( VECTOR_LHS);
      return new_vector *= SCALAR_RHS;
    }

    //! @brief dividevector with scalar
    //! @param VECTOR_LHS lhs vector
    //! @param SCALAR_RHS rhs value to be divided by
    //! @return vector that has the values divided by the scalar
    template< typename t_DataType>
    inline
    Vector< t_DataType>
    operator /
    (
      const VectorConstInterface< t_DataType> &VECTOR_LHS,
      const t_DataType &SCALAR_RHS
    )
    {
      Vector< t_DataType> new_vector( VECTOR_LHS);
      return new_vector /= SCALAR_RHS;
    }

    //! @brief dividescalar by vector
    //! @param SCALAR_LHS lhs value to be divided
    //! @param VECTOR_RHS rhs vector to be used to divide the scalar
    //! @return vector that has the values of scalar divided by each according value of the vector
    template< typename t_DataType>
    inline
    Vector< t_DataType>
    operator /
    (
      const t_DataType &SCALAR_LHS,
      const VectorConstInterface< t_DataType> &VECTOR_RHS
    )
    {
      Vector< t_DataType> new_vector( VECTOR_RHS.GetSize(), SCALAR_LHS);
      return new_vector /= VECTOR_RHS;
    }

    //! operator Value ^ VectorInterface
    template< typename t_DataType>
    inline
    Vector< t_DataType>
    operator ^
    (
      const t_DataType &SCALAR_LHS,
      const VectorConstInterface< t_DataType> &VECTOR_RHS
    )
    {
      Vector< t_DataType> new_vector( VECTOR_RHS.GetSize(), SCALAR_LHS);

      std::transform
      (
        new_vector.Begin(),
        new_vector.End(),
        VECTOR_RHS.Begin(),
        new_vector.Begin(),
        std::ptr_fun( math::Pow< t_DataType>)
      );

      return new_vector;
    }

  //////////////////////
  // Vector functions //
  //////////////////////

    //! @brief distance between two vectors A->B
    //! @param VECTOR_A vector for position A
    //! @param VECTOR_B vector for position B
    //! @return the length of the vector A->B
    template< typename t_DataType>
    inline
    t_DataType
    SquareDistance
    (
      const VectorConstInterface< t_DataType> &VECTOR_A,
      const VectorConstInterface< t_DataType> &VECTOR_B
    )
    {
      BCL_Assert( VECTOR_A.GetSize() == VECTOR_B.GetSize(), "Can only compute distance between vectors of same size");
      t_DataType total( 0);

      // iterate over vector interfaces
      for
      (
        typename VectorConstInterface< t_DataType>::const_iterator
          itr_a( VECTOR_A.Begin()), itr_a_end( VECTOR_A.End()), itr_b( VECTOR_B.Begin());
        itr_a != itr_a_end;
        ++itr_a, ++itr_b
      )
      {
        total += math::Sqr( *itr_a - *itr_b);
      }

      // take the square root of total
      return total;
    }

    //! @brief distance between two vectors A->B
    //! @param VECTOR_A vector for position A
    //! @param VECTOR_B vector for position B
    //! @return the length of the vector A->B
    template< typename t_DataType>
    inline
    t_DataType
    Distance
    (
      const VectorConstInterface< t_DataType> &VECTOR_A,
      const VectorConstInterface< t_DataType> &VECTOR_B
    )
    {
      // take the square root of total
      return math::Sqrt( SquareDistance( VECTOR_A, VECTOR_B));
    }

    //! @brief calculates projection angle between two Vector3D
    //! @param VECTOR_A first vector (point)
    //! @param VECTOR_B second vector (point)
    //! @return projection angle between two Vector3D
    template< typename t_DataType>
    inline
    t_DataType ProjAngle
    (
      const VectorConstInterface< t_DataType> &VECTOR_A,
      const VectorConstInterface< t_DataType> &VECTOR_B
    )
    {
      const double projection_angle_cosinus( ( VECTOR_A * VECTOR_B) / ( VECTOR_A.Norm() * VECTOR_B.Norm()));

      // check that the result is defined; this is more robust than checking if the vectors are defined, since
      if( !util::IsDefined( projection_angle_cosinus))
      {
        return util::GetUndefined< double>();
      }

      // through numerical drift it could be possible that the value is slightly higher or lower than -1 or 1
      // (VECTOR_A * VECTOR_B) / ( Norm( VECTOR_A) * Norm( VECTOR_B)) is the actual cos angle between vectors ab and cd
      if( projection_angle_cosinus >= double( 1.0))
      {
        // acos( 1) == 0.0
        return 0.0;
      }
      else if( projection_angle_cosinus <= double( -1.0))
      {
        // acos( -1) == pi
        return math::g_Pi;
      }

      // -1 < projection_angle_cosinus < 1, so calling acos( projection_angle_cosinus) yields the proper angle
      return std::acos( projection_angle_cosinus);
    }

    //! @brief calculates projection angle from three points (A->B and A->C)
    //! @param VECTOR_A first vector (point)
    //! @param VECTOR_B second vector (point)
    //! @param VECTOR_C third vector (point)
    //! @return projection angle from three points
    template< typename t_DataType>
    inline
    t_DataType ProjAngle
    (
      const VectorConstInterface< t_DataType> &VECTOR_A,
      const VectorConstInterface< t_DataType> &VECTOR_B,
      const VectorConstInterface< t_DataType> &VECTOR_C
    )
    {
      return ProjAngle( VECTOR_B - VECTOR_A, VECTOR_C - VECTOR_A);
    }

    //! @brief calculates projection angle from four points (A->B and C->D)
    //! @param VECTOR_A first vector (point)
    //! @param VECTOR_B second vector (point)
    //! @param VECTOR_C third vector (point)
    //! @param VECTOR_D fourth vector (point)
    //! @return projection angle from four points
    template< typename t_DataType>
    inline
    t_DataType ProjAngle
    (
      const VectorConstInterface< t_DataType> &VECTOR_A,
      const VectorConstInterface< t_DataType> &VECTOR_B,
      const VectorConstInterface< t_DataType> &VECTOR_C,
      const VectorConstInterface< t_DataType> &VECTOR_D
    )
    {
      return ProjAngle( VECTOR_B - VECTOR_A, VECTOR_D - VECTOR_C);
    }

    //! scalar product of two VectorInterfaces
    template< typename t_DataType>
    t_DataType
    ScalarProduct
    (
      const VectorConstInterface< t_DataType> &VECTOR_A,
      const VectorConstInterface< t_DataType> &VECTOR_B
    )
    {
      BCL_Assert( VECTOR_A.GetSize() == VECTOR_B.GetSize(), "ScalarProduct: Vectors must have = size for computation");
      return std::inner_product( VECTOR_A.Begin(), VECTOR_A.End(), VECTOR_B.Begin(), t_DataType( 0));
    }

    //! @brief square the elements of a vector
    //! @param VECTOR the vector to square
    //! @return the squared vector
    template< typename t_DataType>
    Vector< t_DataType> SqrVector( const VectorConstInterface< t_DataType> &VECTOR)
    {
      Vector< t_DataType> result( VECTOR.GetSize());

      t_DataType *itr_result( result.Begin());

      // iterate through all elements and store squared value
      for
      (
        const t_DataType *itr_vec_a( VECTOR.Begin()), *itr_vec_a_end( VECTOR.End());
        itr_vec_a != itr_vec_a_end;
        ++itr_vec_a, ++itr_result
      )
      {
        *itr_result = math::Sqr( *itr_vec_a);
      }

      // return squared vector
      return result;
    }

    //! @brief log of elements of the vector
    //! @param VECTOR the vector whose log is needed
    //! @return log of elements of the vector
    template< typename t_DataType>
    Vector< t_DataType> LogVector( const VectorConstInterface< t_DataType> &VECTOR)
    {
      Vector< t_DataType> result( VECTOR.GetSize());

      t_DataType *itr_result( result.Begin());

      // iterate through all elements and store squared value
      for
      (
        const t_DataType *itr_vec_a( VECTOR.Begin()), *itr_vec_a_end( VECTOR.End());
        itr_vec_a != itr_vec_a_end;
        ++itr_vec_a, ++itr_result
      )
      {
        *itr_result = log( *itr_vec_a);
      }

      // return squared vector
      return result;
    }

    //! @brief Multiply each element of one vector by another, in-place
    //! @param VECTOR_A the vector to transform
    //! @param VECTOR_B the vector whose elements should multiply all elements in vector A
    //! @return the squared vector
    template< typename t_DataType>
    void ElementwiseMultiply( VectorInterface< t_DataType> &VECTOR_A, const VectorConstInterface< t_DataType> &VECTOR_B)
    {
      BCL_Assert( VECTOR_A.GetSize() == VECTOR_B.GetSize(), "Vectors need the same size for elementwise multiply");
      t_DataType *itr_result( VECTOR_A.Begin());

      // iterate through all elements and store squared value
      for
      (
        const t_DataType *itr_b( VECTOR_B.Begin()), *itr_b_end( VECTOR_B.End());
        itr_b != itr_b_end;
        ++itr_b, ++itr_result
      )
      {
        *itr_result *= *itr_b;
      }
    }

    //! @brief Compute the L1 norm (absolute sum) of a vector
    //! @param VECTOR_A the vector to transform
    //! @return the squared vector
    template< typename t_DataType>
    t_DataType AbsSum( const VectorConstInterface< t_DataType> &VECTOR)
    {
      t_DataType sum( 0);

      // iterate through all elements and store squared value
      for
      (
        const t_DataType *itr( VECTOR.Begin()), *itr_end( VECTOR.End());
        itr != itr_end;
        ++itr
      )
      {
        sum += math::Absolute( *itr);
      }
      return sum;
    }

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_VECTOR_OPERATIONS_H_
