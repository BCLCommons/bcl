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

#ifndef BCL_MATH_ROTATION_MATRIX_2D_H_
#define BCL_MATH_ROTATION_MATRIX_2D_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically
#include "coord/bcl_coord.fwd.hh"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix2x2.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RotationMatrix2D
    //! @brief RotationMatrix2D for rotating a Vector in 2D space
    //! @details a 2x2 matrix contains all the euler angle derived coefficient, to rotate a coordinate in space around
    //!          the origin by a simple matrix vector multiplication
    //!
    //! @see @link example_math_rotation_matrix_2d.cpp @endlink
    //! @author mendenjl
    //! @date Jun 12, 2014
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RotationMatrix2D :
      public util::ObjectInterface
    {
    //////////
    // data //
    //////////

    private:

      linal::Matrix2x2< double> m_RotationMatrix2D; //!< RotationMatrix2D

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! construct empty RotationMatrix2D
      RotationMatrix2D();

      //! construct RotationMatrix2D from RotationAngle
      RotationMatrix2D( const double &ANGLE);

      //! construct RotationMatrix2D from a 2x2 matrix
      RotationMatrix2D( const linal::Matrix2x2< double> &ROTATION);

      //! copy constructor
      RotationMatrix2D *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! return the matrix
      const linal::Matrix2x2< double> &GetMatrix() const;

      //! @brief return the angle that this matrix rotates by
      //! @return rotation angle that this matrix represents (in radians)
      double GetAngle() const;

      //! @brief set the angle that this matrix rotates by
      //! @param THETA rotation angle that this matrix represents (in radians)
      void SetAngle( const double &THETA);

      //! @brief set given cos(theta) and sin(theta)
      //! @param COS cosine of theta
      //! @param SIN sine theta
      void SetFromCosSin( const double &COS, const double &SIN);

      //! @brief reverse the angle <-> matrix inversion <-> transposition for 2d rotation matrices
      void Invert()
      {
        m_RotationMatrix2D.Transpose();
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief set random rotation from 0 - MAX_ROTATIONANGLE
      //! @param MAX_ROTATIONANGLE the maximum rotation angle
      RotationMatrix2D &SetRand( const double MAX_ROTATIONANGLE = g_Pi);

      //! @brief set the matrix up to be a givens rotation matrix; which is a numerically stable solution used in tridiagonalization
      //! @param A, B the vector to rotate in the givens rotation matrix problem:
      //! | c -s || A | = | r |
      //! | s  c || B |   | 0 |
      //! @return R (equal to pythag(A,B))
      double MakeGivens( const double &A, const double &B);

    ///////////////
    // operators //
    ///////////////

      //! operator *= RotationMatrix2D
      RotationMatrix2D &operator *=( const RotationMatrix2D &MATRIX);

    //////////////////////
    // input and output //
    //////////////////////

      //! write RotationMatrix2D to std::ostream using the given util::Format
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT, const util::Format &FORMAT) const;

    protected:

      //! write RotationMatrix2D to std::ostream using the given util::Format
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! read RotationMatrix2D from io::IFStream
      std::istream &Read( std::istream &ISTREAM);

    }; //class RotationMatrix2D

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_ROTATION_MATRIX_2D_H_
