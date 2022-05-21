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
#include "math/bcl_math_quaternion.h"

// includes from bcl - sorted alphabetically
#include "random/bcl_random_uniform_distribution.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Quaternion::s_Instance
    (
      GetObjectInstances().AddInstance( new Quaternion())
    );

    Quaternion Conjugate( const Quaternion &QUATERNION)
    {
      Quaternion quat;
      quat[ 0] = QUATERNION[ 0];
      for( size_t i( 1); i < Quaternion::s_Size; ++i)
      {
        quat[ i] = -QUATERNION[ i];
      }
      return quat;
    }

  ///////////////
  // operators //
  ///////////////

    Quaternion &Quaternion::operator *=( const double T)
    {
      std::transform( m_Values, m_Values + s_Size, m_Values, std::bind1st( std::multiplies< double>(), T));

      return *this;
    }

    Quaternion &Quaternion::operator =( const Quaternion &QUATERNION)
    {
      std::copy( QUATERNION.m_Values, QUATERNION.m_Values + s_Size, m_Values);

      return *this;
    }

    Quaternion &Quaternion::operator /=( const double T)
    {
      std::transform( m_Values, m_Values + s_Size, m_Values, std::bind1st( std::multiplies< double>(), double( 1) / T));

      return *this;
    }

  ////////////////
  // operations //
  ////////////////

    // J.J. Kuffner "Effective Sampling and Distance Metrices for 3D Rigid Body Path Planning" 2004 IEEE Int'l Conf. on Robotics and Automation (ICRA 2004)
    //! set random rotation
    Quaternion &Quaternion::SetRand()
    {
      double s( random::GetGlobalRandom().Random< double>( double( 1))),
             sigma1( Sqrt( s - 1)),
             sigma2( Sqrt( s)),
             theta1( 2 * g_Pi * random::GetGlobalRandom().Random< double>( double( 1))),
             theta2( 2 * g_Pi * random::GetGlobalRandom().Random< double>( double( 1)));

      m_Values[ 0] = std::cos( theta2) * sigma2;
      m_Values[ 1] = std::sin( theta1) * sigma1;
      m_Values[ 2] = std::cos( theta1) * sigma1;
      m_Values[ 3] = std::cos( theta2) * sigma2;

      return *this;
    }

    Quaternion &Quaternion::Equal( const double R, const double X, const double Y, const double Z)
    {
      m_Values[ 0] = R;
      m_Values[ 1] = X;
      m_Values[ 2] = Y;
      m_Values[ 3] = Z;

      return *this;
    }

    Quaternion &Quaternion::Equal( const double R, const linal::Vector3D &VECTOR3D)
    {
      m_Values[ 0] = R;
      std::copy( VECTOR3D.Begin(), VECTOR3D.End(), m_Values + 1);

      return *this;
    }

    double Quaternion::Norm() const
    {
      double tmp( 0);
      for( size_t i( 0); i < s_Size; ++i)
      {
        tmp += m_Values[ i] * m_Values[ i];
      }
      return Sqrt( tmp);
    }

  //////////////////////
  // input and output //
  //////////////////////

    std::istream &Quaternion::Read( std::istream &ISTREAM)
    {
      // read all values
      for( double *ptr( m_Values), *ptr_end( m_Values + s_Size); ptr != ptr_end; ++ptr)
      {
        ISTREAM >> *ptr;
      }

      return ISTREAM;
    }

    std::ostream &Quaternion::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // insert indent
      io::Serialize::InsertIndent( OSTREAM, INDENT);

      // write all values
      for( const double *ptr( m_Values), *ptr_end( m_Values + s_Size); ptr != ptr_end; ++ptr)
      {
        OSTREAM << *ptr << '\t';
      }

      // end
      return OSTREAM;
    }

    Quaternion &Quaternion::VectorToQuaternion( const linal::Vector3D &VECTOR3D)
    {
      m_Values[ 0] = 0;
      std::copy( VECTOR3D.Begin(), VECTOR3D.End(), m_Values + 1);

      return *this;
    }

    Quaternion operator *( const Quaternion &QUAT_1, const Quaternion &QUAT_2)
    {
      Quaternion quat;

      quat[ 0] = QUAT_1[ 0] * QUAT_2[ 0] - QUAT_1[ 1] * QUAT_2[ 1] - QUAT_1[ 2] * QUAT_2[ 2] - QUAT_1[ 3] * QUAT_2[ 3];
      quat[ 1] = QUAT_1[ 0] * QUAT_2[ 1] + QUAT_1[ 1] * QUAT_2[ 0] + QUAT_1[ 2] * QUAT_2[ 3] - QUAT_1[ 3] * QUAT_2[ 2];
      quat[ 2] = QUAT_1[ 0] * QUAT_2[ 2] + QUAT_1[ 2] * QUAT_2[ 0] + QUAT_1[ 3] * QUAT_2[ 1] - QUAT_1[ 1] * QUAT_2[ 3];
      quat[ 3] = QUAT_1[ 0] * QUAT_2[ 3] + QUAT_1[ 3] * QUAT_2[ 0] + QUAT_1[ 1] * QUAT_2[ 2] - QUAT_1[ 2] * QUAT_2[ 1];

      // end
      return quat;
    }

    Quaternion operator *( const linal::Vector3D &VECTOR3D, const Quaternion &QUATERNION)
    {
      Quaternion quat;
      quat[ 0] = 0;
      for( size_t i = 0; i < 3; ++i)
      {
        quat[ i + 1] = VECTOR3D( i);
      }

      quat = quat * QUATERNION;

      return quat;
    }

  } // namespace math
} // namespace bcl
