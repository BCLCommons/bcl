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

#ifndef BCL_MATH_MUTATE_TRANSFORMATION_MATRIX_3D_H_
#define BCL_MATH_MUTATE_TRANSFORMATION_MATRIX_3D_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_mutate_interface.h"
#include "bcl_math_mutate_result.h"
#include "bcl_math_rotation_matrix_3d.h"
#include "bcl_math_transformation_matrix_3d.h"
#include "linal/bcl_linal_vector_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateTransformationMatrix3D
    //! @brief This class is a Mutate class for transformation matrices for a Monte Carlo based Minimization
    //! @details When constructing this class the user can decide on a maximum translation and rotation performed
    //! in every mutating step
    //!
    //! @see @link example_math_mutate_transformation_matrix_3d.cpp @endlink
    //! @author woetzen
    //! @date 01.04.2007
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateTransformationMatrix3D :
      public MutateInterface< TransformationMatrix3D>
    {

    private:

    //////////
    // data //
    //////////

      double m_MaxTranslation;   //!< that is the maximum translation to be performed
      double m_MaxRotationAngle; //!< that is the maximum rotation angle in rad

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor allows no translation nor rotation
      MutateTransformationMatrix3D() :
        m_MaxTranslation( 0.0),
        m_MaxRotationAngle( 0.0)
      {
      }

      //! default constructor fro MAX_TRANSALTION and MAX_ROTATION_ANGLE
      MutateTransformationMatrix3D( const double MAX_TRANSALTION, const double MAX_ROTATION_ANGLE) :
        m_MaxTranslation( MAX_TRANSALTION),
        m_MaxRotationAngle( MAX_ROTATION_ANGLE)
      {
      }

      //! virtual copy constructor
      virtual MutateTransformationMatrix3D *Clone() const
      {
        return new MutateTransformationMatrix3D( *this);
      }

      //! returns class name
      //! the class name as const ref std::string
      virtual const std::string &GetClassIdentifier() const
      {
         return GetStaticClassName( *this);
      }

    ///////////////
    // operators //
    ///////////////

      //! virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
      virtual MutateResult< TransformationMatrix3D> operator()
      (
        const TransformationMatrix3D &TRANSFORMATION_MATRIX3D
      ) const
      {
        // make a copy of argument Transformationmatrix
        util::ShPtr< TransformationMatrix3D> new_matrix( TRANSFORMATION_MATRIX3D.Clone());

        // apply random rotation determined with max rotation angle
        new_matrix->operator()( RotationMatrix3D().SetRand( m_MaxRotationAngle));

        // apply random translation not going further than max translation
        new_matrix->operator()( linal::Vector3D().SetRandomTranslation( m_MaxTranslation));

        // return mutated matrix
        return MutateResult< TransformationMatrix3D>( new_matrix, *this);
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      virtual std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    }; //class MutateTransformationMatrix3D

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_MUTATE_TRANSFORMATION_MATRIX_3D_H_
