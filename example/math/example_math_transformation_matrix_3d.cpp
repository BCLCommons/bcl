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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "math/bcl_math_transformation_matrix_3d.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"
#include "math/bcl_math_rotation_matrix_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_transformation_matrix_3d.cpp
  //!
  //! @author woetzen
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathTransformationMatrix3D :
    public ExampleInterface
  {
  public:

    ExampleMathTransformationMatrix3D *Clone() const
    { return new ExampleMathTransformationMatrix3D( *this);}

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

      const double vector6elements[] = { math::g_Pi, math::g_Pi * 0.5, 0.1, 5.5, 6.6, 7.8};

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      math::TransformationMatrix3D transformation_default;

      // construct from vector with 6 elements
      math::TransformationMatrix3D transformation_6elements( linal::Vector< double>( 6, vector6elements));

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

      // apply rotation and transformation
      const math::RotationMatrix3D rotation( &coord::GetAxes().e_X, vector6elements);
      const linal::Vector3D         translation( vector6elements + 3);

      transformation_default( rotation);
      transformation_default( translation);

      // compare the transformation
      BCL_ExampleCheck( SimilarWithinTolerance( transformation_6elements, transformation_default, 0.0001, 0.0001), true);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMathTransformationMatrix3D

  const ExampleClass::EnumType ExampleMathTransformationMatrix3D::s_Instance
  (
    GetExamples().AddEnum( ExampleMathTransformationMatrix3D())
  );

} // namespace bcl
