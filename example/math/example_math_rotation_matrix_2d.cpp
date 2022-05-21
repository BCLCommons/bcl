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
#include "math/bcl_math_rotation_matrix_2d.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_operations.h"
#include "linal/bcl_linal_vector_2d.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_rotation_matrix_2d.cpp
  //!
  //! @author mendenjl
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathRotationMatrix2D :
    public ExampleInterface
  {
  public:

    ExampleMathRotationMatrix2D *Clone() const
    { return new ExampleMathRotationMatrix2D( *this);}

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
    /////////////////
    // constructor //
    /////////////////

      //default constructor
      const math::RotationMatrix2D rotationmatrix2D_default;

      // test the constructor with a few different angles
      BCL_ExampleCheckWithinAbsTolerance( math::RotationMatrix2D( 1.1).GetAngle(), 1.1, 1.0e-6);
      BCL_ExampleCheckWithinAbsTolerance( math::RotationMatrix2D( math::g_Pi).GetAngle(), math::g_Pi, 1.0e-6);
      BCL_ExampleCheckWithinAbsTolerance( math::RotationMatrix2D( 0).GetAngle(), 0.0, 1.0e-6);

      math::RotationMatrix2D rotationmatrix2D_givens;

      // test make givens
      linal::Vector2D arbitrary_vector( 0.1, 1.5);
      BCL_ExampleCheckWithinTolerance
      (
        rotationmatrix2D_givens.MakeGivens( arbitrary_vector( 0), arbitrary_vector( 1)),
        arbitrary_vector.Norm(),
        0.0001
      );
      BCL_ExampleCheckWithinTolerance
      (
        rotationmatrix2D_givens.GetMatrix() * arbitrary_vector,
        linal::Vector2D( arbitrary_vector.Norm(), 0.0),
        0.0001
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMathRotationMatrix2D

  const ExampleClass::EnumType ExampleMathRotationMatrix2D::s_Instance
  (
    GetExamples().AddEnum( ExampleMathRotationMatrix2D())
  );

} // namespace bcl
