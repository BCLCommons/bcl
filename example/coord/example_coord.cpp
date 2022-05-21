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
#include "coord/bcl_coord.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_point_cloud.h"
#include "linal/bcl_linal_vector_3d.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_coord.cpp
  //!
  //! @author mendenjl
  //! @date Mar 12, 2014
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCoord :
    public ExampleInterface
  {
  public:

    ExampleCoord *Clone() const
    {
      return new ExampleCoord( *this);
    }

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

      // create a point cloud for testing various functions
      // just create a 2D box
      coord::PointCloud rectangle_far_from_origin;
      rectangle_far_from_origin.PushBack( linal::Vector3D( 70.0, 30.0, 0.0));
      rectangle_far_from_origin.PushBack( linal::Vector3D( 100.0, 30.0, 0.0));
      rectangle_far_from_origin.PushBack( linal::Vector3D( 70.0, 20.0, 0.0));
      rectangle_far_from_origin.PushBack( linal::Vector3D( 100.0, 20.0, 0.0));

      BCL_ExampleCheckWithinAbsTolerance
      (
        coord::CenterOfMass
        (
          util::SiPtrVector< const linal::Vector3D>
          (
            rectangle_far_from_origin.Begin(),
            rectangle_far_from_origin.End()
          )
        ),
        linal::Vector3D( 85.0, 25.0, 0.0),
        0.001
      );

      // apply a translation to the origin and a rotation; the center of mass should change, but the radius
      // of gyration should not
      coord::PointCloud rectangle_at_origin( rectangle_far_from_origin);
      rectangle_at_origin.Translate( -rectangle_at_origin.GetCenter());
      rectangle_at_origin.Rotate( math::RotationMatrix3D( 1.0, 2.0, 3.0));

      BCL_ExampleCheckWithinAbsTolerance
      (
        coord::CenterOfMass
        (
          util::SiPtrVector< const linal::Vector3D>
          (
            rectangle_at_origin.Begin(),
            rectangle_at_origin.End()
          )
        ),
        linal::Vector3D( 0.0, 0.0, 0.0),
        0.001
      );

      // test that radius of gyration is the same no matter where the rectangle is
      BCL_ExampleCheckWithinAbsTolerance
      (
        coord::RadiusOfGyration
        (
          util::SiPtrVector< const linal::Vector3D>
          (
            rectangle_at_origin.Begin(),
            rectangle_at_origin.End()
          )
        ),
        coord::RadiusOfGyration
        (
          util::SiPtrVector< const linal::Vector3D>
          (
            rectangle_far_from_origin.Begin(),
            rectangle_far_from_origin.End()
          )
        ),
        0.001
      );

      // radius of gyration for a rectangle is width * length ^ 3 / 12
      BCL_ExampleCheckWithinAbsTolerance
      (
        coord::RadiusOfGyration
        (
          util::SiPtrVector< const linal::Vector3D>
          (
            rectangle_at_origin.Begin(),
            rectangle_at_origin.End()
          )
        ),
        15.8114,
        0.001
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleCoord

  const ExampleClass::EnumType ExampleCoord::s_Instance
  (
    GetExamples().AddEnum( ExampleCoord())
  );

} // namespace bcl
