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
#include "coord/bcl_coord_sphere.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_point_cloud.h"
#include "linal/bcl_linal_vector_3d_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_coord_sphere.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCoordSphere :
    public ExampleInterface
  {
  public:

    ExampleCoordSphere *Clone() const
    {
      return new ExampleCoordSphere( *this);
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
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct default sphere
      coord::Sphere sphere_default;

      // construct from position and radius
      const linal::Vector3D position( 1.5, 1.5, 1.5);
      const double radius( 3.0);
      coord::Sphere sphere_constr( position, radius);

      // copy constructor
      coord::Sphere sphere_copy( sphere_constr);

      // clone
      util::ShPtr< util::ObjectInterface> ptr( sphere_constr.Clone());

    /////////////////
    // data access //
    /////////////////

      // check the class identifier
      BCL_ExampleCheck( GetStaticClassName( sphere_default), "bcl::coord::Sphere");

      // class identifier
      BCL_ExampleCheck( ptr->GetClassIdentifier(), GetStaticClassName< coord::Sphere>());

      // position
      BCL_MessageStd( "this is position of constructed sphere: " + util::Format()( sphere_constr.GetPosition()));
      BCL_ExampleCheck( sphere_constr.GetPosition(), position);

      // radius
      BCL_MessageStd( "this is radius of constructed sphere: " + util::Format()( sphere_constr.GetRadius()));
      BCL_ExampleCheck( sphere_constr.GetRadius(), radius);

      // check default sphere
      BCL_ExampleCheck( sphere_default.GetPosition(), linal::Vector3D( 0, 0, 0));
      BCL_ExampleCheck( sphere_default.GetRadius(), double( 0));

      // check copy constructed
      BCL_ExampleCheck( sphere_copy.GetPosition(), sphere_constr.GetPosition());
      BCL_ExampleCheck( sphere_copy.GetRadius(), sphere_constr.GetRadius());

    ////////////////
    // operations //
    ////////////////

      // interpolate surface with points
      const size_t number_interpolation_points( 100);
      const storage::Vector< linal::Vector3D> points_latitudes( sphere_constr.PointsOnSurfaceByLatiudes( number_interpolation_points));
      const storage::Vector< linal::Vector3D> points_random( sphere_constr.PointsOnSurfaceRandom( number_interpolation_points));
      const storage::Vector< linal::Vector3D> points_spiral( sphere_constr.PointsOnSurfaceSpiral( number_interpolation_points));

      {
        const coord::PointCloud pc_latidues( points_latitudes);
        const std::string latidues_filename( AddExampleOutputPathToFilename( sphere_constr, "sphere_latidues_points.pdb"));
        pc_latidues.WriteToPDB( latidues_filename);

        const coord::PointCloud pc_random( points_random);
        const std::string random_filename( AddExampleOutputPathToFilename( sphere_constr, "sphere_random_points.pdb"));
        pc_random.WriteToPDB( random_filename);

        const coord::PointCloud pc_spiral( points_spiral);
        const std::string spiral_filename( AddExampleOutputPathToFilename( sphere_constr, "sphere_spiral_points.pdb"));
        pc_spiral.WriteToPDB( spiral_filename);
      }

      // check the number
      BCL_ExampleCheck( points_latitudes.IsEmpty(), false); // number of points for latitudes is not equal to the desired number
      BCL_ExampleCheck( points_random.GetSize(), number_interpolation_points);
      BCL_ExampleCheck( points_spiral.GetSize(), number_interpolation_points);

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities for coord::Sphere");
      // open output file
      BCL_MessageVrb( "write object");
      WriteBCLObject( sphere_constr);
      BCL_MessageVrb( "read object");
      coord::Sphere sphere_read;
      ReadBCLObject( sphere_read);

      // compare the objects
      BCL_Example_Check
      (
        sphere_constr.GetPosition() == sphere_read.GetPosition() &&
        sphere_constr.GetRadius()   == sphere_read.GetRadius(),
        "read sphere is different from written"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleCoordSphere

  const ExampleClass::EnumType ExampleCoordSphere::s_Instance
  (
    GetExamples().AddEnum( ExampleCoordSphere())
  );

} // namespace bcl
