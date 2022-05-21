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
#include "assemble/bcl_assemble_sse_geometry.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_sse_geometry.cpp
  //!
  //! @author weinerbe
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleSSEGeometry :
    public ExampleInterface
  {
  public:

    ExampleAssembleSSEGeometry *Clone() const
    {
      return new ExampleAssembleSSEGeometry( *this);
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
      // create string "pdb_filename" which has path for example pdb file
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      // create ProteinModel "protein_model" from "pdb"
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 9;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      BCL_MessageStd( "testing default constructor");
      const assemble::SSEGeometry def_construct;
      BCL_Example_Check
      (
        def_construct.GetOrientation() == math::TransformationMatrix3D() &&
        def_construct.GetGeometries().GetSize() == 0,
        "Default constructor failed"
      );

      // test geometry interface constructor
      BCL_MessageStd( "testing construct from GeometryInterface");
      assemble::SSE helix_sse( *protein_model.GetSSEs( biol::GetSSTypes().HELIX).FirstElement());
      assemble::SSEGeometry geometry( helix_sse);
      BCL_Example_Check
      (
        geometry.GetOrientation() == helix_sse.GetOrientation() &&
        geometry.GetType() == helix_sse.GetType() &&
        math::EqualWithinTolerance( geometry.GetMainAxis().GetLength(), helix_sse.GetMainAxis().GetLength()),
        "Constructor from GeometryInterface gave: " + util::Format()( geometry) +
        " instead of " + util::Format()( helix_sse)
      );

      // test copy constructor and == operator implicitly
      BCL_MessageStd( "testing copy constructor");
      assemble::SSEGeometry copy_construct( geometry);
      BCL_Example_Check
      (
        copy_construct == geometry,
        "Copy constructor gave: " + util::Format()( copy_construct) + " instead of " + util::Format()( geometry)
      );

      // test Clone function
      BCL_MessageStd( "testing Clone function");
      util::ShPtr< assemble::SSEGeometry> clone_construct( geometry.Clone());
      BCL_Example_Check
      (
        *clone_construct == geometry,
        "Clone function gave: " + util::Format()( *clone_construct) + " instead of " + util::Format()( geometry)
      );

    /////////////////
    // data access //
    /////////////////

      // check GetStaticClassName
      const std::string correct_static_class_name( "bcl::assemble::SSEGeometry");
      BCL_Example_Check
      (
        GetStaticClassName< assemble::SSEGeometry>() == correct_static_class_name,
        "GetStaticClassName gives " + GetStaticClassName< assemble::SSEGeometry>() + " but should give " +
        correct_static_class_name
      );

      // check GetClassIdentifier
      BCL_Example_Check
      (
        GetStaticClassName< assemble::SSEGeometry>() == clone_construct->GetClassIdentifier(),
        "GetClassIdentifier gives " + clone_construct->GetClassIdentifier() + " but should give " +
        correct_static_class_name
      );

      // check GetGeometries
      const size_t known_number_fragments( 8);
      BCL_Example_Check
      (
        geometry.GetGeometries().GetSize() == known_number_fragments,
        "GetGeometries size is " + util::Format()( geometry.GetGeometries().GetSize()) +
        " but is " + util::Format()( known_number_fragments)
      );

      // check GetMainAxis
      const double known_length( 18.0296);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( geometry.GetMainAxis().GetLength(), known_length),
        "GetMainAxis has a length of : " + util::Format()( geometry.GetMainAxis().GetLength()) + " but should be: " +
        util::Format()( known_length)
      );

      // check GetExtent
      const double known_extent( helix_sse.GetExtent( coord::GetAxes().e_X));
      BCL_Example_Check
      (
        geometry.GetExtent( coord::GetAxes().e_X) == known_extent,
        "GetExtent gives : " + util::Format()( geometry.GetExtent( coord::GetAxes().e_X)) + " but should be: " +
        util::Format()( known_extent)
      );

      // check GetRadialExtent
      BCL_Example_Check
      (
        geometry.GetRadialExtent() == known_extent,
        "GetRadialExtent gives : " + util::Format()( geometry.GetRadialExtent()) + " but should be: " +
        util::Format()( known_extent)
      );

      // check GetOrientation
      const math::TransformationMatrix3D known_orientation( helix_sse.GetOrientation());
      BCL_Example_Check
      (
        geometry.GetOrientation() == known_orientation,
        "GetOrientation gives : " + util::Format()( geometry.GetOrientation()) + " but should be: " +
        util::Format()( known_orientation)
      );

      // check GetCenter
      const linal::Vector3D known_center( helix_sse.GetCenter());
      BCL_Example_Check
      (
        geometry.GetCenter() == known_center,
        "GetCenter gives : " + util::Format()( geometry.GetCenter()) + " but should be: " +
        util::Format()( known_center)
      );

      // check GetAxis
      const linal::Vector3D known_axis( helix_sse.GetAxis( coord::GetAxes().e_Z));
      BCL_Example_Check
      (
        geometry.GetAxis( coord::GetAxes().e_Z) == known_axis,
        "GetAxis gives : " + util::Format()( geometry.GetAxis( coord::GetAxes().e_Z)) + " but should be: " +
        util::Format()( known_axis)
      );

      // check GetSSType
      BCL_Example_Check
      (
        geometry.GetType() == biol::GetSSTypes().HELIX,
        "GetSSType gives : " + util::Format()( geometry.GetType()) + " but should be: " +
        util::Format()( biol::GetSSTypes().HELIX)
      );

    ////////////////
    // operations //
    ////////////////

      // check Translate
      const linal::Vector3D translate( 1.0, 1.0, 1.0);
      geometry.Translate( translate);
      BCL_Example_Check
      (
        geometry.GetCenter() == known_center + translate,
        "GetCenter (after translate) gives : " + util::Format()( geometry.GetCenter()) + " but should be: " +
        util::Format()( known_center + translate)
      );

      // check Transform
      math::RotationMatrix3D rotation( geometry.GetAxis( coord::GetAxes().e_Z), math::g_Pi);
      math::TransformationMatrix3D transform;
      transform( -known_center - translate);
      transform( rotation);
      transform( known_center);
      geometry.Transform( transform);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( geometry.GetCenter(), known_center),
        "GetCenter (after transform) gives : " + util::Format()( geometry.GetCenter()) + " but should be: " +
        util::Format()( known_center)
      );

      // check Rotate
      math::RotationMatrix3D inverse_rotation( geometry.GetAxis( coord::GetAxes().e_Z), math::g_Pi);
      geometry.Translate( -known_center);
      geometry.Rotate( inverse_rotation);
      geometry.Translate( known_center);
      BCL_Example_Check
      (
        math::EqualWithinTolerance
        (
          geometry.GetOrientation().GetAxis( coord::GetAxes().e_X), known_orientation.GetAxis( coord::GetAxes().e_X)
        ),
        "GetOrientation (after rotate) gives : " + util::Format()( geometry.GetOrientation()) + " but should be: " +
        util::Format()( known_orientation)
      );

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities");
      // write the object
      WriteBCLObject( geometry);

      // read the object back in
      assemble::SSEGeometry read_geometry;
      ReadBCLObject( read_geometry);

      // compare the objects
      BCL_MessageStd( "compare written and read objects");
      BCL_Example_Check
      (
        geometry.GetType() == read_geometry.GetType() &&
        math::EqualWithinTolerance
        (
          geometry.GetOrientation().GetAxis( coord::GetAxes().e_Z), read_geometry.GetAxis( coord::GetAxes().e_Z)
        ) &&
        math::EqualWithinTolerance( geometry.GetMainAxis().GetLength(), read_geometry.GetMainAxis().GetLength()) &&
        geometry.GetGeometries().GetSize() == read_geometry.GetGeometries().GetSize(),
        "the written and read FoldTemplate classes differ from each other \n" +
        util::Format()( geometry) + "\nvs\n" + util::Format()( read_geometry)
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleSSEGeometry

  const ExampleClass::EnumType ExampleAssembleSSEGeometry::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleSSEGeometry())
  );

} // namespace bcl
