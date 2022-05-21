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
#include "coord/bcl_coord_point_to_key_spherical_radius.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_coord_point_to_key_spherical_radius.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCoordPointToKeySphericalRadius :
    public ExampleInterface
  {
  public:

    ExampleCoordPointToKeySphericalRadius *Clone() const
    {
      return new ExampleCoordPointToKeySphericalRadius( *this);
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

      const size_t s_quantization_resolution_10( 9);
      const size_t s_quantization_resolution_15( 15);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct with quantization resolution
      coord::PointToKeySphericalRadius ptk_sr_10( s_quantization_resolution_10);

      // copy constructor
      coord::PointToKeySphericalRadius ptk_sr_copy( ptk_sr_10);

      // clone
      util::ShPtr< coord::PointToKeyInterface> sp_ptk_sr( ptk_sr_10.Clone());

      BCL_Example_Check
      (
        ptk_sr_10.GetAngularResolution() == s_quantization_resolution_10
        && ptk_sr_copy.GetAngularResolution() == s_quantization_resolution_10
        && sp_ptk_sr->GetAngularResolution() == s_quantization_resolution_10,
        "problem in constructing and copying"
      );

    /////////////////
    // data access //
    /////////////////

      // class name and identifier
      BCL_MessageStd( "class name: " + sp_ptk_sr->GetClassIdentifier());
      BCL_Example_Check
      (
        GetStaticClassName< coord::PointToKeySphericalRadius>() == "bcl::coord::PointToKeySphericalRadius"
        && sp_ptk_sr->GetClassIdentifier() == GetStaticClassName< coord::PointToKeySphericalRadius>(),
          "incorrect class name: static class name: " + GetStaticClassName< coord::PointToKeySphericalRadius>()
        + " class identifier: " + sp_ptk_sr->GetClassIdentifier()
      );

      // get the resolution
      BCL_MessageStd( "quantization resolution: " + util::Format()( ptk_sr_10.GetAngularResolution()));
      BCL_Example_Check
      (
        ptk_sr_10.GetAngularResolution() == s_quantization_resolution_10,
        "incorrect resolution: " + util::Format()( ptk_sr_10.GetAngularResolution()) +
        " != " + util::Format()( s_quantization_resolution_10)
      );

      // set the resolution
      sp_ptk_sr->SetAngularResolution( s_quantization_resolution_15);
      BCL_Example_Check
      (
        sp_ptk_sr->GetAngularResolution() == s_quantization_resolution_15,
        "incorrect resolution: " + util::Format()( sp_ptk_sr->GetAngularResolution()) +
        " != " + util::Format()( s_quantization_resolution_15)
      );

    ///////////////
    // operators //
    ///////////////

      const linal::Vector3D north_pole( 0, 0, 1);
      const linal::Vector3D equator( 1, 0, 0);
      const linal::Vector3D south_pole( 0, 0, -1);

      const storage::Triplet< int, int, int> north_pole_quantized( ptk_sr_10( north_pole));
      const storage::Triplet< int, int, int> equator_quantized( ptk_sr_10( equator));
      const storage::Triplet< int, int, int> south_pole_quantized( ptk_sr_10( south_pole));

      BCL_MessageStd
      (
        "point on north pole and quantized with resolution: " + util::Format()( ptk_sr_10.GetAngularResolution()) + "\n" +
        util::Format()( north_pole) + "\n->\n" + util::Format()( north_pole_quantized)
      );
      BCL_MessageStd
      (
        "point on equator and quantized with resolution: " + util::Format()( ptk_sr_10.GetAngularResolution()) + "\n" +
        util::Format()( equator) + "\n->\n" + util::Format()( equator_quantized)
      );
      BCL_MessageStd
      (
        "point on south pole and quantized with resolution: " + util::Format()( ptk_sr_10.GetAngularResolution()) + "\n" +
        util::Format()( south_pole) + "\n->\n" + util::Format()( south_pole_quantized)
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // write to file
      WriteBCLObject( ptk_sr_10);
      // read
      coord::PointToKeySphericalRadius ptr_sr_read( 0);
      ReadBCLObject( ptr_sr_read);

      BCL_Example_Check
      (
        ptk_sr_10.GetAngularResolution() == ptr_sr_read.GetAngularResolution(),
        "written and read PointToKeySphericalRadius is different"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleCoordPointToKeySphericalRadius

  const ExampleClass::EnumType ExampleCoordPointToKeySphericalRadius::s_Instance
  (
    GetExamples().AddEnum( ExampleCoordPointToKeySphericalRadius())
  );

} // namespace bcl
