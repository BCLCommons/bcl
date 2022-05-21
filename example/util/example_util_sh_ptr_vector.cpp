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
#include "util/bcl_util_sh_ptr_vector.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //! @example example_util_sh_ptr_vector.cpp
  //!
  //! @author mueller
  //! @date
  //! @remarks status empty
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleUtilShPtrVector :
    public ExampleInterface
  {
  public:

    ExampleUtilShPtrVector *Clone() const
    { return new ExampleUtilShPtrVector( *this);}

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
      // these are two arrays of 3 Vector3D each, like in the bcl_coordinates example
      linal::Vector3D c1[ 3] =
      {
        linal::Vector3D( 40.338,4.848,8.769),
        linal::Vector3D( 42.032,3.315,5.760),
        linal::Vector3D( 41.619,-0.354,6.650),
      };

      linal::Vector3D c2[ 3] =
      {
        linal::Vector3D( 34.816,-3.039,18.802),
        linal::Vector3D( 38.034,-2.978,16.788),
        linal::Vector3D( 39.908,-2.243,20.009),
      };

      BCL_MessageStd( "util::ShPtrVector of 3 coordinates initialized");
      util::ShPtrVector< linal::Vector3D> shpv1( 3, c1);
      BCL_MessageStd( util::Format()( shpv1));

      BCL_MessageStd( util::Format()( shpv1( 0)));
      util::ShPtr< linal::Vector3D> firstelementsoft( shpv1( 0));
      util::ShPtr< linal::Vector3D> firstelementhard( shpv1( 0).HardCopy());
      //  shpv1.PushBack( firstelementsoft.HardCopy());
      BCL_MessageStd( "first element " + util::Format()( shpv1( 0)));
      BCL_MessageStd( "first element in SharedPointVector is set 1.0, 1.0, 1.0");
      *( shpv1( 0)) = linal::Vector3D( 1.0, 1.0, 1.0);
      BCL_MessageStd( util::Format()( shpv1( 0)));
      BCL_MessageStd( "soft copy of first element " + util::Format()( firstelementsoft));
      BCL_MessageStd( "hard copy of first element " + util::Format()( firstelementhard));

      BCL_MessageStd( "second ShPtrVector of coordinates is initialized");
      util::ShPtrVector< linal::Vector3D> shpv2( 3, c2);
      BCL_MessageStd( "third ShPtrVector as soft copy of the second ShPtrVector");
      util::ShPtrVector< linal::Vector3D> shpv3( 3);
      shpv3 = shpv2;
      BCL_MessageStd( util::Format()( shpv2));
      BCL_MessageStd( util::Format()( shpv3));
      BCL_MessageStd( "element 1 of ShPtrVector 3 is set to 3.3");
      *( shpv3( 0)) = linal::Vector3D( 3.3);
      BCL_MessageStd( util::Format()( shpv2));
      BCL_MessageStd( util::Format()( shpv3));

      BCL_MessageStd( "fourth ShPtrVector as hard copy of the second ShPtrVector");
      util::ShPtrVector< linal::Vector3D> shpv4( 3);
      shpv4 = shpv2.HardCopy();
      BCL_MessageStd( util::Format()( shpv2));
      BCL_MessageStd( util::Format()( shpv4));
      BCL_MessageStd( "element 1 of ShPtrVector 3 is set to 1.5");
      *( shpv2( 0)) = linal::Vector3D( 1.5);
      BCL_MessageStd( util::Format()( shpv2));
      BCL_MessageStd( util::Format()( shpv4));
      BCL_MessageStd( "add additional soft copy of element to the end of ShPtrVector 4");
      shpv4.PushBack( firstelementhard);
      BCL_MessageStd( "the origin element is set to 10.0");
      *firstelementhard = linal::Vector3D( 10.0);
      BCL_MessageStd( util::Format()( shpv4));

      BCL_MessageStd( "Add Copy of Object as ShPtr in shpv4");
      linal::Vector3D v1( 39.908, -2.243, 20.009);
      shpv4.PushBack( util::ShPtr< linal::Vector3D>( v1.Clone()));
      BCL_MessageStd( util::Format()( shpv4));

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleUtilShPtrVector

  const ExampleClass::EnumType ExampleUtilShPtrVector::s_Instance
  (
    GetExamples().AddEnum( ExampleUtilShPtrVector())
  );

} // namespace bcl
