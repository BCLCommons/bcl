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
#include "util/bcl_util_si_ptr_vector.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_util_si_ptr_vector.cpp
  //!
  //! @author woetzen
  //! @date
  //! @remarks status empty
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleUtilSiPtrVector :
    public ExampleInterface
  {
  public:

    ExampleUtilSiPtrVector *Clone() const
    { return new ExampleUtilSiPtrVector( *this);}

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
      // these are two arrys of 5 Vector3D each, like in the bcl_coordinates example
      linal::Vector3D c1[ 5] =
      {
        linal::Vector3D( 40.338,4.848,8.769),
        linal::Vector3D( 42.032,3.315,5.760),
        linal::Vector3D( 41.619,-0.354,6.650),
        linal::Vector3D( 37.894,0.039,7.305),
        linal::Vector3D( 37.085,2.778,4.773)
      };

      linal::Vector3D c2[ 5] =
      {
        linal::Vector3D( 34.816,-3.039,18.802),
        linal::Vector3D( 38.034,-2.978,16.788),
        linal::Vector3D( 39.908,-2.243,20.009),
        linal::Vector3D( 38.121,-5.112,21.771),
        linal::Vector3D( 39.082,-7.417,18.900)
      };

      // construct util::SiPtrVector of Vector3D from Length of array and araay
      BCL_MessageStd( "initialize two Simplepointervector from two arrays of Vector3D coordinates");
      util::SiPtrVector< linal::Vector3D> coord_1( 5, c1);
      util::SiPtrVector< linal::Vector3D> coord_2( 5, c2);
      BCL_MessageStd( util::Format()( coord_1));
      BCL_MessageStd( util::Format()( coord_2));
      BCL_MessageStd( "initialize Simplepointervector from util::SiPtrVector of Vector3D coordinates");
      util::SiPtrVector< linal::Vector3D> spv1( coord_1);

      //initialize single Elements and util::SiPtr on Element
      linal::Vector3D v1( 34.816, -3.039, 18.802);
      util::SiPtr< linal::Vector3D> pv1( &v1);
      BCL_MessageStd( "Element and Simplepointer on Element are inizialized and spv1 get them added to the end");
      spv1.PushBack( &v1);
      spv1.PushBack( pv1);
      BCL_MessageStd( util::Format()( spv1));
      BCL_MessageStd( "Insert coord_1 in spv1 at Position 0");
      spv1.InsertElements( 0, coord_2);
      BCL_MessageStd( "Size of spv1: " + util::Format()( spv1.GetSize()));
      BCL_MessageStd( util::Format()( spv1));
      BCL_MessageStd( "Resize spv1 to LENGTH of 4 ");
      spv1.Resize( 4);
      BCL_MessageStd( util::Format()( spv1));
      BCL_MessageStd( "Delete 1st Element of spv1");
      spv1.RemoveElements( 0);
      BCL_MessageStd( util::Format()( spv1));
      BCL_MessageStd( "Delete 1st and 2nd Element of spv1");
      spv1.RemoveElements( 0, 2);
      BCL_MessageStd( util::Format()( spv1));

      BCL_MessageStd( "Initialize spv2 from SubVecttor of coord_1 from POS 2 of  LENGTH 2");
      util::SiPtrVector< linal::Vector3D> spv2( coord_1.SubSiPtrVector( 1, 2));
      BCL_MessageStd( util::Format()( spv2));

      BCL_MessageStd( "swap Position one and two");
      std::swap( spv2( 0), spv2( 1));
      BCL_MessageStd( util::Format()( spv2));

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleUtilSiPtrVector

  const ExampleClass::EnumType ExampleUtilSiPtrVector::s_Instance
  (
    GetExamples().AddEnum( ExampleUtilSiPtrVector())
  );

} // namespace bcl
