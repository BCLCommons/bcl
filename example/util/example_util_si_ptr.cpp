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
#include "util/bcl_util_si_ptr.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //! @example example_util_si_ptr.cpp
  //!
  //! @author riddeljs
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleUtilSiPtr :
    public ExampleInterface
  {
  public:

    ExampleUtilSiPtr *Clone() const
    { return new ExampleUtilSiPtr( *this);}

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

      // default constructor
      util::SiPtr< linal::Vector3D> sp_vec1;
      BCL_Example_Check
      (
        !sp_vec1.IsDefined(), "sp_vec1 is defined but should have pointer NULL"
      );
      // construct from existing object
      linal::Vector3D vector_of_ones( 1.0, 1.0, 1.0);
      util::SiPtr< linal::Vector3D> sp_vec2( vector_of_ones);
      BCL_Example_Check
      (
        sp_vec2.IsDefined(), "sp_vec1 is undefined but should have pointer to existing linal::Vector3D"
      );

      // construct from existing pointer
      linal::Vector3D *vector_of_twos( new linal::Vector3D( 1.0, 1.0, 1.0));
      util::SiPtr< linal::Vector3D> sp_vec3( vector_of_twos);
      BCL_Example_Check
      (
        sp_vec3.IsDefined(), "sp_vec3 is undefined but should have pointer defined based on std pointer"
      );

      // copy constructor
      util::SiPtr< linal::Vector3D> sp_vec4( sp_vec3);
      BCL_Example_Check
      (
        sp_vec4.IsDefined(), "sp_vec4 is undefined but should have pointer copied from sp_vec3"
      );

      // copy constructor from t_DataType to t_OtherDataType
      util::SiPtr< util::ObjectInterface> sp_vec5( sp_vec4);
      BCL_Example_Check
      (
        sp_vec4.IsDefined(), "sp_vec5 is undefined but should have pointer copied from sp_vec4"
      );

      // clone
      util::ShPtr< linal::Vector3D> sp_vec6( sp_vec2->Clone());
      BCL_Example_Check
      (
        sp_vec6.IsDefined(), "sp_vec6 is undefined but should have pointer to new linal::Vector3D"
      );

    ///////////////
    // operators //
    ///////////////

      // assignment from a pointer, we'll use vector_of_twos from above
      util::SiPtr< linal::Vector3D> sp_vec7;
      sp_vec7 = vector_of_twos;
      BCL_Example_Check
      (
        sp_vec6.IsDefined(), "sp_vec7 is undefined but should have pointer from assignment to vector_of_twos"
      );

      // assignment from a pointer to a different data type pointer
      util::SiPtr< linal::Vector3D> sp_vec8;
      util::ObjectInterface *vector_of_threes( new linal::Vector3D( 3.0, 3.0, 3.0));
      sp_vec8 = vector_of_threes;
      BCL_Example_Check
      (
        sp_vec8.IsDefined(),
        "sp_vec8 is undefined but should have pointer from assignment to object interface vector_of_threes"
      );

      // assignment from a SiPtr, we'll just use sp_vec7 from above
      util::SiPtr< linal::Vector3D> sp_vec9;
      sp_vec9 = sp_vec7;
      BCL_Example_Check
      (
        sp_vec9.IsDefined(),
        "sp_vec9 is undefined but should have pointer from assignment to SiPtr sp_vec7"
      );

      // assignment from a SiPtr of a different type, so use sp_vec5 from earlier
      util::SiPtr< linal::Vector3D> sp_vec10;
      sp_vec10 = sp_vec5;
      BCL_Example_Check
      (
        sp_vec10.IsDefined(),
        "sp_vec10 is undefined but should have pointer from assignment to SiPtr< ObjectInterface> sp_vec5"
      );

    /////////////////
    // data access //
    /////////////////

      // get the actual pointer
      BCL_MessageStd
      (
        "pointer pointing to Vector3D in sp_vec7 " + util::Format()( sp_vec7.GetPointer())
      );

      // test ->
      BCL_Example_Check
      (
        sp_vec2->operator()( 0) == 1.0, "operator-> failed"
      );

      // test *
      double x_coord( ( *sp_vec2).X());
      double y_coord( ( *sp_vec2).Y());
      double z_coord( ( *sp_vec2).Z());
      BCL_Example_Check
      (
        x_coord == 1.0 && y_coord == 1.0 && z_coord == 1.0, "operator* failed"
      );

      // test reset
      sp_vec2.Reset();
      BCL_Example_Check
      (
        !sp_vec2.IsDefined(),
        "sp_vec2 is defined but it should have been reset"
      );

      delete vector_of_twos;

      // ensure that we can create a simple pointer to a double
      double favorite_number( math::g_Pi);
      util::SiPtr< double> si_ptr_fav_number( favorite_number);
      BCL_ExampleCheck( *si_ptr_fav_number, favorite_number);

    //////////////////////
    // input and output //
    //////////////////////

      // read and write aren't implemented for SiPtrs so leave this alone.

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleUtilSiPtr

  const ExampleClass::EnumType ExampleUtilSiPtr::s_Instance
  (
    GetExamples().AddEnum( ExampleUtilSiPtr())
  );

} // namespace bcl
