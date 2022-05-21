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
#include "util/bcl_util_sh_ptr.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //! @example example_util_sh_ptr.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleUtilShPtr :
    public ExampleInterface
  {
  public:

    ExampleUtilShPtr *Clone() const
    {
      return new ExampleUtilShPtr( *this);
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

      // empty constructor
      util::ShPtr< util::ObjectInterface> sp_empty;
      BCL_ExampleIndirectCheck( sp_empty.IsDefined(), false, "empty constructor");

      // construct from new object
      util::ShPtr< linal::Vector3D> sp_vec1( new linal::Vector3D( 1.0, 1.0, 1.0));
      BCL_Example_Check
      (
        sp_vec1.IsDefined(),
        "sp_vec1 is undefined but should have pointer to new linal::Vector3D"
      );
      BCL_Example_Check
      (
        sp_vec1.GetSharedCommunitySize() == 1,
        "sp_vec1 should be the only pointer pointing to that object!"
      );

      // construct from new object
      util::ShPtr< util::ObjectInterface> sp_vec2( new linal::Vector3D( 2.0, 2.0, 2.0));
      BCL_Example_Check
      (
        sp_vec2.IsDefined(),
        "sp_vec2 is undefined but should have pointer to new linal::Vector3D"
      );
      BCL_Example_Check
      (
        sp_vec2.GetSharedCommunitySize() == 1,
        "sp_vec2 should be the only pointer pointing to that object!"
      );

      // construct by cloning a second object
      util::ShPtr< linal::VectorInterface< double> > sp_vec3( sp_vec1->Clone());
      BCL_Example_Check
      (
        sp_vec3.IsDefined(),
        "sp_vec3 is undefined but should have pointer to new linal::Vector3D"
      );
      BCL_Example_Check
      (
        sp_vec3.GetSharedCommunitySize() == 1,
        "sp_vec3 should be the only pointer pointing to that object!"
      );

      // default constructor
      util::ShPtr< linal::Vector3D> sp_vec4;
      BCL_Example_Check
      (
        !sp_vec4.IsDefined(),
        "sp_vec4 is defined but should have pointer NULL"
      );
      BCL_ExampleCheck( sp_vec4.GetSharedCommunitySize(), 0);

      // clone
      util::ShPtr< linal::Vector3D> *ptr( sp_vec1.Clone());

      BCL_Example_Check
      (
        *ptr == sp_vec1 && sp_vec1.GetSharedCommunitySize() == 2,
        "clone did not copy sp_vec1 or ring was not increased to size 2"
      );

    /////////////////
    // data access //
    /////////////////

      BCL_Example_Check
      (
        GetStaticClassName< util::ShPtr< linal::Vector3D> >() == ptr->GetClassIdentifier(),
        "incorrect class identifier"
      );

      // get size of ring - number of pointer pointing to the same object
      BCL_MessageStd
      (
        "number of sh_ptr pointing to Vector3D in sp_vec1 ring "
        + util::Format()( sp_vec1.GetSharedCommunitySize())
      );

      // get the actual pointer
      BCL_MessageStd
      (
        "pointer pointing to Vector3D in sp_vec1 " + util::Format()( sp_vec1.GetPointer())
      );

    ///////////////
    // operators //
    ///////////////

      // assignment operator - inserts sp_vec4 into ring of sp_vec3 and lets spvec4 delete its object, if it is
      // the only one left in the ring
      sp_vec4 = sp_vec3;
      BCL_Example_Check
      (
        sp_vec3.GetSharedCommunitySize() == 2 && sp_vec4.GetSharedCommunitySize() == 2,
        "sp_vec3 and sp_vec4 should have ringsize of 2 but are: " +
        util::Format()( sp_vec3.GetSharedCommunitySize()) + " and " + util::Format()( sp_vec4.GetSharedCommunitySize())
      );

      BCL_Example_Check
      (
        sp_vec3.IsDefined() && sp_vec4.IsDefined(),
        "sp_vec3 and sp_vec4 should have ringsize of 2 but are: " +
        util::Format()( sp_vec3.GetSharedCommunitySize()) + " and " + util::Format()( sp_vec4.GetSharedCommunitySize())
      );

      // copy constructor from t_DataType to t_OtherDataType
      // ShPtr of different types point to the same object - dynamic cast assures that pointer are of appropriate type
      util::ShPtr< util::ObjectInterface> sp_vec5( sp_vec4);
      BCL_Example_Check
      (
        sp_vec3.GetSharedCommunitySize() == 3
        && sp_vec4.GetSharedCommunitySize() == 3
        && sp_vec5.GetSharedCommunitySize() == 3,
        "sp_vec3 4 and 5 should be in the ring of size 3"
      );

      // comparison operators
      BCL_Example_Check
      (
        size_t( sp_vec3.GetPointer()) == size_t( sp_vec4.GetPointer()) &&
        size_t( sp_vec3.GetPointer()) == size_t( sp_vec5.GetPointer()),
        "sp_vec3 4 and 5 should point to the same object"
      );

      // -> operator to access object
      // set x of vector3d in pointer
      sp_vec3->operator()( 0) = 2.0;
      BCL_Example_Check
      (
        sp_vec3->operator()( 0) == 2.0 && sp_vec4->X() == 2.0,
        "X of sp_vec 3 4 and 5 should be 2.0"
      );

      BCL_MessageStd( "this is the output when writing sp_vec5: " + util::Format()( sp_vec5));

    //////////////////////
    // input and output //
    //////////////////////

      WriteBCLObject( sp_empty);
      WriteBCLObject( sp_vec3);
      WriteBCLObject( sp_vec4);
      util::ShPtr< linal::VectorInterface< double> > sp_vec6;
      util::ShPtr< linal::Vector3D> sp_vec7;
      util::ShPtr< util::ObjectInterface> sp_empty_read;
      ReadBCLObject( sp_vec6);
      ReadBCLObject( sp_vec7);
      ReadBCLObject( sp_empty_read);

      BCL_MessageStd
      (
        "number of sh_ptr pointing to Vector3D in sp_vec6 ring "
        + util::Format()( sp_vec6.GetSharedCommunitySize())
      );

      BCL_ExampleIndirectCheck( sp_empty_read.IsDefined(), sp_empty.IsDefined(), "reading empty shared pointer");

      // cleanup
      delete ptr;

      // make a sh ptr directly to a double
      util::ShPtr< double> sh_ptr_double( new double( 5.0));
      // perform hard copy
      util::ShPtr< double> cp_sh_ptr_double( sh_ptr_double.HardCopy());
      BCL_ExampleIndirectCheck( *sh_ptr_double, *cp_sh_ptr_double, "HardCopy on pointer to non-bcl type");

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleUtilShPtr

  const ExampleClass::EnumType ExampleUtilShPtr::s_Instance
  (
    GetExamples().AddEnum( ExampleUtilShPtr())
  );

} // namespace bcl
