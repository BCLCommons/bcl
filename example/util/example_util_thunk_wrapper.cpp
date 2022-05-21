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
#include "util/bcl_util_thunk_wrapper.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_util_thunk_wrapper.cpp
  //!
  //! @author mendenjl
  //! @date July 09, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleUtilThunkWrapper :
    public ExampleInterface
  {
  public:

    ExampleUtilThunkWrapper *Clone() const
    {
      return new ExampleUtilThunkWrapper( *this);
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

      //tests for non-const member functions
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      util::ThunkWrapper< linal::Vector3D, double &> call_x( &linal::Vector3D::X);
      util::ThunkWrapper< linal::Vector3D, double &> call_y( &linal::Vector3D::Y);
      util::ThunkWrapper< linal::Vector3D, double &> call_z( &linal::Vector3D::Z);

      linal::Vector3D test_vect( 1.0, 2.0, 3.0);

    ///////////////
    // operators //
    ///////////////
      // test that each of the call_x/y/z functions works properly
      BCL_Example_Check
      (
        call_x( test_vect) == test_vect.X()
        && call_y( test_vect) == test_vect.Y()
        && call_z( test_vect) == test_vect.Z(),
        "util::ThunkWrapper couldn't access X/Y/Z of a linal::Vector3D properly"
      );

      // set test_vect.X() through call_x:
      call_x( test_vect) = 10.0;

      // test that assignment through call X worked properly
      BCL_Example_Check
      (
        call_x( test_vect) == test_vect.X() && test_vect.X() == 10.0,
        "util::ThunkWrapper didn't change X of a linal::Vector3D"
      );

      util::ThunkWrapper< const linal::Vector3D, const double &> call_const_x( &linal::Vector3D::X);
      util::ThunkWrapper< const linal::Vector3D, const double &> call_const_y( &linal::Vector3D::Y);
      util::ThunkWrapper< const linal::Vector3D, const double &> call_const_z( &linal::Vector3D::Z);

    ///////////////
    // operators //
    ///////////////

      // test that each of the call_x/y/z functions works properly
      BCL_Example_Check
      (
        call_const_x( test_vect) == test_vect.X()
        && call_const_y( test_vect) == test_vect.Y()
        && call_const_z( test_vect) == test_vect.Z(),
        "util::ThunkWrapper (Const) couldn't access X/Y/Z of a linal::Vector3D properly"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleUtilThunkWrapper

  const ExampleClass::EnumType ExampleUtilThunkWrapper::s_Instance
  (
    GetExamples().AddEnum( ExampleUtilThunkWrapper())
  );

} // namespace bcl
