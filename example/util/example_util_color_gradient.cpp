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
#include "util/bcl_util_color_gradient.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_operations.h"
#include "util/bcl_util_colors.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_util_color_gradient.cpp
  //!
  //! @author alexanns, karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleUtilColorGradient :
    public ExampleInterface
  {
  public:

    ExampleUtilColorGradient *Clone() const
    {
      return new ExampleUtilColorGradient( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    virtual const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

      // create upper and lower bounds
      const math::Range< double> bounds( 1.0, 9.0);

      // initialize gradient colors to use
      const storage::Vector< util::Color> multiple_gradient_points
      (
        storage::Vector< util::Color>::Create
        (
          util::GetColors().e_Red,
          util::GetColors().e_Yellow,
          util::GetColors().e_White,
          util::GetColors().e_Yellow
        )
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // call default constructor
      util::ColorGradient def_constr;

      // call constructor taking parameters for each member variable
      // create OutputPymolColorGradient "param_constr" and initialize with the parameters
      const util::ColorGradient param_constr( bounds, multiple_gradient_points);
      BCL_ExampleCheck( param_constr( bounds.GetMin()), *multiple_gradient_points.FirstElement());
      BCL_ExampleCheck( param_constr( bounds.GetMax()), *multiple_gradient_points.LastElement());

      // call copy constructor
      // create OutputPymolColorGradient "copy_constr" and initialize with "param_constr"
      const util::ColorGradient copy_constr( param_constr);
      BCL_ExampleIndirectCheck( copy_constr( bounds.GetMin()), *multiple_gradient_points.FirstElement(), "copy constructor");
      BCL_ExampleIndirectCheck( copy_constr( bounds.GetMax()), *multiple_gradient_points.LastElement(), "copy constructor");

      // virtual copy constructor
      // create pointer to OutputPymolColorGradient "clone_constr"
      const util::ShPtr< util::ColorGradient> clone_constr( copy_constr.Clone());
      BCL_ExampleIndirectCheck( clone_constr->operator()( bounds.GetMin()), *multiple_gradient_points.FirstElement(), "clone");
      BCL_ExampleIndirectCheck( clone_constr->operator()( bounds.GetMax()), *multiple_gradient_points.LastElement(), "clone");

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck( clone_constr->GetClassIdentifier(), GetStaticClassName< util::ColorGradient>());

    ///////////////
    // operators //
    ///////////////

      // test operator() with multiple gradients
      // test below range
      {
        const linal::Vector3D &expected_color_points( *multiple_gradient_points.FirstElement());
        const linal::Vector3D calculated_color_points( param_constr( 0.0));
        BCL_ExampleCheck( calculated_color_points, expected_color_points);
      }
      // test above range
      {
        const linal::Vector3D &expected_color_points( *multiple_gradient_points.LastElement());
        const linal::Vector3D calculated_color_points( param_constr( 10.0));
        BCL_ExampleCheck( calculated_color_points, expected_color_points);
      }
      // test in range
      {
        const linal::Vector3D expected_color_points( 1.0, 1.0, 0.5);
        const linal::Vector3D calculated_color_points( param_constr( 5.0));
        BCL_ExampleCheck
        (
          math::EqualWithinTolerance( expected_color_points( 0), calculated_color_points( 0)) &&
          math::EqualWithinTolerance( expected_color_points( 1), calculated_color_points( 1)) &&
          math::EqualWithinTolerance( expected_color_points( 2), calculated_color_points( 2)),
          true
        );
      }
      // test at max
      {
        const linal::Vector3D &expected_color_points( *multiple_gradient_points.LastElement());
        const linal::Vector3D calculated_color_points( param_constr( 9.0));
        BCL_ExampleCheck( calculated_color_points, expected_color_points);
      }
      // test at min
      {
        const linal::Vector3D &expected_color_points( *multiple_gradient_points.FirstElement());
        const linal::Vector3D calculated_color_points( param_constr( 1.0));
        BCL_ExampleCheck( calculated_color_points, expected_color_points);
      }
      // test at exact color point
      {
        const linal::Vector3D expected_color_points( 1.0, 1.0, 1.0);
        const linal::Vector3D calculated_color_points( param_constr( 19.0 / 3.0));
        BCL_ExampleCheck( calculated_color_points, expected_color_points);
      }

    //////////////////////
    // input and output //
    //////////////////////

      // test Write
      BCL_MessageStd( "Testing the write function");
      WriteBCLObject( param_constr);
      // read the file back
      BCL_MessageStd( "Testing the read function");
      util::ColorGradient read_output;
      ReadBCLObject( read_output);
      BCL_ExampleIndirectCheck( read_output( bounds.GetMin()), *multiple_gradient_points.FirstElement(), "reading");
      BCL_ExampleIndirectCheck( read_output( bounds.GetMax()), *multiple_gradient_points.LastElement(), "reading");

    //////////////////////
    // helper functions //
    //////////////////////

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleUtilColorGradient

  const ExampleClass::EnumType ExampleUtilColorGradient::s_Instance
  (
    GetExamples().AddEnum( ExampleUtilColorGradient())
  );

} // namespace bcl
