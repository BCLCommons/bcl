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
#include "opencl/bcl_opencl_transfer_function_gaussian.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opencl_transfer_function_gaussian.cpp
  //!
  //! @author loweew
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOpenclTransferFunctionGaussian :
    public ExampleInterface
  {
  public:

    ExampleOpenclTransferFunctionGaussian *Clone() const
    {
      return new ExampleOpenclTransferFunctionGaussian( *this);
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
      // do not try to run opencl commands if no queue was found
      if( !opencl::GetTools().HasCommandQueues())
      {
        return 1;
      }

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      const opencl::TransferFunctionGaussian< float> tg_default;

      // clone
      const util::ShPtr< opencl::TransferFunctionGaussian< float> > tg_clone( tg_default.Clone());

    /////////////////
    // data access //
    /////////////////

      // check class identifiers
      BCL_ExampleCheck
      (
        tg_default.GetClassIdentifier(),
        "bcl::opencl::TransferFunctionGaussian<float>"
      );

      // check that clone didn't return null
      BCL_ExampleAssert( tg_clone.IsDefined(), true);

      // check that clone returned an equivalent object
      BCL_ExampleCheck( tg_clone->GetClassIdentifier(), tg_default.GetClassIdentifier());

      BCL_ExampleCheck
      (
        tg_default.GetOutputRange().GetMax(),
        1.0
      );

      BCL_ExampleCheck
      (
        tg_default.GetOutputRange().GetMin(),
        0.0
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      float v1[ 3] = { -3.0, 0.0, 3.0};
      float v2[ 3] = { 0.011109, 1.0, 0.011109};
      float m1[ 9] = { -3.0, 0.0, 3.0, -3.0, 0.0, 3.0, -3.0, 0.0, 3.0};
      float m2[ 9] = { 0.011109, 1.0, 0.011109, 0.011109, 1.0, 0.011109, 0.011109, 1.0, 0.011109};
      const linal::Vector< float> features( 3, v1);
      const linal::Vector< float> results( 3, v2);
      const linal::Matrix< float> features_matrix( 3, 3, m1);
      const linal::Matrix< float> results_matrix( 3, 3, m2);

      BCL_ExampleCheck
      (
        tg_default.F( v1[ 0]),
        std::exp( -math::Sqr( v1[ 0]) / float( 2.0))
      );

      BCL_ExampleCheck
      (
        tg_default.F( v1[ 1]),
        std::exp( -math::Sqr( v1[ 1]) / float( 2.0))
      );

      BCL_ExampleCheck
      (
        tg_default.dF( v1[ 0], v2[ 0]),
        v1[ 0] * v2[ 0]
      );

      BCL_ExampleCheck
      (
        tg_default.dF( v1[ 1], v2[ 1]),
        v1[ 1] * v2[ 1]
      );

      BCL_ExampleCheck
      (
        tg_default.dF( v1[ 2], v2[ 2]),
        v1[ 2] * v2[ 2]
      );

      BCL_Example_Check
      (
        math::EqualWithinTolerance( v1[ 0] * v2[ 0], tg_clone->dF( features, results)( 0)) &&
        math::EqualWithinTolerance( v1[ 1] * v2[ 1], tg_clone->dF( features, results)( 1)) &&
        math::EqualWithinTolerance( v1[ 2] * v2[ 2], tg_clone->dF( features, results)( 2)),
        "dF( v1, v2) should return v1 * v2 but returns "
        + util::Format()( tg_clone->dF( features, results))
      );

      BCL_Example_Check
      (
        math::EqualWithinTolerance( v1[ 0] * v2[ 0], tg_clone->dF( features_matrix, results_matrix)( 0,0)) &&
        math::EqualWithinTolerance( v1[ 1] * v2[ 1], tg_clone->dF( features_matrix, results_matrix)( 1,1)) &&
        math::EqualWithinTolerance( v1[ 2] * v2[ 2], tg_clone->dF( features_matrix, results_matrix)( 2,2)),
        "dF( m1, m2) should return m1 * m2 (element-wise) but returns "
        + util::Format()( tg_clone->dF( features_matrix, results_matrix))
      );

      BCL_Example_Check
      (
        math::EqualWithinTolerance( exp( -( -3.0 * -3.0) / 2.0), tg_clone->F( features)( 0), 0.00001, 0.0001) &&
        math::EqualWithinTolerance( exp( -( 0.0 * 0.0) / 2), tg_clone->F( features)( 1), 0.00001, 0.0001) &&
        math::EqualWithinTolerance( exp( -( 3.0 * 3.0) / 2), tg_clone->F( features)( 2), 0.00001, 0.0001) ,
        "F( -3, 0, 3)should return ( 0.011109, 1, 0.011109) but returns "
        + util::Format()( tg_clone->F( features))
      );

      BCL_Example_Check
      (
        math::EqualWithinTolerance( exp( -( -3.0 * -3.0) / 2.0), tg_clone->F( features_matrix)( 0, 0), 0.00001, 0.0001) &&
        math::EqualWithinTolerance( exp( -( 0.0 * 0.0) / 2), tg_clone->F( features_matrix)( 1, 1), 0.00001, 0.0001) &&
        math::EqualWithinTolerance( exp( -( 3.0 * 3.0) / 2), tg_clone->F( features_matrix)( 2, 2), 0.00001, 0.0001) ,
        "F( -3, 0, 3, -3, 0, 3, -3, 0, 3,) should return (  0.011109, 1, 0.011109, 0.011109, 1, 0.011109, 0.011109, 1, 0.011109) but returns "
        + util::Format()( tg_clone->F( features_matrix))
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleOpenclTransferFunctionGaussian

  const ExampleClass::EnumType ExampleOpenclTransferFunctionGaussian::s_Instance
  (
    GetExamples().AddEnum( ExampleOpenclTransferFunctionGaussian())
  );

} // namespace bcl
