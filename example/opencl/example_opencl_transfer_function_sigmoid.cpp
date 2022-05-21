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
#include "opencl/bcl_opencl_transfer_function_sigmoid.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opencl_transfer_function_sigmoid.cpp
  //!
  //! @author loweew
  //! @date Apr 04, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOpenclTransferFunctionSigmoid :
    public ExampleInterface
  {
  public:

    ExampleOpenclTransferFunctionSigmoid *Clone() const
    {
      return new ExampleOpenclTransferFunctionSigmoid( *this);
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

      opencl::CommandQueue queue( opencl::GetTools().GetFirstCommandQueue());
      opencl::Matrix< float> data( 15, 9, queue, 0, 0, 1);
      opencl::Matrix< float> large_data( 5432, 19, queue, 0, 0, 0.5);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      const opencl::TransferFunctionSigmoid< float> ts_default( queue);

      // clone
      const util::ShPtr< opencl::TransferFunctionSigmoid< float> > ts_clone( ts_default.Clone());

    /////////////////
    // data access //
    /////////////////

      // check class identifiers
      BCL_ExampleCheck( ts_default.GetClassIdentifier(), GetStaticClassName< opencl::TransferFunctionSigmoid< float> >());

      // check that clone didn't return null
      BCL_ExampleAssert( ts_clone.IsDefined(), true);

      // check that clone returned an equivalent object
      BCL_ExampleCheck( ts_clone->GetClassIdentifier(), opencl::TransferFunctionSigmoid< float>().GetClassIdentifier());

      BCL_ExampleCheck
      (
        opencl::TransferFunctionSigmoid< float>().GetOutputRange().GetMax(),
        1.0
      );

      BCL_ExampleCheck
      (
        opencl::TransferFunctionSigmoid< float>().GetOutputRange().GetMin(),
        0.0
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      opencl::Matrix< float> ts_output( 15, 9, queue);
      opencl::Matrix< float> dsig_output( 15, 9, queue);

      ts_default.F( data, ts_output);
      queue.finish();
      ts_default.dF( ts_output, dsig_output);
      linal::Matrix< float> sig_results( ts_output.GetHostMatrix());
      linal::Matrix< float> dsig_results( dsig_output.GetHostMatrix());

      opencl::Matrix< float> large_ts_output( 5432, 19, queue, 0, 0, 0.5);
      opencl::Matrix< float> large_dsig_output( 5432, 19, queue, 0, 0, 0.5);

      ts_default.F( large_data, large_ts_output);
      queue.finish();
      ts_default.dF( large_ts_output, large_dsig_output);
      linal::Matrix< float> large_sig_results( large_ts_output.GetHostMatrix());
      linal::Matrix< float> large_dsig_results( large_dsig_output.GetHostMatrix());

      BCL_MessageStd( "data sigmoid: " + util::Format()( sig_results));
      BCL_MessageStd( "data sigmoid derivative: " + util::Format()( dsig_results));

      BCL_ExampleCheckWithinTolerance( double( sig_results( 0, 3)), double( 1 / ( 1 + exp( -1.0))), 0.001);
      BCL_ExampleCheckWithinTolerance( double( dsig_results( 0, 3)), double( sig_results( 0, 3) * ( 1 - sig_results( 0, 3))), 0.001);

      BCL_ExampleCheckWithinTolerance( double( large_sig_results(     0,  3)), double( 1 / ( 1 + exp( -0.5f))), 0.001);
      BCL_ExampleCheckWithinTolerance( double( large_sig_results(   215,  3)), double( 1 / ( 1 + exp( -0.5f))), 0.001);
      BCL_ExampleCheckWithinTolerance( double( large_sig_results(   656, 12)), double( 1 / ( 1 + exp( -0.5f))), 0.001);
      BCL_ExampleCheckWithinTolerance( double( large_sig_results(  1221, 15)), double( 1 / ( 1 + exp( -0.5f))), 0.001);
      BCL_ExampleCheckWithinTolerance( double( large_sig_results(  4512, 18)), double( 1 / ( 1 + exp( -0.5f))), 0.001);
      BCL_ExampleCheckWithinTolerance( double( large_dsig_results(    0,  3)), double( large_sig_results(    0,  3) * ( 1 - large_sig_results(    0,  3))), 0.001);
      BCL_ExampleCheckWithinTolerance( double( large_dsig_results(  215,  3)), double( large_sig_results(  215,  3) * ( 1 - large_sig_results(  215,  3))), 0.001);
      BCL_ExampleCheckWithinTolerance( double( large_dsig_results(  656, 12)), double( large_sig_results(  656, 12) * ( 1 - large_sig_results(  656, 12))), 0.001);
      BCL_ExampleCheckWithinTolerance( double( large_dsig_results( 1221, 15)), double( large_sig_results( 1221, 15) * ( 1 - large_sig_results( 1221, 15))), 0.001);
      BCL_ExampleCheckWithinTolerance( double( large_dsig_results( 4512, 18)), double( large_sig_results( 4512, 18) * ( 1 - large_sig_results( 4512, 18))), 0.001);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleOpenclTransferFunctionSigmoid

  const ExampleClass::EnumType ExampleOpenclTransferFunctionSigmoid::s_Instance
  (
    GetExamples().AddEnum( ExampleOpenclTransferFunctionSigmoid())
  );

} // namespace bcl
