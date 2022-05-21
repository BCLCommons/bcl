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
#include "model/bcl_model_transfer_sigmoid.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_transfer_sigmoid.cpp
  //!
  //! @author mueller, loweew
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelTransferSigmoid :
    public ExampleInterface
  {
  public:

    ExampleModelTransferSigmoid *Clone() const
    {
      return new ExampleModelTransferSigmoid( *this);
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

      // default constructor
      model::TransferSigmoid ts_default;

      // clone
      const util::ShPtr< model::TransferSigmoid> ts_clone( ts_default.Clone());

    /////////////////
    // data access //
    /////////////////

      // check if the class identifiers are working

      BCL_ExampleCheck( model::TransferSigmoid().GetOutputRange().GetString( util::Format()), "[0,1]");

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      BCL_ExampleCheck( ts_default.F( 0.0), float( 0.5));
      BCL_ExampleCheck( ts_default.F( 3.0), float( 1.0 / ( 1.0 + std::exp( -3.0))));
      BCL_ExampleCheck( ts_default.F( -3.0), float( 1.0 / ( 1.0 + std::exp( 3.0))));
      BCL_ExampleCheck( ts_default.dF( 0.0, 0.5), float( 0.25));
      BCL_ExampleCheck( ts_default.dF( 3.0, 0.95257), float( 0.95257) * float( ( 1 - float( 0.95257))));
      BCL_ExampleCheck( ts_default.dF( -3.0, 0.04742), float( 0.04742) * float( ( 1 - float( 0.04742))));

      const float features_fill[ 3] = { -3.0, 0.0, 3.0};
      const float results_fill[ 3] = { 0.04742, 0.5, 0.95257};

      const linal::Vector< float> features( 3, features_fill);
      const linal::Vector< float> results( 3, results_fill);
      linal::Vector< float> multiplied_by_df( 3, double( 2.0));
      model::TransferSigmoid().MultiplyBydF( multiplied_by_df, features, results);
      BCL_ExampleIndirectCheck
      (
        multiplied_by_df( 0),
        float( 2.0) * float( 0.04742) * float( ( 1.0 - float( 0.04742))),
        "MultiplyBydF"
      );
      BCL_ExampleIndirectCheck( multiplied_by_df( 1), float( 0.5), "MultiplyBydF");
      BCL_ExampleIndirectCheck
      (
        multiplied_by_df( 2),
        float( 2.0) * float( 0.95257) * float( ( 1.0 - float( 0.95257))),
        "MultiplyBydF"
      );

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelTransferSigmoid

  const ExampleClass::EnumType ExampleModelTransferSigmoid::s_Instance
  (
    GetExamples().AddEnum( ExampleModelTransferSigmoid())
  );

} // namespace bcl
