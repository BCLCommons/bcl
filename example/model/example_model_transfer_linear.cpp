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
#include "model/bcl_model_transfer_linear.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_transfer_linear.cpp
  //!
  //! @author mendenjl
  //! @date Nov 18, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelTransferLinear :
    public ExampleInterface
  {
  public:

    ExampleModelTransferLinear *Clone() const
    {
      return new ExampleModelTransferLinear( *this);
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
      const model::TransferLinear tg_default;

      // clone
      const util::ShPtr< model::TransferLinear> tg_clone( tg_default.Clone());

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck( model::TransferLinear().GetOutputRange(), math::Range< float>( 0.0, 1.0));

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      BCL_ExampleCheck( model::TransferLinear().F( 0.0), float( 0.0));
      BCL_ExampleCheck( model::TransferLinear().F( 3.0), float( 3.0));
      BCL_ExampleCheck( model::TransferLinear().F( -3.0), float( -3.0));
      BCL_ExampleCheck( model::TransferLinear().dF( 0.0, 1.0), float( 1));
      BCL_ExampleCheck( model::TransferLinear().dF( 3.0, 0.011109), float( 1));

      const float features_fill[ 3] = { -3.0, 0.0, 3.0};
      const float results_fill[ 3] = { 0.011109, 1.0, 0.011109};

      const linal::Vector< float> features( 3, features_fill);
      const linal::Vector< float> results( 3, results_fill);
      linal::Vector< float> multiplied_by_df( 3, float( 2.0));
      model::TransferLinear().MultiplyBydF( multiplied_by_df, features, results);
      BCL_ExampleIndirectCheck( multiplied_by_df( 0), float( 2.0), "MultiplyBydF");
      BCL_ExampleIndirectCheck( multiplied_by_df( 1), float( 2.0), "MultiplyBydF");
      BCL_ExampleIndirectCheck( multiplied_by_df( 2), float( 2.0), "MultiplyBydF");

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelTransferLinear

  const ExampleClass::EnumType ExampleModelTransferLinear::s_Instance
  (
    GetExamples().AddEnum( ExampleModelTransferLinear())
  );

} // namespace bcl
