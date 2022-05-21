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
#include "model/bcl_model_collect_features_top.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_collect_features_top.cpp
  //!
  //! @author mendenjl
  //! @date Apr 03, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelCollectFeaturesTop :
    public ExampleInterface
  {
  public:

    ExampleModelCollectFeaturesTop *Clone() const
    {
      return new ExampleModelCollectFeaturesTop( *this);
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

      // test constructor from label
      util::Implementation< model::CollectFeaturesInterface> top_three( "Top(3)");

      // check that the implementation was defined
      BCL_ExampleIndirectAssert( top_three.IsDefined(), true, "Constructor from label");
      BCL_ExampleIndirectCheck( top_three.GetString(), "Top(3)", "Constructor from label");

    /////////////////
    // data access //
    /////////////////

      // test on an example data set
      const float test_values[] = { 1.0, 2.0, 4.0, 0.5, 0.25, 8.0};
      const linal::Vector< float> test_values_linal( 6, test_values);

      // test collect
      BCL_ExampleCheck( top_three->Collect( test_values_linal), storage::Vector< size_t>::Create( 1, 2, 5));

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelCollectFeaturesTop

  const ExampleClass::EnumType ExampleModelCollectFeaturesTop::s_Instance
  (
    GetExamples().AddEnum( ExampleModelCollectFeaturesTop())
  );

} // namespace bcl
