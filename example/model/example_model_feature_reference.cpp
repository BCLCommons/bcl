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
#include "model/bcl_model_feature_reference.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_feature_reference.cpp
  //!
  //! @author loweew
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelFeatureReference :
    public ExampleInterface
  {
  public:

    ExampleModelFeatureReference *Clone() const
    {
      return new ExampleModelFeatureReference( *this);
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
      // number elements per dimension
      const size_t number_elements( 3);

      // constructing data for feature reference construction
      const linal::Matrix< float> data( number_elements, number_elements, 1.5);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // constructor
      const model::FeatureReference< float> feat_ref_a( number_elements, data.Begin());
      const model::FeatureReference< float> feat_ref_b( number_elements, data.Begin() + number_elements);
      const model::FeatureReference< float> feat_ref_c( number_elements, data.Begin() + 2 * number_elements);

      // clone
      const util::ShPtr< model::FeatureReference< float> > feat_ref_a_clone( feat_ref_a.Clone());

    /////////////////
    // data access //
    /////////////////

      // get size
      BCL_ExampleCheck( feat_ref_a.GetSize(), number_elements);

      // testing begin and end
      BCL_ExampleCheck( size_t( feat_ref_b.End() - feat_ref_b.Begin()), feat_ref_b.GetSize());

      // has ownerships
      BCL_ExampleCheck( feat_ref_c.HasOwnership(), false);

    //////////////////////
    // input and output //
    //////////////////////

      WriteBCLObject( feat_ref_a);
      model::FeatureReference< float> read_feat_ref( 3, linal::Matrix< float>( 3, 0).Begin());
      ReadBCLObject( read_feat_ref);

      // end
      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelFeatureReference

  const ExampleClass::EnumType ExampleModelFeatureReference::s_Instance
  (
    GetExamples().AddEnum( ExampleModelFeatureReference())
  );

} // namespace bcl
