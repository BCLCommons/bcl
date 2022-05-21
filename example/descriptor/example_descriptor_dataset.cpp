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
#include "descriptor/bcl_descriptor_dataset.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_dataset.cpp
  //!
  //! @author mendenjl
  //! @date Dec 13, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorDataset :
    public ExampleInterface
  {
  public:

    ExampleDescriptorDataset *Clone() const
    {
      return new ExampleDescriptorDataset( *this);
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
      // make a dataset with 4 examples, 2 feature cols, 1 result col, 3 id cols
      descriptor::Dataset dataset( 4, 2, 1, 3);

      // check size of dataset
      BCL_ExampleCheck( descriptor::Dataset( 4, 2, 1, 3).GetSize(), 4);
      BCL_ExampleCheck( descriptor::Dataset( 4, 2, 1, 3).GetFeatureSize(), 2);
      BCL_ExampleCheck( descriptor::Dataset( 4, 2, 1, 3).GetResultSize(), 1);
      BCL_ExampleCheck( descriptor::Dataset( 4, 2, 1, 3).GetIdSize(), 3);

      // check that features are undefined
      BCL_ExampleIndirectCheck
      (
        util::IsDefined( dataset.GetFeaturesReference()( 0, 1)),
        false,
        "Uninitialized dataset should have undefined features"
      );
      // check that results are undefined
      BCL_ExampleIndirectCheck
      (
        util::IsDefined( dataset.GetResultsReference()( 0, 1)),
        false,
        "Uninitialized dataset should have undefined results"
      );
      // check that id's are empty
      BCL_ExampleIndirectCheck
      (
        std::string( dataset.GetIdsReference().Begin(), size_t( 12)),
        "            ",
        "Uninitialized dataset should not have ids"
      );

      // ensure that removal of all undefined examples removes everything
      BCL_ExampleCheck( descriptor::Dataset( 4, 2, 1, 3).RemoveUndefinedExamples(), 4);

      // try setting the features to see if all features and results are still removed
      dataset.GetFeaturesReference() = 0.0;

      BCL_ExampleIndirectCheck
      (
        dataset.RemoveUndefinedExamples(),
        4,
        "Undefined results should also be removed"
      );

      // reset the dataset completely
      dataset = descriptor::Dataset( 4, 2, 1, 3);

      // try setting the results
      dataset.GetResultsReference() = 1.0;

      // ensure that the result matrix was changed
      BCL_ExampleIndirectCheck( dataset.GetResultsReference()( 0, 0), 1.0, "GetResultsMatrixReference() = 1.0");

      // see whether the features and results are all still removed
      BCL_ExampleIndirectCheck
      (
        dataset.RemoveUndefinedExamples(),
        4,
        "Undefined features should also be removed"
      );
      BCL_ExampleIndirectCheck
      (
        dataset.GetResultsReference().GetNumberOfElements(),
        0,
        "Removal really removes the rows"
      );

      // reset the dataset completely
      dataset = descriptor::Dataset( 4, 2, 1, 3);

      // set the features and results
      dataset.GetResultsReference() = 1.0;
      dataset.GetFeaturesReference() = 1.0;
      dataset.GetFeaturesReference()( 1, 0) = 2.0;

      // just set one element to undefined to ensure that that row is removed
      dataset.GetFeaturesReference()( 0, 1) = util::GetUndefined< float>();

      // ensure that the first row is removed
      BCL_ExampleCheck( dataset.RemoveUndefinedExamples(), 1);
      BCL_ExampleCheck( dataset.RemoveUndefinedExamples(), 0);
      BCL_ExampleCheck( dataset.GetFeaturesReference()( 0, 0), 2.0);

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorDataset

  const ExampleClass::EnumType ExampleDescriptorDataset::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorDataset())
  );

} // namespace bcl
