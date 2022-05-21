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
#include "model/bcl_model_dtree_binary_partition.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_dtree_binary_partition.cpp
  //!
  //! @author mendenjl
  //! @date May 09, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelDtreeBinaryPartition :
    public ExampleInterface
  {
  public:

    ExampleModelDtreeBinaryPartition *Clone() const
    {
      return new ExampleModelDtreeBinaryPartition( *this);
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

      const size_t feature_index( 5);
      const float split_value( -7.3);
      const float split_rating( 9.4);

      // constructors
      model::DtreeBinaryPartition default_partition;
      model::DtreeBinaryPartition partition( feature_index, split_value, split_rating);

    /////////////////
    // data access //
    /////////////////

      // make sure that the default partition is undefined
      BCL_ExampleCheck( util::IsDefined( model::DtreeBinaryPartition().GetFeatureIndex()), false);
      BCL_ExampleCheck( util::IsDefined( model::DtreeBinaryPartition().GetSplitValue()), false);
      BCL_ExampleCheck( model::DtreeBinaryPartition().GetSplitRating(), -std::numeric_limits< float>::max());

      // make sure that the partition has the correct values
      BCL_ExampleCheck( partition.GetFeatureIndex(), feature_index);
      BCL_ExampleCheck( partition.GetSplitValue(), split_value);
      BCL_ExampleCheck( partition.GetSplitRating(), split_rating);

      // end
      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelDtreeBinaryPartition

  const ExampleClass::EnumType ExampleModelDtreeBinaryPartition::s_Instance
  (
    GetExamples().AddEnum( ExampleModelDtreeBinaryPartition())
  );

} // namespace bcl
