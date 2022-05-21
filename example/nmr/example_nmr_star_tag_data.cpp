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
#include "nmr/bcl_nmr_star_tag_data.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_nmr_star_tag_data.cpp
  //!
  //! @author weinerbe
  //! @date Jul 2, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleNmrStarTagData :
    public ExampleInterface
  {
  public:

    ExampleNmrStarTagData *Clone() const
    {
      return new ExampleNmrStarTagData( *this);
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

      // test default constructor
      nmr::StarTagData def_construct;
      BCL_ExampleIndirectCheck
      (
        def_construct.GetTagCategory() == nmr::GetStarTagCategories().e_Undefined &&
        def_construct.GetDescription() == "" &&
        def_construct.GetDataType() == util::CPPDataTypes::e_Unknown,
        true, "Default constructor"
      );

      // test constructor from save frame and description
      BCL_MessageStd( "testing constructor from category and description");
      nmr::StarTagData tag_data( nmr::GetStarTagCategories().e_DistConstraintValue, "Distance_val", util::CPPDataTypes::e_Double);

      // test clone function
      util::ShPtr< nmr::StarTagData> clone_construct( tag_data.Clone());
      BCL_ExampleIndirectCheck
      (
        tag_data.GetTagCategory() == clone_construct->GetTagCategory() &&
        tag_data.GetDescription() == clone_construct->GetDescription() &&
        tag_data.GetDataType() == clone_construct->GetDataType(),
        true, "Clone construct"
      );

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      const std::string correct_static_class_name( "bcl::nmr::StarTagData");
      BCL_ExampleCheck( GetStaticClassName< nmr::StarTagData>(), clone_construct->GetClassIdentifier());

      // check GetTagCategory
      BCL_ExampleCheck( tag_data.GetTagCategory(), nmr::GetStarTagCategories().e_DistConstraintValue);

      // check GetDescription
      BCL_ExampleCheck( tag_data.GetDescription(), "Distance_val");

      // check GetDataType
      BCL_ExampleCheck( tag_data.GetDataType(), util::CPPDataTypes::e_Double);

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // write the object
      WriteBCLObject( tag_data);

      // read the object back in
      nmr::StarTagData tag_data_read;
      ReadBCLObject( tag_data_read);

      // compare the objects
      BCL_ExampleIndirectCheck
      (
        tag_data.GetTagCategory() == tag_data_read.GetTagCategory() &&
        tag_data.GetDescription() == tag_data_read.GetDescription() &&
        tag_data.GetDataType() == tag_data_read.GetDataType(),
        true, "Read and write"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleNmrStarTagData

  const ExampleClass::EnumType ExampleNmrStarTagData::s_Instance
  (
    GetExamples().AddEnum( ExampleNmrStarTagData())
  );

} // namespace bcl
