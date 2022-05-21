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
#include "nmr/bcl_nmr_star_tag_category_data.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_nmr_star_tag_category_data.cpp
  //!
  //! @author weinerbe
  //! @date Jul 2, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleNmrStarTagCategoryData :
    public ExampleInterface
  {
  public:

    ExampleNmrStarTagCategoryData *Clone() const
    {
      return new ExampleNmrStarTagCategoryData( *this);
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
      nmr::StarTagCategoryData def_construct;
      BCL_ExampleIndirectCheck
      (
        def_construct.GetSaveFrame() == nmr::StarTagCategoryData::s_NumberSaveFrames &&
        def_construct.GetDescription() == "",
        true, "Default Constructor"
      );

      // test constructor from save frame and description
      BCL_MessageStd( "testing constructor from save frame and description");
      nmr::StarTagCategoryData tag_category_data( nmr::StarTagCategoryData::e_DistanceConstraints, "_Dist_constraint");

      // test clone function
      util::ShPtr< nmr::StarTagCategoryData> clone_construct( tag_category_data.Clone());
      BCL_ExampleIndirectCheck
      (
        tag_category_data.GetSaveFrame() == clone_construct->GetSaveFrame() &&
        tag_category_data.GetDescription() == clone_construct->GetDescription(),
        true, "Clone function"
      );

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      const std::string correct_static_class_name( "bcl::nmr::StarTagCategoryData");
      BCL_ExampleCheck( GetStaticClassName< nmr::StarTagCategoryData>(), clone_construct->GetClassIdentifier());

      // check GetSaveFrame
      BCL_ExampleCheck( tag_category_data.GetSaveFrame(), nmr::StarTagCategoryData::e_DistanceConstraints);

      // check GetDescription
      BCL_ExampleCheck( tag_category_data.GetDescription(), "_Dist_constraint");

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities");
      // write the object
      WriteBCLObject( tag_category_data);

      // read the object back in
      nmr::StarTagCategoryData tag_category_read;
      ReadBCLObject( tag_category_read);

      // compare the objects
      BCL_ExampleIndirectCheck
      (
        tag_category_data.GetSaveFrame() == tag_category_read.GetSaveFrame() &&
        tag_category_data.GetDescription() == tag_category_read.GetDescription(),
        true, "Read and write"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleNmrStarTagCategoryData

  const ExampleClass::EnumType ExampleNmrStarTagCategoryData::s_Instance
  (
    GetExamples().AddEnum( ExampleNmrStarTagCategoryData())
  );

} // namespace bcl
