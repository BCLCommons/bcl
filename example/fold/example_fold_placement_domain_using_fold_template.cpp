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
#include "fold/bcl_fold_placement_domain_using_fold_template.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_placement_domain_using_fold_template.cpp
  //!
  //! @author weinerbe
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldPlacementDomainUsingFoldTemplate :
    public ExampleInterface
  {
  public:

    ExampleFoldPlacementDomainUsingFoldTemplate *Clone() const
    {
      return new ExampleFoldPlacementDomainUsingFoldTemplate( *this);
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
      // create string "pdb_filename" which has path for example pdb file
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      // create protein model from pdb
      BCL_MessageStd( "building model");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 9;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      BCL_MessageStd( "testing default constructor");
      fold::PlacementDomainUsingFoldTemplate def_construct;

      // test constructor with a positive body deviation
      BCL_MessageStd( "testing constructor with a positive body deviation");
      const int positive_body_deviation( 2);
      fold::PlacementDomainUsingFoldTemplate positive_construct( positive_body_deviation);

      // test constructor with a negative body deviation
      BCL_MessageStd( "testing constructor with a negative body deviation");
      const int negative_body_deviation( -1);
      fold::PlacementDomainUsingFoldTemplate negative_construct( negative_body_deviation);

      // test Clone function
      BCL_MessageStd( "testing clone function");
      util::ShPtr< fold::PlacementDomainUsingFoldTemplate> clone_construct( def_construct.Clone());

    /////////////////
    // data access //
    /////////////////

      // check GetStaticClassName
      const std::string correct_static_class_name( "bcl::fold::PlacementDomainUsingFoldTemplate");
      BCL_Example_Check
      (
        GetStaticClassName< fold::PlacementDomainUsingFoldTemplate>() == correct_static_class_name,
        "GetStaticClassName gives " + GetStaticClassName< fold::PlacementDomainUsingFoldTemplate>() +
        " but should give " + correct_static_class_name
      );

      // check GetClassIdentifier
      BCL_Example_Check
      (
        GetStaticClassName< fold::PlacementDomainUsingFoldTemplate>() == clone_construct->GetClassIdentifier(),
        "GetClassIdentifier gives " + clone_construct->GetClassIdentifier() + " but should give " +
        correct_static_class_name
      );

    ////////////////
    // operations //
    ////////////////

      // test Place function into a protein model
      BCL_MessageStd( "Testing Place function with protein model");

      assemble::ProteinModel empty_model;
      storage::Pair
      <
        storage::Map< util::SiPtr< const assemble::SSE>, math::TransformationMatrix3D>, bool
      > transformations_a
      (
        def_construct.Place( *( protein_model.GetChains().FirstElement()), empty_model)
      );
      BCL_Example_Check
      (
        transformations_a.Second(),
        "Place function with protein model failed to produce proper transformation matrices:\n" +
        util::Format()( transformations_a.First())
      );

      // test Place function without a protein model
      BCL_MessageStd( "Testing Place function without protein model");
      storage::Pair
      <
        storage::Map< util::SiPtr< const assemble::SSE>, math::TransformationMatrix3D>, bool
      > transformations_b
      (
        def_construct.Place( *( protein_model.GetChains().FirstElement()))
      );
      BCL_Example_Check
      (
        transformations_b.Second(),
        "Place function without a protein model failed to produce proper transformation matrices:\n" +
        util::Format()( transformations_b.First())
      );

      // test Place function into a smaller template
      BCL_MessageStd( "Testing Place function into a smaller template");
      storage::Pair
      <
        storage::Map< util::SiPtr< const assemble::SSE>, math::TransformationMatrix3D>, bool
      > transformations_d
      (
        negative_construct.Place( *( protein_model.GetChains().FirstElement()))
      );
      BCL_Example_Check
      (
        transformations_d.Second(),
        "Place function into a smaller template failed to produce proper transformation matrices:\n" +
        util::Format()( transformations_d.First())
      );

      // test Place function into a larger template
      BCL_MessageStd( "Testing Place function into a larger template");
      storage::Pair
      <
        storage::Map< util::SiPtr< const assemble::SSE>, math::TransformationMatrix3D>, bool
      > transformations_e
      (
        positive_construct.Place( *( protein_model.GetChains().FirstElement()))
      );
      BCL_Example_Check
      (
        transformations_e.Second(),
        "Place function into a larger template failed to produce proper transformation matrices:\n" +
        util::Format()( transformations_e.First())
      );

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities");
      // write the object
      WriteBCLObject( positive_construct);

      // read the object back in
      fold::PlacementDomainUsingFoldTemplate fold_placement_read;
      ReadBCLObject( fold_placement_read);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldPlacementDomainUsingFoldTemplate

  const ExampleClass::EnumType ExampleFoldPlacementDomainUsingFoldTemplate::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldPlacementDomainUsingFoldTemplate())
  );

} // namespace bcl
