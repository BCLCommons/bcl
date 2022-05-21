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
#include "assemble/bcl_assemble_sse_compare_extent.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_sse_compare_extent.cpp
  //!
  //! @author linders
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //! Example to test SSEsMatchExtent class.
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleSSECompareExtent :
    public ExampleInterface
  {
  public:

    ExampleAssembleSSECompareExtent *Clone() const
    {
      return new ExampleAssembleSSECompareExtent( *this);
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

    /////////////////
    // preparation //
    /////////////////

      // create string "pdb_filename" which has path for example pdb file
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi_ideal_model.pdb"));

      // create ProteinModel "protein_model" from "pdb"
      BCL_MessageStd( "building models from pdb chains and sse information");

      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 0;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 0;

      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

      // get the sses from the protein model
      util::SiPtrVector< const assemble::SSE> all_sses( protein_model.GetSSEs());

      // all sses are strands, sse_a is 7 residues, sse_b is 8 residues, sse_c is 6 residues
      const assemble::SSE sse_a( *( all_sses( 0)));
      const assemble::SSE sse_b( *( all_sses( 1)));
      const assemble::SSE sse_c( *( all_sses( 3)));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      assemble::SSECompareExtent extent_matcher_default;

      // construct SSEsMatchExtent from tolerance and axis
      assemble::SSECompareExtent extent_matcher( 3.0, coord::GetAxes().e_Z);

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

      // check operator with two sses whose extents are within the tolerance
      BCL_Example_Check
      (
        extent_matcher( sse_a, sse_b) && extent_matcher( sse_a, sse_c),
        "they should all match within the tolerance"
      );

      // check operator with two sses whose extents are not within the tolerance
      BCL_Example_Check
      (
        !extent_matcher( sse_b, sse_c),
        "sse_b and sse_c should not match within the tolerance"
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities for assemble::SSEsMatchExtent");
      // open output file
      BCL_MessageVrb( "write object");
      WriteBCLObject( extent_matcher);
      BCL_MessageVrb( "read object");
      assemble::SSECompareExtent extent_matcher_read;
      ReadBCLObject( extent_matcher_read);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleSSECompareExtent

  const ExampleClass::EnumType ExampleAssembleSSECompareExtent::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleSSECompareExtent())
  );
  
} // namespace bcl
