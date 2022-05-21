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
#include "assemble/bcl_assemble_sse_compare_type.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_sse_compare_type.cpp
  //!
  //! @author linders
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleSSECompareType :
    public ExampleInterface
  {
  public:

    ExampleAssembleSSECompareType *Clone() const
    {
      return new ExampleAssembleSSECompareType( *this);
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

      // sse_a and sse_b are strands, sse_c is helix
      const assemble::SSE sse_a( *( all_sses( 0)));
      const assemble::SSE sse_b( *( all_sses( 1)));
      const assemble::SSE sse_c( *( all_sses( 2)));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct a SSEsMatchType
      assemble::SSECompareType type_matcher;

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

      // check operator with two sses of same type
      BCL_Example_Check
      (
        type_matcher( sse_a, sse_b), "sse_a and sse_b should both be strands and thus of the same type"
      );

      // check operator with two sses of different type
      BCL_Example_Check
      (
        !type_matcher( sse_a, sse_c) && !type_matcher( sse_b, sse_c),
        "sse_c is a helix and thus should be of the different type as sse_a and sse_b (both strands)"
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities for assemble::SSEsMatchType");
      // open output file
      BCL_MessageVrb( "write object");
      WriteBCLObject( type_matcher);
      BCL_MessageVrb( "read object");
      assemble::SSECompareType type_matcher_read;
      ReadBCLObject( type_matcher_read);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleSSECompareType

  const ExampleClass::EnumType ExampleAssembleSSECompareType::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleSSECompareType())
  );
  
} // namespace bcl
