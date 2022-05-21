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
#include "assemble/bcl_assemble_sse_pair_template.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_sse.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_sse_pair_template.cpp
  //!
  //! @author weinerbe
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleSSEPairTemplate :
    public ExampleInterface
  {
  public:

    ExampleAssembleSSEPairTemplate *Clone() const
    {
      return new ExampleAssembleSSEPairTemplate( *this);
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
      // create string which has path for example pdb file
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      // get SSEs
      util::ShPtr< assemble::SSE> sse_helix( Proteins::GetSSE( pdb_filename, 'A', 23, 34));
      util::ShPtr< assemble::SSE> sse_strand( Proteins::GetSSE( pdb_filename, 'A', 40, 45));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      assemble::SSEPairTemplate def_construct;

      // test construction two SSE geometries and a loop length
      BCL_MessageStd( "test construct from two SSEs");
      assemble::SSEPairTemplate sse_construct
      (
        sse_helix, sse_strand, biol::CalculateSequenceDistance( *sse_helix, *sse_strand)
      );

      // test clone constructor
      util::ShPtr< assemble::SSEPairTemplate> clone_construct( sse_construct.Clone());
      BCL_ExampleIndirectCheck
      (
        clone_construct->GetPacking().GetContactType() == sse_construct.GetPacking().GetContactType() &&
        clone_construct->GetLoopLength() == sse_construct.GetLoopLength() &&
        clone_construct->GetFirstGeometry() == sse_construct.GetFirstGeometry() &&
        clone_construct->GetSecondGeometry() == sse_construct.GetSecondGeometry(),
        true,
        "clone constructor"
      );

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck( GetStaticClassName< assemble::SSEPairTemplate>(), clone_construct->GetClassIdentifier());

      // check GetFirstGeometry
      BCL_ExampleCheck( sse_construct.GetFirstGeometry(), sse_helix);

      // check GetSecondGeometry
      BCL_ExampleCheck( sse_construct.GetSecondGeometry(), sse_strand);

      // check GetLoopLength
      BCL_ExampleCheck( sse_construct.GetLoopLength(), 5);

      // check GetPacking
      BCL_ExampleIndirectCheck
      (
        sse_construct.GetPacking().GetContactType(),
        assemble::SSEGeometryPacking( *sse_helix, *sse_strand).GetContactType(),
        "GetPacking()"
      );

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
      WriteBCLObject( sse_construct);

      // read the object back in
      assemble::SSEPairTemplate sse_pair_template_read;
      ReadBCLObject( sse_pair_template_read);

      // compare the objects
      BCL_MessageStd( "compare written and read objects");
      BCL_ExampleIndirectCheck
      (
        sse_pair_template_read.GetPacking().GetContactType() == sse_construct.GetPacking().GetContactType() &&
        sse_pair_template_read.GetLoopLength() == sse_construct.GetLoopLength() &&
        math::SimilarWithinTolerance
        (
          sse_pair_template_read.GetFirstGeometry()->GetOrientation(),
          sse_construct.GetFirstGeometry()->GetOrientation(), 0.01, 0.01
        ) &&
        math::SimilarWithinTolerance
        (
          sse_pair_template_read.GetSecondGeometry()->GetOrientation(),
          sse_construct.GetSecondGeometry()->GetOrientation(), 0.01, 0.01
        ),
        true,
        "read and write"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleSSEPairTemplate

  const ExampleClass::EnumType ExampleAssembleSSEPairTemplate::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleSSEPairTemplate())
  );

} // namespace bcl
