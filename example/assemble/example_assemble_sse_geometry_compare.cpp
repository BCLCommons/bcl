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
#include "assemble/bcl_assemble_sse_geometry_compare.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_sse.h"
#include "assemble/bcl_assemble_sse_geometry_phi_psi.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_sse_geometry_compare.cpp
  //! @details tests the SSE geometry compare classes
  //!
  //! @author weinerbe
  //! @date Nov 11, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleSSEGeometryCompare :
    public ExampleInterface
  {
  public:

    ExampleAssembleSSEGeometryCompare *Clone() const
    {
      return new ExampleAssembleSSEGeometryCompare( *this);
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

      // get SSEs
      const util::ShPtr< assemble::SSE> sp_sse_helix_23_34( Proteins::GetSSE( pdb_filename, 'A', 23, 34));
      const util::ShPtr< assemble::SSE> sp_sse_strand_40_45( Proteins::GetSSE( pdb_filename, 'A', 40, 45));
      const util::ShPtr< assemble::SSE> sp_sse_strand_64_72( Proteins::GetSSE( pdb_filename, 'A', 64, 72));
      const assemble::SSEGeometryPhiPsi sp_phi_psi_40_45( *sp_sse_strand_40_45);

    ////////////////////////////////////
    // SSEGeometryWithinSizeTolerance //
    ////////////////////////////////////

      // test default constructor
      assemble::SSEGeometryWithinSizeTolerance def_tolerance_construct;

      // test constructor from specific sizes
      assemble::SSEGeometryWithinSizeTolerance tolerance_construct( 4, 4);

      // test Clone
      util::ShPtr< assemble::SSEGeometryWithinSizeTolerance> clone_tolerance_construct
      (
        def_tolerance_construct.Clone()
      );

      // test GetClassIdentifier
      BCL_ExampleCheck
      (
        clone_tolerance_construct->GetClassIdentifier(),
        "bcl::assemble::SSEGeometryWithinSizeTolerance"
      );

      // test () operator
      BCL_ExampleCheck( def_tolerance_construct( *sp_sse_helix_23_34, sp_phi_psi_40_45), false);
      BCL_ExampleCheck( def_tolerance_construct( *sp_sse_strand_64_72, sp_phi_psi_40_45), false);
      BCL_ExampleCheck( tolerance_construct( *sp_sse_strand_64_72, sp_phi_psi_40_45), true);

      // write the object
      WriteBCLObject( tolerance_construct);

      // read the object back in
      assemble::SSEGeometryWithinSizeTolerance tolerance_read;
      ReadBCLObject( tolerance_read);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleSSEGeometryCompare

  const ExampleClass::EnumType ExampleAssembleSSEGeometryCompare::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleSSEGeometryCompare())
  );
  
} // namespace bcl
