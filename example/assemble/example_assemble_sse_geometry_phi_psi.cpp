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
#include "assemble/bcl_assemble_sse_geometry_phi_psi.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_sse.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_sse_geometry_phi_psi.cpp
  //!
  //! @author weinerbe
  //! @date Feb 11, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleSSEGeometryPhiPsi :
    public ExampleInterface
  {
  public:

    ExampleAssembleSSEGeometryPhiPsi *Clone() const
    {
      return new ExampleAssembleSSEGeometryPhiPsi( *this);
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

      // get an SSE
      util::ShPtr< assemble::SSE> sp_sse( Proteins::GetSSE( pdb_filename, 'A', 23, 34));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      const assemble::SSEGeometryPhiPsi def_construct;
      BCL_ExampleIndirectCheck
      (
        !def_construct.GetSSEGeometry().IsDefined() && !def_construct.GetPhiPsi().IsDefined(),
        true,
        "default constructor"
      );

      // test constructor from SSE
      assemble::SSEGeometryPhiPsi sse_construct( *sp_sse);

      // test constructor from AASequencePhiPsi
      const biol::AASequencePhiPsi aa_phi_psi( *sp_sse);
      assemble::SSEGeometryPhiPsi phi_psi_construct( aa_phi_psi, biol::GetSSTypes().HELIX);
      BCL_ExampleIndirectCheck
      (
        phi_psi_construct.GetPhiPsi()->IsDefined() &&
        phi_psi_construct.GetType() == biol::GetSSTypes().HELIX &&
        !phi_psi_construct.GetSSEGeometry()->IsDefined(),
        true,
        "constructor from AASequencePhiPsi"
      );

      // test copy constructor
      const assemble::SSEGeometryPhiPsi copy_construct( sse_construct);
      BCL_ExampleIndirectCheck
      (
        copy_construct.GetPhiPsi()->GetAngles().GetSize() == sse_construct.GetPhiPsi()->GetAngles().GetSize() &&
        copy_construct.GetSSEGeometry()->GetFragments().GetSize() ==
          sse_construct.GetSSEGeometry()->GetFragments().GetSize(),
        true,
        "copy constructor"
      );

    /////////////////
    // data access //
    /////////////////

      // test SetSSEGeometryUsingPhiPsi
      phi_psi_construct.SetSSEGeometryUsingPhiPsi();
      BCL_ExampleIndirectCheck( phi_psi_construct.GetSSEGeometry().IsDefined(), true, "SetSSEGeometryUsingPhiPsi");

    //////////////////////
    // input and output //
    //////////////////////

      // test read and write
      WriteBCLObject( sse_construct);
      assemble::SSEGeometryPhiPsi read_construct;
      ReadBCLObject( read_construct);
      BCL_ExampleIndirectCheck
      (
        read_construct.GetPhiPsi()->GetAngles().GetSize() == sse_construct.GetPhiPsi()->GetAngles().GetSize() &&
        read_construct.GetSSEGeometry()->GetFragments().GetSize() ==
          sse_construct.GetSSEGeometry()->GetFragments().GetSize(),
        true,
        "read and write"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleSSEGeometryPhiPsi

  const ExampleClass::EnumType ExampleAssembleSSEGeometryPhiPsi::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleSSEGeometryPhiPsi())
  );

} // namespace bcl
