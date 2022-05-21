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
#include "assemble/bcl_assemble_pick_sse_short_loops.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_collector_sse.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_pick_sse_short_loops.cpp
  //!
  //! @author weinerbe
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssemblePickSSEShortLoops :
    public ExampleInterface
  {
  public:

    ExampleAssemblePickSSEShortLoops *Clone() const
    {
      return new ExampleAssemblePickSSEShortLoops( *this);
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
      BCL_MessageStd( "testing default constructor");
      assemble::PickSSEShortLoops def_construct;

      // test copy constructor
      BCL_MessageStd( "testing copy constructor");
      assemble::PickSSEShortLoops copy_construct( def_construct);

      // test clone constructor
      BCL_MessageStd( "testing clone constructor");
      util::ShPtr< assemble::PickSSEShortLoops> clone_construct( def_construct.Clone());

    /////////////////
    // data access //
    /////////////////

      // check GetStaticClassName
      const std::string correct_static_class_name( "bcl::assemble::PickSSEShortLoops");
      BCL_Example_Check
      (
        GetStaticClassName< assemble::PickSSEShortLoops>() == correct_static_class_name,
        "GetStaticClassName gives " + GetStaticClassName< assemble::PickSSEShortLoops>() + " but should give " +
        correct_static_class_name
      );

      // check GetClassIdentifier
      BCL_Example_Check
      (
        GetStaticClassName< assemble::PickSSEShortLoops>() == clone_construct->GetClassIdentifier(),
        "GetClassIdentifier gives " + clone_construct->GetClassIdentifier() + " but should give " +
        correct_static_class_name
      );

    ////////////////
    // operations //
    ////////////////

      // create string "pdb_filename" which has path for example pdb file
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      // create protein model from pdb
      BCL_MessageStd( "building model");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 5;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 3;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

      // create collector
      assemble::CollectorSSE collector_sse;

      // create picker w/ a defined max loop length
      assemble::PickSSEShortLoops sse_picker( 20);

      // test pick function
      BCL_MessageStd( "testing Pick function");
      util::SiPtr< const assemble::SSE> sp_sse( sse_picker.Pick( collector_sse.Collect( protein_model), protein_model));
      BCL_ExampleIndirectAssert
      (
        sp_sse->IsDefined(),
        true,
        "sse_picker.Pick( collector_sse.Collect( protein_model), protein_model)"
      );

      BCL_MessageStd( "randomly choosen sse: " + sp_sse->GetIdentification());

    ///////////////
    // operators //
    ///////////////

      // test operator =
      BCL_MessageStd( "testing default constructor");
      def_construct = sse_picker;

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities");
      // write the object
      WriteBCLObject( sse_picker);

      // read the object back in
      assemble::PickSSEShortLoops picker_read;
      ReadBCLObject( picker_read);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssemblePickSSEShortLoops

  const ExampleClass::EnumType ExampleAssemblePickSSEShortLoops::s_Instance
  (
    GetExamples().AddEnum( ExampleAssemblePickSSEShortLoops())
  );

} // namespace bcl
