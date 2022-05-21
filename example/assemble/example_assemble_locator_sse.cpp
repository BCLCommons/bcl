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
#include "assemble/bcl_assemble_locator_sse.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_locator_sse.cpp
  //!
  //! @author loweew
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleLocatorSSE :
    public ExampleInterface
  {
  public:

    ExampleAssembleLocatorSSE *Clone() const
    {
      return new ExampleAssembleLocatorSSE( *this);
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
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1C1D.pdb"));

      // create ProteinModel "protein_model" from "pdb"
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 0;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

      // test default constructor
      assemble::LocatorSSE def_const;

      // test SetSSEID and GetSSEID functions
      BCL_MessageStd( "test SetSSEID and GetSSEID functions");
      def_const.SetSSEID( 48, 70);
      BCL_Example_Check
      (
        def_const.GetSSEID().First() == 48 && def_const.GetSSEID().Second() == 70,
        "SetSSEID or GetSSEID function does not work"
      );

      // test construct from chain ID and SSE identifiers
      assemble::LocatorSSE SSE_known( 'B', 91, 109);
      BCL_MessageStd( "test construct from chain id and SSE identifiers");
      BCL_Example_Check
      (
        SSE_known.GetSSEID().First() == 91 && SSE_known.GetSSEID().Second() == 109 &&
        SSE_known.GetChainID() == 'B', "construct from chain and SSE identifiers does not work"
      );

      // test copy constructor
      assemble::LocatorSSE copy_const( SSE_known);
      BCL_MessageStd( "test copy constructor");
      BCL_Example_Check
      (
        copy_const.GetSSEID().First() == 91 && copy_const.GetSSEID().Second() == 109 &&
        copy_const.GetChainID() == 'B', "copy constructor does not work"
      );

      // test GetChain function
      BCL_MessageStd( "test GetChain function");
      BCL_Example_Check
      (
        copy_const.Locate( protein_model)->GetChainID() == 'B', "GetChain function does not work"
      );

      // test GetSSE function
      BCL_MessageStd( "test GetSSE function");
      BCL_Example_Check
      (
        copy_const.Locate( protein_model)->GetFirstAA()->GetSeqID() == 91 &&
        copy_const.Locate( protein_model)->GetLastAA()->GetSeqID() == 109 &&
        copy_const.Locate( protein_model)->GetChainID() == 'B'
        , "GetSSE function does not work"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // write the object
      WriteBCLObject( def_const);

      // read the object back in
      assemble::LocatorSSE locator_read;
      ReadBCLObject( locator_read);

      BCL_ExampleIndirectCheck( def_const.GetChainID(), locator_read.GetChainID(), "read and write");
      BCL_ExampleIndirectCheck( def_const.GetSSEIDString(), locator_read.GetSSEIDString(), "read and write");

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleLocatorSSE

  const ExampleClass::EnumType ExampleAssembleLocatorSSE::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleLocatorSSE())
  );

} // namespace bcl
