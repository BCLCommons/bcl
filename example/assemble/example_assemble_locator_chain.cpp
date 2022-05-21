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
#include "assemble/bcl_assemble_locator_chain.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_locator_chain.cpp
  //!
  //! @author alexanns
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleLocatorChain :
    public ExampleInterface
  {
  public:

    ExampleAssembleLocatorChain *Clone() const
    {
      return new ExampleAssembleLocatorChain( *this);
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
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 0;
      ssetype_min_size[ biol::GetSSTypes().COIL] = 0;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      assemble::LocatorChain def_const;

      // test construct from chain ID
      assemble::LocatorChain chain_known( 'B');
      BCL_MessageStd( "test construct from chain id");
      BCL_Example_Check
      (
        chain_known.GetChainID() == 'B', "LocatorChain chain ID member should be \"B\" but instead is " +
        util::Format()( chain_known.GetChainID())
      );

      // test copy constructor
      assemble::LocatorChain copy_const( chain_known);
      BCL_MessageStd( "test copy constructor");
      BCL_Example_Check
      (
        copy_const.GetChainID() == chain_known.GetChainID(), "Copy should be " + util::Format()( chain_known) +
        " but instead is " + util::Format()( copy_const)
      );

      // test clone constructor
      util::ShPtr< assemble::LocatorChain> clone_constr( copy_const.Clone());
      BCL_Example_Check
      (
        clone_constr->GetChainID() == copy_const.GetChainID(), "clone should be " + util::Format()( copy_const) +
        " but instead is " + util::Format()( *clone_constr)
      );

    /////////////////
    // data access //
    /////////////////

      // check GetStaticClassName
      const std::string correct_static_class_name( "bcl::assemble::LocatorChain");
      BCL_Example_Check
      (
        GetStaticClassName< assemble::LocatorChain>() == correct_static_class_name,
        "GetStaticClassName gives " + GetStaticClassName< assemble::LocatorChain>() + " but should give " +
        correct_static_class_name
      );

      // check GetClassIdentifier
      BCL_Example_Check
      (
        GetStaticClassName< assemble::LocatorChain>() == clone_constr->GetClassIdentifier(),
        "GetClassIdentifier gives " + clone_constr->GetClassIdentifier() + " but should give " +
        correct_static_class_name
      );

      // test SetChainID and GetChainID functions
      BCL_MessageStd( "test SetChainID and GetChainID functions");
      def_const.SetChainID( 'A');
      BCL_Example_Check
      (
        def_const.GetChainID() == 'A', "SetChainID or GetChainID function does not work"
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // test GetChain function
      BCL_MessageStd( "test GetChain function");
      BCL_Example_Check
      (
        copy_const.Locate( protein_model)->GetChainID() == 'B',
        "GetChain function should have returned chain \"B\" but instead returned chain " +
        util::Format()( copy_const.Locate( protein_model)->GetChainID())
      );

    //////////////////////
    // input and output //
    //////////////////////

      // test Write and read
      WriteBCLObject( copy_const);
      assemble::LocatorChain read_locator_chain;
      ReadBCLObject( read_locator_chain);

      // make sure that the written locator and the read-in locator are the same
      BCL_ExampleIndirectCheck( copy_const.GetChainID(), read_locator_chain.GetChainID(), "I/O");

    //////////////////////
    // helper functions //
    //////////////////////

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleLocatorChain

  const ExampleClass::EnumType ExampleAssembleLocatorChain::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleLocatorChain())
  );

} // namespace bcl
