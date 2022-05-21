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
#include "assemble/bcl_assemble_locator_aa.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_locator_aa.cpp
  //!
  //! @author alexanns
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleLocatorAA :
    public ExampleInterface
  {
  public:

    ExampleAssembleLocatorAA *Clone() const
    {
      return new ExampleAssembleLocatorAA( *this);
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

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      assemble::LocatorAA def_const;

      // test construct from chain ID and AA identifiers
      assemble::LocatorAA AA_known( 'B', 52);
      BCL_MessageStd( "test construct from chain id and AA identifiers");
      BCL_Example_Check
      (
        AA_known.GetAAID() == 52 && AA_known.GetLocatorChain().GetChainID() == 'B',
        "Residue ID should be 52 but is " + util::Format()( AA_known.GetAAID()) +
        " and Chain ID should be \"B\" but is " + util::Format()( AA_known.GetLocatorChain().GetChainID())
      );

      // test copy constructor
      assemble::LocatorAA copy_const( AA_known);
      BCL_MessageStd( "test copy constructor");
      BCL_Example_Check
      (
        copy_const.GetAAID() == AA_known.GetAAID() &&
        copy_const.GetLocatorChain().GetChainID() == AA_known.GetLocatorChain().GetChainID(),
        "copy should be " + util::Format()( AA_known) + " but instead is " + util::Format()( copy_const)
      );

      // test clone constructor
      util::ShPtr< assemble::LocatorAA> clone_constr( copy_const.Clone());
      BCL_Example_Check
      (
        clone_constr->GetAAID() == copy_const.GetAAID() &&
        clone_constr->GetLocatorChain().GetChainID() == copy_const.GetLocatorChain().GetChainID(),
        "clone should be " + util::Format()( copy_const) + " but instead is " + util::Format()( *clone_constr)
      );

    /////////////////
    // data access //
    /////////////////

      // check GetStaticClassName
      const std::string correct_static_class_name( "bcl::assemble::LocatorAA");

      // check GetClassIdentifier
      BCL_ExampleCheck( GetStaticClassName< assemble::LocatorAA>(), clone_constr->GetClassIdentifier());

      // check GetIdentification
      BCL_ExampleCheck( "  B  52  \"UsePDBID\"  0", clone_constr->GetIdentification());

      // test SetAAID and GetAAID functions
      def_const.SetAAID( 15);
      BCL_ExampleIndirectCheck( def_const.GetAAID(), 15, "SetAAID(15)");

      // test GetLocatorChain and SetLocatorChain
      def_const.SetLocatorChain( assemble::LocatorChain( 'C'));
      BCL_ExampleIndirectCheck
      (
        def_const.GetLocatorChain().GetChainID(),
        'C',
        "SetLocatorChain( assemble::LocatorChain( 'C')"
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // test locate function taking an sse
      BCL_MessageStd( "test Locate function taking an sse");

      // get the sses from the protein model
      const util::SiPtrVector< const assemble::SSE> sses( protein_model.GetSSEs());

      // locate a residue from the protein model given an sse
      util::SiPtr< const biol::AABase> located_aa( assemble::LocatorAA( 'A', 7).Locate( *sses.FirstElement()));
      BCL_ExampleCheck( located_aa.IsDefined(), true);
      BCL_ExampleCheck( located_aa->GetSeqID(), 7);

      // test LocateSSE function
      BCL_MessageStd( "test LocateSSE function");
      BCL_Example_Check
      (
        copy_const.LocateSSE( protein_model)->GetFirstAA()->GetSeqID() == 48 &&
        copy_const.LocateSSE( protein_model)->GetLastAA()->GetSeqID() == 70 &&
        copy_const.LocateSSE( protein_model)->GetChainID() == 'B',
        "SeqID of first residue in the located sse should be 48 but is " +
        util::Format()( copy_const.LocateSSE( protein_model)->GetFirstAA()->GetSeqID()) +
        " and SeqID of last residue in the located sse should be 70 but is " +
        util::Format()( copy_const.LocateSSE( protein_model)->GetLastAA()->GetSeqID()) +
        " and chain id should be B but is " +
        util::Format()( copy_const.LocateSSE( protein_model)->GetChainID())
      );

      BCL_ExampleCheck( copy_const.Locate( protein_model)->GetSeqID(), 52);

    //////////////////////
    // input and output //
    //////////////////////

       // test Write and read
       WriteBCLObject( copy_const);
       assemble::LocatorAA read_locator_aa;
       ReadBCLObject( read_locator_aa);

       // make sure that the written locator and the read-in locator are the same
       BCL_ExampleIndirectCheck( copy_const.GetAAID(), read_locator_aa.GetAAID(), "I/O");
       BCL_ExampleIndirectCheck
       (
         copy_const.GetLocatorChain().GetChainID(),
         read_locator_aa.GetLocatorChain().GetChainID(),
         "I/O"
       );

    //////////////////////
    // helper functions //
    //////////////////////

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleLocatorAA

  const ExampleClass::EnumType ExampleAssembleLocatorAA::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleLocatorAA())
  );

} // namespace bcl
