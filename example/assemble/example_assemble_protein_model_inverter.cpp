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
#include "assemble/bcl_assemble_protein_model_inverter.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_protein_model_inverter.cpp
  //!
  //! @author karakam
  //! @date Nov 27, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleProteinModelInverter :
    public ExampleInterface
  {
  public:

    ExampleAssembleProteinModelInverter *Clone() const
    {
      return new ExampleAssembleProteinModelInverter( *this);
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
      // construct min sse sizes map
      storage::Map< biol::SSType, size_t> min_sse_sizes;
      min_sse_sizes[ biol::GetSSTypes().HELIX] = 9;
      min_sse_sizes[ biol::GetSSTypes().STRAND] = 4;
      min_sse_sizes[ biol::GetSSTypes().COIL] = 999;

      // get 1UBI model
      const std::string filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
      assemble::ProteinModel native_model
      (
        Proteins::GetModel( filename, biol::GetAAClasses().e_AABackBone, min_sse_sizes)
      );

      // create reference on the chain
      const assemble::Chain &native_chain( *native_model.GetChains().FirstElement());

      // get individual SSEs
      util::SiPtr< const assemble::SSE> sip_helix_23_34(  assemble::LocatorSSE( 'A', 23, 34).Locate( native_model));
      util::SiPtr< const assemble::SSE> sip_strand_1_7(   assemble::LocatorSSE( 'A',  1,  7).Locate( native_model));
      util::SiPtr< const assemble::SSE> sip_strand_10_17( assemble::LocatorSSE( 'A', 10, 17).Locate( native_model));
      util::SiPtr< const assemble::SSE> sip_strand_40_45( assemble::LocatorSSE( 'A', 40, 45).Locate( native_model));
      util::SiPtr< const assemble::SSE> sip_strand_64_72( assemble::LocatorSSE( 'A', 64, 72).Locate( native_model));

      // change min_sse_sizes to exclude everything
      min_sse_sizes[ biol::GetSSTypes().HELIX] = 999;
      min_sse_sizes[ biol::GetSSTypes().STRAND] = 999;

      // now create a new model where there are no SSEs
      assemble::ProteinModel model( native_model);
      model.FilterByMinSSESizes( min_sse_sizes);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test constructor with no cache
      assemble::ProteinModelInverter inverter_no_cache;
      BCL_ExampleCheck( inverter_no_cache.GetUseCache(), false);
      BCL_ExampleCheck( inverter_no_cache.GetModelCacheSize(), 2);
      BCL_ExampleCheck( inverter_no_cache.GetCoilCacheSize(), 30);

      // test constructo with cache
      assemble::ProteinModelInverter inverter( true, 3, 20);
      BCL_ExampleCheck( inverter.GetUseCache(), true);
      BCL_ExampleCheck( inverter.GetModelCacheSize(), 3);
      BCL_ExampleCheck( inverter.GetCoilCacheSize(), 20);

      // test Clone() constructor
      util::ShPtr< assemble::ProteinModelInverter> sp_inverter( inverter.Clone());

    /////////////////
    // data access //
    /////////////////

      // test the GetStaticClassName
      BCL_ExampleCheck( inverter.GetClassIdentifier(), GetStaticClassName( inverter));

      // test GetDescription( Chain)
      BCL_MessageStd( "testing GetDescription( Chain)");
      const std::string expected_descr( "_1_7_10_17_23_34_40_45_64_72");
      std::string descr_chain( assemble::ProteinModelInverter::GenerateDescription( native_chain));
      const size_t index_a( descr_chain.find_first_of( '_', 0));
      BCL_ExampleIndirectCheck
      (
        index_a != std::string::npos, true,
        "Could not find underscore in the chain description " + descr_chain
      );
      descr_chain = descr_chain.substr( index_a);
      BCL_ExampleCheck( descr_chain, expected_descr);

      // test GetDescription( ProteinModel)
      BCL_MessageStd( "testing GetDescription( ProteinModel)");
      std::string descr_model( assemble::ProteinModelInverter::GenerateDescription( native_model));
      const size_t index_b( descr_model.find_first_of( '_', 0));
      BCL_ExampleIndirectCheck
      (
        index_b != std::string::npos, true,
        "Could not find underscore in the model description " + descr_model
      );
      descr_model = descr_model.substr( index_a);
      BCL_ExampleCheck( descr_model, expected_descr);

    ////////////////
    // operations //
    ////////////////

      // test operator() with empty model
      BCL_MessageStd( "testing GetInvertedModel with empty model");
      {
        // get inverted model
        assemble::ProteinModel inverted_model( *inverter.GetInvertedModel( model));

        // the model should only have one SSE from 1 to 76
        util::SiPtrVector< const assemble::SSE> sses( inverted_model.GetSSEs());

        // check that there is only one SSE
        BCL_ExampleCheck( sses.GetSize(), 1);

        // get the first loop
        const assemble::SSE &loop( *sses.FirstElement());

        // make sure it's correct
        BCL_ExampleCheck( loop.GetFirstAA()->GetSeqID(), 1);
        BCL_ExampleCheck( loop.GetLastAA()->GetSeqID(), 76);
        BCL_ExampleCheck( loop.GetType(), biol::GetSSTypes().COIL);
        BCL_ExampleCheck( loop.GetChainID(), 'A');

        // write the model to PDB
        Proteins::WriteModelToPDB( inverted_model, AddExampleOutputPathToFilename( inverter, "inverted_model_a.pdb"));
      }

      // clone new SSEs
      util::ShPtr< assemble::SSE> sp_helix_23_34( sip_helix_23_34->Clone());
      util::ShPtr< assemble::SSE> sp_strand_1_7( sip_strand_1_7->Clone());
      util::ShPtr< assemble::SSE> sp_strand_10_17( sip_strand_10_17->Clone());
      util::ShPtr< assemble::SSE> sp_strand_40_45( sip_strand_40_45->Clone());
      util::ShPtr< assemble::SSE> sp_strand_64_72( sip_strand_64_72->Clone());

      // construct the sse comparisons for the remaining sections
      assemble::SSECompareByIdentity compare_1_9(    1,  9, 'A', biol::GetSSTypes().COIL);
      assemble::SSECompareByIdentity compare_8_9(    8,  9, 'A', biol::GetSSTypes().COIL);
      assemble::SSECompareByIdentity compare_18_22( 18, 22, 'A', biol::GetSSTypes().COIL);
      assemble::SSECompareByIdentity compare_35_76( 35, 76, 'A', biol::GetSSTypes().COIL);
      assemble::SSECompareByIdentity compare_35_39( 35, 39, 'A', biol::GetSSTypes().COIL);
      assemble::SSECompareByIdentity compare_46_63( 46, 63, 'A', biol::GetSSTypes().COIL);
      assemble::SSECompareByIdentity compare_73_76( 73, 76, 'A', biol::GetSSTypes().COIL);

      // now insert two SSEs into the model
      model.Insert( sp_helix_23_34);
      model.Insert( sp_strand_10_17);

      // test operator() with model with 2 SSEs
      BCL_MessageStd( "testing GetInvertedModel() with model with 2 SSEs");
      {
        // get inverted model
        assemble::ProteinModel inverted_model( *inverter.GetInvertedModel( model));

        // get the SSEs from the model
        util::SiPtrVector< const assemble::SSE> sses( inverted_model.GetSSEs());

        // check that there are three SSEs
        BCL_ExampleCheck( sses.GetSize(), 3);

        // check all three SSEs
        BCL_ExampleCheck( compare_1_9(   sses( 0)), true);
        BCL_ExampleCheck( compare_18_22( sses( 1)), true);
        BCL_ExampleCheck( compare_35_76( sses( 2)), true);

        // now check the inverted to make sure it cached the correct sections
        BCL_ExampleCheck( inverter.DoesContainCoil( native_chain.GetSequence(), 1, 9), true);
        BCL_ExampleCheck( inverter.DoesContainCoil( native_chain.GetSequence(), 18, 22), true);
        BCL_ExampleCheck( inverter.DoesContainCoil( native_chain.GetSequence(), 35, 76), true);
        BCL_ExampleCheck( inverter.DoesContainCoil( native_chain.GetSequence(), 8, 12), false);

        // write the pdb out
        Proteins::WriteModelToPDB( inverted_model, AddExampleOutputPathToFilename( inverter, "inverted_model_b.pdb"));
      }

      // insert the strand 1_7 into the model
      model.Insert( sp_strand_1_7);

      // test operator() with model with 3 SSEs
      BCL_MessageStd( "testing GetInvertedModel() with model with 3 SSEs");
      {
        // get inverted model
        assemble::ProteinModel inverted_model( *inverter.GetInvertedModel( model));

        // get the SSEs from the model
        util::SiPtrVector< const assemble::SSE> sses( inverted_model.GetSSEs());

        // check that there are three SSEs
        BCL_ExampleCheck( sses.GetSize(), 3);

        // check all three SSEs
        BCL_ExampleCheck( compare_8_9(   sses( 0)), true);
        BCL_ExampleCheck( compare_18_22( sses( 1)), true);
        BCL_ExampleCheck( compare_35_76( sses( 2)), true);

        // now check the inverted to make sure it cached the correct sections
        BCL_ExampleCheck( inverter.DoesContainCoil( native_chain.GetSequence(),  1, 9), true);
        BCL_ExampleCheck( inverter.DoesContainCoil( native_chain.GetSequence(), 18, 22), true);
        BCL_ExampleCheck( inverter.DoesContainCoil( native_chain.GetSequence(), 35, 76), true);
        BCL_ExampleCheck( inverter.DoesContainCoil( native_chain.GetSequence(),  8, 9), true);

        // write the pdb out
        Proteins::WriteModelToPDB( inverted_model, AddExampleOutputPathToFilename( inverter, "inverted_model_c.pdb"));
      }

      // insert the remaining two strands
      model.Insert( sp_strand_40_45);
      model.Insert( sp_strand_64_72);

      // test operator() with model with 3 SSEs
      BCL_MessageStd( "testing GetInvertedModel() with model with 5 SSEs");
      {
        // get inverted model
        assemble::ProteinModel inverted_model( *inverter.GetInvertedModel( model));

        // get the SSEs from the model
        util::SiPtrVector< const assemble::SSE> sses( inverted_model.GetSSEs());

        // check that there are three SSEs
        BCL_ExampleCheck( sses.GetSize(), 5);

        // check all three SSEs
        BCL_ExampleCheck( compare_8_9(   sses( 0)), true);
        BCL_ExampleCheck( compare_18_22( sses( 1)), true);
        BCL_ExampleCheck( compare_35_39( sses( 2)), true);
        BCL_ExampleCheck( compare_46_63( sses( 3)), true);
        BCL_ExampleCheck( compare_73_76( sses( 4)), true);

        // now check the inverted to make sure it cached the correct sections
        BCL_ExampleCheck( inverter.DoesContainCoil( native_chain.GetSequence(),  1, 9), true);
        BCL_ExampleCheck( inverter.DoesContainCoil( native_chain.GetSequence(), 18, 22), true);
        BCL_ExampleCheck( inverter.DoesContainCoil( native_chain.GetSequence(), 35, 76), true);
        BCL_ExampleCheck( inverter.DoesContainCoil( native_chain.GetSequence(),  8, 9), true);
        BCL_ExampleCheck( inverter.DoesContainCoil( native_chain.GetSequence(),  35, 39), true);
        BCL_ExampleCheck( inverter.DoesContainCoil( native_chain.GetSequence(),  46, 63), true);
        BCL_ExampleCheck( inverter.DoesContainCoil( native_chain.GetSequence(),  73, 76), true);

        // write the pdb out
        Proteins::WriteModelToPDB( inverted_model, AddExampleOutputPathToFilename( inverter, "inverted_model_d.pdb"));
      }

    //////////////////////
    // input and output //
    //////////////////////

      // write and read
      WriteBCLObject( inverter);
      assemble::ProteinModelInverter inverter_read;
      ReadBCLObject( inverter_read);
      BCL_ExampleCheck( inverter_read.GetUseCache(), inverter.GetUseCache());
      BCL_ExampleCheck( inverter_read.GetModelCacheSize(), inverter.GetModelCacheSize());
      BCL_ExampleCheck( inverter_read.GetCoilCacheSize(), inverter.GetCoilCacheSize());
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleProteinModelInverter

  const ExampleClass::EnumType ExampleAssembleProteinModelInverter::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleProteinModelInverter())
  );

} // namespace bcl
