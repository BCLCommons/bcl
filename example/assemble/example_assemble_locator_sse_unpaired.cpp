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
#include "assemble/bcl_assemble_locator_sse_unpaired.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_pick_sse_furthest_euclidean.h"
#include "assemble/bcl_assemble_protein_model.h"

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_locator_sse_unpaired.cpp
  //!
  //! @author alexanns
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleLocatorSSEUnpaired :
    public ExampleInterface
  {
  public:

    ExampleAssembleLocatorSSEUnpaired *Clone() const
    {
      return new ExampleAssembleLocatorSSEUnpaired( *this);
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
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi_unpaired_strand.pdb"));

      // create ProteinModel "protein_model" from "pdb"
      BCL_MessageStd( "building models from pdb chains and sse information");

      assemble::ProteinModel protein_model( Proteins::GetModel( pdb_filename));

      // create SiPtrVector "sse_vector" and initialize with all strands of "protein_model"
      util::SiPtrVector< const assemble::SSE> sse_vector( protein_model.GetSSEs());

      //util::SiPtr< const assemble::SSE> switch_a( sse_vector( 1));
      //util::SiPtr< const assemble::SSE> switch_b( sse_vector( 3));
      //sse_vector( 1) = switch_b;
      //sse_vector( 3) = switch_a;

      // create SiPtrList "sse_list" and initialize with the contects of "sse_vector"
      util::SiPtrList< const assemble::SSE> sse_list( sse_vector.Begin(), sse_vector.End());

      // iterate through "sse_vector"
      for
      (
        storage::Vector< util::SiPtr< const assemble::SSE> >::const_iterator
          sse_itr( sse_vector.Begin()), sse_itr_end( sse_vector.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // print the distance the SSE currently denoted by "sse_itr" is away from the center of "protein_model"
        BCL_MessageStd
        (
          "Current SSE is " + util::Format()( linal::Distance( ( *sse_itr)->GetCenter(), protein_model.GetCenter()))
          + " away from the center of the protein and begins with seq id "
          + util::Format()( ( *( *sse_itr)->Begin())->GetSeqID())
        );
      }

      // test default constructor
      assemble::LocatorSSEUnpaired locator;

      // test SetCollector, GetCollector, and SetPick functions
      BCL_MessageStd
      (
        "test SetCollector, GetCollector, and SetPick functions"
      );
      // create CollectorUnpairedSSE from Collector of "locator"
      assemble::CollectorSSEUnpaired collector( locator.GetCollector());
      // set "collector" to look for unpaired strand interactions
      collector.SetContactType( contact::GetTypes().STRAND_STRAND);
      // set "collector" to consider distances greater than 20A as unpaired
      collector.SetMaxDistance( 20.0);
      // set the "m_Collector" of "locator" to "collector"
      locator.SetCollector( collector);
      // set "locator" to pick the furthest from a given point in 3D space
      locator.SetPick
      (
        util::ShPtr
        <
          find::PickCriteriaInterface
          <
            util::SiPtr< const assemble::SSE>,
            util::SiPtrList< const assemble::SSE>,
            linal::Vector3D
          >
        >( new assemble::PickSSEFurthestEuclidean())
      );
      // check that the SetData function worked and simultaneously make sure GetData() function works
      BCL_Example_Check
      (
        locator.GetCollector().GetContactType() == contact::GetTypes().STRAND_STRAND
        && locator.GetCollector().GetMaxDistance() == 20.0,
        " set and get data functions are broken"
      );
      // test Locate function
      BCL_MessageStd( "test Locate function");
      util::SiPtr< const assemble::SSE> sp_sse( locator.Locate( protein_model));

      BCL_ExampleAssert( sp_sse.IsDefined(), true);

      // print the distance the located SSE is from the center of the protein
      BCL_MessageStd
      (
        "Furthest unpaired SSE is " +
        util::Format()( linal::Distance( sp_sse->GetCenter(), protein_model.GetCenter()))
        + " away from the center of the protein"
      );
      // make sure that the located SSE is in fact the most distant SSE from the center of the protein
      double correct_distance( 8.03062);
      int correct_first_seqid( 40);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( correct_distance, linal::Distance( sp_sse->GetCenter(), protein_model.GetCenter()))
        && ( *sp_sse->Begin())->GetSeqID() == correct_first_seqid,
        "Furthest unpaired strand located is "
        + util::Format()( linal::Distance( sp_sse->GetCenter(), protein_model.GetCenter()))
        + " away from the center of the protein, but should be "
        + util::Format()( correct_distance)
        + " the beginning SeqID should be "
        + util::Format()( correct_first_seqid)
        + " but is: "
        + util::Format()( ( *sp_sse->Begin())->GetSeqID())
      );

      // test constructor from contact type and maximum distance for two SSEs to be considered paired
      BCL_MessageStd
      (
        "test constructor from contact type and maximum distance for two SSEs to be considered paired"
      );
      assemble::LocatorSSEUnpaired locator_b
      (
        assemble::PickSSEFurthestEuclidean(),
        contact::GetTypes().STRAND_STRAND,
        0.0
      );
      // check that the constructor worked properly
      BCL_Example_Check
      (
        locator_b.GetCollector().GetContactType() == contact::GetTypes().STRAND_STRAND
        && math::EqualWithinTolerance( locator_b.GetCollector().GetMaxDistance(), 0.0),
        " set and get data functions are broken"
      );
      // test Locate function of "locator_b"
      BCL_MessageStd( "test Locate function");
      util::SiPtr< const assemble::SSE> sp_sse_b( locator_b.Locate( protein_model));
      // print the distance the located SSE is from the center of the protein
      BCL_MessageStd
      (
        "Furthest unpaired SSE is " +
        util::Format()( linal::Distance( sp_sse_b->GetCenter(), protein_model.GetCenter()))
        + " away from the center of the protein"
      );
      // make sure that the located SSE is in fact the most distant SSE from the center of the protein
      correct_distance = 8.22513;
      correct_first_seqid = 10;
      BCL_Example_Check
      (
        math::EqualWithinTolerance( correct_distance, linal::Distance( sp_sse_b->GetCenter(), protein_model.GetCenter()))
        && ( *( sp_sse_b->Begin()))->GetSeqID() == correct_first_seqid,
        "Furthest unpaired strand located is "
        + util::Format()( linal::Distance( sp_sse_b->GetCenter(), protein_model.GetCenter()))
        + " away from the center of the protein with SeqID " + util::Format()( ( *( sp_sse_b->Begin()))->GetSeqID())
        + " , but should be "
        + util::Format()( correct_distance)
        + " the beginning SeqID should be "
        + util::Format()( correct_first_seqid)
      );

      // test locating helix without helix-helix interaction
      BCL_MessageStd
      (
        "test locating helix without helix-helix interaction"
      );
      // create LocatorSSEUnpaired "locator_c" constructed as a copy of "locator"
      assemble::LocatorSSEUnpaired locator_c( locator);
      // create CollectorSSEUnpaired "collector_c" and initialize with HELIX-HELIX and 20.0
      assemble::CollectorSSEUnpaired collector_c( contact::GetTypes().HELIX_HELIX, 20.0);
      // set "locator_c" to have "collector_c"
      locator_c.SetCollector( collector_c);
      // check that the constructor worked properly
      BCL_Example_Check
      (
        locator_c.GetCollector().GetContactType() == contact::GetTypes().HELIX_HELIX
        && math::EqualWithinTolerance( locator_c.GetCollector().GetMaxDistance(), 20.0),
        " test locating helix without helix-helix interaction did not construct properly"
      );
      // test Locate function of "locator_c"
      BCL_MessageStd( "test Locate function");
      util::SiPtr< const assemble::SSE> sp_sse_c( locator_c.Locate( protein_model));
      // print the distance the located SSE is from the center of the protein
      BCL_MessageStd
      (
        "Furthest unpaired SSE is " +
        util::Format()( linal::Distance( sp_sse_c->GetCenter(), protein_model.GetCenter()))
        + " away from the center of the protein"
      );
      // make sure that the located SSE is in fact the most distant SSE from the center of the protein
      correct_distance = 7.42538;
      correct_first_seqid = 23;
      BCL_Example_Check
      (
        math::EqualWithinTolerance( correct_distance, linal::Distance( sp_sse_c->GetCenter(), protein_model.GetCenter()))
        && ( *( sp_sse_c->Begin()))->GetSeqID() == correct_first_seqid,
        "Furthest unpaired strand located is "
        + util::Format()( linal::Distance( sp_sse_c->GetCenter(), protein_model.GetCenter()))
        + " away from the center of the protein with SeqID " + util::Format()( ( *( sp_sse_c->Begin()))->GetSeqID())
        + " , but should be "
        + util::Format()( correct_distance)
        + " the beginning SeqID should be "
        + util::Format()( correct_first_seqid)
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleLocatorSSEUnpairedStrand

  const ExampleClass::EnumType ExampleAssembleLocatorSSEUnpaired::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleLocatorSSEUnpaired())
  );

} // namespace bcl
