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
#include "assemble/bcl_assemble_collector_sse_unpaired.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_collector_sse_unpaired.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleCollectorSSEUnpaired :
    public ExampleInterface
  {
  public:

    ExampleAssembleCollectorSSEUnpaired *Clone() const
    {
      return new ExampleAssembleCollectorSSEUnpaired( *this);
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

      // create SiPtrList "sse_list" and initialize with the contencts of "sse_vector"
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
      assemble::CollectorSSEUnpaired collector;

      // test SetContactType function, GetContactType, SetMaxDistance, and GetMaxDistance functions
      BCL_MessageStd
      (
        "test SetContactType function, GetContactType, SetMaxDistance, and GetMaxDistance functions"
      );
      // set "collector" to look for unpaired strand interactions
      collector.SetContactType( contact::GetTypes().STRAND_STRAND);
      // set "collector" to consider distances greater than 20A as unpaired
      collector.SetMaxDistance( 20.0);
      // check that the SetData function worked and simultaneously make sure GetData() function works
      BCL_Example_Check
      (
        collector.GetContactType() == contact::GetTypes().STRAND_STRAND
        && collector.GetMaxDistance() == 20.0,
        " set and get data functions are broken"
      );
      // test Collect function
      BCL_MessageStd( "test Collect function");
      util::SiPtrList< const assemble::SSE> sse_collection( collector.Collect( protein_model));
      // make sure that the correct number of SSEs are collected
      BCL_Example_Check
      (
        sse_collection.GetSize() == 1,
        "number of SSEs collected is " + util::Format()( sse_collection.GetSize()) + "but should be 1"
      );

      // test constructor from contact type and maximum distance for two SSEs to be considered paired
      BCL_MessageStd
      (
        "test constructor from contact type and maximum distance for two SSEs to be considered paired"
      );
      assemble::CollectorSSEUnpaired collector_b( contact::GetTypes().STRAND_STRAND, 0.0);
      // check that the constructor worked properly
      BCL_Example_Check
      (
        collector_b.GetContactType() == contact::GetTypes().STRAND_STRAND
        && math::EqualWithinTolerance( collector_b.GetMaxDistance(), 0.0),
        " constructor from contact type and maximum distance for two SSEs to be considered paired is broken"
      );
      // test Collect function of "collector_b"
      BCL_MessageStd( "test Collect function");
      util::SiPtrList< const assemble::SSE> sse_collection_b( collector_b.Collect( protein_model));
      // make sure that the correct number of SSEs are collected
      BCL_Example_Check
      (
        sse_collection_b.GetSize() == 4,
        "number of SSEs collected is " + util::Format()( sse_collection_b.GetSize()) + "but should be 4"
      );

      // test locating helix without helix-helix interaction
      BCL_MessageStd
      (
        "test locating helix without helix-helix interaction"
      );
      assemble::CollectorSSEUnpaired collector_c( contact::GetTypes().HELIX_HELIX, 20.0);
      // check that the constructor worked properly
      BCL_Example_Check
      (
        collector_c.GetContactType() == contact::GetTypes().HELIX_HELIX
        && math::EqualWithinTolerance( collector_c.GetMaxDistance(), 20.0),
        " test locating helix without helix-helix interaction did not construct properly"
      );
      // test Collect function of "collector_c"
      BCL_MessageStd( "test Collect function");
      util::SiPtrList< const assemble::SSE> sse_collection_c( collector_c.Collect( protein_model));
      // make sure that the correct number of SSEs are collected
      BCL_Example_Check
      (
        sse_collection_c.GetSize() == 1,
        "number of SSEs collected is " + util::Format()( sse_collection_c.GetSize()) + "but should be 1"
      );

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleCollectorSSEUnpaired

  const ExampleClass::EnumType ExampleAssembleCollectorSSEUnpaired::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleCollectorSSEUnpaired())
  );

} // namespace bcl
