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
#include "assemble/bcl_assemble_collector_sse_paired.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_geometry_packing.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_collector_sse_paired.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleCollectorSSEPaired :
    public ExampleInterface
  {
  public:

    ExampleAssembleCollectorSSEPaired *Clone() const
    {
      return new ExampleAssembleCollectorSSEPaired( *this);
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

      // initialize contact types set
      storage::Set< contact::Type> contact_types_a;
      contact_types_a.InsertElement( contact::GetTypes().HELIX_HELIX);
      contact_types_a.InsertElement( contact::GetTypes().HELIX_SHEET);
      contact_types_a.InsertElement( contact::GetTypes().SHEET_HELIX);

      storage::Set< contact::Type> contact_types_b;
      contact_types_b.InsertElement( contact::GetTypes().HELIX_SHEET);
      contact_types_b.InsertElement( contact::GetTypes().SHEET_HELIX);

      storage::Set< contact::Type> contact_types_c;
      contact_types_c.InsertElement( contact::GetTypes().STRAND_STRAND);
      contact_types_c.InsertElement( contact::GetTypes().SHEET_SHEET);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      BCL_MessageStd( "Testing default constructor");

      // call default constructor
      assemble::CollectorSSEPaired collector_a;

      // test constructor
      BCL_MessageStd( "Testing constructor");

      // call constructor
      assemble::CollectorSSEPaired collector_b
      (
        assemble::GetSSEGeometryPackingPickers().e_BestInteractionWeight, contact_types_a, 10.0, true
      );

      // test copy constructor
      BCL_MessageStd( "Testing copy constructor");

      // call copy constructor
      assemble::CollectorSSEPaired collector_c( collector_b);

      // check the copy constructor worked correctly
      BCL_Example_Check
      (
        collector_c.GetContactTypes().InternalData() == collector_b.GetContactTypes().InternalData() &&
        collector_c.GetMaxDistance() == collector_b.GetMaxDistance() &&
        collector_c.GetOrthogonalConnection() == collector_b.GetOrthogonalConnection(),
        "The copy constructor failed!"
      );

    /////////////////
    // data access //
    /////////////////

      // test GetContactTypes()
      BCL_MessageStd( "Testing GetContactTypes()");

      // check the contact types are set correctly initially
      BCL_Example_Check
      (
        collector_b.GetContactTypes().InternalData() == contact_types_a.InternalData(),
        "The collector was unable to set contact_types correctly" + util::Format()( collector_b.GetContactTypes())
      );

      // test SetContactTypes()
      BCL_MessageStd( "Testing SetContactTypes()");

      // set the contact types
      collector_b.SetContactTypes( contact_types_b);

      // check the contact types are set correctly after change
      BCL_Example_Check
      (
        collector_b.GetContactTypes().InternalData() == contact_types_b.InternalData(),
        "The collector was unable to set contact_types correctly" + util::Format()( collector_b.GetContactTypes())
      );

      // test GetMaxDistance()
      BCL_MessageStd( "Testing GetMaxDistance()");

      // check the contact types are set correctly initially
      BCL_Example_Check
      (
        collector_b.GetMaxDistance() == 10.0,
        "The collector was unable to set MaxDistance correctly" + util::Format()( collector_b.GetMaxDistance())
      );

      // test SetMaxDistance()
      BCL_MessageStd( "Testing SetMaxDistance()");

      // set the MaxDistance
      collector_b.SetMaxDistance( 25.0);

      // check the MaxDistance is set correctly after change
      BCL_Example_Check
      (
        collector_b.GetMaxDistance() == 25.0,
        "The collector was unable to set MaxDistance correctly" + util::Format()( collector_b.GetMaxDistance())
      );

      // test GetOrthogonalConnection()
      BCL_MessageStd( "Testing GetOrthogonalConnection()");

      // check the OrthogonalConnection is set correctly initially
      BCL_Example_Check
      (
        collector_b.GetOrthogonalConnection(),
        "The collector was unable to set OrthogonalConnection correctly" +
          util::Format()( collector_b.GetOrthogonalConnection())
      );

      // test SetOrthogonalConnection()
      BCL_MessageStd( "Testing SetOrthogonalConnection()");

      // set OrthogonalConnection
      collector_b.SetOrthogonalConnection( false);

      // check the OrthogonalConnection is set correctly after change
      BCL_Example_Check
      (
        !collector_b.GetOrthogonalConnection(),
        "The collector was unable to set OrthogonalConnection correctly" +
          util::Format()( collector_b.GetOrthogonalConnection())
      );

    ////////////////
    // operations //
    ////////////////

      // create ProteinModel "protein_model" from pdb file
      BCL_MessageStd( "building models from pdb chains and sse information");
      assemble::ProteinModel protein_model( Proteins::GetModel( pdb_filename));

      // create SiPtrVector "sse_vector" and initialize with all strands of "protein_model"
      util::SiPtrVector< const assemble::SSE> sse_vector( protein_model.GetSSEs());

      // iterate over sses in the protein model
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator sse_itr_a( sse_vector.Begin()),
          sse_itr_end( sse_vector.End());
        sse_itr_a != sse_itr_end; ++sse_itr_a
      )
      {
        // iterate over every other sse
        for
        (
          util::SiPtrVector< const assemble::SSE>::const_iterator sse_itr_b( sse_itr_a + 1);
          sse_itr_b != sse_itr_end; ++sse_itr_b
        )
        {
          // get the packing for this pair
          const assemble::SSEGeometryPacking this_packing
          (
            **sse_itr_a, **sse_itr_b, assemble::GetSSEGeometryPackingPickers().e_BestInteractionWeight
          );

          // print out information
          BCL_MessageVrb
          (
            ( *sse_itr_a)->GetIdentification() + " vs " +
            ( *sse_itr_b)->GetIdentification() + " ==> " +
            util::Format()( this_packing.GetContactType().GetName()) + " with distance " +
            util::Format()( this_packing.GetShortestConnection().GetLength()) + " orthogonal: " +
            util::Format()( this_packing.GetOrthogonalConnection())
          );
        }
      }

      // create SiPtrList "sse_list" and initialize with the contencts of "sse_vector"
      util::SiPtrList< const assemble::SSE> sse_list( sse_vector.Begin(), sse_vector.End());

      // test Collect function
      BCL_MessageStd( "test Collect function for Helix-Strand, Strand-Helix");
      storage::List< storage::VectorND< 2, util::SiPtr< const assemble::SSE> > > sse_pairs
      (
        collector_b.Collect( protein_model)
      );

      // make sure that the correct number of SSEs are collected
      BCL_Example_Check
      (
        sse_pairs.GetSize() == 2,
        "number of SSEs collected is " + util::Format()( sse_pairs.GetSize()) + " but should be 2"
      );

      // test for strand-strand sheet-sheet
      BCL_MessageStd( "test Collect function for Strand-Strand, Sheet-Sheet");

      // create a collector with strand strand and sheet
      assemble::CollectorSSEPaired collector_d
      (
        assemble::GetSSEGeometryPackingPickers().e_BestInteractionWeight, contact_types_c, 15.0, false
      );

      storage::List< storage::VectorND< 2, util::SiPtr< const assemble::SSE> > > sse_pairs_b
      (
        collector_d.Collect( protein_model)
      );

      // make sure that the correct number of SSEs are collected
      BCL_Example_Check
      (
        sse_pairs_b.GetSize() == 6,
        "number of SSEs pairs collected is " + util::Format()( sse_pairs_b.GetSize()) + " but should be 6"
      );

    ///////////////
    // operators //
    ///////////////

      // test assignment operator
      BCL_MessageStd( "Testing assignment operator");

      // initialize default collector
      assemble::CollectorSSEPaired collector_assigned;

      // assign it collector_b
      collector_assigned = collector_b;

      // check that assignment was done correctly
      BCL_Example_Check
      (
        collector_assigned.GetContactTypes().InternalData() == collector_b.GetContactTypes().InternalData() &&
        collector_assigned.GetMaxDistance() == collector_b.GetMaxDistance() &&
        collector_assigned.GetOrthogonalConnection() == collector_b.GetOrthogonalConnection(),
        "The assignment operator failed!"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // test Write
      BCL_MessageStd( "Testing the write function");
      WriteBCLObject( collector_b);
      // read the file back
      BCL_MessageStd( "Testing the read function");
      assemble::CollectorSSEPaired collector_read;
      ReadBCLObject( collector_read);

      // check both collectors are same
      BCL_MessageStd( "Checking the written and read objects are the same");
      BCL_Example_Check
      (
        collector_read.GetContactTypes().InternalData() == collector_b.GetContactTypes().InternalData() &&
        collector_read.GetMaxDistance() == collector_b.GetMaxDistance() &&
        collector_read.GetOrthogonalConnection() == collector_b.GetOrthogonalConnection(),
        "The written and read objects are different " + util::Format()( collector_b) +
          "\nvs\n" + util::Format()( collector_read)
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleCollectorSSEPaired

  const ExampleClass::EnumType ExampleAssembleCollectorSSEPaired::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleCollectorSSEPaired())
  );

} // namespace bcl
