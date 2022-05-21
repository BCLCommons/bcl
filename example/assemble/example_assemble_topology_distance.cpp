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
#include "assemble/bcl_assemble_topology_distance.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_topology_distance.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleTopologyDistance :
    public ExampleInterface
  {
  public:

    ExampleAssembleTopologyDistance *Clone() const
    {
      return new ExampleAssembleTopologyDistance( *this);
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
      BCL_MessageStd( "Reading native structure and perturbed models");

      // initialize stream and pdb factory and minsse_sizes
      storage::Map< biol::SSType, size_t> sse_min_sizes;
      sse_min_sizes[ biol::GetSSTypes().HELIX] = 9;
      sse_min_sizes[ biol::GetSSTypes().STRAND] = 5;
      sse_min_sizes[ biol::GetSSTypes().COIL] = 999;

      // initialize native structure pdb filename
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
      //build models from pdb
      BCL_MessageStd( "building model from pdb file");
      assemble::ProteinModel native_structure
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, sse_min_sizes)
      );

      // initialize perturbed model filename
      const std::string perturbed_model_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb.perturbed"));
      //build models from pdb
      BCL_MessageStd( "building model from pdb file");
      assemble::ProteinModel perturbed_model
      (
        Proteins::GetModel( perturbed_model_filename, biol::GetAAClasses().e_AABackBone, sse_min_sizes)
      );
      // initialize default distance and angle cutoffs
      const double default_distance_cutoff( 5.0);
      const double default_angle_cutoff( math::g_Pi / 2.0);

      // initialize distance and angle cutoffs
      const double distance_cutoff( 8.0);
      const double angle_cutoff( math::g_Pi / 3.0);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      assemble::TopologyDistance distance_default;

      // test construction
      BCL_MessageStd( "test construct from distance and angle cutoff");
      assemble::TopologyDistance distance( distance_cutoff, angle_cutoff);

      // test copy constructor
      assemble::TopologyDistance distance_copy( distance);
      BCL_MessageStd( "test copy constructor");
      BCL_Example_Check
      (
        distance_copy.GetDistanceCutoff() == distance.GetDistanceCutoff() &&
        distance_copy.GetAngleCutoff() == distance.GetAngleCutoff(),
        "copy should be " + util::Format()( distance) + " but instead is " + util::Format()( distance_copy)
      );

      // test clone constructor
      util::ShPtr< assemble::TopologyDistance> distance_clone( distance.Clone());
      BCL_MessageStd( "test clone constructor");
      BCL_Example_Check
      (
        distance_clone->GetDistanceCutoff() == distance.GetDistanceCutoff() &&
        distance_clone->GetAngleCutoff() == distance.GetAngleCutoff(),
        "clone  should be " + util::Format()( distance) + " but instead is " + util::Format()( *distance_clone)
      );

    /////////////////
    // data access //
    /////////////////

      // check GetStaticClassName
      BCL_MessageStd( "test GetStatisClassName");
      const std::string correct_static_class_name( "bcl::assemble::TopologyDistance");
      BCL_Example_Check
      (
        GetStaticClassName< assemble::TopologyDistance>() == correct_static_class_name,
        "GetStaticClassName gives " + GetStaticClassName< assemble::TopologyDistance>() + " but should give " +
        correct_static_class_name
      );

      // check GetClassIdentifier
      BCL_MessageStd( "test GetClassIdentifier");
      BCL_Example_Check
      (
        GetStaticClassName< assemble::TopologyDistance>() == distance.GetClassIdentifier(),
        "GetClassIdentifier gives " + distance.GetClassIdentifier() + " but should give " +
        correct_static_class_name
      );

      // check GetDefaultDistanceCutoff()
      BCL_MessageStd( "test GetDefaultDistanceCutoff");
      BCL_Example_Check
      (
        assemble::TopologyDistance::GetDefaultDistanceCutoff() == default_distance_cutoff,
        "GetDefaultDistanceCutoff gives " + util::Format()( assemble::TopologyDistance::GetDefaultDistanceCutoff()) +
        " instead of " + util::Format()( default_distance_cutoff)
      );

      // check GetDefaultAngleCutoff()
      BCL_MessageStd( "test GetDefaultAngleCutoff");
      BCL_Example_Check
      (
        assemble::TopologyDistance::GetDefaultAngleCutoff() == default_angle_cutoff,
        "GetDefaultAngleCutoff gives " + util::Format()( distance.GetDefaultAngleCutoff()) +
        " instead of " + util::Format()( default_angle_cutoff)
      );

      // check GetDistanceCutoff()
      BCL_MessageStd( "test GetDistanceCutoff");
      BCL_Example_Check
      (
        distance.GetDistanceCutoff() == distance_cutoff,
        "GetDistanceCutoff gives " + util::Format()( distance.GetDistanceCutoff()) +
        " instead of " + util::Format()( distance_cutoff)
      );

      // check GetAngleCutoff()
      BCL_MessageStd( "test GetAngleCutoff");
      BCL_Example_Check
      (
        distance.GetAngleCutoff() == angle_cutoff,
        "GetAngleCutoff gives " + util::Format()( distance.GetAngleCutoff()) +
        " instead of " + util::Format()( angle_cutoff)
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////////////////////////////////
    // test operator() with ProteinModel pair //
    ////////////////////////////////////////////
      BCL_MessageStd( "test operator() with ProteinModel pair");

      // compare native structure to itself
      BCL_MessageStd( "comparing the native structure with itself");
      const double distance_to_itself( distance_default( native_structure, native_structure));
      // check the distance
      BCL_Example_Check
      (
        distance_to_itself == double( 0.0),
        "The topology distance of native structure to itself should be 0 not " + util::Format()( distance_to_itself)
      );

      // compare native structure to empty model
      BCL_MessageStd( "comparing the native structure with an empty model");
      // initialize empty model
      const assemble::ProteinModel empty_model;
      const double distance_to_empty( distance_default( empty_model, native_structure));
      // check the distance
      BCL_Example_Check

      (
        distance_to_empty == double( 100.0),
        "The topology distance of native structure to empty model should be 10 not " + util::Format()( distance_to_empty)
      );

      // compare native structure to perturbed models
      BCL_MessageStd( "calculating distance to native structure");

      // calculate the distance to native_structure
      const double model_distance_expected( 50.0);

      // calculate the distance to native_structure
      const double model_distance( distance_default( perturbed_model, native_structure));

      // print out and check the calculated distance
      BCL_MessageStd
      (
        "The topology distance to native for model " + util::Format()( model_distance)
      );
      BCL_Example_Check
      (
        model_distance == model_distance_expected,
        "The topology distance of perturbed model to native structure should be: " +
        util::Format()( model_distance_expected) + " but instead it is: " + util::Format()( model_distance)
      );

    /////////////////////////////////////
    // test operator() with Chain pair //
    /////////////////////////////////////
      BCL_MessageStd( "test operator() with Chain pair");

      // calculate the distance to native_structure
      const double chain_distance( distance_default( perturbed_model.GetChain( 'A'), native_structure.GetChain( 'A')));

      // print out and check the calculated distance
      BCL_MessageStd
      (
        "The topology distance to native for model for chain A " + util::Format()( chain_distance)
      );
      BCL_Example_Check
      (
        chain_distance == model_distance_expected,
        "The topology distance of perturbed model to native structure for Chain A: " +
        util::Format()( model_distance_expected) + " but instead it is: " + util::Format()( chain_distance)
      );

    ///////////////////////////////////
    // test operator() with SSE pair //
    ///////////////////////////////////
      BCL_MessageStd( "test operator() with SSE pair");

      // calculate the distance to native_structure
      const double sse_distance_expected( 1.0);

      // create reference to sses
      assemble::LocatorSSE sse_locator( 'A', 23, 34);
      util::SiPtr< const assemble::SSE> native_sse( sse_locator.Locate( native_structure));
      util::SiPtr< const assemble::SSE> perturbed_sse( sse_locator.Locate( perturbed_model));

      // make sure SSEs were located correctly
      BCL_ExampleIndirectAssert
      (
        native_sse.IsDefined(), true,
        "Helix 23-34  was not located in the native structure"
      );
      BCL_ExampleIndirectAssert
      (
        perturbed_sse.IsDefined(), true,
        "Helix 23-34  was not located in the perturbed model"
      );

      // calculate the distance to native_structure
      const double sse_distance( distance_default( *perturbed_sse, *native_sse));

      // print out and check the calculated distance
      BCL_MessageStd
      (
        "The topology distance to native for SSE " + native_sse->GetIdentification() + " : " +
        util::Format()( sse_distance)
      );
      BCL_Example_Check
      (
        sse_distance == sse_distance_expected,
        "The topology distance of perturbed model to native structure for SSE " + native_sse->GetIdentification()
        + " should be: " + util::Format()( sse_distance_expected) + " but instead it is: " + util::Format()( sse_distance)
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // initialize output stream
      BCL_MessageStd( "testing read and write functionalities");
      // write object
      WriteBCLObject( distance);
      // read object
      assemble::TopologyDistance distance_read;
      ReadBCLObject( distance_read);

      // compare the objects
      BCL_MessageStd( "compare written and read objects");
      BCL_Example_Check
      (
        distance_read.GetDistanceCutoff() == distance.GetDistanceCutoff() &&
        math::EqualWithinTolerance( distance_read.GetAngleCutoff(), distance.GetAngleCutoff()),
        "the written and read TopologyDistance classes differ from each other" +
        util::Format()( distance_read) + "\nvs\n" + util::Format()( distance)
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleTopologyDistance

  const ExampleClass::EnumType ExampleAssembleTopologyDistance::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleTopologyDistance())
  );

} // namespace bcl

