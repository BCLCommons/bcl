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
#include "fold/bcl_fold_loop_parameters.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"

// external includes - sorted alphabetically

namespace bcl
{

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_loop_parameters.cpp
  //! @brief this example tests the implementation of fold::LoopParameters
  //!
  //! @author fischea
  //! @date Dec 16, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ExampleFoldLoopParameters :
     public ExampleInterface
  {

  //////////
  // data //
  //////////

  public:

    //! single instance of this class
    static const ExampleClass::EnumType s_Instance;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to a new ExampleFoldLoopParameters
    ExampleFoldLoopParameters *Clone() const
    {
      return new ExampleFoldLoopParameters( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief run routine
    //! this is performing the execution of the example
    int Run() const
    {
      // create the members of the loop parameters class
      const linal::Vector3D translation( 1.0, 2.0, 3.0);
      const linal::Vector3D euler_angles( 4.0, 5.0, 6.0);
      storage::Vector< int> anchors;
      anchors.PushBack( 5);
      anchors.PushBack( 8);
      const size_t sequence_distance( 2);
      const storage::Vector< double> dihedral_angles( storage::Vector< double>::Create( 7.0, 8.0, 9.0, 10.0, 11.0, 12.0));

      // prepare protein model input files for testing
      const std::string test_pdb( AddExampleInputPathToFilename( e_Biology, "2AP3_no_loops.pdb"));
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, test_pdb);
      const pdb::Handler pdb( read);
      io::File::CloseClearFStream( read);
      const assemble::ProteinModel model( pdb::Factory().ProteinModelFromPDB( pdb));
      const biol::AABase &anchor_1( *model.GetAminoAcids()( 50));
      const biol::AABase &anchor_2( *model.GetAminoAcids()( 59));
      const char chain_id( anchor_1.GetChainID());

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      fold::LoopParameters param_default;

      // construct from members
      fold::LoopParameters param( translation, euler_angles, anchors, chain_id, dihedral_angles);

        // check create function
        // the following checks will assume that GetLocalCoordinateSystem() and GetEulerAngles() are working
        // correctly. those two functions are checked independently below
        const util::ShPtr< fold::LoopParameters> sp_param( fold::LoopParameters::Create( anchor_1, anchor_2));
        BCL_ExampleCheck( sp_param->GetSequenceDistance(), anchor_2.GetSeqID() - anchor_1.GetSeqID() - 1);
        const linal::Vector3D &trans( sp_param->GetTranslation());
        const linal::Vector3D exp_trans( -3.37596, 4.37533, 15.5913);
        BCL_ExampleCheck( linal::EqualWithinTolerance( trans, exp_trans, 0.05), true);

    /////////////////
    // data access //
    /////////////////

      // check class name
      BCL_ExampleCheck( param.GetClassIdentifier(), ( GetStaticClassName< fold::LoopParameters>()));

      // check getter functions
      BCL_ExampleCheck( param.GetTranslation(), translation);
      BCL_ExampleCheck( param.GetRotation(), euler_angles);
      BCL_ExampleCheck( param.GetAnchors(), anchors);
      BCL_ExampleCheck( param.GetChainID(), anchor_1.GetChainID());
      BCL_ExampleCheck( param.GetSequenceDistance(), sequence_distance);
      BCL_ExampleCheck( param.GetSequenceDirection(), biol::AASequenceFlexibility::e_CTerminal);
      BCL_ExampleCheck( param.GetAngles(), dihedral_angles);

    ////////////////
    // operations //
    ////////////////

      // test calculation of local coordinate system
      const biol::AABase &aa( *model.GetAminoAcids()( 96));
      const storage::VectorND< 3, linal::Vector3D> coord_system
      (
        fold::LoopParameters::GetLocalCoordinateSystem( aa)
      );
      const storage::VectorND< 3, linal::Vector3D> correct_coord_system
      (
        linal::Vector3D( 0.5322, 0.5669, 0.6287),   // x-axis
        linal::Vector3D( -0.3535, 0.8237, -0.4434), // y-axis
        linal::Vector3D( 0.7654, -0.0137, -0.6388)  // z-axis
      );
      BCL_ExampleCheck
      (
        linal::EqualWithinTolerance( coord_system.First(), correct_coord_system.First(), 0.05), true
      );
      BCL_ExampleCheck
      (
        linal::EqualWithinTolerance( coord_system.Second(), correct_coord_system.Second(), 0.05), true
      );
      BCL_ExampleCheck
      (
        linal::EqualWithinTolerance( coord_system.Third(), correct_coord_system.Third(), 0.05), true
      );

      // create coordinate systems to test the calculation of the Euler angles
      const storage::VectorND< 3, linal::Vector3D> frame_a
      (
        linal::Vector3D( 1.0, 0.0, 0.0),  // x-axis
        linal::Vector3D( 0.0, 1.0, 0.0),  // y-axis
        linal::Vector3D( 0.0, 0.0, 1.0)   // z-axis
      );
      const storage::VectorND< 3, linal::Vector3D> frame_b
      (
        linal::Vector3D( 0.0, 1.0, 0.0),  // x-axis
        linal::Vector3D( -1.0, 0.0, 0.0), // y-axis
        linal::Vector3D( 0.0, 0.0, 1.0)   // z-axis
      );
      const storage::VectorND< 3, linal::Vector3D> frame_c
      (
        linal::Vector3D( 1.0, 0.0, 0.0),  // x-axis
        linal::Vector3D( 0.0, 0.0, 1.0),  // y-axis
        linal::Vector3D( 0.0, -1.0, 0.0) // z-axis
      );
      const storage::VectorND< 3, linal::Vector3D> frame_d
      (
        linal::Vector3D( 0.0, 1.0, 0.0),  // x-axis
        linal::Vector3D( 0.0, 0.0, 1.0),  // y-axis
        linal::Vector3D( 1.0, 0.0, 0.0)   // z-axis
      );

      // check the computed Euler angles
      const linal::Vector3D euler_a_b( fold::LoopParameters::GetEulerAngles( frame_a, frame_b));
      BCL_ExampleCheckWithinAbsTolerance( euler_a_b( 0), 0.0, 0.05);
      BCL_ExampleCheckWithinAbsTolerance( euler_a_b( 1), 0.0, 0.05);
      BCL_ExampleCheckWithinAbsTolerance( euler_a_b( 2), 1.5 * math::g_Pi, 0.05);
      const linal::Vector3D euler_a_c( fold::LoopParameters::GetEulerAngles( frame_a, frame_c));
      BCL_ExampleCheckWithinAbsTolerance( euler_a_c( 0), 1.5 * math::g_Pi, 0.05);
      BCL_ExampleCheckWithinAbsTolerance( euler_a_c( 1), 0.0, 0.05);
      BCL_ExampleCheckWithinAbsTolerance( euler_a_c( 2), 0.0, 0.05);
      const linal::Vector3D euler_a_d( fold::LoopParameters::GetEulerAngles( frame_a, frame_d));
      BCL_ExampleCheckWithinAbsTolerance( euler_a_d( 0), 1.5 * math::g_Pi, 0.05);
      BCL_ExampleCheckWithinAbsTolerance( euler_a_d( 1), 0.0, 0.05);
      BCL_ExampleCheckWithinAbsTolerance( euler_a_d( 2), 1.5 * math::g_Pi, 0.05);

    //////////////////////
    // input and output //
    //////////////////////

      // read object
      fold::LoopParameters param_read;
      param_read.AssertRead( param.GetString());

      // compare read in object to written object
      BCL_ExampleCheck( param.GetTranslation(), param_read.GetTranslation());
      BCL_ExampleCheck( param.GetRotation(), param_read.GetRotation());
      BCL_ExampleCheck( param.GetAnchors(), param_read.GetAnchors());
      BCL_ExampleCheck( param.GetChainID(), param_read.GetChainID());
      BCL_ExampleCheck( param.GetSequenceDistance(), param_read.GetSequenceDistance());
      BCL_ExampleCheck( param.GetSequenceDirection(), param_read.GetSequenceDirection());
      BCL_ExampleCheck( param.GetAngles(), param_read.GetAngles());

      return 0;
    }

  }; // class ExampleFoldLoopParameters

  //! single instance of this class
  const ExampleClass::EnumType ExampleFoldLoopParameters::s_Instance
  (
     GetExamples().AddEnum( ExampleFoldLoopParameters())
  );

} // namespace bcl
