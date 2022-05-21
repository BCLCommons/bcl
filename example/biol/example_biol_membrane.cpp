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
#include "example_proteins.h"
// include the header of the class which this example is for
#include "biol/bcl_biol_membrane.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_protein_model_multiplier.h"
#include "io/bcl_io_file.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_membrane.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolMembrane :
    public ExampleInterface
  {
  public:

    ExampleBiolMembrane *Clone() const
    { return new ExampleBiolMembrane( *this);}

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

      // default constructor
      biol::Membrane membrane_default;

      storage::Vector< double> thicknesses;
      thicknesses.PushBack(   10.0);
      thicknesses.PushBack(    2.5);
      thicknesses.PushBack(   10.0);
      thicknesses.PushBack(    2.5);
      thicknesses.PushBack( 1000.0);
      thicknesses.PushBack(    0.0);
      thicknesses.PushBack(    0.0);
      thicknesses.PushBack(    0.0);
      thicknesses.PushBack(    0.0);

      // normal
      const linal::Vector3D normal( 0.0, 0.0, 1.0);
      // construct from thicknesses and normal used by the db connector
      biol::Membrane membrane_t_n( thicknesses, normal);

      // construct from thickness core, transition, gap
      const double t_core( 11.0);
      const double t_trans( 10.5);
      const double t_gap( 2.4);
      biol::Membrane membrane_ctg_n( t_core, t_trans, t_gap, normal);

      // copy
      biol::Membrane membrane_copy( membrane_t_n);

      // clone
      util::ShPtr< util::ObjectInterface> ptr( membrane_t_n.Clone());

    /////////////////
    // data access //
    /////////////////

      // check the class identifier
      BCL_ExampleCheck( GetStaticClassName( membrane_default), "bcl::biol::Membrane");

      // class identifier
      BCL_ExampleCheck( ptr->GetClassIdentifier(), GetStaticClassName< biol::Membrane>());

      // membrane normal
      BCL_ExampleCheck( membrane_t_n.GetNormal(), normal);

      // thicknesses
      BCL_ExampleCheck( membrane_t_n.GetThicknesses(), thicknesses);

      // limits
      storage::Vector< double> limits( thicknesses);
      for( size_t i( 1); i < thicknesses.GetSize(); ++i)
      {
        limits( i) = limits( i - 1) + thicknesses( i);
      }
      BCL_ExampleCheck( membrane_t_n.GetLimits(), limits);

      // thickness
      BCL_ExampleIndirectCheck
      (
        membrane_ctg_n.GetThickness( biol::GetEnvironmentTypes().e_MembraneCore) == t_core &&
        membrane_ctg_n.GetThickness( biol::GetEnvironmentTypes().e_GapCoreTransition) == t_gap &&
        membrane_ctg_n.GetThickness( biol::GetEnvironmentTypes().e_Transition) == t_trans &&
        membrane_ctg_n.GetThickness( biol::GetEnvironmentTypes().e_GapTransitionSolution) == t_gap,
        true,
        "thicknesses"
      );

      // limit
      BCL_ExampleIndirectCheck
      (
        membrane_ctg_n.GetLimit( biol::GetEnvironmentTypes().e_MembraneCore) == t_core &&
        membrane_ctg_n.GetLimit( biol::GetEnvironmentTypes().e_GapCoreTransition) == t_core + t_gap &&
        membrane_ctg_n.GetLimit( biol::GetEnvironmentTypes().e_Transition) == t_core + t_gap + t_trans &&
        membrane_ctg_n.GetLimit( biol::GetEnvironmentTypes().e_GapTransitionSolution) == t_core + t_gap + t_trans + t_gap,
        true,
        "limits"
      );

      // IsDefined
      BCL_ExampleCheck( membrane_ctg_n.IsDefined(), true);
      const biol::Membrane &membrane_undefined( biol::Membrane::GetUndefinedMembrane());
      BCL_ExampleCheck( membrane_undefined.IsDefined(), false);

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // determine environment type from zcoordinate
      linal::Vector3D coords_a( linal::Vector3D( 0.0, 0.0, 5.5));
      const linal::Vector3D coords_b( linal::Vector3D( 0.0, 0.0, 11.0));
      const linal::Vector3D coords_c( linal::Vector3D( 0.0, 0.0, 12.2));
      const linal::Vector3D coords_d( linal::Vector3D( 0.0, 0.0, 18.65));
      const linal::Vector3D coords_e( linal::Vector3D( 0.0, 0.0, 25.1));
      const linal::Vector3D coords_f( linal::Vector3D( 0.0, 0.0, 27.0));
      const biol::EnvironmentType env1( membrane_ctg_n.DetermineEnvironmentType( coords_a));
      const biol::EnvironmentType env2( membrane_ctg_n.DetermineEnvironmentType( coords_b));
      const biol::EnvironmentType env3( membrane_ctg_n.DetermineEnvironmentType( coords_c));
      const biol::EnvironmentType env4( membrane_ctg_n.DetermineEnvironmentType( coords_d));
      const biol::EnvironmentType env5( membrane_ctg_n.DetermineEnvironmentType( coords_e));
      const biol::EnvironmentType env6( membrane_ctg_n.DetermineEnvironmentType( coords_f));
      BCL_ExampleIndirectCheck
      (
        env1 == biol::GetEnvironmentTypes().e_MembraneCore &&
        env2 == biol::GetEnvironmentTypes().e_GapCoreTransition &&
        env3 == biol::GetEnvironmentTypes().e_GapCoreTransition &&
        env4 == biol::GetEnvironmentTypes().e_Transition &&
        env5 == biol::GetEnvironmentTypes().e_GapTransitionSolution &&
        env6 == biol::GetEnvironmentTypes().e_Solution,
        true,
        "environment types for different z coordinates"
      );

      // determine environment and weight
      const storage::Pair< biol::EnvironmentType, double> env_w1( membrane_ctg_n.DetermineEnvironmentTypeAndWeight( coords_a));
      const storage::Pair< biol::EnvironmentType, double> env_w2( membrane_ctg_n.DetermineEnvironmentTypeAndWeight( coords_b));
      const storage::Pair< biol::EnvironmentType, double> env_w3( membrane_ctg_n.DetermineEnvironmentTypeAndWeight( coords_c));
      const storage::Pair< biol::EnvironmentType, double> env_w4( membrane_ctg_n.DetermineEnvironmentTypeAndWeight( coords_d));
      const storage::Pair< biol::EnvironmentType, double> env_w5( membrane_ctg_n.DetermineEnvironmentTypeAndWeight( coords_e));
      const storage::Pair< biol::EnvironmentType, double> env_w6( membrane_ctg_n.DetermineEnvironmentTypeAndWeight( coords_f));
      BCL_ExampleIndirectCheck
      (
        env1 == env_w1.First() && math::EqualWithinTolerance( env_w1.Second(), 1.0) &&
        env2 == env_w2.First() && math::EqualWithinTolerance( env_w2.Second(), 1.0) &&
        env3 == env_w3.First() && math::EqualWithinTolerance( env_w3.Second(), 0.5) &&
        env4 == env_w4.First() && math::EqualWithinTolerance( env_w4.Second(), 1.0) &&
        env5 == env_w5.First() && math::EqualWithinTolerance( env_w5.Second(), 0.5) &&
        env6 == env_w6.First() && math::EqualWithinTolerance( env_w6.Second(), 1.0),
        true,
        "environment types for different z coordinates and their weights"
      );

      // calculate solvation energies for different positions in the membrane
      // free energy for core, interface and solution
      const linal::Vector3D sol_energies( 2, -1, 1);

      const double energy1( membrane_ctg_n.CalculateSolvationEnergy( linal::Vector3D( 0.0, 0.0, -30.0), sol_energies));
      const double energy2( membrane_ctg_n.CalculateSolvationEnergy( linal::Vector3D( 0.0, 0.0, -25.4), sol_energies));
      const double energy3( membrane_ctg_n.CalculateSolvationEnergy( linal::Vector3D( 0.0, 0.0, -18.5), sol_energies));
      const double energy4( membrane_ctg_n.CalculateSolvationEnergy( linal::Vector3D( 0.0, 0.0, -11.6), sol_energies));
      const double energy5( membrane_ctg_n.CalculateSolvationEnergy( linal::Vector3D( 0.0, 0.0, -00.1), sol_energies));
      const double energy6( membrane_ctg_n.CalculateSolvationEnergy( linal::Vector3D( 0.0, 0.0, 11.4), sol_energies));
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance( energy1,  1.0) &&
        math::EqualWithinTolerance( energy2,  0.382683) &&
        math::EqualWithinTolerance( energy3, -1.0) &&
        math::EqualWithinTolerance( energy4,  1.56066) &&
        math::EqualWithinTolerance( energy5,  2.0) &&
        math::EqualWithinTolerance( energy6,  1.79904),
        true,
        "solvation energies"
      );

      // test environment after transformation
      math::TransformationMatrix3D random_transform
      (
        coord::OrientationInterface::GenerateRandomTransformationAroundCenter
        (
          10.0,
          math::g_Pi / 2.0,
          linal::Vector3D( 0.0, 0.0, 0.0)
        )
      );
      membrane_ctg_n.Transform( random_transform);
      coords_a.Transform( random_transform);
      BCL_ExampleIndirectCheck
      (
        membrane_ctg_n.DetermineEnvironmentType( coords_a),
        biol::GetEnvironmentTypes().e_MembraneCore,
        "environment type after transform"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      // transformation matrix and membrane (core thickness) from a pdbtm xml file
      // open the file
      io::IFStream read_xml;
      BCL_ExampleMustOpenInputFile( read_xml, AddExampleInputPathToFilename( e_Biology, "3b5x.xml"));
      const storage::Pair< biol::Membrane, math::TransformationMatrix3D> membrane_trans
      (
        biol::Membrane::MembraneAndTransformationFromPDBTMXML( read_xml, 10.0, 2.5)
      );
      io::File::CloseClearFStream( read_xml);

      // check the translation
      const linal::Vector3D expected_translation( -13.52850080, -39.09541154, -21.92984962);
      BCL_ExampleCheck( membrane_trans.Second().GetTranslation(), expected_translation);

      linal::Vector3D expected_normal( -1.01088989, 2.11328077, 10.23533821);
      expected_normal.Normalize();
      BCL_ExampleCheckWithinAbsTolerance( membrane_trans.First().GetNormal(), expected_normal, 0.00001);

      // use transformation from xml file for the pdb
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "3b5x.pdb"));
      assemble::ProteinModel model( Proteins::GetModel( pdb_filename));
      const linal::Vector3D center_before( model.GetCenter());
      model.Transform( membrane_trans.Second());
      const linal::Vector3D center_after( model.GetCenter());

      // check that center is close to the membrane core
      BCL_ExampleCheck( math::Absolute( center_before.Z()) >= math::Absolute( center_after.Z()), true);

      // write membrane model
      const std::string
        output_filename( AddExampleOutputPathToFilename( membrane_default, "3b5x_membrane.pdb"));
      Proteins::WriteModelToPDB( model, output_filename);

      const std::string bio_pdb_filename( AddExampleInputPathToFilename( e_Biology, "1lgh.pdb"));
      assemble::ProteinModel model_bio( Proteins::GetModel( bio_pdb_filename));

      // transformation matrix to create the biomolecule
      // open the file
      BCL_ExampleMustOpenInputFile( read_xml, AddExampleInputPathToFilename( e_Biology, "1lgh.xml"));
      const storage::Vector< storage::Triplet< char, char, math::TransformationMatrix3D> > bio_trans( biol::Membrane::BioTransformationMatricesFromPDBTMXML( read_xml, model_bio.GetChainIDs()));
      io::File::CloseClearFStream( read_xml);

      BCL_ExampleCheck( !bio_trans.IsEmpty(), true);
      const assemble::ProteinModelMultiplier multiplier( bio_trans, model_bio);
      model_bio = multiplier( model_bio);
      BCL_ExampleMustOpenInputFile( read_xml, AddExampleInputPathToFilename( e_Biology, "1lgh.xml"));
      const storage::Pair< biol::Membrane, math::TransformationMatrix3D> bio_membrane_trans( biol::Membrane::MembraneAndTransformationFromPDBTMXML( read_xml, 10.0, 2.5));
      io::File::CloseClearFStream( read_xml);
      model_bio.Transform( bio_membrane_trans.Second());

      // write bio model
      const std::string
        bio_output_filename( AddExampleOutputPathToFilename( membrane_default, "1lgh_bio.pdb"));
      Proteins::WriteModelToPDB( model_bio, bio_output_filename);

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities for biol::Membrane");
      // open output file
      BCL_MessageVrb( "write object");
      WriteBCLObject( membrane_default);
      BCL_MessageVrb( "read object");
      biol::Membrane membrane_read;
      ReadBCLObject( membrane_read);

      // compare written and read object
      BCL_ExampleIndirectCheck
      (
        membrane_default.GetNormal() == membrane_read.GetNormal() &&
        membrane_default.GetLimits() == membrane_read.GetLimits() &&
        membrane_default.GetThicknesses() == membrane_read.GetThicknesses(),
        true,
        "read and write"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleBiolMembrane

  const ExampleClass::EnumType ExampleBiolMembrane::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolMembrane())
  );

} // namespace bcl
