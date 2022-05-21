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
#include "quality/bcl_quality_maxsub.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_quality_maxsub.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleQualityMaxSub :
    public ExampleInterface
  {
  public:

    ExampleQualityMaxSub *Clone() const
    {
      return new ExampleQualityMaxSub( *this);
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

      // get the protein model 1ubi
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
      assemble::ProteinModel
        native( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone));

      // get the model for 1ubi
      const std::string model_filename( AddExampleInputPathToFilename( e_Biology, "1ubi_rosetta_model.pdb"));
      assemble::ProteinModel
      model( Proteins::GetModel( model_filename, biol::GetAAClasses().e_AABackBone));

      // extract the CA atom coordinates from both native and model coordinates
      storage::Set< biol::AtomType> atom_types;
      atom_types.Insert( biol::GetAtomTypes().CA);
      util::SiPtrVector< const linal::Vector3D> native_coords( native.GetAtomCoordinates( atom_types));
      util::SiPtrVector< const linal::Vector3D> model_coords( model.GetAtomCoordinates( atom_types));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // initialize values
      const double default_cutoff( 3.5);
      const double cutoff( 4.5);
      const size_t default_seed_length( 3);
      const size_t seed_length( 3);
      const size_t default_nr_iterations( 4);
      const size_t nr_iterations( 2);

      // default constructor
      quality::MaxSub maxsub;
      BCL_ExampleCheck( maxsub.GetRMSDCutoff(), default_cutoff);
      BCL_ExampleCheck( maxsub.GetSeedLength(), default_seed_length);
      BCL_ExampleCheck( maxsub.GetNumberIterations(), default_nr_iterations);

      // construct with values
      quality::MaxSub maxsub_b( cutoff, seed_length, nr_iterations);
      BCL_ExampleCheck( maxsub_b.GetRMSDCutoff(), cutoff);
      BCL_ExampleCheck( maxsub_b.GetSeedLength(), seed_length);
      BCL_ExampleCheck( maxsub_b.GetNumberIterations(), nr_iterations);

      // clone constructor
      util::ShPtr< quality::MaxSub> sp_maxsub( maxsub.Clone());
      BCL_ExampleCheck( sp_maxsub->GetRMSDCutoff(), maxsub.GetRMSDCutoff());
      BCL_ExampleCheck( sp_maxsub->GetSeedLength(), maxsub.GetSeedLength());
      BCL_ExampleCheck( sp_maxsub->GetNumberIterations(), maxsub.GetNumberIterations());

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_MessageStd( "test GetClassIdentifier");
      BCL_ExampleCheck( GetStaticClassName< quality::MaxSub>(), maxsub.GetClassIdentifier());

    ////////////////
    // operations //
    ////////////////

      // expected values
      const double expected_maxsub_0_5_s3( 15.7895);
      const double expected_maxsub_0_5_s4( 15.7895);
      const double expected_maxsub_0_5_s5( 15.7895);
      const double expected_maxsub_0_5_s6( 18.4211);
      const double expected_maxsub_0_5_s7( 23.6842);

      // calculate the maxsub values with differing seed lengths
      BCL_MessageStd( "Test CalculateMeasure() with different seed lengths");
      const double maxsub_0_5_s3( quality::MaxSub( 0.5, 3).CalculateMeasure( native_coords, model_coords));
      BCL_MessageStd( "Maxsub 0.5 seed 3): " + util::Format()( maxsub_0_5_s3));
      const double maxsub_0_5_s4( quality::MaxSub( 0.5, 4).CalculateMeasure( native_coords, model_coords));
      BCL_MessageStd( "Maxsub 0.5 seed 4): " + util::Format()( maxsub_0_5_s4));
      const double maxsub_0_5_s5( quality::MaxSub( 0.5, 5).CalculateMeasure( native_coords, model_coords));
      BCL_MessageStd( "Maxsub 0.5 seed 5): " + util::Format()( maxsub_0_5_s5));
      const double maxsub_0_5_s6( quality::MaxSub( 0.5, 6).CalculateMeasure( native_coords, model_coords));
      BCL_MessageStd( "Maxsub 0.5 seed 6): " + util::Format()( maxsub_0_5_s6));
      const double maxsub_0_5_s7( quality::MaxSub( 0.5, 7).CalculateMeasure( native_coords, model_coords));
      BCL_MessageStd( "Maxsub 0.5 seed 7): " + util::Format()( maxsub_0_5_s7));

      // check the results
      BCL_ExampleCheckWithinTolerance( maxsub_0_5_s3, expected_maxsub_0_5_s3, 0.001);
      BCL_ExampleCheckWithinTolerance( maxsub_0_5_s4, expected_maxsub_0_5_s4, 0.001);
      BCL_ExampleCheckWithinTolerance( maxsub_0_5_s5, expected_maxsub_0_5_s5, 0.001);
      BCL_ExampleCheckWithinTolerance( maxsub_0_5_s6, expected_maxsub_0_5_s6, 0.001);
      BCL_ExampleCheckWithinTolerance( maxsub_0_5_s7, expected_maxsub_0_5_s7, 0.001);

      // expected values
      const double expected_maxsub_1( 34.2105);
      const double expected_maxsub_2( 51.3158);
      const double expected_maxsub_4( 84.2105);
      const double expected_maxsub_8( 97.3684);

      // calculate the maxsub values with differing RMSD cutoffs
      BCL_MessageStd( "Test CalculateMeasure() with different RMSD cutoffs");

      const double maxsub_1( quality::MaxSub( 1).CalculateMeasure( native_coords, model_coords));
      BCL_MessageStd( "Maxsub 1.0): " + util::Format()( maxsub_1));
      const double maxsub_2( quality::MaxSub( 2).CalculateMeasure( native_coords, model_coords));
      BCL_MessageStd( "Maxsub 2.0): " + util::Format()( maxsub_2));
      const double maxsub_4( quality::MaxSub( 4).CalculateMeasure( native_coords, model_coords));
      BCL_MessageStd( "Maxsub 4.0): " + util::Format()( maxsub_4));
      const double maxsub_8( quality::MaxSub( 8).CalculateMeasure( native_coords, model_coords));
      BCL_MessageStd( "Maxsub 8.0): " + util::Format()( maxsub_8));

      // check the results
      BCL_ExampleCheckWithinTolerance( maxsub_1, expected_maxsub_1, 0.001);
      BCL_ExampleCheckWithinTolerance( maxsub_2, expected_maxsub_2, 0.001);
      BCL_ExampleCheckWithinTolerance( maxsub_4, expected_maxsub_4, 0.001);
      BCL_ExampleCheckWithinTolerance( maxsub_8, expected_maxsub_8, 0.001);

      // expected superimposition matrix
      const double expected_superimposition_array[ 16] =
      {
        0.373727, -0.914483, 0.155076, 0,
        -0.0304841, 0.154991, 0.987445, 0,
        -0.927037, -0.373763, 0.0300472, 0,
        -0.650747, 20.8936, -30.9863, 1
      };
      const linal::Vector< double> expected_superimposition_vector( 16, expected_superimposition_array);

      BCL_MessageStd( "testing CalculateSuperimposition");
      const linal::Matrix< double> superimposition_matrix
      (
        quality::MaxSub( 2.0).CalculateSuperimposition( native_coords, model_coords).GetMatrix()
      );
      BCL_MessageStd( "the superimposition matrix\n" + util::Format()( superimposition_matrix));
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance
        (
          linal::Vector< double>( 16, superimposition_matrix.Begin()), expected_superimposition_vector, 0.001, 0.001
        ), true, "comparison failed for transformation matrix calculated"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // test read write functionalities
      WriteBCLObject( maxsub);
      quality::MaxSub maxsub_read;
      ReadBCLObject( maxsub_read);
      BCL_ExampleCheck( maxsub_read.GetRMSDCutoff(), maxsub.GetRMSDCutoff());
      BCL_ExampleCheck( maxsub_read.GetSeedLength(), maxsub.GetSeedLength());
      BCL_ExampleCheck( maxsub_read.GetNumberIterations(), maxsub.GetNumberIterations());

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleQualityMaxSub

  const ExampleClass::EnumType ExampleQualityMaxSub::s_Instance
  (
    GetExamples().AddEnum( ExampleQualityMaxSub())
  );

} // namespace bcl
