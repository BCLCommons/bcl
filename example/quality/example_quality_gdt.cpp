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
#include "quality/bcl_quality_gdt.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "quality/bcl_quality_average.h"
#include "quality/bcl_quality_measures.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_quality_gdt.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleQualityGDT :
    public ExampleInterface
  {
  public:

    ExampleQualityGDT *Clone() const
    {
      return new ExampleQualityGDT( *this);
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
      storage::Set< biol::AtomType> atom_types( biol::GetAtomTypes().CA);
      util::SiPtrVector< const linal::Vector3D> native_coords( native.GetAtomCoordinates( atom_types));
      util::SiPtrVector< const linal::Vector3D> model_coords( model.GetAtomCoordinates( atom_types));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // initialize variables
      const storage::Set< double> cutoff_set( storage::Set< double>::Create( 2.0, 3.0, 4.0, 5.0));
      const double cutoff( 8.0);
      const size_t default_seed_length( 3);
      const size_t seed_length( 8);

      // construct with a single cutoff and seed length
      quality::GDT gdt_8( cutoff);
      BCL_ExampleCheck( gdt_8.GetDistanceCutoff(), 8.0);
      BCL_ExampleCheck( gdt_8.GetSeedLength(), quality::GDT::GetDefaultSeedLength());

      // construct with cutoffs vector and seed length
      quality::Average gdt_avg( quality::GDT::CreateAverageGDT( cutoff_set, seed_length));

      // construct gdt_ts total score
      quality::Average gdt_ts( quality::GDT::CreateAverageGDT( quality::Measures::GetDistanceCutoffsTS()));

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck( GetStaticClassName< quality::GDT>(), gdt_8.GetClassIdentifier());

      // check GetDefaultSeedLength()
      BCL_ExampleCheck( quality::GDT::GetDefaultSeedLength(), default_seed_length);

    ////////////////
    // operations //
    ////////////////

      // initialize expected values
      const double expected_gdt_8_result( 98.6842);
      const double expected_gdt_ts_result( 70.7237);
      const double expected_gdt_ha_result( 51.3158);

      // if the current message level is verbose
      if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Verbose))
      {
        // iterate over different distances
        for( double cutoff( 0.5), end_cutoff( 10.0); cutoff <= end_cutoff; cutoff += 0.50)
        {
          // calculate the corresponding gdt_avg with this cutoff
          const double this_gdt( quality::GDT( cutoff).CalculateMeasure( native_coords, model_coords));
          // print out the result
          BCL_MessageVrb( "GDT " + util::Format()( cutoff) + " : " + util::Format()( this_gdt));
        }
      }

      // test operator with a specified cutoff
      BCL_MessageStd( "Checking CalculateMeasure() with a single cutoff");
      const double gdt_8_result( gdt_8.CalculateMeasure( native_coords, model_coords));
      BCL_MessageStd( "GDT_8: " + util::Format()( gdt_8_result));
      BCL_ExampleCheckWithinTolerance( gdt_8_result, expected_gdt_8_result, 0.001);

      // calculate GDT_TS which is ( GDT_1 + GDT_2 + GDT_4 + GDT_8) / 4
      BCL_MessageStd( "Checking CalculateMeasure() fort GDT_TS");
      const double gdt_ts_result( ( *quality::GetMeasures().e_GDT_TS)->CalculateMeasure( native_coords, model_coords));
      BCL_MessageStd( "GDT_TS: " + util::Format()( gdt_ts_result));
      BCL_ExampleCheckWithinTolerance( gdt_ts_result, expected_gdt_ts_result, 0.001);

      // calculate GDT_HA which is ( GDT_0.5 + GDT_1 + GDT_2 + GDT_4) / 4
      BCL_MessageStd( "Checking CalculateMeasure() fort GDT_HA");
      const double gdt_ha_result( ( *quality::GetMeasures().e_GDT_HA)->CalculateMeasure( native_coords, model_coords));
      BCL_MessageStd( "GDT_HA: " + util::Format()( gdt_ha_result));
      BCL_ExampleCheckWithinTolerance( gdt_ha_result, expected_gdt_ha_result, 0.001);

      // expected superimposition matrix
      const double expected_superimposition_array[ 16] =
      {
         0.367862 ,  -0.919359,   0.139487 , 0,
         0.0407068,   0.165783,   0.985322 , 0,
        -0.928989 ,  -0.356784,   0.0984093, 0,
        -2.66104  ,  20.3012  , -31.5503   , 1
      };

      const linal::Vector< double> expected_superimposition_vector( 16, expected_superimposition_array);

      BCL_MessageStd( "testing CalculateSuperimposition gdt 0.5");
      const storage::Pair< double, math::TransformationMatrix3D> superimposition_matrix
      (
        quality::GDT( 0.5).CalculateGDTAndSuperimposition( native_coords, model_coords)
      );
      BCL_MessageStd( "GDT 0.5: " + util::Format()( superimposition_matrix.First()));
      BCL_MessageStd( "the superimposition matrix\n" + util::Format()( superimposition_matrix));
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance
        (
          linal::Vector< double>( 16, superimposition_matrix.Second().GetMatrix().Begin()), expected_superimposition_vector, 0.001, 0.001
        ), true, "comparison failed for transformation matrix calculated"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // test read write
      WriteBCLObject( gdt_8);
      quality::GDT gdt_read( 0.0, 20);
      ReadBCLObject( gdt_read);
      BCL_ExampleCheck( gdt_read.GetDistanceCutoff(), gdt_8.GetDistanceCutoff());
      BCL_ExampleCheck( gdt_read.GetSeedLength(), gdt_8.GetSeedLength());

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleQualityGDT

  const ExampleClass::EnumType ExampleQualityGDT::s_Instance
  (
    GetExamples().AddEnum( ExampleQualityGDT())
  );

} // namespace bcl
