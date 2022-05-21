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
#include "opencl/bcl_opencl_quality_gdt.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "command/bcl_command_flag_interface.h"
#include "opencl/bcl_opencl_tools.h"
#include "quality/bcl_quality_measures.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opencl_quality_gdt.cpp
  //!
  //! @author woetzen
  //! @date Jan 15, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOpenclQualityGDT :
    public ExampleInterface
  {
  public:

    ExampleOpenclQualityGDT *Clone() const
    {
      return new ExampleOpenclQualityGDT( *this);
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
      // currently only operational if -opencl flag was given
      if( !opencl::Platform::GetPlatformFlag().IsDefined() || !opencl::Platform::GetPlatformFlag()->GetFlag())
      {
        return 0;
      }

      // do not try to run opencl commands if no queue was found
      if( !opencl::GetTools().HasCommandQueues())
      {
        return 0;
      }

      if( !opencl::GetTools().GetFirstCommandQueue().GetDevice( NULL).Extensions( NULL).Contains( opencl::GetExtensions().e_khr_fp64))
      {
        return 0;
      }

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
      const storage::Set< double> default_cutoff_set( storage::Set< double>::Create( 1.0, 2.0, 4.0, 8.0));
      const storage::Set< double> cutoff_set( storage::Set< double>::Create( 2.0, 3.0, 4.0, 5.0));
      const storage::Set< double> cutoff_set_8( 8.0);
      const double cutoff( 8.0);
//      const size_t seed_length( 8);

      // construct with a single cutoff and seed length
      opencl::QualityGDT gdt_8( opencl::GetTools().GetFirstCommandQueue(), cutoff);
//      BCL_ExampleCheck( gdt_8.GetDistanceCutoffs(), cutoff_set_8);
      BCL_ExampleCheck( gdt_8.GetSeedLength(), quality::GDT::GetDefaultSeedLength());

//      // construct with cutoffs vector and seed length
//      opencl::QualityGDT gdt( opencl::GetTools().GetFirstCommandQueue(), cutoff_set, seed_length);
//      BCL_ExampleCheck( gdt.GetDistanceCutoffs(), cutoff_set);
//      BCL_ExampleCheck( gdt.GetSeedLength(), seed_length);
//
//      // clone construct
//      util::ShPtr< opencl::QualityGDT> sp_gdt( gdt.Clone());
//      BCL_ExampleCheck( sp_gdt->GetDistanceCutoffs(), gdt.GetDistanceCutoffs());
//      BCL_ExampleCheck( sp_gdt->GetSeedLength(), gdt.GetSeedLength());

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck( GetStaticClassName< opencl::QualityGDT>(), gdt_8.GetClassIdentifier());

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
          // calculate the corresponding gdt with this cutoff
          const double this_gdt( opencl::QualityGDT( opencl::GetTools().GetFirstCommandQueue(), cutoff).CalculateMeasure( native_coords, model_coords));
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
      const double gdt_ts_result( ( *quality::Measure( "OpenclGDT_TS"))->CalculateMeasure( native_coords, model_coords));
      BCL_MessageStd( "GDT_TS: " + util::Format()( gdt_ts_result));
      BCL_ExampleCheckWithinTolerance( gdt_ts_result, expected_gdt_ts_result, 0.001);

      // calculate GDT_HA which is ( GDT_0.5 + GDT_1 + GDT_2 + GDT_4) / 4
      BCL_MessageStd( "Checking CalculateMeasure() fort GDT_HA");
      const double gdt_ha_result( ( *quality::Measure( "OpenclGDT_HA"))->CalculateMeasure( native_coords, model_coords));
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
      storage::Pair< double, math::TransformationMatrix3D> superimposition_matrix
      (
        opencl::QualityGDT( opencl::GetTools().GetFirstCommandQueue(), 0.5).CalculateGDTAndSuperimposition( native_coords, model_coords)
      );
      BCL_MessageStd( "GDT 0.5: " + util::Format()( superimposition_matrix.First()));
      BCL_MessageStd( "the superimposition matrix\n" + util::Format()( superimposition_matrix));
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance
        (
          linal::Vector< double>( 16, superimposition_matrix.Second().GetMatrix().Begin()), expected_superimposition_vector, 0.001, 0.001
        ), true, "check transformation matrix calculated"
      );

//      // time test between opencl implementation and cpu implementation
//      {
//        storage::Set< double> high_res_cutoffs( quality::Measures::GetDistanceCutoffsHA());
//
//        {
//          quality::GDT calculator( high_res_cutoffs);
//          util::Stopwatch gdt_watch( "gdt cpu", util::Message::e_Standard, true);
//          const double high_res_gdt( calculator.CalculateMeasure( native_coords, model_coords));
//          BCL_MessageStd( "gdt high res cpu: " + util::Format()( high_res_gdt));
//        }
//        {
//          opencl::QualityGDT calculator( opencl::GetTools().GetFirstCommandQueue(), high_res_cutoffs);
//          util::Stopwatch gdt_watch( "gdt gpu", util::Message::e_Standard, true);
//          const double high_res_gdt( calculator.CalculateMeasure( native_coords, model_coords));
//          BCL_MessageStd( "gdt high res gpu: " + util::Format()( high_res_gdt));
//        }
//      }

    //////////////////////
    // input and output //
    //////////////////////

      // test read write
      WriteBCLObject( gdt_8);
      opencl::QualityGDT gdt_read( 0.0, 7);
      ReadBCLObject( gdt_read);
      BCL_ExampleCheck( gdt_read.GetDistanceCutoff(), gdt_8.GetDistanceCutoff());
      BCL_ExampleCheck( gdt_read.GetSeedLength(), gdt_8.GetSeedLength());

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleQualityGDT

  const ExampleClass::EnumType ExampleOpenclQualityGDT::s_Instance
  (
    GetExamples().AddEnum( ExampleOpenclQualityGDT())
  );

} // namespace bcl
