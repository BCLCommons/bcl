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
#include "opencl/bcl_opencl_quality_rmsd.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_collector_common_aa.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "command/bcl_command_flag_interface.h"
#include "opencl/bcl_opencl_tools.h"
#include "quality/bcl_quality_rmsd.h"
#include "util/bcl_util_si_ptr_list.h"
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opencl_quality_rmsd.cpp
  //!
  //! @author woetzen
  //! @date Jan 16, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOpenclQualityRMSD :
    public ExampleInterface
  {
  public:

    ExampleOpenclQualityRMSD *Clone() const
    {
      return new ExampleOpenclQualityRMSD( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
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

      {
        const linal::Vector3D c1[ 5] =
        {
          linal::Vector3D( 40.338, 4.848, 8.769),
          linal::Vector3D( 42.032, 3.315, 5.760),
          linal::Vector3D( 41.619, -0.354, 6.650),
          linal::Vector3D( 37.894, 0.039, 7.305),
          linal::Vector3D( 37.085, 2.778, 4.773)
        };

        linal::Vector3D c2[ 5] =
        {
          linal::Vector3D( 34.816, -3.039, 18.802),
          linal::Vector3D( 38.034, -2.978, 16.788),
          linal::Vector3D( 39.908, -2.243, 20.009),
          linal::Vector3D( 38.121, -5.112, 21.771),
          linal::Vector3D( 39.082, -7.417, 18.900)
        };

        const util::SiPtrVector< const linal::Vector3D> coord_const_1( 5, c1);
        const util::SiPtrVector< const linal::Vector3D> coord_const_2( 5, c2);
        util::SiPtrVector< linal::Vector3D> coord_2( 5, c2);

        const double expected_orig_rmsd( 14.8248);
        const double expected_suim_rmsd( 0.495305);

        // Calculate RMSD before superimposition
        const double calculated_orig_rmsd( opencl::QualityRMSD( opencl::GetTools().GetFirstCommandQueue(), false).CalculateMeasure( coord_const_1, coord_const_2));
        BCL_MessageStd( "rmsd of two sets of 5 coordinates each before superimposing: ");
        BCL_MessageStd( util::Format()( calculated_orig_rmsd));
        BCL_ExampleIndirectCheckWithinTolerance( calculated_orig_rmsd, expected_orig_rmsd, 0.001, "CalculateMeasure");

        // transformation for super imposition
        const math::TransformationMatrix3D superimpose_t( opencl::QualityRMSD( opencl::GetTools().GetFirstCommandQueue(), true).CalculateSuperimposition( coord_const_2, coord_const_1));
        BCL_MessageStd( "this is the calculated transformationmatrix:\n " + util::Format()( superimpose_t));

        const double calculated_suim_rmsd_virtual( opencl::QualityRMSD( opencl::GetTools().GetFirstCommandQueue(), true).CalculateMeasure( coord_const_1, coord_const_2));
        coord::TransformCoordinates( coord_2, superimpose_t);
        const double calculated_suim_rmsd_direct( opencl::QualityRMSD( opencl::GetTools().GetFirstCommandQueue(), false).CalculateMeasure( coord_const_1, coord_const_2));
        BCL_MessageStd( "rmsd of two sets of 5 coordinates each after superimposing: virtual=" + util::Format()( calculated_suim_rmsd_virtual) + " direct=" + util::Format()( calculated_suim_rmsd_direct));
        BCL_ExampleIndirectCheckWithinTolerance
        (
          calculated_suim_rmsd_virtual,
          expected_suim_rmsd,
          0.001,
          "RMSD of non-superimposed coordinates"
        );
        BCL_ExampleIndirectCheckWithinTolerance
        (
          calculated_suim_rmsd_direct,
          expected_suim_rmsd,
          0.001,
          "RMSD of superimposed coordinates"
        );
      }
      // create string "pdb_filename_a" which has path for example pdb file
      const std::string pdb_filename_a( AddExampleInputPathToFilename( e_Biology, "1ubi_ideal_model.pdb"));
      // create ProteinModel "protein_model_a" from "pdb_a"
      BCL_MessageStd( "building model_a from pdb_a chains and sse information");
      assemble::ProteinModel protein_model_a( Proteins::GetModel( pdb_filename_a));

      // create string "pdb_filename_b" which has path for example pdb file
      const std::string pdb_filename_b( AddExampleInputPathToFilename( e_Biology, "1ubi_unpaired_strand.pdb"));
      // create ProteinModel "protein_model_b" from "pdb_a"
      BCL_MessageStd( "building model_b from pdb_a chains and sse information");
      assemble::ProteinModel protein_model_b( Proteins::GetModel( pdb_filename_b));

      // create assemble::CollectorCommonAA to get the amino acids common to "protein_model_a" and "protein_model_b"
      assemble::CollectorCommonAA collector;

      // collect the common amino acids
      storage::VectorND< 2, util::SiPtrList< const biol::AABase> > common_amino_acids
      (
        collector.Collect( storage::VectorND< 2, assemble::ProteinModel>( protein_model_a, protein_model_b))
      );

      // create Set of atoms of interest "bb_atoms" with CA atom type
      storage::Set< biol::AtomType> bb_atoms( biol::GetAtomTypes().CA);

      // create SiPtrVectors "atom_coordinates_a" and "atom_coordinates_b" to hold the coordinates those in "bb_atoms"
      util::SiPtrVector< const linal::Vector3D> atom_coordinates_a;
      util::SiPtrVector< const linal::Vector3D> atom_coordinates_b;

      // get the coordinates of atoms in "bb_atoms"
      for
      (
        util::SiPtrList< const biol::AABase>::const_iterator
          itr_a( common_amino_acids.First().Begin()), itr_a_end( common_amino_acids.First().End()),
          itr_b( common_amino_acids.Second().Begin()), itr_b_end( common_amino_acids.Second().End());
        itr_a != itr_a_end && itr_b != itr_b_end;
        ++itr_a, ++itr_b
      )
      {
        atom_coordinates_a.Append( ( *itr_a)->GetAtomCoordinates( bb_atoms));
        atom_coordinates_b.Append( ( *itr_b)->GetAtomCoordinates( bb_atoms));
      }

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      opencl::QualityRMSD rmsd( opencl::GetTools().GetFirstCommandQueue(), true);
      BCL_ExampleCheck( rmsd.GetSuperimposeCoordinates(), true);

      // construct with boolean
      opencl::QualityRMSD rmsd_no_superimpose( opencl::GetTools().GetFirstCommandQueue(), false);
      BCL_ExampleCheck( rmsd_no_superimpose.GetSuperimposeCoordinates(), false);

    /////////////////
    // data access //
    /////////////////

      // check GetStaticClassName
      BCL_MessageStd( "test GetStatisClassName");
      const std::string correct_static_class_name( "bcl::opencl::QualityRMSD");
      BCL_ExampleCheck( GetStaticClassName< opencl::QualityRMSD>(), correct_static_class_name);

      // check GetClassIdentifier
      BCL_MessageStd( "test GetClassIdentifier");
      BCL_ExampleCheck( GetStaticClassName< opencl::QualityRMSD>(), rmsd.GetClassIdentifier());

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // expected values
      const double correct_rmsd( 10.2721);
      const double correct_rmsd_wo_superimpose( 39.4201);

      // test superimposed RMSD
      BCL_MessageStd( "test CalculateMeasure() with superimposition");
      const double calculated_rmsd( rmsd.CalculateMeasure( atom_coordinates_a, atom_coordinates_b));
      BCL_MessageStd( "rmsd: " + util::Format()( calculated_rmsd));

      // test RMSD without superimposition
      BCL_MessageStd( "test CalculateMeasure() without superimposition");
      const double calculated_rmsd_wo_superimpose
      (
        rmsd_no_superimpose.CalculateMeasure( atom_coordinates_a, atom_coordinates_b)
      );
      BCL_MessageStd( "rmsd_no_superimpose: " + util::Format()( calculated_rmsd_wo_superimpose));

      // check results
      BCL_ExampleCheckWithinTolerance( correct_rmsd, calculated_rmsd, 0.001);
      BCL_ExampleCheckWithinTolerance( correct_rmsd_wo_superimpose, calculated_rmsd_wo_superimpose, 0.001);

      // expected superimposition matrix
      const double expected_superimposition_array[ 16] =
      {
          0.937071,  0.233191,  -0.259845, 0,
          0.200901,  0.24855 ,   0.947556, 0,
          0.285546, -0.94013 ,   0.186061, 0,
        -33.6756  ,  5.3619  , -24.1805  , 1
      };
      const linal::Vector< double> expected_superimposition_vector( 16, expected_superimposition_array);

      BCL_MessageStd( "testing CalculateSuperimposition");
      const linal::Matrix< double> superimposition_matrix
      (
        rmsd.CalculateSuperimposition( atom_coordinates_a, atom_coordinates_b).GetMatrix()
      );
      BCL_MessageStd( "the superimposition matrix\n" + util::Format()( superimposition_matrix));
      BCL_ExampleIndirectCheckWithinTolerance
      (
        linal::Vector< double>( 16, superimposition_matrix.Begin()),
        expected_superimposition_vector,
        0.001,
        "check transformation matrix calculated"
      );

      // time test between opencl implementation and cpu implementation
      {
        storage::Set< double> high_res_cutoffs;
        for( double cut( 1.0); cut < 8.0; cut += 0.1)
        {
          high_res_cutoffs.Insert( cut);
        }

        {
          quality::RMSD calculator( true);
          util::Stopwatch rmsd_watch( "rmsd cpu", util::Message::e_Standard, true);
          const double rmsd_time( calculator.CalculateMeasure( atom_coordinates_a, atom_coordinates_b));
          BCL_MessageStd( "rmsd cpu: " + util::Format()( rmsd_time));
        }
        {
          opencl::QualityRMSD calculator( opencl::GetTools().GetFirstCommandQueue(), true);
          util::Stopwatch rmsd_watch( "rmsd gpu", util::Message::e_Standard, true);
          const double rmsd_time( calculator.CalculateMeasure( atom_coordinates_a, atom_coordinates_b));
          BCL_MessageStd( "rmsd gpu: " + util::Format()( rmsd_time));
        }
      }

    //////////////////////
    // input and output //
    //////////////////////

      // test read write functionality
      WriteBCLObject( rmsd);
      opencl::QualityRMSD rmsd_read( opencl::GetTools().GetFirstCommandQueue());
      ReadBCLObject( rmsd_read);
      BCL_ExampleCheck( rmsd_read.GetSuperimposeCoordinates(), rmsd.GetSuperimposeCoordinates());

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleOpenclQualityRMSD

  const ExampleClass::EnumType ExampleOpenclQualityRMSD::s_Instance
  (
    GetExamples().AddEnum( ExampleOpenclQualityRMSD())
  );

} // namespace bcl
