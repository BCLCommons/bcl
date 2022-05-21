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
#include "quality/bcl_quality_rmsd.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_collector_common_aa.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "quality/bcl_quality_measures.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_quality_rmsd.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleQualityRMSD :
    public ExampleInterface
  {
  public:

    ExampleQualityRMSD *Clone() const
    {
      return new ExampleQualityRMSD( *this);
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
        const double calculated_orig_rmsd( quality::RMSD::RealSpaceRMSD( coord_const_1, coord_const_2));
        BCL_MessageStd( "rmsd of two sets of 5 coordinates each before superimposing: ");
        BCL_MessageStd( util::Format()( calculated_orig_rmsd));
        BCL_Example_Check
        (
          math::EqualWithinTolerance( calculated_orig_rmsd, expected_orig_rmsd, 0.001),
          "RMSDCoordinates does not return expected rmsd of " + util::Format()( expected_orig_rmsd)
        );

        // transformation for super imposition
        const math::TransformationMatrix3D superimpose_t( quality::RMSD::SuperimposeCoordinates( coord_const_1, coord_const_2));
        BCL_MessageStd( "this is the calculated transformationmatrix:\n " + util::Format()( superimpose_t));

        const double calculated_suim_rmsd_virtual( quality::RMSD::SuperimposedRMSD( coord_const_1, coord_const_2));
        coord::TransformCoordinates( coord_2, superimpose_t);
        const double calculated_suim_rmsd_direct( quality::RMSD::RealSpaceRMSD( coord_const_1, coord_const_2));
        BCL_MessageStd( "rmsd of two sets of 5 coordinates each after superimposing: virtual=" + util::Format()( calculated_suim_rmsd_virtual) + " direct=" + util::Format()( calculated_suim_rmsd_direct));
        BCL_Example_Check
        (
          math::EqualWithinTolerance( calculated_suim_rmsd_virtual, expected_suim_rmsd, 0.001)
          && math::EqualWithinTolerance( calculated_suim_rmsd_direct, expected_suim_rmsd, 0.001),
          "RMSDCoordinates does not return expected rmsd of " + util::Format()( expected_suim_rmsd)
        );

        // superimpose coordinates on them self
        const math::TransformationMatrix3D self_impose( quality::RMSD::SuperimposeCoordinates( coord_const_1, coord_const_1));
        BCL_MessageStd( "superimpose coordinates onto them self: " + util::Format()( self_impose));
        BCL_ExampleCheckWithinAbsTolerance( self_impose.GetRotation().EffectiveRotationAngle(), 0.0, 0.00001);
        BCL_ExampleCheckWithinAbsTolerance( self_impose.GetTranslation().Norm(), 0.0, 0.00001);
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
      quality::RMSD rmsd;
      BCL_ExampleCheck( rmsd.GetSuperimposeCoordinates(), true);
      BCL_ExampleCheck( rmsd.GetIgnoreZCoordinates(), false);

      // construct with boolean
      quality::RMSD rmsd_no_superimpose( false, false);
      BCL_ExampleCheck( rmsd_no_superimpose.GetSuperimposeCoordinates(), false);
      BCL_ExampleCheck( rmsd_no_superimpose.GetIgnoreZCoordinates(), false);

    /////////////////
    // data access //
    /////////////////

      // check GetStaticClassName
      BCL_MessageStd( "test GetStatisClassName");
      const std::string correct_static_class_name( "bcl::quality::RMSD");
      BCL_ExampleCheck( GetStaticClassName< quality::RMSD>(), correct_static_class_name);

      // check GetClassIdentifier
      BCL_MessageStd( "test GetClassIdentifier");
      BCL_ExampleCheck( GetStaticClassName< quality::RMSD>(), rmsd.GetClassIdentifier());

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // expected values
      const double correct_rmsd( 10.2721);
      const double correct_rmsd_wo_superimpose( 39.4201);
      const double correct_rmsd_xy_superimpose( 18.4646);

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

      // test RMSD with superimposition with only XY coordinates
      BCL_MessageStd( "test CalculateMeasure() with superimposition without z coordinates");
      const double calculated_rmsd_xy_superimpose
      (
        ( *quality::GetMeasures().e_RMSD_XYSuperimposition)->CalculateMeasure( atom_coordinates_a, atom_coordinates_b)
      );
      BCL_MessageStd( "rmsd_xy: " + util::Format()( calculated_rmsd_xy_superimpose));

      // check results
      BCL_ExampleCheckWithinTolerance( correct_rmsd, calculated_rmsd, 0.001);
      BCL_ExampleCheckWithinTolerance( correct_rmsd_wo_superimpose, calculated_rmsd_wo_superimpose, 0.001);
      BCL_ExampleCheckWithinTolerance( correct_rmsd_xy_superimpose, calculated_rmsd_xy_superimpose, 0.001);

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
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance
        (
          linal::Vector< double>( 16, superimposition_matrix.Begin()), expected_superimposition_vector, 0.001, 0.001
        ), true, "comparison failed for transformation matrix calculated"
      );

//      {
//        // copy coordinates a
//        storage::Vector< linal::Vector3D> coord_a_copy( util::ConvertToStorageVector< linal::Vector3D>( atom_coordinates_a));
//
//        const quality::RMSD rmsd_xy_no_superimpose( false, true);
//        const quality::RMSD rmsd_xy_superimpose( true, true);
//
//        const double rmsd_no_superimposition( rmsd_xy_no_superimpose.CalculateMeasure( util::ConvertToConstSiPtrVector( coord_a_copy), atom_coordinates_b));
//        const double rmsd_superimposition( rmsd_xy_superimpose.CalculateMeasure( util::ConvertToConstSiPtrVector( coord_a_copy), atom_coordinates_b));
//
//        math::TransformationMatrix3D transform( rmsd_xy_superimpose.SuperimposeCoordinates( util::ConvertToConstSiPtrVector( coord_a_copy), atom_coordinates_b));
//        transform.Invert();
//        for( storage::Vector< linal::Vector3D>::iterator itr( coord_a_copy.Begin()), itr_end( coord_a_copy.End()); itr != itr_end; ++itr)
//        {
//          itr->Transform( transform);
//        }
//
//        const double rmsd_no_superimposition_transform( rmsd_xy_no_superimpose.CalculateMeasure( util::ConvertToConstSiPtrVector( coord_a_copy), atom_coordinates_b));
//        const double rmsd_superimposition_transform( rmsd_xy_superimpose.CalculateMeasure( util::ConvertToConstSiPtrVector( coord_a_copy), atom_coordinates_b));
//        BCL_Message
//        (
//          util::Message::e_Standard,
//          "testing xy superimposition:\n"
//          "rmsd xy no superimpose: " + util::Format()( rmsd_no_superimposition) + "\n" +
//          "rmsd xy    superimpose: " + util::Format()( rmsd_superimposition) + "\n" +
//          "rmsd xy no superimpose after transform: " + util::Format()( rmsd_no_superimposition_transform) + "\n" +
//          "rmsd xy    superimpose after transform: " + util::Format()( rmsd_superimposition_transform)
//        );
//
//        BCL_ExampleIndirectCheck
//        (
//          math::EqualWithinAbsoluteTolerance( rmsd_no_superimposition_transform, rmsd_superimposition, 0.25), true,
//          "checking that rmsd xy of transformed coordinates are equal to the rmsd xy superimposed without explicit transformation"
//        );
//      }

    //////////////////////
    // input and output //
    //////////////////////

      // test read write functionality
      WriteBCLObject( rmsd);
      quality::RMSD rmsd_read;
      ReadBCLObject( rmsd_read);
      BCL_ExampleCheck( rmsd_read.GetSuperimposeCoordinates(), rmsd.GetSuperimposeCoordinates());
      BCL_ExampleCheck( rmsd_read.GetIgnoreZCoordinates(), rmsd.GetIgnoreZCoordinates());

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleQualityRMSD

  const ExampleClass::EnumType ExampleQualityRMSD::s_Instance
  (
    GetExamples().AddEnum( ExampleQualityRMSD())
  );

} // namespace bcl
