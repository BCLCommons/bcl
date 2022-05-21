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
#include "fold/bcl_fold_phi_psi_generator_ccd.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_aa.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_phi_psi_generator_ccd.cpp
  //!
  //! @author alexanns
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldPhiPsiGeneratorCCD :
    public ExampleInterface
  {
  public:

    ExampleFoldPhiPsiGeneratorCCD *Clone() const
    {
      return new ExampleFoldPhiPsiGeneratorCCD( *this);
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

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConstantPhiPsi
    //!
    //! @author alexanns
    //!
    //! @date Sep 5, 2010
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class RandomNumberGenerator :
      public random::UniformDistribution
    {

    private:

    //////////
    // data //
    //////////

      bool m_RotatePhi;
      double m_RandomDouble;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      RandomNumberGenerator( const bool ROTATE_PHI, const double RANDOM_DOUBLE) :
        m_RotatePhi( ROTATE_PHI),
        m_RandomDouble( RANDOM_DOUBLE)
      {
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief random double
      //! @return random number in double range
      double Double() const
      {
        return m_RandomDouble;
      }

      //! @brief generate a random boolean
      //! @return random true or false
      bool Boolean() const
      {
        return m_RotatePhi;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    }; // class RandomNumberGenerator

    int Run() const
    {
      {
        BCL_MessageDbg( "test rotate phi");
        const RandomNumberGenerator rotate_phi_half( true, 0.5);

        // moving points
        //ATOM   1305  N   SER A  38      41.567  13.649  18.682  1.00 41.80           N
        //ATOM   1306  CA  SER A  38      41.761  14.279  20.012  1.00 24.25           C
        //ATOM   1307  C   SER A  38      42.506  13.445  20.987  1.00 29.88           C

        // target points
        //ATOM    304  N   SER A  38      47.398  10.001  19.166  1.00 41.80           N
        //ATOM    305  CA  SER A  38      48.246   9.790  20.366  1.00 24.25           C
        //ATOM    306  C   SER A  38      48.288   8.389  20.853  1.00 29.88           C

        const linal::Vector3D target_n_atom_coordinates(  47.398, 10.001, 19.166);
        const linal::Vector3D target_ca_atom_coordinates( 48.246,  9.790, 20.366);
        const linal::Vector3D target_c_atom_coordinates(  48.288,  8.389, 20.853);
        const linal::Vector3D moving_n_atom_coordinates(  41.567, 13.649, 18.682);
        const linal::Vector3D moving_ca_atom_coordinates( 41.761, 14.279, 20.012);
        const linal::Vector3D moving_c_atom_coordinates(  42.506, 13.445, 20.987);

        const double starting_distance_sum
        (
          math::Sqr( linal::Distance( target_n_atom_coordinates,  moving_n_atom_coordinates)) +
          math::Sqr( linal::Distance( target_ca_atom_coordinates, moving_ca_atom_coordinates)) +
          math::Sqr( linal::Distance( target_c_atom_coordinates,  moving_c_atom_coordinates))
        );

        BCL_MessageDbg( "starting distance sum is " + util::Format()( starting_distance_sum));

        storage::List< coord::CyclicCoordinateDescent::TargetAndMovingPointPair> target_and_moving_points;
        target_and_moving_points.PushBack( coord::CyclicCoordinateDescent::TargetAndMovingPointPair( target_n_atom_coordinates,  moving_n_atom_coordinates));
        target_and_moving_points.PushBack( coord::CyclicCoordinateDescent::TargetAndMovingPointPair( target_ca_atom_coordinates, moving_ca_atom_coordinates));
        target_and_moving_points.PushBack( coord::CyclicCoordinateDescent::TargetAndMovingPointPair( target_c_atom_coordinates,  moving_c_atom_coordinates));

        const fold::PhiPsiGeneratorCCD ccd_generator( target_and_moving_points, rotate_phi_half, biol::AASequenceFlexibility::e_CTerminal, math::Range< double>( 0.0, 1.0));

        const assemble::ProteinModel model( Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "2LZM_rotated_35_38.pdb")));

        const assemble::LocatorAA mutate_residue_locator( 'A', 35);
        const assemble::LocatorAA previous_residue_locator( 'A', 34);
        const assemble::LocatorAA following_residue_locator( 'A', 36);

        const util::ShPtr< biol::AABase> mutate_residue( mutate_residue_locator.Locate( model)->Clone());
        const util::ShPtr< biol::AABase> previous_residue( previous_residue_locator.Locate( model)->Clone());
        const util::ShPtr< biol::AABase> following_residue( following_residue_locator.Locate( model)->Clone());

        const fold::MutationResidue mutation_residue( mutate_residue, previous_residue, following_residue);

        const storage::VectorND< 2, double> phi_psi_undefined( ccd_generator( mutation_residue));
        BCL_ExampleCheck
        (
          !util::IsDefined( phi_psi_undefined.Second()), true
        );

        const double ccd_generated_phi( phi_psi_undefined.First());
        BCL_MessageDbg( "ccd generated phi is " + util::Format()( ccd_generated_phi));
        const double expected_ccd_generated_phi( -2.06423);
        BCL_ExampleCheckWithinTolerance( expected_ccd_generated_phi, ccd_generated_phi, 0.00001);

        // this should be deleted since PhiPsiGeneratorCCD does nothing with the distance sum
//        const double expected_post_rotation_distance_sum( 44.5621);

//        BCL_ExampleCheck
//        (
//          math::EqualWithinTolerance( expected_post_rotation_distance_sum, ccd_calculated_post_rotation_distance_sum, 0.00001), true
//        );

        // post ccd proposed rotation coordinates
        //ATOM    304  N   SER A  38      44.833  12.446  19.311  1.00 41.80           N
        //ATOM    305  CA  SER A  38      45.446  12.719  20.636  1.00 24.25           C
        //ATOM    306  C   SER A  38      45.920  11.512  21.356  1.00 29.88           C
        const linal::Vector3D post_rotation_n_atom_coordinates(  44.833, 12.446, 19.311);
        const linal::Vector3D post_rotation_ca_atom_coordinates( 45.446, 12.719, 20.636);
        const linal::Vector3D post_rotation_c_atom_coordinates(  45.920, 11.512, 21.356);

        // not as precise as ccd calculated distance sum, so this is a little bit larger than the ccd calculated sum
        const double coordinate_calculated_post_rotation_distance_sum
        (
          math::Sqr( linal::Distance( target_n_atom_coordinates,  post_rotation_n_atom_coordinates)) +
          math::Sqr( linal::Distance( target_ca_atom_coordinates, post_rotation_ca_atom_coordinates)) +
          math::Sqr( linal::Distance( target_c_atom_coordinates,  post_rotation_c_atom_coordinates))
        );

        BCL_MessageDbg( "post rotation distance sum is " + util::Format()( coordinate_calculated_post_rotation_distance_sum));

      }

      // test rotating psi
      {
        BCL_MessageDbg( "test rotate psi");
        const RandomNumberGenerator rotate_phi_half( false, 0.25);

        // moving points
        //ATOM    304  N   SER A  38      42.642  14.282  14.442  1.00 41.80           N
        //ATOM    305  CA  SER A  38      41.999  15.462  13.813  1.00 24.25           C
        //ATOM    306  C   SER A  38      40.550  15.301  13.536  1.00 29.88           C

        // target points
        //ATOM    304  N   SER A  38      47.398  10.001  19.166  1.00 41.80           N
        //ATOM    305  CA  SER A  38      48.246   9.790  20.366  1.00 24.25           C
        //ATOM    306  C   SER A  38      48.288   8.389  20.853  1.00 29.88           C

        const linal::Vector3D target_n_atom_coordinates(  47.398, 10.001, 19.166);
        const linal::Vector3D target_ca_atom_coordinates( 48.246,  9.790, 20.366);
        const linal::Vector3D target_c_atom_coordinates(  48.288,  8.389, 20.853);
        const linal::Vector3D moving_n_atom_coordinates(  42.642, 14.282, 14.442);
        const linal::Vector3D moving_ca_atom_coordinates( 41.999, 15.462, 13.813);
        const linal::Vector3D moving_c_atom_coordinates(  40.550, 15.301, 13.536);

        const double starting_distance_sum
        (
          math::Sqr( linal::Distance( target_n_atom_coordinates,  moving_n_atom_coordinates)) +
          math::Sqr( linal::Distance( target_ca_atom_coordinates, moving_ca_atom_coordinates)) +
          math::Sqr( linal::Distance( target_c_atom_coordinates,  moving_c_atom_coordinates))
        );

        BCL_MessageDbg( "starting distance sum is " + util::Format()( starting_distance_sum));

        storage::List< coord::CyclicCoordinateDescent::TargetAndMovingPointPair> target_and_moving_points;
        target_and_moving_points.PushBack( coord::CyclicCoordinateDescent::TargetAndMovingPointPair( target_n_atom_coordinates,  moving_n_atom_coordinates));
        target_and_moving_points.PushBack( coord::CyclicCoordinateDescent::TargetAndMovingPointPair( target_ca_atom_coordinates, moving_ca_atom_coordinates));
        target_and_moving_points.PushBack( coord::CyclicCoordinateDescent::TargetAndMovingPointPair( target_c_atom_coordinates,  moving_c_atom_coordinates));

        const fold::PhiPsiGeneratorCCD ccd_generator( target_and_moving_points, rotate_phi_half, biol::AASequenceFlexibility::e_CTerminal, math::Range< double>( 0.0, 1.0));

        const assemble::ProteinModel model( Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "2LZM_rotated_35_38_psi.pdb")));

        const assemble::LocatorAA mutate_residue_locator( 'A', 35);
        const assemble::LocatorAA previous_residue_locator( 'A', 34);
        const assemble::LocatorAA following_residue_locator( 'A', 36);

        const util::SiPtr< const biol::AABase> sp_mutate_residue( mutate_residue_locator.Locate( model));
        const util::SiPtr< const biol::AABase> sp_previous_residue( previous_residue_locator.Locate( model));
        const util::SiPtr< const biol::AABase> sp_following_residue( following_residue_locator.Locate( model));

        BCL_ExampleAssert( sp_mutate_residue.IsDefined() && sp_previous_residue.IsDefined() && sp_following_residue.IsDefined(), true);

        const util::ShPtr< biol::AABase> mutate_residue( sp_mutate_residue->Clone());
        const util::ShPtr< biol::AABase> previous_residue( sp_previous_residue->Clone());
        const util::ShPtr< biol::AABase> following_residue( sp_following_residue->Clone());

        const fold::MutationResidue mutation_residue( mutate_residue, previous_residue, following_residue);

        const storage::VectorND< 2, double> phi_psi( ccd_generator( mutation_residue));
        BCL_ExampleCheck
        (
          !util::IsDefined( phi_psi.First()), true
        );

        const double ccd_generated_psi( phi_psi.Second());
        BCL_MessageDbg( "ccd generated phi is " + util::Format()( ccd_generated_psi));
        const double expected_ccd_generated_psi( -2.01493);
        BCL_ExampleCheckWithinTolerance( expected_ccd_generated_psi, ccd_generated_psi, 0.00001);

        // this should be deleted since PhiPsiGeneratorCCD does nothing with the distance sum
//        const double expected_post_rotation_distance_sum( 255.233);

//        BCL_ExampleCheck
//        (
//          math::EqualWithinTolerance( expected_post_rotation_distance_sum, ccd_calculated_post_rotation_distance_sum, 0.00001), true
//        );

        // post ccd proposed rotation coordinates
        //ATOM    304  N   SER A  38      42.290  14.170  17.105  1.00 41.80           N
        //ATOM    305  CA  SER A  38      41.456  15.366  17.382  1.00 24.25           C
        //ATOM    306  C   SER A  38      40.062  15.065  17.793  1.00 29.88           C

        const linal::Vector3D post_rotation_n_atom_coordinates(  42.290, 14.170, 17.105);
        const linal::Vector3D post_rotation_ca_atom_coordinates( 41.456, 15.366, 17.382);
        const linal::Vector3D post_rotation_c_atom_coordinates(  40.062, 15.065, 17.793);

        // not as precise as ccd calculated distance sum, so this is a little bit larger than the ccd calculated sum
        const double coordinate_calculated_post_rotation_distance_sum
        (
          math::Sqr( linal::Distance( target_n_atom_coordinates,  post_rotation_n_atom_coordinates)) +
          math::Sqr( linal::Distance( target_ca_atom_coordinates, post_rotation_ca_atom_coordinates)) +
          math::Sqr( linal::Distance( target_c_atom_coordinates,  post_rotation_c_atom_coordinates))
        );

        BCL_MessageDbg( "post rotation distance sum is " + util::Format()( coordinate_calculated_post_rotation_distance_sum));

      }
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldPhiPsiGeneratorCCD

  const ExampleClass::EnumType ExampleFoldPhiPsiGeneratorCCD::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldPhiPsiGeneratorCCD())
  );

} // namespace bcl
