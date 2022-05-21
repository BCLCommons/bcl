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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "fold/bcl_fold_phi_psi_generator_ccd.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> PhiPsiGeneratorCCD::s_Instance
    (
      GetObjectInstances().AddInstance( new PhiPsiGeneratorCCD())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PhiPsiGeneratorCCD::PhiPsiGeneratorCCD() :
      m_TargetAndMovingPoints(),
      m_RandomNumberGenerator( random::GetGlobalRandom()),
      m_Direction( biol::AASequenceFlexibility::e_CTerminal),
      m_RandomFraction( 1.0, 1.0)
    {
    }

    //! @brief constructor taking member variable parameters
    //! @param TARGET_AND_MOVING_POINTS the list of target and moving points whose distances will be minimized
    //! @param RANDOM_NUMBER_GENERATOR the random number generator that should be used
    //! @param DIRECTION indicates if the rotation needs to occur in the C to N direction
    //! @param RANDOM_FRACTION_RANGE range a random fraction is drawn from and multiplied with suggested rotation
    PhiPsiGeneratorCCD::PhiPsiGeneratorCCD
    (
      const storage::List< coord::CyclicCoordinateDescent::TargetAndMovingPointPair> &TARGET_AND_MOVING_POINTS,
      const random::DistributionInterface &RANDOM_NUMBER_GENERATOR,
      const biol::AASequenceFlexibility::SequenceDirection &DIRECTION,
      const math::Range< double> &RANDOM_FRACTION_RANGE
    ) :
      m_TargetAndMovingPoints( TARGET_AND_MOVING_POINTS),
      m_RandomNumberGenerator( RANDOM_NUMBER_GENERATOR),
      m_Direction( DIRECTION),
      m_RandomFraction( RANDOM_FRACTION_RANGE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new PhiPsiGeneratorCCD
    PhiPsiGeneratorCCD *PhiPsiGeneratorCCD::Clone() const
    {
      return new PhiPsiGeneratorCCD( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &PhiPsiGeneratorCCD::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator taking MutationResidue and returning storage::VectorND< 2, double> with phi psi, respectively
    //! @param MUTATION_RESIDUE MutationResidue whose phi or psi will be changed
    //! @return phi or psi value, respectively, which will minimize the distance between target and moving points.
    //!         the other value will be undefined
    storage::VectorND< 2, double> PhiPsiGeneratorCCD::operator()( const MutationResidue &MUTATION_RESIDUE) const
    {
      // result initialized as undefined
      storage::VectorND< 2, double> phi_psi( util::GetUndefined< double>(), util::GetUndefined< double>());

      // get a random boolean indicating whether to rotate around phi or psi
      // "rotate_phi" is true  if rotation around phi will occur
      // "rotate_phi" is false if rotation around psi will occur
      const bool rotate_phi( m_RandomNumberGenerator.Boolean());

      // message if rotation around phi will occur or not
      BCL_MessageDbg( "rotate phi is " + util::Format()( rotate_phi));

      // create a CyclicCoordinateDescent which will be used to do the ccd calculation
      const coord::CyclicCoordinateDescent ccd;

      // define the rotation axis depending on whether phi or psi is being rotated around
      // also calculate the current value of whichever dihedral angle is being changed
      coord::LineSegment3D rotation_axis;
      double current_dihedral_angle;

      // true if rotation around phi will occur
      if( rotate_phi)
      {
        // set the origin and end-point of the rotation axis
        rotation_axis.SetStartPoint
        (
          MUTATION_RESIDUE.GetMutationResidue()->GetCA().GetCoordinates()
        );
        rotation_axis.SetEndPoint( MUTATION_RESIDUE.GetMutationResidue()->GetAtom( biol::GetAtomTypes().N).GetCoordinates());

        // calculate phi and set "current_dihedral_angle" to it
        current_dihedral_angle = MUTATION_RESIDUE.GetMutationResidue()->Phi();

        // does amino acid have internal phi
        if( !util::IsDefined( current_dihedral_angle))
        {
          // if it does not have internal phi, previous aa is required
          if( !MUTATION_RESIDUE.GetPreviousResidue().IsDefined())
          {
            return phi_psi;
          }
          current_dihedral_angle = MUTATION_RESIDUE.GetMutationResidue()->CalculatePhi
          (
            MUTATION_RESIDUE.GetPreviousResidue()->GetAtom( biol::GetAtomTypes().C)
          );
        }
      }
      else // rotation around psi will occur
      {
        // set the origin and end-point of the rotation axis
        rotation_axis.SetStartPoint( MUTATION_RESIDUE.GetMutationResidue()->GetAtom( biol::GetAtomTypes().C).GetCoordinates());
        rotation_axis.SetEndPoint
        (
          MUTATION_RESIDUE.GetMutationResidue()->GetCA().GetCoordinates()
        );

        // calculate psi and set "current_dihedral_angle" to it
        current_dihedral_angle = MUTATION_RESIDUE.GetMutationResidue()->Psi();

        // does amino acid have internal psi
        if( !util::IsDefined( current_dihedral_angle))
        {
          // if it does not have internal phi, following aa is required
          if( !MUTATION_RESIDUE.GetFollowingResidue().IsDefined())
          {
            return phi_psi;
          }
          current_dihedral_angle = MUTATION_RESIDUE.GetMutationResidue()->CalculatePsi
          (
            MUTATION_RESIDUE.GetFollowingResidue()->GetAtom( biol::GetAtomTypes().N)
          );
        }
      }

      // dihedral could not be determined
      if( !util::IsDefined( current_dihedral_angle))
      {
        return phi_psi;
      }

      // invert the rotation axis, if direction is c terminal
      if( m_Direction == biol::AASequenceFlexibility::e_CTerminal)
      {
        rotation_axis = rotation_axis.GetReverse();
      }

      // message how many target and moving points are left for the ccd calculation
      BCL_MessageDbg
      (
        "number of points being used is " + util::Format()( m_TargetAndMovingPoints.GetSize())
      );

//        // calculate the starting distance so that it can be output and checked against to make sure the ccd is working
//        // i.e. the ending sum distance should be less than this starting distance
//        double starting_sum_distance( fold::LocatorLoopDomain::CalculateDistanceSum( m_TargetAndMovingPoints));
//
//        // message the starting distance
//        BCL_MessageDbg( "starting sum difference is " + util::Format()( starting_sum_distance));
//
//        // message the starting dihedral angle
//        BCL_MessageDbg( "current dihedral angle is " + util::Format()( current_dihedral_angle));

      // calculate optimal rotation to minimize distance and put the results in
      const coord::CyclicCoordinateDescent::ResultsAndCoefficients rotation_result
      (
        ccd.GetOptimalRotationandDistance( rotation_axis, m_TargetAndMovingPoints)
      );

      // get the optimal rotation value from "rotation_result"
      const double optimal_rotation( rotation_result.GetOptimalRotation());

      // message the optimal rotation
      BCL_MessageDbg( "optimal rotation is " + util::Format()( optimal_rotation));

      // get a random number between 0 and 1 which will be used as the fraction of optimal rotation
      const double random_number( m_RandomNumberGenerator.Random( m_RandomFraction));

      // message the random fraction
      BCL_MessageDbg( "random fraction is " + util::Format()( random_number));

      // use "optimal_rotation" and "random_number" to determine the rotation that will be used
      double rotation( optimal_rotation * random_number);

      // message the effective rotation
      BCL_MessageDbg( "effective rotation will be " + util::Format()( rotation));

//        // calculate how far away from the target the moving points will be after "rotation"
//        const double distance_difference_sum
//        (
//          ccd.CalculateDistanceSum
//          (
//            rotation,
//            rotation_result.GetCoefficientA(), rotation_result.GetCoefficientB(), rotation_result.GetCoefficientC()
//          )
//        );
//
//        // message the resulting distance sum
//        BCL_MessageDbg( "resulting sum difference is " + util::Format()( distance_difference_sum));
//
//        // make sure that the resulting distance sum is less or equal within tolerance to the starting distance sum
//        BCL_Assert
//        (
//          starting_sum_distance > distance_difference_sum ||
//          math::EqualWithinTolerance( starting_sum_distance, distance_difference_sum),
//          "starting distance is " + util::Format()( starting_sum_distance) +
//          " while ending distance sum would be " + util::Format()( distance_difference_sum)
//        );

      // calculate what the dihedral angle needs to be. Negate rotation so that it happens in the correct direction
      const double new_dihedral_angle( current_dihedral_angle - rotation);

      // message what the new dihedral angle will be
      BCL_MessageDbg( "new dihedral angle will be " + util::Format()( new_dihedral_angle));

      // put the "new_dihedral_angle" into the correct place in vectorND "phi_psi"
      // "new_dihedral_angle" will be first if phi was rotated around, otherwise "new_dihedral_angle" will be second
      if( rotate_phi) //< rotated around phi
      {
        // set phi to "new_dihedral_angle"
        phi_psi.First() = new_dihedral_angle;
      }
      else //< rotated around psi
      {
        // set psi as "new_dihedral_angle"
        phi_psi.Second() = new_dihedral_angle;
      }

      // return the phi and psi values
      return phi_psi;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PhiPsiGeneratorCCD::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_TargetAndMovingPoints, ISTREAM);
      io::Serialize::Read( m_Direction, ISTREAM);
      io::Serialize::Read( m_RandomFraction, ISTREAM);

      // return "ISTREAM"
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PhiPsiGeneratorCCD::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_TargetAndMovingPoints, OSTREAM, INDENT);
      io::Serialize::Write( m_Direction, OSTREAM, INDENT);
      io::Serialize::Write( m_RandomFraction, OSTREAM, INDENT);

      // return "OSTREAM"
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
