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
#include "chemistry/bcl_chemistry_mutate_dihedral_bond.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"
#include "chemistry/bcl_chemistry_bond_dihedral_angles.h"
#include "chemistry/bcl_chemistry_bond_isometry_handler.h"
#include "chemistry/bcl_chemistry_hydrogens_handler.h"
#include "chemistry/bcl_chemistry_mutate_fragment.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
  //////////
  // data //
  //////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor taking the member variable parameters
    //! @param FRAGMENT fragment which the mutate mutates
    //! @param GRAPH constitution graph of the molecule of interest
    //! @param MOLECULE molecule of interest
    //! @param WIGGLE_ONLY -- if true, do not change dihedral bin, just wiggle the dihedral within its bin
    MutateDihedralBond::MutateDihedralBond
    (
      const AtomConformationalInterface &ATOM_B,
      const AtomConformationalInterface &ATOM_C,
      const graph::ConstGraph< size_t, size_t> &GRAPH,
      const FragmentComplete &MOLECULE,
      const bool &FLIP_ONLY,
      const bool &WIGGLE_ONLY
    ) :
      m_AtomIndexA( 0),
      m_AtomIndexB( MOLECULE.GetAtomIndex( ATOM_B)),
      m_AtomIndexC( MOLECULE.GetAtomIndex( ATOM_C)),
      m_AtomIndexD( 0),
      m_FlipOnly( FLIP_ONLY),
      m_WiggleOnly( WIGGLE_ONLY),
      m_ConnectedAtoms(),
      m_MoleculeDihedralBondIndex
      (
        size_t( 1),
        PriorityDihedralAngles::GetDihedralEdges( MOLECULE).Find
        (
          graph::UndirectedEdge< ConfigurationalBondType>( m_AtomIndexB, m_AtomIndexC, ATOM_B.FindBondTo( ATOM_C)->GetBondType())
        )
      ),
      m_HasAmide( false),
      m_IsAmide( false),
      m_IsDoubleBond( ATOM_B.FindBondTo( ATOM_C)->GetBondType()->GetNumberOfElectrons() == size_t( 4))
    {
      BCL_Assert( ATOM_B.GetBonds().GetSize() > size_t( 1), "Atom B must have at least two bonds for dihedral rotations");
      BCL_Assert( ATOM_C.GetBonds().GetSize() > size_t( 1), "Atom C must have at least two bonds for dihedral rotations");
      m_AtomIndexA = MOLECULE.GetAtomIndex( ATOM_B.GetBonds().Begin()->GetTargetAtom());
      if( m_AtomIndexA == m_AtomIndexC)
      {
        m_AtomIndexA = MOLECULE.GetAtomIndex( ( ++ATOM_B.GetBonds().Begin())->GetTargetAtom());
      }
      m_AtomIndexD = MOLECULE.GetAtomIndex( ATOM_C.GetBonds().Begin()->GetTargetAtom());
      if( m_AtomIndexD == m_AtomIndexB)
      {
        m_AtomIndexD = MOLECULE.GetAtomIndex( ( ++ATOM_C.GetBonds().Begin())->GetTargetAtom());
      }

      m_ConnectedAtoms =
        MutateFragment::GetConnectedAtoms
        (
          graph::UndirectedEdge< ConfigurationalBondType>( m_AtomIndexB, m_AtomIndexC, ConfigurationalBondType()),
          GRAPH
        );

      // get the existing angle
      double existing_angle
      (
        math::Angle::Degree
        (
          linal::Dihedral
          (
            MOLECULE.GetAtomVector()( m_AtomIndexA).GetPosition(),
            MOLECULE.GetAtomVector()( m_AtomIndexB).GetPosition(),
            MOLECULE.GetAtomVector()( m_AtomIndexC).GetPosition(),
            MOLECULE.GetAtomVector()( m_AtomIndexD).GetPosition()
          )
        )
      );

      // take care of the wrap around
      if( existing_angle < 15.0)
      {
        existing_angle = 360.0 + existing_angle;
      }
      m_OriginalAngle = existing_angle;

      if
      (
        MOLECULE.GetAtomVector()( m_AtomIndexB).CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsAmide, 1)
        + MOLECULE.GetAtomVector()( m_AtomIndexC).CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsAmide, 1)
      )
      {
        m_HasAmide = true;
      }
      if( ATOM_B.FindBondTo( ATOM_C)->GetBondType()->GetConjugation() == ConstitutionalBondTypeData::e_Amide)
      {
        m_IsAmide = true;
      }
    }

    //! @brief Clone function
    //! @return pointer to new MutateDihedralBond
    MutateDihedralBond *MutateDihedralBond::Clone() const
    {
      return new MutateDihedralBond( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateDihedralBond::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns atom indices that are connected to different bonds
    //! @return atom indices that are connected to different rotatable bonds
    const storage::Vector< size_t> &MutateDihedralBond::GetConnectedAtoms() const
    {
      return m_ConnectedAtoms;
    }

    //! @brief get the dihedral bond indices referenced by this mutate
    const storage::Vector< size_t> &MutateDihedralBond::GetMoleculeDihedralBondIndices() const
    {
      return m_MoleculeDihedralBondIndex;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator taking a conformation and returns a mutated conformation
    //! @param MOLECULE conformation of interest
    //! @return MutateResult that results from mutating to the argument
    math::MutateResult< FragmentComplete> MutateDihedralBond::Mutate
    (
      const FragmentComplete &MOLECULE
    ) const
    {
//      BCL_Debug( "Mutating dihedral bond");
      // get atom info of the molecule of interest
      storage::Vector< sdf::AtomInfo> atom_info( MOLECULE.GetAtomInfo());

      // get the existing angle
      double existing_angle
      (
        math::Angle::Degree
        (
          linal::Dihedral
          (
            atom_info( m_AtomIndexA).GetCoordinates(),
            atom_info( m_AtomIndexB).GetCoordinates(),
            atom_info( m_AtomIndexC).GetCoordinates(),
            atom_info( m_AtomIndexD).GetCoordinates()
          )
        )
      );

      // take care of the wrap around
      if( existing_angle < 15.0)
      {
        existing_angle = 360.0 + existing_angle;
      }

      const double effective_std
      (
        m_HasAmide
        ? 4.5
        : BondDihedralAngles::GetEstimatedStdForDihedralBondAngleBin
          (
            MOLECULE.GetAtomVector()( m_AtomIndexB),
            MOLECULE.GetAtomVector()( m_AtomIndexC)
          )
      );
      double new_angle
      (
        m_WiggleOnly
        ? existing_angle // m_OriginalAngle
        :
          (
            m_FlipOnly
            ? random::GetGlobalRandom().Random( 1, 144) * 180.0
            : random::GetGlobalRandom().Random( 1, 144) * 30.0
          )
      );
      double new_angle_w_wiggle
      (
        effective_std < 7.5
        ? random::GetGlobalRandom().RandomGaussian( new_angle, effective_std)
        : random::GetGlobalRandom().Random( new_angle - 2.0 * effective_std, new_angle + 2.0 * effective_std)
      );
      if( m_WiggleOnly)
      {
        if( new_angle_w_wiggle < 15.0)
        {
          new_angle_w_wiggle += 360.0;
        }
        while( rint( existing_angle / 30.0) != rint( new_angle_w_wiggle / 30.0))
        {
          new_angle_w_wiggle =
           effective_std < 7.5
           ? random::GetGlobalRandom().RandomGaussian( new_angle, effective_std)
           : random::GetGlobalRandom().Random( new_angle - 2.0 * effective_std, new_angle + 2.0 * effective_std);
          if( new_angle_w_wiggle < 15.0)
          {
            new_angle_w_wiggle += 360.0;
          }
        }
      }
      MutateFragment::RotateBond
      (
        atom_info,
        m_ConnectedAtoms,
        new_angle_w_wiggle,
        existing_angle
      );
      // create new fragment from the updated coordinaes
      util::ShPtr< FragmentComplete> new_molecule
      (
        new FragmentComplete( AtomVector< AtomComplete>( atom_info, MOLECULE.GetBondInfo()), MOLECULE.GetName())
      );
      if( new_molecule->GetTotalAmideBondNonPlanarity( 20.0) > MOLECULE.GetTotalAmideBondNonPlanarity( 20.0) + 0.01)
      {
        return math::MutateResult< FragmentComplete>();
      }

      return math::MutateResult< FragmentComplete>( new_molecule, *this);
    }

    //! @brief remove the clash of two particular atoms by the most conservative movement of the given atoms if possible
    math::MutateResult< FragmentComplete> MutateDihedralBond::RemoveClash
    (
      const FragmentComplete &MOLECULE,
      const size_t &ATOM_A,
      const size_t &ATOM_B,
      const double &MIN_DISTANCE,
      const double &LOWER_BOUND,
      const double &UPPER_BOUND
    ) const
    {
      // get atom info of the molecule of interest
      storage::Vector< sdf::AtomInfo> atom_info( MOLECULE.GetAtomInfo());

      // get the existing angle
      double existing_angle
      (
        math::Angle::Degree
        (
          linal::Dihedral
          (
            atom_info( m_AtomIndexA).GetCoordinates(),
            atom_info( m_AtomIndexB).GetCoordinates(),
            atom_info( m_AtomIndexC).GetCoordinates(),
            atom_info( m_AtomIndexD).GetCoordinates()
          )
        )
      );

      const double wrapping_angle( m_IsAmide ? 20.0 : 15.0);
      // take care of the wrap around
      if( existing_angle <= wrapping_angle)
      {
        existing_angle = 360.0 + existing_angle;
      }
      const double original_angle( existing_angle);

      // get lower and upper bounds for the angle
      double lower_bound( m_WiggleOnly ? rint( existing_angle / 30.0) * 30.0 - wrapping_angle + 0.1 : wrapping_angle);
      double upper_bound( m_WiggleOnly ? rint( existing_angle / 30.0) * 30.0 + wrapping_angle - 0.1 : 360.0 + wrapping_angle);
      if( m_IsAmide)
      {
        const double tolerance( 20.0);
        lower_bound = existing_angle < 270.0 && existing_angle > 90.0 ? 180.0 - tolerance + 0.1 : 360.0 - tolerance + 0.1;
        upper_bound = existing_angle < 270.0 && existing_angle > 90.0 ? 180.0 + tolerance - 0.1 : 360.0 + tolerance - 0.1 ;
        lower_bound = std::min( lower_bound, original_angle);
        upper_bound = std::max( upper_bound, original_angle);
      }
      bool continue_expanding_range( true);
      lower_bound = std::max( lower_bound, LOWER_BOUND);
      upper_bound = std::min( upper_bound, UPPER_BOUND);
      double original_lower_bound( lower_bound);
      double original_upper_bound( upper_bound);
      const double dist_original
      (
        linal::Distance( atom_info( ATOM_A).GetCoordinates(), atom_info( ATOM_B).GetCoordinates())
      );
      const double new_angle( std::min( upper_bound, existing_angle + 0.5));
      double new_dist( dist_original);
      if( new_angle != existing_angle)
      {
        //BCL_MessageStd( "Track0.5");
        MutateFragment::RotateBond
        (
          atom_info,
          m_ConnectedAtoms,
          new_angle,
          existing_angle
        );
        new_dist = linal::Distance( atom_info( ATOM_A).GetCoordinates(), atom_info( ATOM_B).GetCoordinates());
        existing_angle = new_angle;
      }

      std::ostringstream oss;
      while( upper_bound - lower_bound < 360.0 && continue_expanding_range)
      {
        // test the distances at upper and lower bounds.
        MutateFragment::RotateBond
        (
          atom_info,
          m_ConnectedAtoms,
          upper_bound,
          existing_angle
        );
        const double dist_at_upper_bound
        (
          linal::Distance( atom_info( ATOM_A).GetCoordinates(), atom_info( ATOM_B).GetCoordinates())
        );
        existing_angle = upper_bound;
        MutateFragment::RotateBond
        (
          atom_info,
          m_ConnectedAtoms,
          lower_bound,
          existing_angle
        );
        const double dist_at_lower_bound
        (
          linal::Distance( atom_info( ATOM_A).GetCoordinates(), atom_info( ATOM_B).GetCoordinates())
        );
        existing_angle = lower_bound;

        double last_good( original_angle);

        if
        (
          ( m_HasAmide || m_IsDoubleBond)
          &&
          ( dist_at_lower_bound <= dist_original && dist_at_upper_bound <= dist_original && new_dist <= dist_original)
        )
        {
          //        BCL_MessageStd( "Nothing can be done");
          // nothing to do
          return math::MutateResult< FragmentComplete>();
        }

        // decide whether to go up or down
        bool can_go_up( dist_at_upper_bound >= MIN_DISTANCE), can_go_down( dist_at_lower_bound >= MIN_DISTANCE);
        bool just_max( false);
        if( !can_go_up && !can_go_down)
        {
          if( m_HasAmide || m_IsDoubleBond)
          {
            just_max = true;
            if( dist_at_upper_bound > dist_at_lower_bound)
            {
              can_go_up = true;
              last_good = upper_bound;
            }
            else
            {
              can_go_down = true;
              last_good = lower_bound;
            }
          }
          else
          {
            lower_bound -= 30.0;
            upper_bound += 30.0;
            continue;
          }
        }
        else if( can_go_up && can_go_down)
        {
          can_go_up = random::GetGlobalRandom().Double() < 0.5;
          can_go_down = !can_go_up;
        }
        continue_expanding_range = false;
        oss << "Refinement: "
          << lower_bound << ' ' << dist_at_lower_bound << ' '
          << upper_bound << ' ' << dist_at_upper_bound << ' '
          << original_angle << ' ' << dist_original << ' ' << MIN_DISTANCE << ' ';

        if( !just_max)
        {
          if( can_go_down)
          {
            upper_bound = original_angle;
            while( upper_bound - lower_bound > 0.5)
            {
              const double new_angle( ( upper_bound + lower_bound) * 0.5);
              MutateFragment::RotateBond
              (
                atom_info,
                m_ConnectedAtoms,
                new_angle,
                existing_angle
              );
              const double new_dist_in
              (
                linal::Distance( atom_info( ATOM_A).GetCoordinates(), atom_info( ATOM_B).GetCoordinates())
              );
              oss << '(' << new_angle << ',' << new_dist_in << ") ";
              if( new_dist_in >= MIN_DISTANCE)
              {
                lower_bound = new_angle;
                last_good = lower_bound;
              }
              else
              {
                upper_bound = new_angle;
              }
              existing_angle = new_angle;
            }
          }
          else if( can_go_up)
          {
            lower_bound = new_angle;
            while( upper_bound - lower_bound > 0.5)
            {
              const double new_angle( ( upper_bound + lower_bound) * 0.5);
              MutateFragment::RotateBond
              (
                atom_info,
                m_ConnectedAtoms,
                new_angle,
                existing_angle
              );
              const double new_dist_in
              (
                linal::Distance( atom_info( ATOM_A).GetCoordinates(), atom_info( ATOM_B).GetCoordinates())
              );
              oss << '(' << new_angle << ',' << new_dist_in << ") ";
              if( new_dist_in >= MIN_DISTANCE)
              {
                upper_bound = new_angle;
                last_good = upper_bound;
              }
              else
              {
                lower_bound = new_angle;
              }
              existing_angle = new_angle;
            }
          }
        }
        if( existing_angle != last_good)
        {
          oss << " rotated back to last good: " << last_good << " at distance: ";
          MutateFragment::RotateBond
          (
            atom_info,
            m_ConnectedAtoms,
            last_good,
            existing_angle
          );
          oss << linal::Distance( atom_info( ATOM_A).GetCoordinates(), atom_info( ATOM_B).GetCoordinates());
        }
      }
      BCL_MessageVrb( "Clash Removal trace: " + oss.str());
      if( ( linal::Distance( atom_info( ATOM_A).GetCoordinates(), atom_info( ATOM_B).GetCoordinates()) - dist_original) < 0.01)
      {
        return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
      }

      // create new fragment from the updated coordinaes
      util::ShPtr< FragmentComplete> new_molecule
      (
        new FragmentComplete( AtomVector< AtomComplete>( atom_info, MOLECULE.GetBondInfo()), MOLECULE.GetName())
      );

      // check whether amide bonds were messed up. This will generally only happen if m_HasAmide is true, so skip this
      // somewhat expensive test if we don't need it.
      const double amide_change
      (
        m_HasAmide
        ? new_molecule->GetTotalAmideBondNonPlanarity( 20.0) - MOLECULE.GetTotalAmideBondNonPlanarity( 20.0)
        : 0.0
      );
      if( amide_change > 0.01)
      {
        return
          RemoveClash
          (
            MOLECULE,
            ATOM_A,
            ATOM_B,
            MIN_DISTANCE,
            original_lower_bound + amide_change + 0.01,
            original_upper_bound - amide_change - 0.01
          );
      }

      return math::MutateResult< FragmentComplete>( new_molecule, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateDihedralBond::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_AtomIndexA, ISTREAM);
      io::Serialize::Read( m_AtomIndexB, ISTREAM);
      io::Serialize::Read( m_FlipOnly, ISTREAM);
      io::Serialize::Read( m_ConnectedAtoms, ISTREAM);

      // return "ISTREAM"
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateDihedralBond::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_AtomIndexA, OSTREAM, INDENT) << '\t' << m_AtomIndexB << '\t' << m_FlipOnly << '\n';
      io::Serialize::Write( m_ConnectedAtoms, OSTREAM, INDENT);

      // return "OSTREAM"
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
