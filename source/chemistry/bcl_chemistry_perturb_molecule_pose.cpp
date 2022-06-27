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
#include "linal/bcl_linal_vector_3d.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_result.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "math/bcl_math_transformation_matrix_3d.h"
#include "quality/bcl_quality_rmsd.h"
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "chemistry/bcl_chemistry_perturb_molecule_pose.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_fragment_map_conformer.h"
#include "chemistry/bcl_chemistry_fragment_split_largest_component.h"
#include "chemistry/bcl_chemistry_fragment_track_mutable_atoms.h"
#include "chemistry/bcl_chemistry_voxel_grid_atom.h"
#include "find/bcl_find_collector_interface.h"
#include "random/bcl_random_uniform_distribution.h"
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

    //! @brief default constructor
    PerturbMoleculePose::PerturbMoleculePose() :
        m_Ensemble( FragmentEnsemble()),
        m_RotAmount( float( 0.2 * math::g_Pi)),
        m_TransAmount( float( 0.3))
    {
    }

    //! @brief construct with ensemble
    PerturbMoleculePose::PerturbMoleculePose
    (
      const FragmentEnsemble &ENSEMBLE
    ) :
      m_Ensemble( ENSEMBLE),
      m_RotAmount( float( 0.2 * math::g_Pi)),
      m_TransAmount( float( 0.3))
    {
    }

    //! @brief construct with ensemble and rot/trans amounts
    PerturbMoleculePose::PerturbMoleculePose
    (
      const FragmentEnsemble &ENSEMBLE,
      const float &ROT,
      const float &TRANS
    ) :
      m_Ensemble( ENSEMBLE),
      m_RotAmount( ROT),
      m_TransAmount( TRANS)
    {
    }

    //! @brief clone constructor
    PerturbMoleculePose *PerturbMoleculePose::Clone() const
    {
      return new PerturbMoleculePose( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PerturbMoleculePose::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an fragment and generating a new fragment by growing on a valence
    //! @param FRAGMENT small molecule of interest
    //! @return MutateResult with Constitution after the mutate
    math::MutateResult< FragmentComplete> PerturbMoleculePose::operator()
    (
      const FragmentComplete &MOLECULE
    ) const
    {
      // copy the molecule
      FragmentComplete molecule( MOLECULE);

        // perform a random rotation
        math::RotationMatrix3D rot;
        rot.SetRand( m_RotAmount); // in radians
        linal::Vector3D mol_center_coords( molecule.GetCenter());
        molecule.Translate( -mol_center_coords);
        molecule.Rotate( rot);
        molecule.Translate( mol_center_coords);

        // perform a random translation
        const double mag( random::GetGlobalRandom().Random( m_TransAmount));
        linal::Vector3D trans( mag, 0.0, 0.0);
        molecule.Translate( trans.Rotate( math::RotationMatrix3D().SetRand()));

      // swap one conformer for another
      if( m_Ensemble.GetSize())
      {
        // choose a random conformer
        FragmentEnsemble::const_iterator mol_itr( m_Ensemble.Begin());
        size_t pos( random::GetGlobalRandom().Random< double>( m_Ensemble.GetSize()));
        std::advance( mol_itr, pos);

        // superimpose with original conformer
        util::SiPtrVector< const linal::Vector3D> scaffold_coords( molecule.GetAtomCoordinates());
        util::SiPtrVector< const linal::Vector3D> molecule_coords( mol_itr->GetAtomCoordinates());
        math::TransformationMatrix3D transform( quality::RMSD::SuperimposeCoordinates( scaffold_coords, molecule_coords));
        molecule = *mol_itr;
        molecule.Transform( transform);
      }

      // return mutate result
      util::ShPtr< FragmentComplete> mutate_result( new FragmentComplete( molecule));
      return math::MutateResult< FragmentComplete>( mutate_result, *this);
    }

  ////////////////
  // operations //
  ////////////////

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace chemistry
} // namespace bcl
