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

#ifndef BCL_CHEMISTRY_PERTURB_MOLECULE_POSE_H_
#define BCL_CHEMISTRY_PERTURB_MOLECULE_POSE_H_

// include the namespace header
#include "bcl_chemistry.h"
#include <descriptor/bcl_descriptor.fwd.hh>

// include other forward headers - sorted alphabetically
#include "find/bcl_find.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_conformational_interface.h"
#include "bcl_chemistry_collector_valence.h"
#include "bcl_chemistry_fragment_complete.h"
#include "bcl_chemistry_fragment_constitution_shared.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "descriptor/bcl_descriptor_base.h"
#include "find/bcl_find_pick_interface.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_result.h"
#include "util/bcl_util_function_interface.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_si_ptr_list.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PerturbMoleculePose
    //! @brief Used to transform atom types inside a molecule
    //!
    //! @see @link example_chemistry_perturb_molecule_pose.cpp @endlink
    //! @author brownbp1
    //! @date May 22, 2020
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PerturbMoleculePose :
      public math::MutateInterface< FragmentComplete>
    {

    /////////////
    // friends //
    /////////////

    private:

    //////////
    // data //
    //////////

      // the conformational ensemble of the target small molecule
      FragmentEnsemble m_Ensemble;

      // rotation amount (radians)
      float m_RotAmount;

      // translation amount (angstroms)
      float m_TransAmount;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PerturbMoleculePose();

      //! @brief construct with ensemble
      PerturbMoleculePose
      (
        const FragmentEnsemble &ENSEMBLE
      );

      //! @brief construct with ensemble and rot/trans amounts
      PerturbMoleculePose
      (
        const FragmentEnsemble &ENSEMBLE,
        const float &ROT,
        const float &TRANS
      );

      //! @brief clone constructor
      PerturbMoleculePose *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking an SmallMolecule and returning a grown SmallMolecule
      //! @param FRAGMENT small molecule of interest
      //! @return Constitution after the mutate
      math::MutateResult< FragmentComplete> operator()
      (
        const FragmentComplete &MOLECULE
      ) const;

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class PerturbMoleculePose

  } // namespace chemistry
} // namespace bcl

#endif //BCL_CHEMISTRY_PERTURB_MOLECULE_POSE_H_
