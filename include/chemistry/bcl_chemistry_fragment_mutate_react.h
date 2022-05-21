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

#ifndef BCL_CHEMISTRY_FRAGMENT_MUTATE_REACT_H_

#define BCL_CHEMISTRY_FRAGMENT_MUTATE_REACT_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically
#include "find/bcl_find.fwd.hh"

// include parent class headers
#include "bcl_chemistry_fragment_react.h"
#include "math/bcl_math_mutate_interface.h"

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_conformational_interface.h"
#include "bcl_chemistry_fragment_complete.h"
#include "bcl_chemistry_fragment_constitution_shared.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "bcl_chemistry_fragment_mutate_interface.h"
#include "bcl_chemistry_reaction_ensemble.h"
#include "bcl_chemistry_reaction_search.h"
#include "bcl_chemistry_reaction_worker.h"
#include "find/bcl_find_pick_interface.h"
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
    //! @class FragmentMutateReact
    //! @brief Perform reactions on molecules within the FragmentMutateInterface framework. Multiple inheritance
    //! from both the abstract base / interface class FragmentMutateInterface and the generic class FragmentReact.
    //! FragmentMutateInterface holds a bunch of common data, declares virtual functions necessary for the
    //! design framework, and is never constructed apart from a derived class. Derived classes from this one should
    //! not encounter diamond inheritance problems, but please avoid deriving additional classes from here.
    //!
    //! @see @link example_chemistry_fragment_mutate_react.cpp @endlink
    //! @author brownbp1
    //! @date Oct 10, 2021
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentMutateReact :
      public FragmentMutateInterface, public FragmentReact
    {

    /////////////
    // friends //
    /////////////

    private:

    //////////
    // data //
    //////////

      //! overrides 3D conformer settings to just produce an arbitrary conformer without preserving spatial information
      bool m_LigandBased;

      //! pose-dependent; if 3D conformer matters, fix atoms with bad geometry even if they are in reference structure
      bool m_CorrectGeometry;

      //! pose-dependent; if 3D conformer matters, add all ring atoms from non-reference scaffolds to mobile selection
      bool m_CorrectNonReferenceRingGeometry;

      //! pose-dependent; if 3D conformer matters, fix atoms this many bonds out from any other mobile atom
      size_t m_AdditionalAdjacentAtoms;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FragmentMutateReact();

      //! @brief constructor with fragment react object
      FragmentMutateReact
      (
        const FragmentReact &REACT
      );

      //! @brief clone constructor
      FragmentMutateReact *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking an SmallMolecule and returning a grown SmallMolecule
      //! @param FRAGMENT small molecule of interest
      //! @return Constitution after the mutate
      math::MutateResult< FragmentComplete> operator()( const FragmentComplete &FRAGMENT) const;

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERROR_STREAM the stream to write errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    }; // class FragmentMutateReact

  } // namespace chemistry
} // namespace bcl
#endif //BCL_CHEMISTRY_FRAGMENT_MUTATE_REACT_H_
