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

#ifndef BCL_CHEMISTRY_FRAGMENT_GROW_H_
#define BCL_CHEMISTRY_FRAGMENT_GROW_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically
#include "find/bcl_find.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_conformational_interface.h"
#include "bcl_chemistry_fragment_complete.h"
#include "bcl_chemistry_fragment_constitution_shared.h"
#include "bcl_chemistry_fragment_ensemble.h"
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
    //! @class FragmentGrow
    //! @brief Used for attaching a grow fragment to a base fragment
    //!
    //! @see @link example_chemistry_fragment_grow.cpp @endlink
    //! @author loweew, geanesar
    //! @date  03/30/2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentGrow :
      public math::MutateInterface< FragmentComplete>
    {

    /////////////
    // friends //
    /////////////

    private:

    //////////
    // data //
    //////////

      //! pool of fragments to be picked from
      util::ShPtr< FragmentEnsemble> m_FragmentPool;

      //! collector that collects possible atoms to which fragments can be added
      util::ShPtr< find::CollectorInterface< util::SiPtrList< const AtomConformationalInterface>, FragmentComplete> >
        m_AtomCollector;

      //! picker that decides at which atom in the molecule to add the fragment to
      util::ShPtr
      <
        find::PickInterface
        <
          util::SiPtr< const AtomConformationalInterface>,
          util::SiPtrList< const AtomConformationalInterface>
        >
      > m_AtomPicker;

      //! picker that picks the fragments to add to the molecule
      util::ShPtr< find::PickInterface< const FragmentComplete &, FragmentEnsemble> > m_FragmentPicker;

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
      FragmentGrow();

      //! @brief constructor from provided arguments
      //! @param FRAGMENT_POOL
      //! @param ATOM_COLLECTOR
      //! @param ATOM_PICKER
      //! @param FRAGMENT_PICKER
      FragmentGrow
      (
        const util::ShPtr< FragmentEnsemble> &FRAGMENT_POOL,
        const find::CollectorInterface< util::SiPtrList< const AtomConformationalInterface>, FragmentComplete> &ATOM_COLLECTOR,
        const find::PickInterface< util::SiPtr< const AtomConformationalInterface>, util::SiPtrList< const AtomConformationalInterface> > &ATOM_PICKER,
        const find::PickInterface< const FragmentComplete &, FragmentEnsemble> &FRAGMENT_PICKER
      );

      //! @brief clone constructor
      FragmentGrow *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief retrieves the fragment pool that this class uses by default
      //! @return a ShPtr to the fragment list
      util::ShPtr< FragmentEnsemble> GetFragmentPool() const;

      //! @brief Set the fragment pool
      //! @param FRAGMENT_POOL the fragment pool to use
      void SetFragmentPool( const util::ShPtr< FragmentEnsemble> &FragmentPool);

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

      //! @brief a function that adds a FragmentComplete to another FragmentComplete
      //! @param FRAGMENT the small molecule of interest
      //! @return a pointer to the original fragment with another piece added to it, or undefined if there's an error
      util::ShPtr< FragmentComplete> GrowFragment( const FragmentComplete &FRAGMENT) const;

      //! @brief create a copy of a molecule with a single random H removed
      //! @param MOLECULE the molecule to modify
      //! @return a new molecule with a single hydrogen removed
      FragmentComplete RemoveSingleH( const FragmentComplete &MOLECULE) const;

      //! @brief a function that adds a FragmentComplete to another FragmentComplete
      //! @param FRAGMENT the small molecule of interest
      //! @param FRAGMENT_LIST a list of fragments to use
      //! @return FragmentComplete with the newly added fragment
      util::ShPtr< FragmentComplete> AddFragmentFromList
      (
        const FragmentComplete &FRAGMENT,
        const FragmentEnsemble &FRAGMENT_LIST
      ) const;

      //! @brief adds a fragment to a base molecule at a specific atom
      //! @param BASE the base molecule to add a fragment to
      //! @param FRAGMENT the piece to add to BASE
      //! @param BASE_ATOM the atom on BASE to connect to FRAGMENT
      //! @param FRAGMENT_ATOM the atom on FRAGMENT to connect to BASE
      //! @param MAX_BOND_ORDER the maximum multiplicity of the formed bond (0=max allowed)
      static util::ShPtr< FragmentComplete> AddFragment
      (
        const FragmentComplete &BASE,
        const FragmentComplete &FRAGMENT,
        const util::SiPtr< const AtomConformationalInterface> &BASE_ATOM,
        const util::SiPtr< const AtomConformationalInterface> &FRAGMENT_ATOM,
        const size_t &MAX_BOND_ORDER = 0
      );

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream
      //! @param INDENT number of indentations
      //! @return ostream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class FragmentGrow

  } // namespace chemistry
} // namespace bcl

#endif //BCL_CHEMISTRY_FRAGMENT_GROW_H_
