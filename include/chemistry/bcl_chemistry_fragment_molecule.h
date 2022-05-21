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

#ifndef BCL_CHEMISTRY_FRAGMENT_MOLECULE_H_
#define BCL_CHEMISTRY_FRAGMENT_MOLECULE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_set.h"
#include "bcl_chemistry_fragment_constitution_shared.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "bcl_chemistry_hydrogens_handler.h"
#include "bcl_chemistry_sub_fragment.h"
#include "math/bcl_math_range.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FragmentMolecule
    //! @brief Fragments a given molecule into all possible fragments.
    //!
    //! @see @link example_chemistry_fragment_molecule.cpp @endlink
    //! @author kothiwsk
    //! @date May 03, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentMolecule :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      //! share pointer to configuration set (may be undefined, if finding constitutions)
      util::ShPtr< ConfigurationSet> m_ConfigurationSet;

      //! share pointer to constitution set (may be undefined, if finding configurations)
      util::ShPtr< ConstitutionSet> m_ConstitutionSet;

      //! maximum number of rotatable bonds that a molecule can contain to consider for fragmentation
      size_t                          m_MaxRot;

      //! whether to consider configurations
      bool                            m_ConsiderConfigurations;

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
      FragmentMolecule()
      {}

      //! @brief construct from list of FragmentConstitutionShared
      //! @param CONFORMATIONS shptr to conformation statistics
      FragmentMolecule
      (
        util::ShPtr< ConfigurationSet> &CONFIGURATION, const size_t MAX_ROT
      );

      //! @brief construct from list of FragmentConstitutionShared
      //! @param CONSTITUTION_SET shptr to constitution set
      FragmentMolecule
      (
        util::ShPtr< ConstitutionSet> &CONSTITUTION_SET, const size_t MAX_ROT
      );

      //! @brief Clone function
      //! @return pointer to new FragmentMolecule
      FragmentMolecule *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return a reference to conformation statistics
      //! @return a
      const util::ShPtr< ConfigurationSet> &GetConfigurations() const
      {
        return m_ConfigurationSet;
      }

      //! @brief return a reference to conformation statistics
      //! @return a
      const util::ShPtr< ConstitutionSet> &GetConstitutions() const
      {
        return m_ConstitutionSet;
      }

      //! @brief returns the total number of configurations
      //! @return the total number of configurations
      size_t GetSize() const
      {
        return m_ConfigurationSet->GetSize();
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief fragment a molecule/fragment that is passed
      //! @param MOLECULE molecule that needs to be fragmented
      void operator()( const FragmentComplete &MOLECULE);

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief the function that calls the fragmentation algorithm
      //! @brief MOLECULE_GRAPH graph containing atom information for molecule that needs fragmenting
      //! @param FRAGMENTS list containing sub fragments
      //! @param UNIQUE_FRAGMENTS list containing unique sub fragments
      //! @param NODE_ISO a set containing atom indices of original molecule being fragmented, from which unique FragmentComplete are constructed
      //! @return
      bool Initialize
      (
        const graph::ConstGraph< util::SiPtr< const AtomConformationalInterface>, size_t> &MOLECULE_GRAPH,
        storage::List< SubFragment> &FRAGMENTS,
        storage::List< SubFragment> &UNIQUE_FRAGMENTS,
        storage::Set< storage::Set< size_t> > &NODE_ISO
      );

      //! @brief fragment a molecule/fragment that is passed.
      //! @brief MOLECULE_GRAPH graph containing atom information for molecule that needs fragmenting
      //! @param MOLECULE part of original molecule which needs to be further fragmented
      //! @param EDGES_TO_BREAK a list of edges which are not in a ring and can be broken
      //! @param FRAGMENTS list containing sub fragments
      //! @param UNIQUE_FRAGMENTS list containing unique FragmentComplete
      //! @param NODE_ISO a set containing atom indices of original molecule being fragmented, from which unique FragmentComplete are constructed
      void Fragmentation
      (
        const graph::ConstGraph< util::SiPtr< const AtomConformationalInterface>, size_t> &MOLECULE_GRAPH,
        const FragmentComplete &MOLECULE,
        const SubFragment &SUB_FRAGMENT,
        storage::List< storage::Vector< sdf::BondInfo> > &EDGES_TO_BREAK,
        storage::List< SubFragment> &FRAGMENTS,
        storage::List< SubFragment> &UNIQUE_FRAGMENTS,
        storage::Set< storage::Set< size_t> > &NODE_ISO
      );

      //! @brief get a list of bonds that need to be broken i.e. all non ring bonds and a pair of ring bonds( usually)
      //! @brief that need to broken to get bonds sticking out of a ring.
      //! @param MOLECULE the molecule to be fragmented
      //! @return list of bonds that need be broken
      storage::List< storage::Vector< sdf::BondInfo> > GetBreakableConnectivityPairs
      (
        const FragmentComplete &MOLECULE
      ) const;

      //! @brief inserts fragments from list of sub fragments into list of unique_fragments and update fragment vertices
      //! that have been seen
      //! @param FRAGMENTS list of sub fragments that need to be fragmented
      //! @param UNIQUE_FRAGMENTS list to store unique fragment in FragmentComplete form
      //! @param NODE_ISO set to store set of vertices(isomorphic to fragment) which have been seen in fragments
      void UniqueSetInsert
      (
        storage::List< SubFragment> &FRAGMENTS,
        storage::List< SubFragment> &UNIQUE_FRAGMENTS,
        storage::Set< storage::Set< size_t> > &NODE_ISO
      );

      //! @brief function performs checks to see if subfragment passed in contained in existing sub fragments
      //! @param FRAGMENT_LIST list containing sub fragments
      //! @param FRAGMENT sub fragment which need to be checked for uniqueness for inserting into list
      void FragmentElimination
      (
        const graph::ConstGraph< util::SiPtr< const AtomConformationalInterface>, size_t> &MOLECULE_GRAPH,
        storage::List< SubFragment> &FRAGMENT_LIST,
        SubFragment &FRAGMENT
      );

    }; // class FragmentMolecule

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_FRAGMENT_MOLECULE_H_
