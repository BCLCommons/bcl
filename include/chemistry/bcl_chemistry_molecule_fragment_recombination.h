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

#ifndef BCL_CHEMISTRY_MOLECULE_FRAGMENT_RECOMBINATION_H_
#define BCL_CHEMISTRY_MOLECULE_FRAGMENT_RECOMBINATION_H_

// forward headers from bcl - sorted alphabetically

// headers from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_graph_converter.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "bcl_chemistry_fragment_split_interface.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_triplet.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_object_data_label.h"
#include "util/bcl_util_serializable_interface.h"
#include "util/bcl_util_sh_ptr.h"

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeFragmentRecombination
    //! @brief Recombines differing substructures from two molecules into an ensemble of new molecules
    //!
    //! @see @link example_chemistry_molecule_fragment_recombination.cpp @endlink
    //! @author brownbp1
    //! @date Mar 04, 2021
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeFragmentRecombination :
      public FragmentSplitInterface
    {

    //////////
    // data //
    //////////

      //! the filename to read molecules from
      std::string m_File;

      //! whether molecules have been read or not
      mutable bool m_FileWasRead;

      //! the molecules that are read from the file, used to compare incoming molecules
      mutable FragmentEnsemble m_Molecules;

      //! the type of atom comparison to perform
      ConformationGraphConverter::AtomComparisonTypeEnum m_AtomComparison;

      //! the type of bond comparison to perform
      ConfigurationalBondTypeData::DataEnum m_BondComparison;

      //! the conformation graph converter to use
      ConformationGraphConverter m_Converter;

      //! use pre-aligned mutually-matching atoms instead of substructure comparison
      bool m_MutuallyMatchingAtoms;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MoleculeFragmentRecombination
      (
        const std::string &FILENAME = "",
        const ConformationGraphConverter::AtomComparisonType &ATOM_COMPARISON = ConformationGraphConverter::e_ElementType,
        const ConfigurationalBondTypeData::Data &BOND_COMPARISON = ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness,
        const bool &MUTUALLY_MATCHING_ATOMS = false
      );

      //! @brief Clone function
      //! @return pointer to new MoleculeFragmentRecombination
      MoleculeFragmentRecombination *Clone() const
      {
        return new MoleculeFragmentRecombination( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the name of this class
      //! @return the name of this class
      const std::string &GetAlias() const;

      //! @brief Get a description for what this class does (used when writing help)
      //! @return a description for what this class does (used when writing help)
      const std::string &GetClassDescription() const;

      //! @brief gets the minimum size of fragments
      //! @return the minimum size of fragments
      const size_t GetMinSize() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief returns an ensemble of fragments of a molecule
      //! @param CONFORMATION molecule of interest
      //! @return an ensemble of common substructures relative to those in a file
      FragmentEnsemble operator()( const ConformationInterface &CONFORMATION) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief compare a molecule to other input molecules and assign differences in score to differences in substructure
      //! @details all atoms shared by a fragment are assigned the score, unless normalization by atoms per fragment is requested
      //! @details if two or more fragments share atoms, then that atom will receive the sum of contributions
      //! @param MOLECULE parent molecule to be perturbed
      //! @param ENSEMBLE ensemble of molecules against which parent molecule will be compared
      //! @param MOL_INDEX index of parent molecule in ENSEMBLE; leave undefined if parent molecule not in ENSEMBLE
      //! @return indices for each atom in parent molecule and their corresponding score contributions
      FragmentEnsemble CompareSubstructures
      (
        const FragmentComplete &MOL_A,
        const FragmentComplete &MOL_B,
        const bool &MUTUALLY_MATCHING_ATOMS = false,
        const bool &OUTPUT_INTERMEDIATES = false
      ) const;

      FragmentComplete RecombineComponents
      (
        const FragmentComplete &PARENT_A,
        const FragmentComplete &PARENT_B,
        const FragmentComplete &BASE_MOL_A,
        const FragmentComplete &BASE_MOL_B,
        const storage::Map< size_t, size_t> &ISO,
        const FragmentComplete &FRAG_A,
        const storage::Vector< size_t> &FRAG_A_COMPONENT_INDICES,
        const FragmentComplete &FRAG_B,
        const storage::Vector< size_t> &FRAG_B_COMPONENT_INDICES
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief reads in molecules from a given file if it is necessary
      void ReadFile() const;

    protected:

      //! @brief Set the members with LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook
      (
        const util::ObjectDataLabel &LABEL,
        std::ostream &ERR_STREAM
      );

      io::Serializer GetSerializer() const;

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class MoleculeFragmentRecombination

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_MOLECULE_FRAGMENT_RECOMBINATION_H_
