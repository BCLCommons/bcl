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

#ifndef BCL_CHEMISTRY_MOLECULE_FEATURE_MAPPER_H_
#define BCL_CHEMISTRY_MOLECULE_FEATURE_MAPPER_H_

// forward headers from bcl - sorted alphabetically

// headers from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_graph_converter.h"
#include "bcl_chemistry_element_types.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_possible_atom_types_for_atom.h"
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
    //! @class MoleculeFeatureMapper
    //! @brief Maps pharmacophores onto a scaffold
    //!
    //! @see @link example_chemistry_molecule_feature_mapper.cpp @endlink
    //! @author geanesar, brownbp1
    //! @date Jul 11, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeFeatureMapper :
      public util::SerializableInterface
    {

    //////////
    // data //
    //////////

      //! descriptor that should be used for analysis
      descriptor::CheminfoProperty m_Descriptor;

      //! the type of atom comparison to perform
      ConformationGraphConverter::AtomComparisonTypeEnum m_AtomComparison;

      //! the type of bond comparison to perform
      ConfigurationalBondTypeData::DataEnum m_BondComparison;

    ////////////////////
    // nested structs //
    ////////////////////

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MoleculeFeatureMapper();

      //! @brief constructor
      //! @param SCORE_LABEL the score descriptor to use
      MoleculeFeatureMapper
      (
        const util::ObjectDataLabel &SCORE_LABEL
      );

      //! @brief Clone function
      //! @return pointer to new MoleculeFeatureMapper
      MoleculeFeatureMapper *Clone() const
      {
        return new MoleculeFeatureMapper( *this);
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

      //! @brief get the name of this class
      //! @return the name of this class
      const std::string &GetAlias() const
      {
        static const std::string s_name( "MoleculeFeatureMapper");
        return s_name;
      }

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief perturb a chemical structure and measure the score changes it makes to a model
      //! @details currently only removes atoms from a structure
      std::vector< std::map< size_t, float> > Perturb
      (
        const FragmentComplete &MOLECULE,
        const bool &IGNORE_H = false,
        const bool &SPLIT_LARGEST = false, // false to be compatible with what Alex originally coded
        storage::Vector< size_t> ATOM_INDICES = storage::Vector< size_t>(),
        storage::Vector< ElementType> ELEMENTS = storage::Vector< ElementType>()
      ) const;

      //! @brief perturb a chemical structure by splitting out fragments and measure the score changes it makes to a model
      //! @details all atoms shared by a fragment are assigned the score, unless normalization by atoms per fragment is requested
      //! @details if two or more fragments share atoms, then that atom will receive the sum of contributions
      //! @param MOLECULE parent molecule to be perturbed
      //! @param SPLITTER split interface determining how parent molecule will be perturbed
      //! @param NORMALIZE normalize per atom score by the number of atoms in the fragment
      //! @param IGNORE_H pass this to exclude hydrogen atoms from normalization; no effect if NORMALIZE is false
      //! @return indices for each atom in parent molecule and their corresponding score contributions
      std::vector< std::map< size_t, float> > PerturbByFragment
      (
        const FragmentComplete &MOLECULE,
        const util::Implementation< FragmentSplitInterface> &SPLITTER,
        const ConformationGraphConverter::AtomComparisonType &ATOM_TYPE = ConformationGraphConverter::e_ElementType,
        const ConfigurationalBondTypeData::Data &BOND_TYPE = ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness,
        const bool &MUTUALLY_MATCHING_ATOMS = false,
        const bool &NORMALIZE = false,
        const bool &IGNORE_H = false,
        const bool &AVERAGE = false,
        const size_t &MOL_INDEX = util::GetUndefinedSize_t(),
        const bool &OUTPUT_INTERMEDIATES = false
      ) const;

      //! @brief compare a molecule to other input molecules and assign differences in score to differences in substructure
      //! @details all atoms shared by a fragment are assigned the score, unless normalization by atoms per fragment is requested
      //! @details if two or more fragments share atoms, then that atom will receive the sum of contributions
      //! @param MOLECULE parent molecule to be perturbed
      //! @param ENSEMBLE ensemble of molecules against which parent molecule will be compared
      //! @param NORMALIZE normalize per atom score by the number of atoms in the fragment
      //! @param IGNORE_H pass this to exclude hydrogen atoms from normalization; no effect if NORMALIZE is false
      //! @param AVERAGE if true, compute average value per atom across all comparisons; otherwise take cumulative sums
      //! @param MOL_INDEX index of parent molecule in ENSEMBLE; leave undefined if parent molecule not in ENSEMBLE
      //! @return indices for each atom in parent molecule and their corresponding score contributions
      std::vector< std::map< size_t, float> > CompareSubstructuresNaive
      (
        const FragmentComplete &MOLECULE,
        const FragmentEnsemble &ENSEMBLE,
        const ConformationGraphConverter::AtomComparisonType &ATOM_TYPE = ConformationGraphConverter::e_ElementType,
        const ConfigurationalBondTypeData::Data &BOND_TYPE = ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness,
        const bool &NORMALIZE = false,
        const bool &IGNORE_H = false,
        const bool &AVERAGE = false,
        const size_t &MOL_INDEX = util::GetUndefinedSize_t(),
        const bool &OUTPUT_INTERMEDIATES = false
      ) const;

      //! @brief compare a molecule to other input molecules and assign differences in score to differences in substructure
      //! @details all atoms shared by a fragment are assigned the score, unless normalization by atoms per fragment is requested
      //! @details if two or more fragments share atoms, then that atom will receive the sum of contributions
      //! @param MOLECULE parent molecule to be perturbed
      //! @param ENSEMBLE ensemble of molecules against which parent molecule will be compared
      //! @param NORMALIZE normalize per atom score by the number of atoms in the fragment
      //! @param IGNORE_H pass this to exclude hydrogen atoms from normalization; no effect if NORMALIZE is false
      //! @param AVERAGE if true, compute average value per atom across all comparisons; otherwise take cumulative sums
      //! @param MOL_INDEX index of parent molecule in ENSEMBLE; leave undefined if parent molecule not in ENSEMBLE
      //! @return indices for each atom in parent molecule and their corresponding score contributions
      std::vector< std::map< size_t, float> > CompareSubstructuresRigorous
      (
        const FragmentComplete &MOLECULE,
        const FragmentEnsemble &ENSEMBLE,
        const ConformationGraphConverter::AtomComparisonType &ATOM_TYPE = ConformationGraphConverter::e_ElementType,
        const ConfigurationalBondTypeData::Data &BOND_TYPE = ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness,
        const bool &MUTUALLY_MATCHING_ATOMS = false,
        const bool &NORMALIZE = false,
        const bool &AVERAGE = false,
        const size_t &MOL_INDEX = util::GetUndefinedSize_t(),
        const bool &OUTPUT_INTERMEDIATES = false
      ) const;

      //! @brief calculates pharmacophore maps by iteratively removing atoms and determining the QSAR score
      //! @brief calculates average features over a common scaffold
      std::vector< std::map< size_t, float> > AverageFeatures
      (
        const FragmentComplete &SCAFFOLD,
        const FragmentEnsemble &MOLECULES,
        const std::vector< std::vector< std::map< size_t, float> > > &MAPS = std::vector< std::vector< std::map< size_t, float> > >(),
        const ConformationGraphConverter::AtomComparisonType &ATOM_COMPARISON_TYPE = ConformationGraphConverter::e_ElementType,
        const ConfigurationalBondTypeData::Data &BOND_TYPE_INFO = ConfigurationalBondTypeData::e_BondOrderOrAromatic
      ) const;

      storage::Vector< storage::Triplet< size_t, size_t, float> > GetFeatureRMSDOnScaffold
      (
        const FragmentComplete &SCAFFOLD,
        const FragmentEnsemble &MOLS,
        const std::vector< std::vector< std::map< size_t, float> > > &MAPS,
        const ConformationGraphConverter::AtomComparisonType &ATOM_COMPARISON_TYPE = ConformationGraphConverter::e_ElementType,
        const ConfigurationalBondTypeData::Data &BOND_TYPE_INFO = ConfigurationalBondTypeData::e_BondOrderOrAromatic
      ) const;

      float CompareCommonStructs
      (
        const storage::Vector< size_t> &SCAFF_ISO_MOL_1,
        const std::map< size_t, float> &PERTURBS_MOL_1,
        const storage::Vector< size_t> &SCAFF_ISO_MOL_2,
        const std::map< size_t, float> &PERTURBS_MOL_2
      ) const;

      //! @brief extracts fragments from a molecule; keeps atoms whose scores are outside of a given confidence interval
      //! @param MOLECULE the molecule to extract from
      //! @param ATOM_SCORES scores of each atom
      //! @param CONFIDENCE_INTERVAL any scores inside this range are considered statistically insignificant
      //! @return an ensemble of "significant" fragments
      FragmentEnsemble ExtractFragments
      (
        const FragmentComplete &MOLECULE,
        const std::map< size_t, float> &ATOM_SCORES,
        const math::Range< float> &CONFIDENCE_INTERVAL
      ) const;

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

    }; // class MoleculeFeatureMapper

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_MOLECULE_FEATURE_MAPPER_H_
