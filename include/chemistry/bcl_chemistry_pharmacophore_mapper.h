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

#ifndef BCL_CHEMISTRY_PHARMACOPHORE_MAPPER_H_
#define BCL_CHEMISTRY_PHARMACOPHORE_MAPPER_H_

// forward headers from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor.fwd.hh"

// headers from bcl - sorted alphabetically
#include "bcl_chemistry_fragment_complete.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "descriptor/bcl_descriptor_combine.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_running_average_sd.h"
#include "opti/bcl_opti_print_interface.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_triplet.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_object_data_label.h"
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_si_ptr_vector.h"

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PharmacophoreMapper
    //! @brief Maps pharmacophores onto a scaffold
    //!
    //! @see @link example_chemistry_pharmacophore_mapper.cpp @endlink
    //! @author geanesar, loweew
    //! @date 05/05/2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PharmacophoreMapper :
      public util::ObjectInterface
    {

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class PharmMapMolData
      //! @brief Data for an individual pharm-map
      //!
      //! @see @link example_chemistry_pharmacophore_mapper.cpp @endlink
      //! @author geanesar, loweew
      //! @date 05/05/2014
      //!
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      struct PharmMapMolData
      {
        FragmentComplete m_Molecule; //!< the molecule
        graph::ConstGraph< size_t, size_t> m_MolGraph; //!< molecule graph
        storage::Vector< size_t> m_ScaffoldIsomorphism; //!< isomorphism to scaffold
        storage::Vector< size_t> m_Substitution; //! the substitution pattern of this molecule
        linal::Vector< float> m_SummedProperties;
        float m_Score;
      };

    //////////
    // data //
    //////////

      //! molecule to graph converter
      ConformationGraphConverter m_GraphMaker;

      //! Name of the property that should be used for predicted activities
      descriptor::Combine< AtomConformationalInterface, float> m_ScoreDescriptor;

      //! The properties to calculate
      descriptor::Combine< AtomConformationalInterface, float> m_Properties;

      //! Property names
      storage::Vector< std::string> m_PropertyStrings;

      //! the m_MolData for the scaffold
      PharmMapMolData m_ScaffoldData;

      //std::vector< MoleculeInfo> m_MoleculeInfo;
      std::vector< PharmMapMolData> m_MolData;

      //! an association of scaffold atoms with which pairs of indices in m_MolData change at a single substitition point
      storage::Map< size_t, storage::Vector< storage::Pair< size_t, size_t> > > m_SinglePointChanges;

      //! map from grow points to solutions of multiple least squares regression across all property differences at this grow point
      storage::Map< size_t, storage::Pair< linal::Vector< float>, float> > m_MLSSolutions;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor
      //! @param SCORE_LABEL the score descriptor to use
      //! @param PROPERTY_LABEL the label of which properties to use
      //! @param ATOM_COMPARISON_TYPE how to compare the scaffold to the molecules (atom info)
      //! @param BOND_TYPE_INFO how to compare the scaffold to molecules (bond info)
      PharmacophoreMapper
      (
//        const FragmentEnsemble &MOLECULES,
//        const FragmentComplete &SCAFFOLD,
        const util::ObjectDataLabel &SCORE_LABEL,
        const util::ObjectDataLabel &PROPERTY_LABEL,
        const ConformationGraphConverter::AtomComparisonType &ATOM_COMPARISON_TYPE = ConformationGraphConverter::e_ElementType,
        const ConfigurationalBondTypeData::Data &BOND_TYPE_INFO = ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness
      );

      //! @brief Clone function
      //! @return pointer to new PharmacophoreMapper
      PharmacophoreMapper *Clone() const
      {
        return new PharmacophoreMapper( *this);
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

      //! @brief set atom properties map
      //! @param MAP atom properties map
      void SetProperties( const util::ObjectDataLabel &PROPERTIES);

      const descriptor::Combine< AtomConformationalInterface, float> &GetProperties() const;

      const FragmentEnsemble &GetMolecules() const
      {
        static FragmentEnsemble dummy;
        return dummy;
        //return m_Molecules;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief calculate the pharmacophore map
      void Calculate( const std::string &DETAILS_FILENAME);

      util::ShPtr< storage::Pair< linal::Matrix< float>, linal::Vector< float> > > CalculatePropertyDifferences
      (
        const std::vector< PharmMapMolData> &DATA,
        const storage::Vector< storage::Pair< size_t, size_t> > &CHANGED_PAIRS
      ) const;

      void SetupMolData( const FragmentComplete &SCAFFOLD, const FragmentEnsemble &MOLECULES);

      void CalculateMaps
      (
        const FragmentComplete &SCAFFOLD,
        const FragmentEnsemble &INPUT_MOLECULES,
        const size_t &USE // TODO remove this
      );

      //! @brief calculates pharmacophore maps given a scaffold and a set of molecules
      //! @param SCAFFOLD the scaffold to use
      //! @param INPUT_MOLECULES the nominal set of molecules to use
      //! @param INCLUDE_INTERCEPT whether the least-squares regression should include an intercept column
      //! @param REPLACE_ZERO_PROPERTIES if true, this will replace zero-difference properties with constant values
      //!                                this is useful to do an MLS regression when some of the properties are nonsensical for
      //!                                certain grow points
      void CalculateMaps
      (
        const FragmentComplete &SCAFFOLD,
        const FragmentEnsemble &INPUT_MOLECULES,
        const bool &REPLACE_ZERO_PROPERTIES = false,
        const float &SCORE_TOLERANCE = 0.0,
        const bool &BINARY = false
      );

      //! @brief Scores a set of molecules using a pharmacophore map
      //! @param SCAFFOLD the scaffold to use
      //! @param MOLECULES the molecules to score
      //! @param COEFF_MAP a list of property coefficients for changes at each grow point
      //!                  (calculated by CalculateMaps)
      //! @return a map from molecule indices (keys) to pairs (values) of compared molecule indices (First) and
      //!         values relevant to score calculation.  The first value will always be the changed grow point (scaffold atom)
      //!         and the last value is the activity change going from m1<->m2. Intermediate values are property changes
      storage::Map< size_t, storage::Vector< storage::Pair< size_t, linal::Vector< float> > > > ScoreMolecules
      (
        const FragmentComplete &SCAFFOLD,
        const FragmentEnsemble &MOLECULES,
        const storage::Map< size_t, linal::Vector< float> > &COEFF_MAP
      ) const;

      //!@brief find common scaffolds for molecules
      //! @param MOLECULES the molecules to look through
      //! @param SAMPLING the percentage of pairwise molecule comparisons that should be done
      //! @return a mapping from scaffolds to which molecules match that scaffold
      storage::Vector< storage::Pair< graph::ConstGraph< size_t, size_t>, storage::Vector< size_t> > >
        FindScaffolds
      (
        const FragmentEnsemble &MOLECULES,
        const float &SAMPLING,
        const size_t &MIN_SCAFF_SIZE,
        const bool &PRINT_STATUS = false
      ) const;

      void MolecularFeatureMap
      (
        const FragmentComplete &SCAFFOLD,
        const FragmentEnsemble &MOLECULES,
        const util::ObjectDataLabel &MODEL,
        const std::string &OUT_PREFIX
      ) const;

      //! @brief get human-readable names of properties that are mapped
      const storage::Vector< std::string> &GetPropertyStrings() const;

      storage::Vector< size_t> GetChangedGrowPoints
      (
        const storage::Map< size_t, storage::Vector< size_t> > &MOL_1_FRAGS,
        const storage::Map< size_t, storage::Vector< size_t> > &MOL_2_FRAGS,
        const graph::ConstGraph< size_t, size_t> &MOL_1_GRAPH,
        const graph::ConstGraph< size_t, size_t> &MOL_2_GRAPH,
        const storage::Vector< size_t> &GROW_POINT_INDICES
      ) const;

      storage::Vector< math::RunningAverageSD< float> > CalculatePropertyDiffStatistics
      (
        const std::vector< PharmMapMolData> &DATA,
        const storage::Vector< storage::Pair< size_t, size_t> > &CHANGED_PAIRS
      ) const;

      math::RunningAverageSD< float> CalculateScoreDiffStatistics
      (
        const std::vector< PharmMapMolData> &DATA,
        const storage::Vector< storage::Pair< size_t, size_t> > &CHANGED_PAIRS
      ) const;

      //! @brief writes maps to a PML file
      void WriteMapsPML( const std::string &SCAFFOLD_FILENAME, const std::string &OUTPUT_FILENAME) const;

      //! @param DETAILS_PREFIX the prefix to use for files
      void WriteDetailsMols( const std::string &DETAILS_PREFIX) const;

      //! @brief write details of the calculations
      //! @param DETAILS_PREFIX the prefix to use for files
      void WriteDetailsCSV( const std::string &DETAILS_PREFIX) const;

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

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief collects fragments from a substructure given a substructure
      //! @param GRAPH the graph to use
      //! @param SUBSTRUCT vector of indices representing a substructure to match
      //! @return a map from vertices (keys) to the fragment associated with them (values)
      storage::Map< size_t, storage::Vector< size_t> > CollectFragmentsFromSubstructure
      (
        const graph::ConstGraph< size_t, size_t> &GRAPH,
        const storage::Vector< size_t> &SUBSTRUCT
      ) const;

      //! @brief gets pairs of molecules that differ at a single point of a common scaffold
      //! @return a vector of pairs of molecules
      static storage::Vector< storage::Triplet< size_t, size_t, size_t> > GetChangedPairs( const std::vector< PharmMapMolData> &DATA);

      static storage::Vector< size_t> GetScaffoldSubstitution
      (
        const graph::ConstGraph< size_t, size_t> GRAPH,
        const storage::Vector< size_t> &ISOMORPHISM
      );

      static storage::Vector< size_t> GetSubstitutionDifferences
      (
        graph::ConstGraph< size_t, size_t> GRAPH_MOL_1,
        const storage::Vector< size_t> &ISOMORPHISM_MOL_1,
        graph::ConstGraph< size_t, size_t> GRAPH_MOL_2,
        const storage::Vector< size_t> &ISOMORPHISM_MOL_2
      );

      //! @brief determine which atom of the scaffold differs between these two molecules
      //! @param GRAPH_MOL_1 the first molecule graph (used in calculation, must be mutable)
      //! @param ISOMORPHISM_MOL_1 the isomorphism to a scaffold for mol 1
      //! @param GRAPH_MOL_2 the second molecule graph (used in calculation, must be mutable)
      //! @param ISOMORPHISM_MOL_2 the isomorphism to a scaffold for mol 2
      //! @return the index of the scaffold at which these two graphs differ, or the size of the scaffold if none or multiple
      size_t GetSingleSubstitutionDifference
      (
        graph::ConstGraph< size_t, size_t> &GRAPH_MOL_1,
        const storage::Vector< size_t> &ISOMORPHISM_MOL_1,
        graph::ConstGraph< size_t, size_t> &GRAPH_MOL_2,
        const storage::Vector< size_t> &ISOMORPHISM_MOL_2
      ) const
      {
        BCL_Assert
        (
          ISOMORPHISM_MOL_1.GetSize() == ISOMORPHISM_MOL_2.GetSize(),
          "Trying to get a substitution difference of molecules with different scaffolds"
        );

        return 0;
      }

    }; // class PharmacophoreMapper

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_PHARMACOPHORE_MAPPER_H_
