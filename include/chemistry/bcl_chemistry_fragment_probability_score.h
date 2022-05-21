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

#ifndef BCL_CHEMISTRY_FRAGMENT_PROBABILITY_SCORE_H_
#define BCL_CHEMISTRY_FRAGMENT_PROBABILITY_SCORE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_one_four_interaction_score.h"
#include "bcl_chemistry_bond_angle_assignment.h"
#include "bcl_chemistry_conformation_comparison_by_dihedral_bins.h"
#include "bcl_chemistry_conformation_interface.h"
#include "bcl_chemistry_rotamer_dihedral_bond_data.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "storage/bcl_storage_hash_map.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FragmentProbabilityScore
    //! @brief Scores the conformation based on propensity of observing rotamers of constituent fragments
    //! @details for a given conformation, score is calculated based on which rotamers have been used for sampling
    //!           conformations
    //!
    //! @see @link example_chemistry_fragment_probability_score.cpp @endlink
    //! @author kothiwsk, mendenjl, brownbp1
    //! @date Sep 26, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentProbabilityScore :
      public math::FunctionInterfaceSerializable< FragmentComplete, double>
    {
    private:

      //! rotamer data of fragments that are part of the molecule of interest
      util::ShPtrVector< RotamerDihedralBondData>                        m_RotamerData;

      //! Number of dihedral bonds
      size_t                                                             m_NumberDihedralBonds;

      // Conformer comparer for dihedrals
      ConformationComparisonByDihedralBins                               m_DihedralBinComparer;

      //! whether to consider changes in isometry
      bool                                                               m_ConsiderIsometryChange;

      //! size of each component associated with each dihedral bond; 1 for chain bonds,
      //! !N-atoms-in-complete-ring-structure for rings
      storage::Vector< size_t>                                           m_DihedralComponentSize;

      //! internal 1-4 interaction score
      AtomOneFourInteractionScore                                        m_14Score;

    public:

// MinGW and apple (newer releases of clang) do not support storage::HashMap/unordered_map, so drop back to storage::Map
#if defined(__MINGW32__) || defined(__APPLE__)
      typedef storage::Map< uint64_t, std::pair< int, double> > t_Map;
#else
      typedef storage::HashMap< uint64_t, std::pair< int, double> > t_Map;
#endif

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      FragmentProbabilityScore();

      //! @brief constructor
      //! @param FRAGMENT_DATA vector of fragemnts that are part of the molecule of interest
      //! @param PROBABILITY probability of seeing a fragment in the CSD
      //! @param PRIORITY_ANGLES priority dihedral angle object for molecule whose conformations are being sampled
      //! @param FRAGMENT_SIZE_WT weight to be given to fragment size
      //! @param FRAGMENT_COUNTS_WT weight to be given to fragment counts
      //! @param ROTAMER_PROPENSITY_WT weight to be given to rotmaer propensity
      FragmentProbabilityScore
      (
        const util::ShPtrVector< RotamerDihedralBondData> &FRAGMENT_DATA,
        const bool &CONSIDER_ISOMETRY_CHANGE = true
      );

      //! @brief Clone function
      //! @return pointer to new FragmentProbabilityScore
      FragmentProbabilityScore *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief get the rotamer data used for the score
      const util::ShPtrVector< RotamerDihedralBondData> &GetRotamerData() const;

      //! @brief get the dihedral component sizes used to normalize per-dihedral scores
      const storage::Vector< size_t> &GetDihedralComponentSize() const;

      //! @brief get the number of dihedrals
      const size_t GetNumberDihedralBonds() const;

      //! @brief evaluate clashes for given atom pair
      //! @param MOLECULE molecule that needs to scored
      //! @return propensity score for observing rotamers that exist in conformation
      double operator()
      (
        const FragmentComplete &MOLECULE
      ) const;

      //! @brief evaluate 1-4 interaction score
      //! @param MOLECULE molecule that needs to scored
      //! @return 1-4 interaction score
      double Get14InteractionScore
      (
        const FragmentComplete &MOLECULE
      ) const;

      //! @brief get back each component of the score as a separate molecule with full details/energies/etc
      FragmentEnsemble GetScoreComponents
      (
        const FragmentComplete &FRAG_COMP
      ) const;

      //! @brief get back the weighted average dihedral scores prior to aggregation
      storage::Vector< double> GetDihedralScores
      (
        const FragmentComplete &MOLECULE
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief returns the returns the bin signature for molecule derived from dihedral angles of fragments
      //! @param MOLECULE molecule of interest
      //! @param FRAGMENT_ISOMORPHISM fragment that will be used to calculate the bin signature
      //! @oaran BIN_SIZE bin size being used for determining dihedral keys
      //! @return all possible dihedral bin signatures for part of molecule corresponding to the fragment being provided
      std::pair< linal::Vector< int>, linal::Vector< double> > GetMolecularBinFromFragment
      (
        const FragmentComplete &MOLECULE,
        const storage::Vector< storage::VectorND< 4, size_t> > &FRAGMENT_ISOMORPHISM,
        t_Map &BIN_HASH_TO_BIN
      ) const;

      //! @brief returns the dihedral bin in which dihedral bond belongs
      //! @param MOLECULE molecule of interest
      //! @param DIHEDRAL_BOND atom indices that make up the dihedral bond
      //! @oaran BIN_SIZE bin size being used for determining dihedral keys
      //! @return dihedral bin in which dihedral bond belongs
      std::pair< int, double> GetDihedralBin
      (
        const FragmentComplete &MOLECULE,
        const storage::VectorND< 4, size_t> &DIHEDRAL_BOND,
        t_Map &BIN_HASH_TO_BIN
      ) const;

    }; // class FragmentProbabilityScore

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_FRAGMENT_PROBABILITY_SCORE_H_
