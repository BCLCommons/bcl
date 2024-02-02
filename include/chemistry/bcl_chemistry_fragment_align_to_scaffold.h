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

#ifndef BCL_CHEMISTRY_FRAGMENT_ALIGN_TO_SCAFFOLD_H_
#define BCL_CHEMISTRY_FRAGMENT_ALIGN_TO_SCAFFOLD_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically
#include "descriptor/bcl_descriptor.fwd.hh"
#include "find/bcl_find.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_conformational_interface.h"
#include "bcl_chemistry_configurational_bond_type_data.h"
#include "bcl_chemistry_conformation_comparison_property_field_correlation.h"
#include "bcl_chemistry_conformation_graph_converter.h"
#include "bcl_chemistry_fragment_complete.h"
#include "bcl_chemistry_fragment_constitution_shared.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "find/bcl_find_pick_interface.h"
#include "graph/bcl_graph_common_subgraph_isomorphism_base.h"
#include "graph/bcl_graph_subgraph_isomorphism.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_result.h"
#include "storage/bcl_storage_vector.h"
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
    //! @class FragmentAlignToScaffold
    //! @brief Align molecules based on substructure.
    //!
    //! @see @link example_chemistry_fragment_align_to_scaffold.cpp @endlink
    //! @author brownbp1
    //! @date May 18, 2020
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentAlignToScaffold :
      public util::ObjectInterface
    {

    private:

      //////////
      // data //
      //////////

      // atom comparison type
      ConformationGraphConverter::AtomComparisonType m_AtomType;

      // bond comparison type
      ConfigurationalBondTypeData::Data m_BondType;

      // minimum isomorphism size for alignment to occur
      size_t m_MinIsoSize;

      // solution type for the isomorphism search
      graph::CommonSubgraphIsomorphismBase::SolutionType m_SolutionType;

      // common subgraph isomorphism
      mutable graph::CommonSubgraphIsomorphism< size_t, size_t> m_CommonSubgraphIsomorphism;

      // subgraph isomorphism
      mutable graph::SubgraphIsomorphism< size_t, size_t> m_SubgraphIsomorphism;

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
      FragmentAlignToScaffold();

      //! @brief constructor with graph isomorphism settings
      //! @param ATOM_TYPE the atom comparison type to be used for substructure comparison
      //! @param BOND_TYPE the bond comparison type to be used for substructure comparison
      //! @param MIN_ISO_SIZE the minimum size a substructure can be and be a solution
      //! @param SOLUTION_TYPE the solution type for the graph isomorphism
      FragmentAlignToScaffold
      (
        const ConformationGraphConverter::AtomComparisonType &ATOM_TYPE,
        const ConfigurationalBondTypeData::Data &BOND_TYPE,
        const size_t &MIN_ISO_SIZE,
        const graph::CommonSubgraphIsomorphismBase::SolutionType &SOLUTION_TYPE = graph::CommonSubgraphIsomorphismBase::e_GreedyUnconnected
      );

      //! @brief clone constructor
      FragmentAlignToScaffold *Clone() const;

      /////////////////
      // data access //
      /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the atom type used for comparison
      const ConformationGraphConverter::AtomComparisonType &GetAtomType() const;

      //! @brief return the bond type used for comparison
      const ConfigurationalBondTypeData::Data &GetBondType() const;

      //! @brief return the size of the minimum allowed isomorphism
      const size_t GetMinIsoSize() const;

      //! @brief return the solution type used for comparison
      const graph::CommonSubgraphIsomorphismBase::SolutionType &GetSolutionType() const;

      //! @brief return the common subgraph isomorphism object
      const graph::CommonSubgraphIsomorphism< size_t, size_t> GetCommonSubgraphIsomorphism() const;

      //! @brief return the subgraph isomorphism object
      const graph::SubgraphIsomorphism< size_t, size_t> GetSubgraphIsomorphism() const;

      //! @brief set the atom type used for comparison
      void SetAtomType( const ConformationGraphConverter::AtomComparisonType &ATOM_TYPE);

      //! @brief set the bond type used for comparison
      void SetBondType( const ConfigurationalBondTypeData::Data &BOND_TYPE);

      // @brief set the minimum isomorphism size
      void SetMinIsoSize( const size_t MIN_ISO_SIZE);

      //! @brief set the solution type used for comparison
      void SetSolutionType( const graph::CommonSubgraphIsomorphismBase::SolutionType &SOLUTION_TYPE);

      ////////////////
      // operations //
      ////////////////

      //! @brief maximum common substructure alignment of small molecules
      //! @param TARGET_MOL the molecule to be aligned
      //! @param SCAFFOLD_MOL the molecule against which the target is aligned
      //! @param TARGET_MOL_INDICES atoms in TARGET_MOL that will be searched for a common substructure
      //! @param SCAFFOLD_MOL_INDICES atoms in SCAFFOLD_MOL that will be searched for a common substructure
      //! @return true if alignment occurs, false if the isomorphism size
      //! is below the minimum size allowed
      bool AlignToScaffold
      (
        FragmentComplete &TARGET_MOL,
        const FragmentComplete &SCAFFOLD_MOL,
        const storage::Vector< size_t> &TARGET_MOL_INDICES = storage::Vector< size_t>(),
        const storage::Vector< size_t> &SCAFFOLD_MOL_INDICES = storage::Vector< size_t>()
      ) const;

      //! @brief maximum common substructure alignment of small molecules with pose-dependent scoring
      //! @param TARGET_MOL the molecule to be aligned
      //! @param SCAFFOLD_MOL the molecule against which the target is aligned
      //! @param MDL the SDF file MDL property specifying the binding pocket filename
      //! @param BINDING_POCKET_FILENAME name of protein/protein binding pocket PDB file
      storage::Pair< bool, float> PoseSensitiveAlignToScaffold
      (
        FragmentComplete &TARGET_MOL,
        const FragmentComplete &SCAFFOLD_MOL,
        const descriptor::CheminfoProperty &SCORE,
        const std::string &MDL,
        const std::string &BINDING_POCKET_FILENAME
      );

      //! @brief maximum common substructure alignment of small molecule conformer ensembles
      //! @param TARGET_ENS the molecule ensemble to be aligned
      //! @param SCAFFOLD_MOL the molecule against which the targets are aligned
      //! @param COMPARER the metric to be used to compare alignments
      //! @return true if alignment occurs, false if the isomorphism size
      //! is below the minimum size allowed or the ensemble is empty
      storage::Vector< storage::Pair< bool, float> > AlignEnsembleToScaffold
      (
        FragmentEnsemble &TARGET_ENS,
        const FragmentComplete &SCAFFOLD_MOL,
        const util::Implementation< ConformationComparisonInterface> &COMPARER
          = util::Implementation< ConformationComparisonInterface>( "PropertyFieldDistance")
      ) const;

      //! @brief build a new conformer of the target molecule starting from the conformation of the largest shared substructure with a scaffold molecule
      //! @param TARGET_MOL the molecule for which a new conformer will be generated
      //! @param SCAFFOLD_MOL the molecule whose MCS with TARGET_MOL will be the core of the new conformation for TARGET_MOL
      //! @param TARGET_MOL_INDICES atoms in TARGET_MOL that will be searched for a common substructure
      //! @param SCAFFOLD_MOL_INDICES atoms in SCAFFOLD_MOL that will be searched for a common substructure
      //! @param COMPARER if provided, a conformer will be selected that has the minimum value of the specified property
      //! @param UPPER_BOUND upper bound on isomorphism search; if 0 then use EstimateUpperBounds; lower bound is default 1
      //! @return true if the conformer is generated successfully, false otherwise
      bool ConformerFromScaffoldMCS
      (
        FragmentComplete &TARGET_MOL,
        const FragmentComplete &SCAFFOLD_MOL,
        const storage::Vector< size_t> &TARGET_MOL_INDICES = storage::Vector< size_t>(),
        const storage::Vector< size_t> &SCAFFOLD_MOL_INDICES = storage::Vector< size_t>(),
        const descriptor::CheminfoProperty &COMPARER = descriptor::CheminfoProperty(),
        const size_t UPPER_BOUND = 0
      ) const;

      //! @brief build an ensemble of conformers of the target molecule starting from the conformation of a shared substructure with a scaffold molecule
      //! @details this function will generate one 3D conformer for each subgraph isomorphism up to a specified limit. Each conformer will be saved
      //! with the corresponding SampleByParts atom indices as an MDL property so that additional conformers can be generated from this starting point
      //! while keeping the common subgraph more-or-less fixed in space (barring any potential lever-arm effects). substructures are chosen in descending
      //! order based on size.
      //! @param TARGET_MOL the molecule for which a new conformer will be generated
      //! @param SCAFFOLD_MOL the molecule whose MCS with TARGET_MOL will be the core of the new conformation for TARGET_MOL
      //! @param TARGET_MOL_INDICES atoms in TARGET_MOL that will be searched for a common substructure
      //! @param SCAFFOLD_MOL_INDICES atoms in SCAFFOLD_MOL that will be searched for a common substructure
      //! @param COMPARER if provided, a conformer will be selected that has the minimum value of the specified property
      //! @param N_MAX_SOLUTIONS maximum number of largest subgraphs to consider when generating conformers
      //! @return an ensemble of conformers
      FragmentEnsemble ConformersFromScaffoldCS
      (
        const FragmentComplete &TARGET_MOL,
        const FragmentComplete &SCAFFOLD_MOL,
        const storage::Vector< size_t> &TARGET_MOL_INDICES = storage::Vector< size_t>(),
        const storage::Vector< size_t> &SCAFFOLD_MOL_INDICES = storage::Vector< size_t>(),
        const descriptor::CheminfoProperty &COMPARER = descriptor::CheminfoProperty(),
        const size_t N_MAX_SOLUTIONS = 100
      ) const;

      //! @brief build an ensemble of conformers of the target molecule starting from the conformation of a shared substructure with a scaffold molecule;
      //! iterates between ConformersFromScaffoldCS and ConformerFromScaffoldMCS reducing the size of the largest allowed subgraph progressively
      //! @param TARGET_MOL the molecule for which a new conformer will be generated
      //! @param SCAFFOLD_MOL the molecule whose MCS with TARGET_MOL will be the core of the new conformation for TARGET_MOL
      //! @param TARGET_MOL_INDICES atoms in TARGET_MOL that will be searched for a common substructure
      //! @param SCAFFOLD_MOL_INDICES atoms in SCAFFOLD_MOL that will be searched for a common substructure
      //! @param COMPARER if provided, a conformer will be selected that has the minimum value of the specified property
      //! @param N_MAX_SOLUTIONS maximum number of largest subgraphs to consider when generating conformers with ConformersFromScaffoldCS
      //! @param N_ITERATIONS the number of iterations
      //! @return an ensemble of conformers
      FragmentEnsemble ConformersFromScaffoldIterative
      (
        const FragmentComplete &TARGET_MOL,
        const FragmentComplete &SCAFFOLD_MOL,
        const storage::Vector< size_t> &TARGET_MOL_INDICES = storage::Vector< size_t>(),
        const storage::Vector< size_t> &SCAFFOLD_MOL_INDICES = storage::Vector< size_t>(),
        const descriptor::CheminfoProperty &COMPARER = descriptor::CheminfoProperty(),
        const size_t N_MAX_SOLUTIONS = 100,
        const size_t N_ITERATIONS = 1
      ) const;

      //////////////////////
      // helper functions //
      //////////////////////

      //! @brief extract a fragment from a molecule based on its indices
      static FragmentComplete ExtractFragmentByIndices( const FragmentComplete &MOLECULE, const storage::Vector< size_t> &INDICES);

      //! @brief get the maximum common substructure between two molecules
      graph::CommonSubgraphIsomorphism< size_t, size_t> FindCommonSubgraphIsomorphism( const FragmentComplete &MOL_A, const FragmentComplete &MOL_B) const;

      //! @brief get the common subgraphs between two molecules
      graph::SubgraphIsomorphism< size_t, size_t> FindSubgraphIsomorphism( const FragmentComplete &MOL_A, const FragmentComplete &MOL_B) const;

      //! @brief get a mutex
      static sched::Mutex &GetMutex();

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

    }; // class FragmentAlignToScaffold

  } // namespace chemistry
} // namespace bcl

#endif //BCL_CHEMISTRY_FRAGMENT_ALIGN_TO_SCAFFOLD_H_
