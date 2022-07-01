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

#ifndef BCL_CHEMISTRY_MOLECULE_EVOLUTIONARY_OPTIMIZER_H_
#define BCL_CHEMISTRY_MOLECULE_EVOLUTIONARY_OPTIMIZER_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// headers from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_graph_converter.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "bcl_chemistry_fragment_evolve_implementations.h"
#include "bcl_chemistry_fragment_mutate_interface.h"
#include "bcl_chemistry_fragment_react.h"
#include "bcl_chemistry_molecule_evolution_info.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_comparisons.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_triplet.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_object_data_label.h"
#include "util/bcl_util_serializable_interface.h"
#include "util/bcl_util_sh_ptr.h"

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeEvolutionaryOptimizer
    //! @brief Evolves a molecule according to a fitness function and a tournament-style selection process.
    //!
    //! @see @link example_chemistry_molecule_evolutionary_optimizer.cpp @endlink
    //! @author geanesar,brownbp1
    //! @date Jun 18, 2022
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeEvolutionaryOptimizer :
      public util::SerializableInterface
      {

    private:

      //! @brief methods for how molecules should be selected
      enum SelectionType
      {
        e_Top,                  //!< Select the top members
        e_Tournament,           //!< Select members by tournament
        s_NumberSelectionTypes
      };

      //! @brief how models are read in (BCL internal or external script)
      enum ModelType
      {
        e_Internal,         //!< BCL models only
        e_External,         //!< Launch an external script
        s_NumberModelTypes
      };

      //! @brief how parent molecules are treated
      enum ParentRetirementType
      {
        e_DoNotRetire,              //!< Retirement policy "none", keep all parents
        e_RetireProbabilisticByAge, //!< Retirement policy "probabilistic", retire parents by age
        e_RetireAll,                //!< Retirement policy "all", do not keep parents
        s_NumberRetirementTypes
      };

      //! @brief how the molecule evolution is balanced
      //! InsertMols always steady at 15%
      enum EvolutionBalance
      {
        e_AlchemicalMutate,
        e_ReactionInsertionOnly,          //!< 100% MDL-style reaction
        e_ReactionDominant,               //!< 75% MDL Reaction, 10% Recombination
        e_RecombinationDominant,          //!< 10% MDL Reaction, 75% Recombination
        e_RecombinationInsertionOnly,     //!< 100% recombination
        e_ReactionRecombinationBalanced,  //!< 42.5% Reaction, 42.5% Recombination
        s_NumberEvolutionBalances
      };

      //! array for storing population information
      //! should contain m_FinalPopSize MoleculeEvolutionInfos per element
      std::vector< std::vector< MoleculeEvolutionInfo> > m_Populations;

      //! the maximum number of molecules to generate, per population
      size_t m_GenerateMax;

      //! the number of molecules to keep, per population
      size_t m_FinalPopSize;

      //! the number of times to try to generate a molecule before terminating early
      size_t m_MaxFailedAttempts;

      //! The reaction operation class, used for modifying structures
      FragmentReact m_ReactOp;

      //! The alchemical mutate to use to combine molecules
      util::Implementation< FragmentMutateInterface> m_Mutate; //!< obtains a implementation

      //! vector of molecules for addition/insertion operations
      util::ShPtrVector< FragmentComplete> m_InsertMols;

      //! molecule scorer, if model type uses BCL-only models
      descriptor::CheminfoProperty m_Scorer;

      //! coordinates for the center of the binding pocket
      mutable linal::Vector3D m_ReferenceCoordinates;

      //! tracker for how many times different operations are used
      mutable std::vector< size_t> m_Operations;

      //! how molecules should be selected
      //! @see SelectionType
      SelectionType m_SelectionType;

      //! size of tournaments for the replacement step, if using tournament selection
      float m_ReplacementTournSizeFactor;

      //! size of the tournaments for each modification step, if using tournament selection
      float m_ModifyTournSizeFactor;

      //! the type of models that will be used
      //! @see ModelType
      ModelType m_ModelType;

      //! the shell command that should be run if using external scoring functions
      std::string m_ModelCmd;

      //! parent retirement policy
      //! @see ParentRetirementType
      ParentRetirementType m_RetirementType;

      //! how molecule evolution processes are balanced
      //! @see EvolutionBalance
      EvolutionBalance m_EvolutionBalanceType;

      //! vector of molecules for recombination operations
      mutable FragmentEvolveImplementations m_RecombineOp;

      //! components to filter druglikeness by comparing two descriptors
      storage::Triplet
      <
        descriptor::CheminfoProperty,
        math::Comparisons< float>::Comparison,
        descriptor::CheminfoProperty
      > m_DruglikenessFilter;
      std::string m_DruglikenessFilterStr;
//      descriptor::CheminfoProperty m_DruglikenessFilterLHS;
//      descriptor::CheminfoProperty m_DruglikenessFilterRHS;
//      math::Comparisons< float>::Comparison m_DruglikenessComparisonType;

      //! filename for where EvoGen activity should be logged
      std::string m_LogFile;

      //! output stream for EvoGen logging information
//      io::OFStream m_LogStream;

    //////////
    // data //
    //////////

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MoleculeEvolutionaryOptimizer();

      //! @brief Clone function
      //! @return pointer to new MoleculeEvolutionaryOptimizer
      MoleculeEvolutionaryOptimizer *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the name of this class
      //! @return the name of this class
      const std::string &GetAlias() const;

      //! @brief Get the internal molecule data
      //! @return a vector containing population data for each iteration so far
      const std::vector< std::vector< MoleculeEvolutionInfo> > &GetMoleculeEvolutionInfos() const;

      //! @brief Get filename for EvoGen log file
      const std::string &GetLogFile() const;

      //! @brief Get the final size of each population.
      const size_t GetFinalPopSize() const;

      //! @brief Get the maximum number of molecules to generate during each iteration.
      const size_t GetMaxToGenerate() const;

      //! @brief Get molecule selection type to keep highest-scoring molecules
      const SelectionType &GetSelectionType() const;

      //! @brief Get molecule replacement method to use tournament selection
      const float GetReplacementTypeTournament() const;

      //! @brief Get molecule modification method to use tournament selection
      const float GetModifyTypeTournament() const;

      //! @brief Get retirement policy so that all parents are discarded
      const ParentRetirementType &GetRetirementType() const;

      //! @brief Get molecule evolution type to primarily perform balanced processes
      const EvolutionBalance &GetEvolutionBalanceType() const;

      //! @brief Use BCL-internal models for molecule scoring
      const ModelType &GetModelType() const;

      //! @brief Get the reaction operation structure
      const std::string &GetReactionOperationLabel() const;

      //! @brief Get the one-shot reaction operation structure
      const std::string &GetAlchemicalOperationLabel() const;

      //! @brief Get the descriptor to use for scoring, if internal models are used
      const std::string &GetModelDescriptorLabel() const;

      //! @brief Get the descriptor to use for scoring, if internal models are used
      const descriptor::CheminfoProperty &GetModelDescriptor() const;

      //! @brief Get the external script path if external scoring is used
      const std::string &GetModelCmd();

    ////////////////
    // operations //
    ////////////////

      //! @brief set filename for EvoGen log file
      //! @details this file is a json-formatted log file containing information about the generated molecules
      void SetLogFile( const std::string &FILENAME);

      //! @brief set up the initial population from an ensemble of molecules
      //! @details molecules must have >= 1 atoms or they will be discarded.
      //! @param MOLS the structures to use, will be copied
      void SetInitialPopulation( const FragmentEnsemble &MOLS);

      //! @brief set the final size of each population.
      //! @param POP_SIZE the final size of populations
      void SetFinalPopSize( const size_t &POP_SIZE);

      //! @brief set the maximum number of molecules to generate during each iteration.  This will be pruned to
      //!  m_FinalPopSize afterwards
      //! @param MAX maximum number of molecules to generate
      void SetMaxToGenerate( const size_t &MAX);

      //! @brief set molecule selection type to keep highest-scoring molecules
      void SetSelectionTop();

      //! @brief Set molecule replacement method to use tournament selection
      //! @param FACTOR the percentage of the available data that will be used in a single tournament round
      //! @details will warn if FACTOR < 0 or FACTOR > 1 and set FACTOR to 0.0 or 1.0, respectively.  The number of
      //!  molecules selected will be 1 <= NUMBER <= data size regardless of FACTOR
      void SetReplacementTypeTournament( const float &FACTOR);

      //! @brief Set molecule modification method to use tournament selection
      //! @param FACTOR the percentage of the available data that will be used in a single tournament round
      //! @details will warn if FACTOR < 0 or FACTOR > 1 and set FACTOR to 0.0 or 1.0, respectively.  The number of
      //!  molecules selected will be 1 <= NUMBER <= data size regardless of FACTOR
      void SetModifyTypeTournament( const float &FACTOR);

      //! @brief set retirement policy so that no parents are discarded
      void SetRetirementTypeNone();

      //! @brief set retirement policy so that parents are discarded probabilistically based on age
      void SetRetirementTypeProbabilistic();

      //! @brief set retirement policy so that all parents are discarded
      void SetRetirementTypeAll();

      //! @brief set molecule evolution type to primarily perform reactions
      void SetEvolutionBalanceAlchemicalMutate();

      void SetEvolutionBalanceReactionDominant();

      //! @brief set molecule evolution type to primarily perform reactions and insertions
      void SetEvolutionBalanceReactionInsertionOnly();

      //! @brief set molecule evolution type to primarily perform recombinations
      void SetEvolutionBalanceRecombinationDominant();

      //! @brief set molecule evolution type to primarily perform recombinations and insertions
      void SetEvolutionBalanceRecombinationInsertionOnly();

      //! @brief set molecule evolution type to primarily perform balanced processes
      void SetEvolutionBalancedBalanced();

      //! @brief use BCL-internal models for molecule scoring
      void SetModelTypeInternal();

      //! @brief use external script for molecule scoring
      void SetModelTypeExternal();

      //! @brief initialize the reaction operation structure
      //! @brief REACTANT_FILENAME filename from which to read reactant molecules
      //! @brief REACTION_DIRNAME directory in which RXN files should be found
      void SetupReactOperation( const std::string &REACTANT_FILENAME, const std::string &REACTION_DIRNAME);

      //! @brief initialize the one-shot reaction operation structure
      //! @brief IMPLEMENTATION alchemical mutate implementation
      void SetupAlchemicalMutate( const std::string &IMPLEMENTATION);

      //! @brief set the descriptor to use for scoring, if internal models are used
      //! @param DESCRIPTOR a string which should property encode a descriptor::CheminfoProperty
      //! @details this asserts that internal models are used to prevent programming errors
      void SetModelDescriptor( const std::string &DESCRIPTOR);

      //! @brief set the external script path if external scoring is used
      //! @param CMD the command to use.  Should be executed with the format: CMD <sdf from BCL> <output SDF>
      //! @details this asserts that internal models are used to prevent programming errors
      void SetModelCmd( const std::string &CMD);

      //! @brief set up the molecule insertion object
      //! @param INSERT_MOL_FILENAME sdf files which should be added to populations during a run
      void SetupInsertOperation( const std::string &INSERT_MOL_FILENAME);

      //! @brief select a molecule by tournament selection
      //! @param MOLS the population to select from
      //! @param TOURN_SIZE the tournament size (0 <= TOURN_SIZE <= MOLS.size())
      //! @param IGNORE_INDICES the indices to ignore when selecting
      //! @return an index in the MOLS vector, or MOLS.size() if an error occurs
      //! @details if TOURN_SIZE == 0 then the highest-scoring molecule will be used
      size_t SelectMoleculeTournament
      (
        const std::vector< MoleculeEvolutionInfo> &MOLS,
        size_t TOURN_SIZE,
        const std::set< int> &IGNORE_INDICES = std::set< int>()
      ) const;

      //! @brief score a vector of MoleculeEvolutionInfos
      //! @param MOL_INFOS the MoleculeEvolutionInfos containing molecules to score
      //! @details this updates the MoleculeEvolutionInfo.m_Fitness function of all molecules
      void ScoreMolecules( std::vector< MoleculeEvolutionInfo> &MOL_INFOS);

      //! @brief helper function for scoring a single molecule using internal BCL models
      //! @param MOL the molecule to score
      //! @return the score of the molecule
      //! @details does not do any preprocessing of the molecule, this is up to the caller
      float ScoreMoleculeInternal( const FragmentComplete &MOL) const;

      //! @brief helper function for scoring a single molecule using an external script
      //! @param MOL_INFOS the molecules to score
      //! @details this will update the MOL_INFOS structure directly with the new fitness scores as read from the
      //!  output of the called script
      void ScoreMoleculesExternal( std::vector< MoleculeEvolutionInfo> &MOL_INFOS) const;

      //! @brief copy the highest scoring molecules from one population to a number
      //! @param FROM the population to copy from
      //! @param TO the population to copy to
      //! @param NUM the number to copy
      void CopyHighestScoring( std::vector< MoleculeEvolutionInfo> &FROM, std::vector< MoleculeEvolutionInfo> &TO, size_t NUM);

      //! @brief picks molecules from FROM and copies them into TO
      //! @param FROM the population to copy molecules from
      //! @param TO the population to copy molecules to
      //! @param NUM the number to select
      //! @param TOURN_SIZE the tournament size to use
      //! @param REPLACEMENT whether to replace selected molecules (i.e. allow multiple picking)
      void CopyByTournament
      (
        const std::vector< MoleculeEvolutionInfo> &FROM,
        std::vector< MoleculeEvolutionInfo> &TO,
        size_t NUM,
        size_t TOURN_SIZE,
        const bool &REPLACEMENT = false
      ) const;

      //! @brief evaluates a descriptor to determine whether a molecule
      //! passes or fails the druglikeness filter
      bool EvaluateDruglikeness( const MoleculeEvolutionInfo &MOL) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief execute a reaction/addition operation to generate new molecules
      //! @brief MOLS the vector to add the new molecules to
      //! @details this only generates structures and sets history values appropriately, but does nothing for ensuring
      //!  molecules are suitable for scoring or that their
      void GenerateMolecules( std::vector< MoleculeEvolutionInfo> &MOLS) const;

      //! @brief execute an iteration of molecule generation/scoring/replacement
      //! @return 0 on success, negative value on error
      int Next();

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief open log file for writing; continues if file cannot be opened
//      void StartLogging();

      //! @brief close/flush logging file stream
//      void StopLogging();

      //! @brief remove whitespace (via isspace) from a string
      //! @param STR the string to remove whitespace from
      //! @return STR without any whitespace
      std::string RemoveWhitespace( const std::string &STR) const;

      //! @brief prepare a string for writing to CSV by escaping quotes
      std::string PrepareForCSV( const std::string &STR) const;

      //! @brief escape quotes with '\' in a string
      //! @param STR the string to escape
      //! @return a copy of STR with escaped quotes
      std::string EscapeQuotes( const std::string &STR) const;

      //! @brief write data to the json log file
      //! @param STR the string to write
//      void WriteLog( const std::string &STR);

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
    // input and output //
    //////////////////////

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream
      //! @param INDENT number of indentations
      //! @return ostream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      }; // class MoleculeEvolutionaryOptimizer

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_MOLECULE_EVOLUTIONARY_OPTIMIZER_H_
