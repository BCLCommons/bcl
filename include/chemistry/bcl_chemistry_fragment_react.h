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

#ifndef BCL_CHEMISTRY_FRAGMENT_REACT_H_
#define BCL_CHEMISTRY_FRAGMENT_REACT_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically
#include "find/bcl_find.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_conformational_interface.h"
#include "bcl_chemistry_fragment_complete.h"
#include "bcl_chemistry_fragment_constitution_shared.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "bcl_chemistry_reaction_search.h"
#include "bcl_chemistry_reaction_worker.h"
#include "find/bcl_find_pick_interface.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_result.h"
#include "sched/bcl_sched_mutex.h"
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
    //! @class FragmentReact
    //! @brief Perform reactions on molecules
    //!
    //! @see @link example_chemistry_fragment_react.cpp @endlink
    //! @author geanesar, brownbp1
    //! @date Jun 18, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentReact :
      public util::ObjectInterface
    {

    /////////////
    // friends //
    /////////////

    protected:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! debug output fragments as they match reactant positions
      bool m_OutputMatchedReagents;

      //! object that does the heavy lifting of fragment recombination during reaction
      mutable ReactionWorker m_ReactionWorker;

      //! reaction catalog to store reactants/reaction relationships
      mutable ReactionSearch m_ReactionSearch;

      //! the directory containing the reaction files
      std::string m_ReactionsDirectory;
      mutable ReactionEnsemble m_Reactions;

      //! the molecules to be reacted with the input molecules
      std::string m_ReagentsFilename;
      mutable FragmentEnsemble m_Reagents;

      //! only allow the input molecule to try matching with specific reactant positions
      std::string m_TargetReactantPositionsStr;
      mutable storage::Vector< size_t> m_TargetReactantPositions;

    public:

    //////////
    // data //
    //////////

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FragmentReact
      (
        const ReactionSearch &CATALOG = ReactionSearch()
      );

      //! @brief reaction worker constructor
      FragmentReact
      (
        const ReactionWorker &WORKER
      );

      //! @brief full constructor
      FragmentReact
      (
        const ReactionWorker &WORKER,
        const ReactionSearch &CATALOG
      );

      //! @brief clone constructor
      FragmentReact *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief returns the reaction search object
      ReactionSearch &GetReactionSearch() const;

      //! @brief returns the reaction worker object
      ReactionWorker &GetReactionWorker() const;

      //! @brief returns the reactions
      const ReactionEnsemble &GetReactions() const;

      //! @brief returns the reagents
      const FragmentEnsemble &GetReagents() const;

      //! @brief returns the allowed reactant positions for the target molecule
      const storage::Vector< size_t> &GetTargetReactantPositions() const;

      //! @brief return a static mutex object for optional debug output
      static sched::Mutex &GetMutex();

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      //! @brief set the reaction search object
      void SetReactionSearch( const ReactionSearch &REACTION_SEARCH);

      //! @brief set the reaction worker
      void SetReactionWorker( const ReactionWorker &REACTION_WORKER);

      //! @brief set the reactions
      void SetReactions( const ReactionEnsemble &REACTIONS);

      //! @brief set the reagents
      void SetReagents( const FragmentEnsemble &REAGENTS);

      //! @brief set the allowed reactant positions for the target molecule
      void SetTargetReactantPositions( const storage::Vector< size_t> &TARGET_REACTANT_POS);

      //! @brief return true if the reaction search is not empty
      bool IsGood() const
      {
        return !m_ReactionSearch.IsEmpty();
      }

      //! @brief react a molecule at a random reactive center
      storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble>
      ReactRandom
      ( 
        const FragmentComplete &MOL
      ) const;

      //! @brief react a molecule at all reactive centers with all reagents for a single specified reaction
      storage::Pair
      <
        storage::Vector< storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble> >,
            storage::Vector< storage::Vector< std::string> >
      >
      ReactExhaustiveOneReaction
      (
        const FragmentComplete &MOL,
        const ReactionComplete &RXN
              ) const;

      //! @brief react a molecule at all reactive centers with all reagents
      storage::Pair
      <
        storage::Vector< storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble> >,
        storage::Vector< storage::Vector< std::string> >
      >
      ReactExhaustive
      (
        const FragmentComplete &MOL
              ) const;

    protected:

      //! @brief worker function for three reactant reactions
      //! @param RXN the reaction to be performed
      //! @param MOL the molecule on which to perform the reaction
      //! @return returns a pair where the first element is a pair
      //! containing the full ReactionComplete and the ensemble of products,
      //! and the second element is a vector of strings corresponding to
      //! each product in the output ensemble indicating the reaction and reactants
      //! used to produce said product
      storage::Pair
      <
        storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble>,
        storage::Vector< std::string>
      >
      ReactExhaustiveOneReactantHelper
      (
        const util::SiPtr< const ReactionComplete> &RXN,
        const FragmentComplete &MOL
      ) const;

      //! @brief worker function for three reactant reactions
      //! @param RXN the reaction to be performed
      //! @param MOL the molecule on which to perform the reaction
      //! @return returns a pair where the first element is a pair
      //! containing the full ReactionComplete and the ensemble of products,
      //! and the second element is a vector of strings corresponding to
      //! each product in the output ensemble indicating the reaction and reactants
      //! used to produce said product
      storage::Pair
      <
        storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble>,
        storage::Vector< std::string>
      >
      ReactExhaustiveTwoReactantHelper
      (
        const util::SiPtr< const ReactionComplete> &RXN,
        const FragmentComplete &MOL
      ) const;

      //! @brief worker function for three reactant reactions
      //! @param RXN the reaction to be performed
      //! @param MOL the molecule on which to perform the reaction
      //! @return returns a pair where the first element is a pair
      //! containing the full ReactionComplete and the ensemble of products,
      //! and the second element is a vector of strings corresponding to
      //! each product in the output ensemble indicating the reaction and reactants
      //! used to produce said product
      storage::Pair
      <
        storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble>,
        storage::Vector< std::string>
      >
      ReactExhaustiveThreeReactantHelper
      (
        const util::SiPtr< const ReactionComplete> &RXN,
        const FragmentComplete &MOL
      ) const;

      //! @brief worker function for four reactant reactions
      //! @param RXN the reaction to be performed
      //! @param MOL the molecule on which to perform the reaction
      //! @return returns a pair where the first element is a pair
      //! containing the full ReactionComplete and the ensemble of products,
      //! and the second element is a vector of strings corresponding to
      //! each product in the output ensemble indicating the reaction and reactants
      //! used to produce said product
      storage::Pair
      <
        storage::Pair< util::SiPtr< const ReactionComplete>, FragmentEnsemble>,
        storage::Vector< std::string>
      >
      ReactExhaustiveFourReactantHelper
      (
        const util::SiPtr< const ReactionComplete> &RXN,
        const FragmentComplete &MOL
      ) const;

      //! @brief react a molecule with a center containing a specified set of atoms
      //! @param REQUIRE_ALL whether all atom indices must be in a reactive center (true), or
      //!        if only a single atom must be in a reactive center (false)
      FragmentEnsemble ReactCenterContaining
      ( 
        const FragmentComplete &MOL, 
        storage::Vector< size_t> &ATOM_INDICES,
        const bool &REQUIRE_ALL = true
      ) const;

      //! @brief remove whitespace (via isspace) from a string
      //! @param STR the string to remove whitespace from
      //! @return STR without any whitespace
      std::string RemoveWhitespace( const std::string &STR) const;

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

    }; // class FragmentReact

  } // namespace chemistry
} // namespace bcl
#endif //BCL_CHEMISTRY_FRAGMENT_REACT_H_
