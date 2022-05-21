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

#ifndef BCL_FOLD_MUTATE_TREE_H_
#define BCL_FOLD_MUTATE_TREE_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_serializable_interface.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateTree
    //! @brief Convenience class to build up two layer tree with different mutates and corresponding
    //! probabilities
    //! @details This class allows building up a tree with mutates and associated probabilities which can later
    //! be converted to a MutateDecisionNode in order to be used in a minimization process such as BCL::Fold.
    //! It allows inserting new mutates, setting probabilities for mutates or resetting all of them. It stores
    //! the mutates categorized by the mutate types, e_Add, e_Remove, e_Swap and e_Move. It also has two
    //! convenience functions: CreateTable() for creating a weight table and ConstructMutates() for building up
    //! the actual mutate tree with MutateDecisionNode.
    //!
    //! @see @link example_fold_mutate_tree.cpp @endlink
    //! @author karakam, fischea
    //! @date Nov 18, 2010
    //!
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API MutateTree :
      public virtual util::SerializableInterface
    {

    public:

      //! @enum MutateTypes
      //! enumeration of types of mutates
      enum MutateTypes
      {
        e_Add,
        e_Remove,
        e_Swap,
        e_SSE,
        e_Helix,
        e_Strand,
        e_SSEPair,
        e_HelixPair,
        e_HelixDomain,
        e_Sheet,
        e_Domain,
        e_Model,
        e_Ensemble,
        s_NumberMutateTypes
      };

      //! @brief conversion to a string from a MutateTypes
      //! @param MUTATE_TYPE the mutate type to get a string for
      //! @return a string representing that mutate type
      static const std::string &GetMutateTypeName( const MutateTypes &MUTATE_TYPE);

      //! @brief enum class wrapper for MutateTypes
      typedef util::WrapperEnum< MutateTypes, &GetMutateTypeName, s_NumberMutateTypes> MutateTypesEnum;

    private:

      //! map that stores the probabilities for each mutate type
      storage::Map< MutateTypesEnum, double> m_MutateTypeProbabilities;

      //! map that stores whether or not each mutate type has been set
      storage::Map< MutateTypesEnum, bool> m_MutateTypeSet;

      //! map that stores for each mutate type, list of corresponding mutates and the associated probabilities
      storage::Map< MutateTypesEnum, storage::Map< Mutate, double> > m_MutateProbabilities;

      //! map that stores whether or not each mutate has been set
      storage::Map< Mutate, bool> m_MutateSet;

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutateTree();

      //! @brief constructor from a weight table
      //! @param WEIGHT_TABLE table with mutate weights
      MutateTree( const storage::Table< double> &WEIGHT_TABLE);

      //! @brief Clone function
      //! @return pointer to new MutateTree
      MutateTree *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief get the map that stores the probabilities for each mutate type
      //! @return map that stores the probabilities for each mutate type
      const storage::Map< MutateTypesEnum, double> &GetMutateTypeProbabilities() const
      {
        return m_MutateTypeProbabilities;
      }

      //! @brief get the map that stores for each mutate type, list of corresponding mutates and the associated
      //! probabilities
      //! @return map that stores for each mutate type, list of corresponding mutates and the associated
      //! probabilities
      const storage::Map< MutateTypesEnum, storage::Map< Mutate, double> > &GetMutateProbabilities() const
      {
        return m_MutateProbabilities;
      }

      //! @brief create a table from mutate tree
      //! @return table created from the mutate tree
      storage::Table< double> CreateTable() const;

      //! @brief creates the mutate function using all the mutates in the mutate tree
      //! @return the mutate function constructed from all the mutates in the mutate tree
      util::ShPtr< math::MutateInterface< assemble::ProteinModel> > ConstructMutates() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief reset all the information in the tree
      void Reset();

      //! @brief sets the probability for a given mutate enum
      //! @param MUTATE_ENUM enum of the mutate of interest
      //! @param PROBABILITY probability of the mutate to be assigned
      void SetMutateProbability( const Mutate &MUTATE_ENUM, const double PROBABILITY);

      //! @brief sets the default probability for a given mutate enum
      //! @param MUTATE_ENUM enum of the mutate of interest
      //! @param PROBABILITY probability of the mutate to be assigned
      void SetDefaultMutateProbability( const Mutate &MUTATE_ENUM, const double PROBABILITY);

      //! @brief gets the probability for a given mutate enum
      //! @param MUTATE_ENUM enum of the mutate of interest
      //! @return probability of the mutate
      double GetMutateProbability( const Mutate &MUTATE_ENUM) const;

      //! @brief sets the probability for a given mutate
      //! @param MUTATE_NAME Name of the mutate of interest
      //! @param PROBABILITY probability of the mutate to be assigned
      void SetMutateProbability( const std::string &MUTATE_NAME, const double PROBABILITY);

      //! @brief get the probability of a given mutate type
      //! @param MUTATE_TYPE mutate type of interest
      //! return probability of mutate type
      double GetMutateTypeProbability( const MutateTypes &MUTATE_TYPE) const;

      //! @brief set the probabilities for a specific mutate type
      //! @param MUTATE_TYPE MutateTypes of interest
      //! @param PROBABILITY Probability to set to
      void SetMutateTypeProbability( const MutateTypes &MUTATE_TYPE, const double PROBABILITY);

      //! @brief set the default probabilities for a specific mutate type
      //! @param MUTATE_TYPE MutateTypes of interest
      //! @param PROBABILITY Probability to set to
      void SetDefaultMutateTypeProbability( const MutateTypes &MUTATE_TYPE, const double PROBABILITY);

      //! @brief get whether or not a given mutate type has been set
      //! @param MUTATE_TYPE mutate type of interest
      //! @return whether or not this mutate type has been set
      bool GetMutateTypeSet( const MutateTypes &MUTATE_TYPE) const;

      //! @brief set the mutate type to having been set
      //! @param MUTATE_TYPE mutate type of interest
      //! @param SET whether or not it has been set
      void SetMutateTypeSet( const MutateTypes &MUTATE_TYPE, const bool SET);

      //! @brief get whether or not the probability of the given mutate type is zero
      //! @param MUTATE_TYPE mutate type of interest
      //! @return whether or not the probability is zero
      bool GetMutateTypeZero( const MutateTypes &MUTATE_TYPE) const;

      //! @brief get whether or not a given mutate has been set
      //! @param MUTATE mutate of interest
      //! @return whether or not this mutate has been set
      bool GetMutateSet( const Mutate &MUTATE) const;

      //! @brief set the mutate to having been set
      //! @param MUTATE mutate of interest
      //! @param SET whether or not it has been set
      void SetMutateSet( const Mutate &MUTATE, const bool SET);

      //! @brief get whether or not the probability of the given mutate is zero
      //! @param MUTATE mutate of interest
      //! @return whether or not the probability is zero
      bool GetMutateZero( const Mutate &MUTATE) const;

      //! @brief parses the given mutate name and returns the mutate type
      //! @param MUTATE_NAME name of the mutate
      //! @return corresponding mutate type
      static MutateTypes ParseMutateType( const std::string &MUTATE_NAME);

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief merge two mutate trees
      void Merge( const MutateTree &MUTATE_TREE);

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief initializes the mutate tree using a given table
      //! @param WEIGHT_TABLE table that contains the weights for mutates
      void InitializeFromTable( const storage::Table< double> &WEIGHT_TABLE);

    }; // class MutateTree

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_TREE_H_
