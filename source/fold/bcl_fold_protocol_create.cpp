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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "fold/bcl_fold_protocol_create.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_pool.h"
#include "fold/bcl_fold_default_flags.h"
#include "fold/bcl_fold_default_mutates.h"
#include "fold/bcl_fold_mutate_tree.h"
#include "fold/bcl_fold_setup.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProtocolCreate::ProtocolCreate()
    {
    }

    //! @brief Clone function
    //! @return pointer to new ProtocolCreate
    ProtocolCreate *ProtocolCreate::Clone() const
    {
      return new ProtocolCreate( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief access to only instance
    //! @return reference to only instance
    ProtocolCreate &ProtocolCreate::GetInstance()
    {
      static ProtocolCreate s_protocol_instance;
      return s_protocol_instance;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProtocolCreate::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ProtocolCreate::GetAlias() const
    {
      static const std::string s_name( "ProtocolCreate");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProtocolCreate::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Initial protocol for de-novo folding. Gets all SSEs into the model w/o clashes or loop-closure issues");

      return serializer;
    }

  ///////////
  // flags //
  ///////////

    //! @brief returns all flags that are specialized for this protocol
    //! @return all flags that are specialized for this protocol
    const util::ShPtrVector< command::FlagInterface> &ProtocolCreate::GetAllFlags() const
    {
      // construct static vector to hold the flags
      static const util::ShPtrVector< command::FlagInterface> s_flags;

      // end
      return s_flags;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief modifies the start model
    //! @param START_MODEL Start model to be modified
    void ProtocolCreate::ModifyStartModel( assemble::ProteinModel &START_MODEL) const
    {
    }

    //! @brief initialize the scores and add them to Scores enumerator
    void ProtocolCreate::InitializeScores()
    {
    }

    //! @brief modify the score weight set
    //! @param SCORE_WEIGHT_SET Score weight set
    void ProtocolCreate::ModifyScoreWeightSet( ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
    }

    //! @brief modify the terminate object
    //! @param CRITERION which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolCreate::ModifyCriterion
    (
      opti::CriterionCombine< assemble::ProteinModel, double> &CRITERION,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief modify the printer object
    //! @param PRINTER which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolCreate::ModifyPrinter
    (
      mc::PrinterCombined< assemble::ProteinModel, double> &PRINTER,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief modify the pdb factory object
    //! @param FACTORY pdb factory to be modified
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolCreate::ModifyFactory
    (
      util::ShPtr< pdb::Factory> &FACTORY,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief initialize the mutates and add them to Mutates enumerator
    void ProtocolCreate::InitializeMutates()
    {
    }

    //! @brief modify the mutate tree used
    //! @param MUTATE_TREE MutateTree to be modified
    void ProtocolCreate::ModifyMutateTree( MutateTree &MUTATE_TREE) const
    {
      // reset MutateTree
      MUTATE_TREE.Reset();

      // get the pool using the empty model
      util::ShPtr< assemble::SSEPool> sp_pool
      (
        GetSetup().GetEmptyModel()->GetProteinModelData()->GetData( assemble::ProteinModelData::e_Pool)
      );
      BCL_Assert( sp_pool.IsDefined(), "pool is not initialized!");

      // get the average number of helices and strands, and calculate the ratios
      const storage::Pair< double, double> avg_sse_counts( sp_pool->CalculateAverageHelixStrandCounts());
      const double avg_helix_count( avg_sse_counts.First());
      const double avg_strand_count( avg_sse_counts.Second());
      const double avg_sse_count( avg_helix_count + avg_strand_count);
      const double avg_strand_ratio( avg_strand_count / avg_sse_count);

      // first choose probabilities for add, remove, swap, move
      MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_Add,         0.75);
      MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_Remove,      0.25);
      MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_Swap,        0.25);
      MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_SSE,         0.01);
      MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_Helix,       0.01);
      MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_Strand,      0.01);
      MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_SSEPair,     0.01);
      MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_HelixPair,   0.01);
      MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_HelixDomain, 0.01);
      MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_Sheet,       0.01);

      // adds
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_AddSSENextToSSE     , 0.375);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_AddSSEShortLoop     , 0.5);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_AddStrandNextToSheet, 1.00 * avg_strand_ratio);

      // remove
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_RemoveRandom   , 1.00);

      // swap
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SwapSSEs       , 1.00);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SwapSSEWithPool, 1.00);
    }
    //! @brief get the mutate tree associated with this protocol
    //! @return the mutate tree associated with this protocol
    util::ShPtr< MutateTree> ProtocolCreate::GetMutateTree() const
    {
      util::ShPtr< MutateTree> sp_mutate_tree( new MutateTree());
      ModifyMutateTree( *sp_mutate_tree);
      return sp_mutate_tree;
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void ProtocolCreate::MergeAndModifyMutateTree( MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }

  ////////////
  // readme //
  ////////////

    //! @brief returns string containing short description of the protocol
    //! @return string containing short description of the protocol
    const std::string &ProtocolCreate::GetDescription() const
    {
      // initialize static const string variable to hold the description
      static const std::string s_description
      (
        "assembly protocol for de novo folding"
      );

      // end
      return s_description;
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ProtocolCreate::GetReadMe() const
    {
      // initialize static const string variable to hold the readme
      static const std::string s_readme
      (
        "This protocol consists solely of adding sses and swapping them them those already in the pool."
      );

      // end
      return s_readme;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProtocolCreate::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProtocolCreate::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
