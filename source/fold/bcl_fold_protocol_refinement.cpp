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
#include "fold/bcl_fold_protocol_refinement.h"

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
    ProtocolRefinement::ProtocolRefinement()
    {
    }

    //! @brief Clone function
    //! @return pointer to new ProtocolRefinement
    ProtocolRefinement *ProtocolRefinement::Clone() const
    {
      return new ProtocolRefinement( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief access to only instance
    //! @return reference to only instance
    ProtocolRefinement &ProtocolRefinement::GetInstance()
    {
      static ProtocolRefinement s_protocol_instance;
      return s_protocol_instance;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProtocolRefinement::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

     //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ProtocolRefinement::GetAlias() const
    {
      static const std::string s_name( "ProtocolRefinement");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProtocolRefinement::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "refinement protocol that uses only small-scale moves");

      return serializer;
    }

  ///////////
  // flags //
  ///////////

    //! @brief returns all flags that are specialized for this protocol
    //! @return all flags that are specialized for this protocol
    const util::ShPtrVector< command::FlagInterface> &ProtocolRefinement::GetAllFlags() const
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
    void ProtocolRefinement::ModifyStartModel( assemble::ProteinModel &START_MODEL) const
    {
    }

    //! @brief initialize the scores and add them to Scores enumerator
    void ProtocolRefinement::InitializeScores()
    {
    }

    //! @brief modify the score weight set
    //! @param SCORE_WEIGHT_SET Score weight set
    void ProtocolRefinement::ModifyScoreWeightSet( ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
    }

    //! @brief modify the terminate object
    //! @param CRITERION which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolRefinement::ModifyCriterion
    (
      opti::CriterionCombine< assemble::ProteinModel, double> &CRITERION,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief modify the printer object
    //! @param PRINTER which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolRefinement::ModifyPrinter
    (
      mc::PrinterCombined< assemble::ProteinModel, double> &PRINTER,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief modify the pdb factory object
    //! @param FACTORY pdb factory to be modified
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolRefinement::ModifyFactory
    (
      util::ShPtr< pdb::Factory> &FACTORY,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief initialize the mutates and add them to Mutates enumerator
    void ProtocolRefinement::InitializeMutates()
    {
    }

    //! @brief modify the mutate tree used
    //! @param MUTATE_TREE MutateTree to be modified
    void ProtocolRefinement::ModifyMutateTree( MutateTree &MUTATE_TREE) const
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
      const double avg_helix_ratio( avg_helix_count / avg_sse_count);
      const double avg_strand_ratio( avg_strand_count / avg_sse_count);

      BCL_MessageStd
      (
        "average helix/strand/sse counts from pool :" + util::Format()( avg_helix_count) + " / " +
         util::Format()( avg_strand_count) + " / " + util::Format()( avg_sse_count) + " / "
      );

      // initialize variables that are based on the counts
      const bool has_helix( avg_helix_count > 0);
      const bool has_strand( avg_strand_count > 0);
      const bool has_helixdomain( avg_helix_count >= 1.5);
      const bool has_helixpair( avg_helix_count >= 1.5);
      const bool has_sheet( avg_strand_count >= 1.5);

      // first choose probabilities for add, remove, swap, move
      MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_Add,         0.00);
      MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_Remove,      0.00);
      MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_Swap,        0.00);
      MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_SSE,         0.10);
      MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_Helix,       0.20 * avg_helix_ratio);
      MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_Strand,      0.20 * avg_strand_ratio);
      const double ssepair_weight(                                     0.20);
      MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_SSEPair,     ssepair_weight);
      const double domain_weight(                                      0.30);

      // if has a helix pair
      if( has_helixpair)
      {
        // update helix pair weight
        MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_HelixPair, ssepair_weight * avg_helix_ratio);
        MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_SSEPair,   ssepair_weight * avg_strand_ratio);
      }

      // update domain weights
      if( has_helixdomain)
      {
        MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_HelixDomain, domain_weight * avg_helix_ratio);
      }
      if( has_sheet)
      {
        MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_Sheet,       domain_weight * avg_strand_ratio);
      }

      // single SSE moves
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSEBendRamachandran,    1.0);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSEBendRandomSmall,    1.0);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSEBendRandomLarge,    1.0);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSETranslateSmall,      3.0);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSETranslateXSmall,    1.0);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSETranslateYSmall,    1.0);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSETranslateZSmall,    1.0);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSERotateSmall,         3.0);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSERotateXSmall,       1.0);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSERotateYSmall,       1.0);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSERotateZSmall,       1.0);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSETransformSmall,      6.0);

      // if secondary structure prediction were provided
      if
      (
        sspred::Methods::GetFlagReadSSPredictions()->GetFlag() &&
        DefaultFlags::GetFlagEnableSSEResize()->GetFlag()
      )
      {
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSEResizeCTerm, 1.5);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSEResizeNTerm, 1.5);
      }

      // single helix moves
      if( has_helix)
      {
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixTranslateXYSmall,     1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixTranslateZSmall,      1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixRotateXYSmall,        1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixRotateZSmall,         1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixTransformXYSmall,     1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixTransformZSmall,      1.0);
      }

      // single strand moves
      if( has_strand)
      {
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_StrandTranslateZSmall,     1.0);
      }

      // sse pair moves
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSEPairTranslateNoHingeSmall,    1.0);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSEPairTranslateSmall,             1.0);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSEPairRotateSmall,                1.0);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSEPairTransformSmall,             1.0);

      // helix pair moves
      if( has_helixpair)
      {
        // helix pair moves
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixPairRotateZSmallNoHinge,    1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixPairRotateZSmallHinge,       1.0);
      }

      // helix domain moves
      if( has_helixdomain)
      {
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixDomainTranslateSmall,  1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixDomainRotateSmall,     1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixDomainTransformSmall,  1.0);
      }

      // sheet moves
      if( has_sheet)
      {
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SheetRotateSmall,           1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SheetTranslateSmall,        1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SheetTransformSmall,        1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_Sheet_twist_small,            1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_Sheet_register_fix,           1.0);
      }
    }

    //! @brief get the mutate tree associated with this protocol
    //! @return the mutate tree associated with this protocol
    util::ShPtr< MutateTree> ProtocolRefinement::GetMutateTree() const
    {
      util::ShPtr< MutateTree> sp_mutate_tree( new MutateTree());
      ModifyMutateTree( *sp_mutate_tree);
      return sp_mutate_tree;
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void ProtocolRefinement::MergeAndModifyMutateTree( MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }

  ////////////
  // readme //
  ////////////

    //! @brief returns string containing short description of the protocol
    //! @return string containing short description of the protocol
    const std::string &ProtocolRefinement::GetDescription() const
    {
      // initialize static const string variable to hold the description
      static const std::string s_description
      (
        "refinement protocol for de novo folding"
      );

      // end
      return s_description;
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ProtocolRefinement::GetReadMe() const
    {
      // initialize static const string variable to hold the readme
      static const std::string s_readme
      (
        "The refinement protocol focuses on small amplitude moves that maintain the current topology but optimize "
        "interactions between SSEs and introduce bends into SSEs."
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
    std::istream &ProtocolRefinement::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProtocolRefinement::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
