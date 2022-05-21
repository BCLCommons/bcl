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
#include "fold/bcl_fold_protocol_assembly.h"

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
    ProtocolAssembly::ProtocolAssembly()
    {
    }

    //! @brief Clone function
    //! @return pointer to new ProtocolAssembly
    ProtocolAssembly *ProtocolAssembly::Clone() const
    {
      return new ProtocolAssembly( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief access to only instance
    //! @return reference to only instance
    ProtocolAssembly &ProtocolAssembly::GetInstance()
    {
      static ProtocolAssembly s_protocol_instance;
      return s_protocol_instance;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProtocolAssembly::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ProtocolAssembly::GetAlias() const
    {
      static const std::string s_name( "ProtocolAssembly");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProtocolAssembly::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Protocol that modifies the default protocol to use only large-scale mutates");

      return serializer;
    }

  ///////////
  // flags //
  ///////////

    //! @brief returns all flags that are specialized for this protocol
    //! @return all flags that are specialized for this protocol
    const util::ShPtrVector< command::FlagInterface> &ProtocolAssembly::GetAllFlags() const
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
    void ProtocolAssembly::ModifyStartModel( assemble::ProteinModel &START_MODEL) const
    {
    }

    //! @brief initialize the scores and add them to Scores enumerator
    void ProtocolAssembly::InitializeScores()
    {
    }

    //! @brief modify the score weight set
    //! @param SCORE_WEIGHT_SET Score weight set
    void ProtocolAssembly::ModifyScoreWeightSet( ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
    }

    //! @brief modify the terminate object
    //! @param CRITERION which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolAssembly::ModifyCriterion
    (
      opti::CriterionCombine< assemble::ProteinModel, double> &CRITERION,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief modify the printer object
    //! @param PRINTER which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolAssembly::ModifyPrinter
    (
      mc::PrinterCombined< assemble::ProteinModel, double> &PRINTER,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief modify the pdb factory object
    //! @param FACTORY pdb factory to be modified
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolAssembly::ModifyFactory
    (
      util::ShPtr< pdb::Factory> &FACTORY,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief initialize the mutates and add them to Mutates enumerator
    void ProtocolAssembly::InitializeMutates()
    {
    }

    //! @brief modify the mutate tree used
    //! @param MUTATE_TREE MutateTree to be modified
    void ProtocolAssembly::ModifyMutateTree( MutateTree &MUTATE_TREE) const
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
      MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_Add,         0.075);
      MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_Remove,      0.025);
      MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_Swap,        0.10);
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

      // adds
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_AddSSENextToSSE     , 0.75);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_AddSSEShortLoop     , 0.25);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_AddStrandNextToSheet, 1.00 * avg_strand_ratio);
      // remove
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_RemoveRandom        , 0.5);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_RemoveUnpairedStrand, 0.5);

      // swap
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SwapSSEs       , 0.8);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SwapSSEWithPool, 0.2);
      // if the pool is overlapping
      if( sp_pool->IsOverlapping())
      {
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SwapSSEWithPoolOverlap, 0.2);
      }

      // single SSE moves
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSEMoveNext        , 3.0);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSEMoveShortLoop   , 3.0);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSEFurthestMoveNext, 3.0);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSETranslateLarge  , 3.0);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSETranslateXLarge , 1.0);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSETranslateYLarge , 1.0);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSETranslateZLarge , 1.0);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSERotateLarge     , 3.0);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSERotateXLarge    , 1.0);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSERotateYLarge    , 1.0);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSERotateZLarge    , 1.0);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSETransformLarge  , 6.0);

      // if secondary structure prediction were provided and resizing was not disabled
      if
      (
        sspred::Methods::GetFlagReadSSPredictions()->GetFlag() &&
        DefaultFlags::GetFlagEnableSSEResize()->GetFlag()
      )
      {
        // average number of non overlapping sses
        const double number_non_overlapping_long_sses( sp_pool->CountLongNonOverlappingSSEs( 30, 20));
        const double split_probability( ( 2.0 + 4.0 * number_non_overlapping_long_sses) / DefaultMutates::GetInstance().m_SSESplit.GetSize());

        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSEResizeCTerm, 1.5);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSEResizeNTerm, 1.5);

        // iterate over ss prediction methods
        for
        (
          storage::Map< sspred::Method, Mutate>::const_iterator
            sspred_itr( DefaultMutates::GetInstance().m_SSESplit.Begin()),
            sspred_itr_end( DefaultMutates::GetInstance().m_SSESplit.End());
          sspred_itr != sspred_itr_end;
          ++sspred_itr
        )
        {
          MUTATE_TREE.SetMutateProbability( sspred_itr->second, split_probability);
        }
      }

      // single helix moves
      if( has_helix)
      {
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixMoveNext,         2.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixMoveShortLoop,    2.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixFurthestMoveNext, 2.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixTranslateXYLarge, 2.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixTranslateZLarge,  2.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixRotateXYLarge,    2.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixRotateZLarge,     2.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixTransformXYLarge, 2.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixTransformZLarge,  2.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixFlipXY,           2.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixFlipZ,            2.0);
      }

      // single strand moves
      if( has_strand)
      {
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_StrandMoveNext,          2.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_StrandFurthestMoveNext,  2.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_StrandMoveSheet,         2.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_StrandFurthestMoveSheet, 2.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_StrandTranslateZLarge,   2.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_StrandFlipX,             2.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_StrandFlipY,             2.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_StrandFlipZ,             2.0);
      }

      // sse pair moves
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSEPairTranslateNoHingeLarge, 2.0);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSEPairTranslateLarge,        2.0);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSEPairRotateLarge,           2.0);
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSEPairTransformLarge,        2.0);

      // helix pair moves
      if( has_helixpair)
      {
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixPairRotateZLargeNoHinge, 1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixPairRotateZLargeHinge,   1.0);
      }

      // helix domain moves
      if( has_helixdomain)
      {
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixDomainShuffle,         8.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixDomainTranslateLarge, 2.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixDomainRotateLarge,    2.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixDomainTransformLarge, 2.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixDomainFlipExt,        2.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_HelixDomainFlipInt,        2.0);
      }

      // sheet moves
      if( has_sheet)
      {
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SheetRotateLarge,        1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SheetTranslateLarge,     1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SheetTransformLarge,     1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SheetPairStrands,        1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SheetSwitchStrand,       1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SheetFlipExt,            1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SheetFlipInt,            1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_Sheet_flip_int_sub,        1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_Sheet_flip_int_sub_diff,   1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_Sheet_divide_sandwich,     1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_Sheet_divide,              1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_Sheet_twist_large,         1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_Sheet_shuffle,             4.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_Sheet_cycle,               1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_Sheet_cycle_intact,        1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_Sheet_cycle_subset,        1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_Sheet_cycle_subset_intact, 1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_Sheet_register_shift,      1.0);
        MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_Sheet_register_shift_flip, 1.0);
      }
    }
    //! @brief get the mutate tree associated with this protocol
    //! @return the mutate tree associated with this protocol
    util::ShPtr< MutateTree> ProtocolAssembly::GetMutateTree() const
    {
      util::ShPtr< MutateTree> sp_mutate_tree( new MutateTree());
      ModifyMutateTree( *sp_mutate_tree);
      return sp_mutate_tree;
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void ProtocolAssembly::MergeAndModifyMutateTree( MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }

  ////////////
  // readme //
  ////////////

    //! @brief returns string containing short description of the protocol
    //! @return string containing short description of the protocol
    const std::string &ProtocolAssembly::GetDescription() const
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
    const std::string &ProtocolAssembly::GetReadMe() const
    {
      // initialize static const string variable to hold the readme
      static const std::string s_readme
      (
        "The assembly protocol consists of large amplitude translation or rotations and moves that add or remove SSEs. "
        "Other moves central to this phase shuffle beta-strands within beta-sheets or break large beta-sheets to "
        "create beta-sandwiches."
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
    std::istream &ProtocolAssembly::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProtocolAssembly::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
