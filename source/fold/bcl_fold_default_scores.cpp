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
#include "fold/bcl_fold_default_scores.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_count.h"
#include "fold/bcl_fold_score_weight_set.h"
#include "score/bcl_score_aa_neighborhood_exposure.h"
#include "score/bcl_score_aa_pair_clash.h"
#include "score/bcl_score_aa_pair_contact_energy.h"
#include "score/bcl_score_aa_pair_distance.h"
#include "score/bcl_score_aa_pair_hi_res_clash.h"
#include "score/bcl_score_aa_sequence_pair.h"
#include "score/bcl_score_contact_order.h"
#include "score/bcl_score_loop.h"
#include "score/bcl_score_loop_angle.h"
#include "score/bcl_score_loop_closure.h"
#include "score/bcl_score_phi_psi.h"
#include "score/bcl_score_protein_model_aa_neighborhood.h"
#include "score/bcl_score_protein_model_inverted.h"
#include "score/bcl_score_protein_model_sse.h"
#include "score/bcl_score_protein_model_sse_chirality.h"
#include "score/bcl_score_protein_model_sse_completeness.h"
#include "score/bcl_score_protein_model_sse_linear_loop_proximity.h"
#include "score/bcl_score_protein_model_sse_neighbors.h"
#include "score/bcl_score_protein_model_sse_packing.h"
#include "score/bcl_score_protein_model_sse_pairs.h"
#include "score/bcl_score_protein_model_wrapper.h"
#include "score/bcl_score_radius_of_gyration.h"
#include "score/bcl_score_sse_pair_clash.h"
#include "score/bcl_score_sse_pair_gap.h"
#include "score/bcl_score_sse_pair_packing.h"
#include "score/bcl_score_sse_pairs_fragments.h"
#include "score/bcl_score_sse_predictions.h"
#include "score/bcl_score_strand_pairing.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> DefaultScores::s_Instance
    (
      GetObjectInstances().AddInstance( new DefaultScores())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    DefaultScores::DefaultScores() :
      m_ProteinInverter( new assemble::ProteinModelInverter( true))
    {
    }

    //! @brief Clone function
    //! @return pointer to new DefaultScores
    DefaultScores *DefaultScores::Clone() const
    {
      return new DefaultScores( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief acces to static instance of this class
    DefaultScores &DefaultScores::GetInstance()
    {
      static DefaultScores s_default_scores_instance;
      return s_default_scores_instance;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DefaultScores::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the protein inverter
    //! @return the protein inverter
    const util::ShPtr< assemble::ProteinModelInverter> &DefaultScores::GetProteinInverter() const
    {
      return m_ProteinInverter;
    }

    //! @brief get the sspred scores
    //! the sspred scores
    const storage::Map< sspred::Method, storage::VectorND< 2, Score> > &DefaultScores::GetSSPredScores() const
    {
      return m_SSPredScores;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initialize the scores and add them to Scores enumerator
    void DefaultScores::InitializeScores()
    {
      // was already initialized?
      if( e_ScoreAAPairClash.IsDefined())
      {
        return;
      }

      // aa distance clash
      e_ScoreAAPairClash = GetScores().AddScore
      (
//        util::ShPtr< score::ProteinModel>
//        (
//          new score::ProteinModelSSEPairs
//          (
//            Scores::WrapCacheSSEPairScore
//            (
//              util::CloneToShPtr( score::AASequencePair( score::AAPairClash(), false)),
//              false // symmetric
//            ),
//            false, // no normalization
//            score::ProteinModel::e_Sequence,
//            "AA pair distance"
//          )
//        )
        util::ShPtr< score::ProteinModel>
        (
          new score::AAPairClash()
        )
//        util::ShPtr< score::ProteinModel>
//        (
//          new score::ProteinModelSSEPairs
//          (
//            fold::Scores::WrapCacheSSEPairScore( util::CloneToShPtr( score::AAPairClash()), false),
//            false, // no normalization
//            score::ProteinModel::e_Sequence,
//            "AA pair distance"
//          )
//        )
      );

      // aa distance score
      e_ScoreAAPairDistance = GetScores().AddScore
      (
        util::ShPtr< score::ProteinModel>
        (
          new score::AAPairDistance()
        )
      );

      // aa_clash_hires
      e_ScoreAAPairHiResClash = GetScores().AddScore
      (
        Scores::WrapCacheProteinModelScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::AAPairHiResClash()
          )
        )
      );

      // aa_pair_interaction
      e_ScoreAAPairSCInteraction = GetScores().AddScore
      (
        util::ShPtr< score::ProteinModel>
        (
          new score::AAPairContactEnergy()
        )
      );

      // aaneigh score
      e_ScoreAANeighborCount = GetScores().AddScore
      (
        util::ShPtr< score::ProteinModel>
        (
            new score::ProteinModelAANeighborhood
            (
              util::CloneToShPtr( score::AANeighborhoodExposure( assemble::AANeighborCount())),
              score::ProteinModelAANeighborhood::e_None,
              true,
              score::ProteinModel::e_Sequence,
              "AA environment"
            )
        )
      );

      // aaneigh entropy score
      e_ScoreAANeighborCountEntropy = GetScores().AddScore
      (
        util::ShPtr< score::ProteinModel>
        (
          new score::ProteinModelInverted
          (
            Scores::WrapCacheProteinModelScore
            (
              util::ShPtr< score::ProteinModel>
              (
                  new score::ProteinModelAANeighborhood
                  (
                    util::CloneToShPtr( score::AANeighborhoodExposure( assemble::AANeighborCount())),
                    score::ProteinModelAANeighborhood::e_None, // normalize
                    false // consider different chain
                  )
              )
            ),
            m_ProteinInverter,
            assemble::AANeighborCount::GetDefaultScheme() + "_ent",
            score::ProteinModel::e_Sequence,
            "AA environment"
          )
        )
      );

      // loop score; faster to compute than to cache
      e_ScoreLoop = GetScores().AddScore
      (
        util::ShPtr< score::ProteinModel>
        (
          new score::ProteinModelSSEPairs
          (
            Scores::WrapCacheSSEPairScore( util::CloneToShPtr( score::Loop()), false),
            true,
            score::ProteinModel::e_Structure,
            "Loop closure"
          )
        )
      );

      // loop angle score
      e_ScoreLoopAngle = GetScores().AddScore
      (
        util::ShPtr< score::ProteinModel>( new score::LoopAngle())
      );

      // loop closure score; computation is faster than caching, so no wrapping with cache
      e_ScoreLoopClosure = GetScores().AddScore
      (
        util::ShPtr< score::ProteinModel>
        (
          new score::ProteinModelWrapper
          (
            util::CloneToShPtr
            (
              score::ProteinModelSSENeighbors
              (
                util::CloneToShPtr( score::LoopClosure( 1, 1.0, 1.0)),
                false // no normalization
              )
            ),
            score::ProteinModel::e_Structure,
            "Loop closure 1A bin",
            "loop_closure"
          )
        )
      );

      // loop closure score with wider transition region, calculation is faster than caching
      e_ScoreLoopClosureGradient = GetScores().AddScore
      (
        util::ShPtr< score::ProteinModel>
        (
          new score::ProteinModelWrapper
          (
            util::CloneToShPtr
            (
              score::ProteinModelSSENeighbors
              (
                util::CloneToShPtr( score::LoopClosure( 1, 20, 1.0)),
                false // no normalization
              )
            ),
            score::ProteinModel::e_Structure,
            "Loop closure 20A bin",
            "loop_closure_gradient"
          )
        )
      );

      // phi psi
      e_ScorePhiPsi = GetScores().AddScore
      (
        util::ShPtr< score::ProteinModel>
        (
          new score::ProteinModelSSE
          (
            Scores::WrapCacheSSEScore
            (
              util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> > >
              (
                new score::PhiPsi( score::PhiPsi::GetDefaultScheme())
              )
            ),
            false
          )
        )
      );

      e_ScoreSSELoopClash = GetScores().AddScore
      (
        Scores::WrapCacheProteinModelScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::ProteinModelSSELinearLoopProximity( false, true, 0.125)
          )
        )
      );

      // sse confidence
      // iterate over sse methods
      storage::Set< sspred::Method> ss_pred_methods( sspred::Methods::GetCommandLineMethods());
      for
      (
        storage::Set< sspred::Method>::const_iterator
          sspred_itr( ss_pred_methods.Begin()), sspred_itr_end( ss_pred_methods.End());
        sspred_itr != sspred_itr_end;
        ++sspred_itr
      )
      {
        const std::string scheme( "ss_" + sspred_itr->GetName());
        util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> > >
        sp_sse_preds
        (
          util::CloneToShPtr( score::SSEPredictions( *sspred_itr, 0.0, scheme))
        );

        // add the wrapped ProteinModel score
        const util::ShPtr< score::ProteinModel> sp_score
        (
          util::CloneToShPtr
          (
            score::ProteinModelSSE
            (
              sp_sse_preds,
              false,
              score::ProteinModel::e_Sequence,
              "SS prediction"
            )
          )
        );
        m_SSPredScores[ *sspred_itr].First() = GetScores().AddScore( sp_score);

        // construct the entropy term form the same score
        const std::string scheme_ent( scheme + "_ent");
        util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> > >
        sp_sse_preds_ent
        (
          Scores::WrapCacheSSEScore( util::CloneToShPtr( score::SSEPredictions( *sspred_itr, 0.0, "")))
        );
        m_SSPredScores[ *sspred_itr].Second() = GetScores().AddScore
        (
          Scores::WrapCacheProteinModelScore
          (
            util::CloneToShPtr
            (
              score::ProteinModelInverted
              (
                Scores::WrapCacheProteinModelScore( sp_score),
                m_ProteinInverter,
                scheme_ent,
                score::ProteinModel::e_Sequence,
                "SS prediction"
              )
            )
          )
        );
      }

      // radius of gyration
      e_ScoreRadiusOfGyration = GetScores().AddScore
      (
        util::CloneToShPtr( score::RadiusOfGyration())
      );

      // relative contact order score
      e_ScoreContactOrder = GetScores().AddScore
      (
        util::CloneToShPtr
        (
          score::ContactOrder
          (
            contact::Order::e_RelativeAAsUsed,
            false, // normalize
            score::ContactOrder::GetDefaultScheme(),
            false
          )
        )
      );

      // sse packing
      e_ScoreSSEFragmentPairPacking = GetScores().AddScore
      (
        util::ShPtr< score::ProteinModel>
        (
          new score::ProteinModelSSEPairs
          (
            Scores::WrapCacheSSEPairScore
            (
              util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, assemble::SSE, double> >
              (
                new score::SSEPairsFragments
                (
                  assemble::GetSSEGeometryPackingListPickers().e_BestInteractionWeight,
                  score::SSEPairPacking
                  (
                    score::SSEPairPacking::GetDefaultScheme() + "_fr", "sse_fragment_angle_distance.histograms2D"
                  ),
                  false
                )
              ),
              false // symmetric
            ),
            false,
            score::ProteinModel::e_Structure,
            "SSE packing"
          )
        )
      );

      // strand pairing. Very fast to compute, so caching is not beneficial
      e_ScoreStrandFragmentPairing = GetScores().AddScore
      (
        util::ShPtr< score::ProteinModel>
        (
          new score::ProteinModelSSEPairs
          (
            Scores::WrapCacheSSEPairScore
            (
              util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, assemble::SSE, double> >
              (
                new score::SSEPairsFragments
                (
                  assemble::GetSSEGeometryPackingListPickers().e_BestInteractionWeight,
                  score::StrandPairing
                  (
                    score::StrandPairing::GetDefaultScheme() + "_fr", "strand_fragment_angle_distance.histograms2D"
                  ), false
                )
              ),
              false
            ),
            false,
            score::ProteinModel::e_Structure,
            "Strand pairing"
          )
        )
      );

       // between sses aa clash checking all atoms of each aa
      e_ScoreSSECompleteness = GetScores().AddScore
      (
        Scores::WrapCacheProteinModelScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::ProteinModelSSECompleteness()
          )
        )
      );
      e_ScoreSSETripletChirality = GetScores().AddScore
      (
        Scores::WrapCacheProteinModelScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::ProteinModelSSEChirality( 120.0, false)
          )
        )
      );
      e_ScoreSSEDerivedTripletChirality = GetScores().AddScore
      (
        Scores::WrapCacheProteinModelScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::ProteinModelSSEChirality( 120.0, true)
          )
        )
      );
      e_ScoreSSEContactType = GetScores().AddScore
      (
        Scores::WrapCacheProteinModelScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::ProteinModelSSEPacking( score::ProteinModelSSEPacking::e_ContactType)
          )
        )
      );
      e_ScoreSSEAdjacentContact = GetScores().AddScore
      (
        Scores::WrapCacheProteinModelScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::ProteinModelSSEPacking( score::ProteinModelSSEPacking::e_AdjacentSSEContactPropensity)
          )
        )
      );
      e_ScoreSSEOrientation = GetScores().AddScore
      (
        Scores::WrapCacheProteinModelScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::ProteinModelSSEPacking( score::ProteinModelSSEPacking::e_Orientation)
          )
        )
      );
      e_ScoreSSEInteractionWeight = GetScores().AddScore
      (
        Scores::WrapCacheProteinModelScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::ProteinModelSSEPacking( score::ProteinModelSSEPacking::e_InteractionWeight)
          )
        )
      );
    }

    //! @brief modify the score weight set
    //! @param SCORE_WEIGHT_SET Score weight set
    void DefaultScores::ModifyScoreWeightSet( ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
      // reset the map
      SCORE_WEIGHT_SET.Reset();

      // adjust the weights
      SCORE_WEIGHT_SET.SetWeight( e_ScoreAAPairClash           ,     0.00);
      SCORE_WEIGHT_SET.SetWeight( e_ScoreAAPairDistance        ,     0.00);
      SCORE_WEIGHT_SET.SetWeight( e_ScoreAAPairHiResClash      ,    50.00);
      SCORE_WEIGHT_SET.SetWeight( e_ScoreAAPairSCInteraction   ,     0.00);
      SCORE_WEIGHT_SET.SetWeight( e_ScoreAANeighborCount       ,    50.00);
      SCORE_WEIGHT_SET.SetWeight( e_ScoreAANeighborCountEntropy,    50.00);
      SCORE_WEIGHT_SET.SetWeight( e_ScoreContactOrder          ,     0.00);
      SCORE_WEIGHT_SET.SetWeight( e_ScoreLoop                  ,    10.00);
      SCORE_WEIGHT_SET.SetWeight( e_ScoreLoopAngle             ,     0.00);
      SCORE_WEIGHT_SET.SetWeight( e_ScoreLoopClosureGradient   , 50000.00);
      SCORE_WEIGHT_SET.SetWeight( e_ScoreRadiusOfGyration      ,     5.00);
      SCORE_WEIGHT_SET.SetWeight( e_ScoreSSELoopClash          ,    10.00);
      SCORE_WEIGHT_SET.SetWeight( e_ScoreSSEFragmentPairPacking,     8.00);
      SCORE_WEIGHT_SET.SetWeight( e_ScoreStrandFragmentPairing ,    20.00);
      SCORE_WEIGHT_SET.SetWeight( e_ScoreSSECompleteness       ,     0.00);
      SCORE_WEIGHT_SET.SetWeight( e_ScoreSSETripletChirality   ,     0.00);
      SCORE_WEIGHT_SET.SetWeight( e_ScoreSSEDerivedTripletChirality, 0.00);
      SCORE_WEIGHT_SET.SetWeight( e_ScoreSSEContactType,             0.00);
      SCORE_WEIGHT_SET.SetWeight( e_ScoreSSEAdjacentContact,         0.00);
      SCORE_WEIGHT_SET.SetWeight( e_ScoreSSEOrientation,             0.00);
      SCORE_WEIGHT_SET.SetWeight( e_ScoreSSEInteractionWeight,       0.00);

      for
      (
        storage::Map< sspred::Method, storage::VectorND< 2, Score> >::const_iterator sspred_itr( m_SSPredScores.Begin()),
          sspred_itr_end( m_SSPredScores.End());
        sspred_itr != sspred_itr_end;
        ++sspred_itr
      )
      {
        double weight( 1.0);
        if( sspred_itr->first == sspred::GetMethods().e_PSIPRED)
        {
          weight = 20.0;
        }
        else if( sspred_itr->first == sspred::GetMethods().e_JUFO)
        {
          weight = 5.0;
        }

        SCORE_WEIGHT_SET.SetWeight( sspred_itr->second.First(), weight);  // normal score
        SCORE_WEIGHT_SET.SetWeight( sspred_itr->second.Second(), weight); // entropy score
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DefaultScores::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DefaultScores::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return the stream
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
