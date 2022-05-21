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
#include "fold/bcl_fold_protocol_membrane.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_domain_sse_pool_overlapping.h"
#include "assemble/bcl_assemble_pick_sse_short_loops.h"
#include "biol/bcl_biol_exposure_prediction.h"
#include "biol/bcl_biol_membrane.h"
#include "coord/bcl_coord_move_rotate_random_external_reference.h"
#include "coord/bcl_coord_move_translate_random.h"
#include "fold/bcl_fold_default_flags.h"
#include "fold/bcl_fold_default_mutates.h"
#include "fold/bcl_fold_default_scores.h"
#include "fold/bcl_fold_mutate_protein_model.h"
#include "fold/bcl_fold_mutate_protein_model_sse_add.h"
#include "fold/bcl_fold_mutate_tree.h"
#include "fold/bcl_fold_placement_sse_short_loop.h"
#include "fold/bcl_fold_score_weight_set.h"
#include "fold/bcl_fold_setup.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_printer_membrane.h"
#include "pdb/bcl_pdb_printer_quality_membrane.h"
#include "score/bcl_score_aa_neighborhood_exposure_prediction.h"
#include "score/bcl_score_environment_predictions.h"
#include "score/bcl_score_protein_model_aa_neighborhood.h"
#include "score/bcl_score_protein_model_membrane_topology.h"
#include "score/bcl_score_protein_model_sse.h"
#include "score/bcl_score_sse_membrane_alignment.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProtocolMembrane::ProtocolMembrane() :
      m_EnvironmentScores()
    {
    }

    //! @brief Clone function
    //! @return pointer to new ProtocolMembrane
    ProtocolMembrane *ProtocolMembrane::Clone() const
    {
      return new ProtocolMembrane( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief access to only instance
    //! @return reference to only instance
    ProtocolMembrane &ProtocolMembrane::GetInstance()
    {
      static ProtocolMembrane s_protocol_instance;
      return s_protocol_instance;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProtocolMembrane::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

     //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ProtocolMembrane::GetAlias() const
    {
      static const std::string s_name( "ProtocolMembrane");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProtocolMembrane::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "protocol for folding membrane proteins");
      serializer.AddInitializer
        (
         "mutate add mem SSE short loop",
         "add membrane",
         io::Serialization::GetAgent( &e_MutateAddMemSSEShortLoop)
         );
      serializer.AddInitializer
        (
         "mutate model global z translate",
         "translate the model along its z axis",
         io::Serialization::GetAgent( &e_MutateModelGlobalZTranslate)
         );
      serializer.AddInitializer
        (
         "mutate model global z rotate",
         "rotate the model",
         io::Serialization::GetAgent( &e_MutateModelGlobalRotate)
         );
      serializer.AddInitializer
      (
        "score SSE membrane alignment",
        "scores alignment of SSEs with the membrane",
        io::Serialization::GetAgent( &e_ScoreSSEMembraneAlignment)
       );
      serializer.AddInitializer
      (
        "score helix topology",
        "scores liekliness of helix topology",
        io::Serialization::GetAgent( &e_ScoreHelixTopology)
       );
      serializer.AddInitializer
      (
        "score exposure",
        "scores model exposure outside of the membrane",
        io::Serialization::GetAgent( &e_ScoreExposure)
       );
      serializer.AddInitializer
      (
        "environment scores",
        "sspred environment score",
        io::Serialization::GetAgent( &m_EnvironmentScores)
       );

      return serializer;
    }

  ///////////
  // flags //
  ///////////

    //! @brief returns all flags that are specialized for this protocol
    //! @return all flags that are specialized for this protocol
    const util::ShPtrVector< command::FlagInterface> &ProtocolMembrane::GetAllFlags() const
    {
      // initialize static ShPtrVector of FlagInterfaces to form the comment line
      static util::ShPtrVector< command::FlagInterface> s_all_flags_vector;

      // if the flag vector is initialize for the first time
      if( s_all_flags_vector.IsEmpty())
      {
        // insert all the flags in the vector
        s_all_flags_vector.PushBack( biol::Membrane::GetFlagMembrane());
        s_all_flags_vector.PushBack
        (
          score::ProteinModelMembraneTopology::GetFlagExpectedTransmembraneHelicesPoolFile()
        );
      }

      // end
      return s_all_flags_vector;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief modifies the start model
    //! @param START_MODEL Start model to be modified
    void ProtocolMembrane::ModifyStartModel( assemble::ProteinModel &START_MODEL) const
    {
      // add the membrane to the protein model data
      util::ShPtr< assemble::ProteinModelData> sp_data( START_MODEL.GetProteinModelData());
      sp_data->Insert( assemble::ProteinModelData::e_Membrane, util::CloneToShPtr( biol::Membrane::GetCommandLineMembrane()));
      START_MODEL.SetProteinModelData( sp_data);

      // read environment predictions
      if( score::AANeighborhoodExposurePrediction::GetFlagScoreExposure()->GetFlag())
      {
        biol::ExposurePrediction::ReadPredictions
        (
          START_MODEL,
          DefaultFlags::GetFlagReadSequenceDataPath()->GetParameterList()( 1)->GetValue(),
          DefaultFlags::GetFlagReadSequenceDataPath()->GetFirstParameter()->GetValue()
        );
        util::ShPtr< assemble::ProteinModel> sp_native_model
        (
          START_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_NativeModel).HardCopy()
        );
        if( sp_native_model.IsDefined())
        {
          biol::ExposurePrediction::ReadPredictions
          (
            *sp_native_model,
            DefaultFlags::GetFlagReadSequenceDataPath()->GetParameterList()( 1)->GetValue(),
            DefaultFlags::GetFlagReadSequenceDataPath()->GetFirstParameter()->GetValue()
          );
        }
        util::ShPtr< assemble::ProteinModel> sp_native_filtered_model
        (
          START_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_NativeFilteredModel).HardCopy()
        );
        if( sp_native_filtered_model.IsDefined())
        {
          biol::ExposurePrediction::ReadPredictions
          (
            *sp_native_filtered_model,
            DefaultFlags::GetFlagReadSequenceDataPath()->GetParameterList()( 1)->GetValue(),
            DefaultFlags::GetFlagReadSequenceDataPath()->GetFirstParameter()->GetValue()
          );
        }
      }
    }

    //! @brief initialize the scores and add them to Scores enumerator
    void ProtocolMembrane::InitializeScores()
    {
      if( e_ScoreSSEMembraneAlignment.IsDefined())
      {
        return;
      }

      // radius of gyration for membrane
      // membrane alignment
      e_ScoreSSEMembraneAlignment = GetScores().AddScore
      (
        Scores::WrapCacheProteinModelScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::ProteinModelSSE
            (
              Scores::WrapCacheSSEScore
              (
                util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> > >
                (
                  new score::SSEMembraneAlignment()
                )
              ),
              false
            )
          )
        )
      );

      // make sure no loops will need to go through the membrane to connect consecutive sses
      score::ProteinModelMembraneTopology topology;
      if( topology.InitializeFromFlag())
      {
        e_ScoreHelixTopology = GetScores().AddScore
        (
          Scores::WrapCacheProteinModelScore( util::ShPtr< score::ProteinModel>( topology.Clone()))
        );
      }

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
        // add environment score
        m_EnvironmentScores[ *sspred_itr] = GetScores().AddScore
        (
          Scores::WrapCacheProteinModelScore
          (
            util::ShPtr< score::ProteinModel>
            (
              new score::ProteinModelSSE
              (
                Scores::WrapCacheSSEScore
                (
                  util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> > >
                  (
                    new score::EnvironmentPredictions
                    (
                      *sspred_itr,
                      "ss_" + sspred_itr->GetName() + "_env"
                    )
                  )
                ),
                false
              )
            )
          )
        );
      }

      // if scoring exposure
      if( score::AANeighborhoodExposurePrediction::GetFlagScoreExposure()->GetFlag())
      {
        // construct the score
        e_ScoreExposure = GetScores().AddScore
        (
          Scores::WrapCacheProteinModelScore
          (
            util::ShPtr< score::ProteinModel>
            (
              new score::ProteinModelAANeighborhood
              (
                util::CloneToShPtr
                (
                  score::AANeighborhoodExposurePrediction()
                ),
                score::ProteinModelAANeighborhood::e_RMSD
              )
            )
          )
        );

        // construct the score
        e_ScoreExposureExplained = GetScores().AddScore
        (
          Scores::WrapCacheProteinModelScore
          (
            util::ShPtr< score::ProteinModel>
            (
              new score::ProteinModelAANeighborhood
              (
                util::CloneToShPtr
                (
                  score::AANeighborhoodExposurePrediction()
                ),
                score::ProteinModelAANeighborhood::e_FractionExplained
              )
            )
          )
        );
      }
    }

    //! @brief modify the score weight set
    //! @param SCORE_WEIGHT_SET Score weight set
    void ProtocolMembrane::ModifyScoreWeightSet( ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
      // set the weights
      SCORE_WEIGHT_SET.SetWeight( e_ScoreSSEMembraneAlignment, 8.0);
      SCORE_WEIGHT_SET.SetWeight( e_ScoreHelixTopology, 500);

      // iterate over sspred scores
      for
      (
        storage::Map< sspred::Method, storage::VectorND< 2, Score> >::const_iterator
          sspred_itr( DefaultScores::GetInstance().GetSSPredScores().Begin()),
          sspred_itr_end( DefaultScores::GetInstance().GetSSPredScores().End());
        sspred_itr != sspred_itr_end; ++sspred_itr
      )
      {
        // set weight for octopus and jufo9d
        if( sspred_itr->first == sspred::GetMethods().e_OCTOPUS)
        {
          const double weight( 20.0);
          SCORE_WEIGHT_SET.SetWeight( sspred_itr->second.First(), weight);  // normal score
          SCORE_WEIGHT_SET.SetWeight( sspred_itr->second.Second(), weight); // entropy score
          SCORE_WEIGHT_SET.SetWeight( m_EnvironmentScores.Find( sspred_itr->first)->second, weight); // environment score
        }
        else if( sspred_itr->first == sspred::GetMethods().e_JUFO9D)
        {
          const double weight( 5.0);
          SCORE_WEIGHT_SET.SetWeight( sspred_itr->second.First(), weight);  // normal score
          SCORE_WEIGHT_SET.SetWeight( sspred_itr->second.Second(), weight); // entropy score
          SCORE_WEIGHT_SET.SetWeight( m_EnvironmentScores.Find( sspred_itr->first)->second, weight); // environment score
        }
      }
    }

    //! @brief modify the terminate object
    //! @param CRITERION which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolMembrane::ModifyCriterion
    (
      opti::CriterionCombine< assemble::ProteinModel, double> &CRITERION,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief modify the printer object
    //! @param PRINTER which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolMembrane::ModifyPrinter
    (
      mc::PrinterCombined< assemble::ProteinModel, double> &PRINTER,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief modify the pdb factory object
    //! @param FACTORY pdb factory to be modified
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolMembrane::ModifyFactory
    (
      util::ShPtr< pdb::Factory> &FACTORY,
      const mc::Stage &STAGE
    ) const
    {
      FACTORY->AppendPrinter
      (
        util::ShPtr< util::FunctionInterface< assemble::ProteinModel, util::ShPtrList< pdb::Line> > >
        (
          new pdb::PrinterMembrane()
        )
      );
      FACTORY->AppendPrinter
      (
        util::ShPtr< util::FunctionInterface< assemble::ProteinModel, util::ShPtrList< pdb::Line> > >
        (
          new pdb::PrinterQualityMembrane
          (
            GetSetup().GetQualityMeasures(),
            storage::Set< biol::EnvironmentType>( biol::GetEnvironmentTypes().e_MembraneCore),
            GetSetup().GetNativeModel()
          )
        )
      );
    }

    //! @brief initialize the mutates and add them to Mutates enumerator
    void ProtocolMembrane::InitializeMutates()
    {
      if( e_MutateAddMemSSEShortLoop.IsDefined())
      {
        return;
      }

      // add sse short loop using a wider angle than for soluble additions
      const size_t short_loop_length( 7);
      e_MutateAddMemSSEShortLoop = GetMutates().AddMutate
      (
        MutateProteinModelSSEAdd
        (
          assemble::PickSSEShortLoops( short_loop_length),
          PlacementSSEShortLoop( short_loop_length, 0.25, math::g_Pi / 2.0),
          "add_mem_sse_next_to_sse"
        )
      );

      // global z translation
      e_MutateModelGlobalZTranslate = GetMutates().AddMutate
      (
        util::ShPtr< math::MutateInterface< assemble::ProteinModel> >
        (
          new MutateProteinModel
          (
            coord::MoveTranslateRandom( linal::Vector3D( 0.0, 0.0, 2.0), linal::Vector3D( 0.0, 0.0, 10.0), false),
            "model_global_z_translate"
          )
        )
      );

      // global rotation
      e_MutateModelGlobalRotate = GetMutates().AddMutate
      (
        util::ShPtr< math::MutateInterface< assemble::ProteinModel> >
        (
          new MutateProteinModel
          (
            coord::MoveRotateRandomExternalReference
            (
              linal::Vector3D( 0.0, 0.0, 0.0),
              linal::Vector3D( math::g_Pi, math::g_Pi, math::g_Pi),
              math::TransformationMatrix3D(),
              coord::MoveRotateRandomExternalReference::e_InternalTranslate
            ),
            "model_global_rotate"
          )
        )
      );
    }

    //! @brief modify the mutate tree used
    //! @param MUTATE_TREE MutateTree to be modified
    void ProtocolMembrane::ModifyMutateTree( MutateTree &MUTATE_TREE) const
    {
      DefaultMutates::GetInstance().InitializeMutates();
      MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_AddSSEShortLoop, 0.0);
      MUTATE_TREE.SetMutateProbability( e_MutateAddMemSSEShortLoop, 0.75);

      // reduce split probability
      // if secondary structure prediction were provided and sse resizing is not disabled
      if
      (
        sspred::Methods::GetFlagReadSSPredictions()->GetFlag() &&
        DefaultFlags::GetFlagEnableSSEResize()->GetFlag()
      )
      {
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
          const double probability( 0);
          BCL_MessageStd
          (
            "setting split probability for " + sspred_itr->first->GetName() +
            " to " + util::Format()( probability)
          );
          MUTATE_TREE.SetMutateProbability( sspred_itr->second, probability);
        }
      }

      MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_Model, 0.05);
      MUTATE_TREE.SetMutateProbability( e_MutateModelGlobalZTranslate, 1.5);
      MUTATE_TREE.SetMutateProbability( e_MutateModelGlobalRotate, 1.5);
    }

    //! @brief get the mutate tree associated with this protocol
    //! @return the mutate tree associated with this protocol
    util::ShPtr< MutateTree> ProtocolMembrane::GetMutateTree() const
    {
      util::ShPtr< MutateTree> sp_mutate_tree( new MutateTree());
      ModifyMutateTree( *sp_mutate_tree);
      return sp_mutate_tree;
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void ProtocolMembrane::MergeAndModifyMutateTree( MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }

  ////////////
  // readme //
  ////////////

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ProtocolMembrane::GetDescription() const
    {
      // initialize string to store the description
      static const std::string s_description
      (
        "protocol to be applied to membrane proteins with specialized moves and scores"
      );

      // end
      return s_description;
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ProtocolMembrane::GetReadMe() const
    {
      // initialize string to store the readme information
      static const std::string s_readme
      (
        "BCL::MP-Fold adapts the BCL::Fold method to function with membrane proteins. KB scores are adjusted to "
        "account for the membrane environment, which is represented statically. Additional MC moves allow for "
        "efficient sampling of the membrane environment.\n"
        "When using BCL::MP-Fold in a publication, and for more detailed information regarding the method, please "
        "cite the publication describing the application's development which is currently in preparation. Refer "
        "to www.meilerlab.org for future details.\n\n"
        "Speficific flags:\n"
        "-membrane #turns on membrane scoring, allows the thickness to be set\n"
        "-tm_helices [tm_helices.pool] #define TM spanning helices in order to score the MP topology\n"
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
    std::istream &ProtocolMembrane::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProtocolMembrane::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return the stream
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
