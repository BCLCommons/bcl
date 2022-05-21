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
#include "fold/bcl_fold_protocol_em.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_sse_random.h"
#include "assemble/bcl_assemble_pick_sse_random.h"
#include "assemble/bcl_assemble_sse_compare_type.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_flag_static_and_dynamic.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "command/bcl_command_parameter_check_extension.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "coord/bcl_coord_move_rotate_defined.h"
#include "coord/bcl_coord_move_rotate_random.h"
#include "coord/bcl_coord_move_translate_random.h"
#include "density/bcl_density_protein_agreements.h"
#include "density/bcl_density_simulators.h"
#include "find/bcl_find_collector_criteria_combined.h"
#include "find/bcl_find_locator_criteria.h"
#include "find/bcl_find_pick_body_extent.h"
#include "find/bcl_find_pick_body_random.h"
#include "find/bcl_find_pick_criteria_wrapper.h"
#include "fold/bcl_fold_default_mutates.h"
#include "fold/bcl_fold_default_scores.h"
#include "fold/bcl_fold_mutate_protein_model_sse.h"
#include "fold/bcl_fold_mutate_protein_model_sse_add.h"
#include "fold/bcl_fold_mutate_protein_model_sse_remove.h"
#include "fold/bcl_fold_mutate_protein_model_sse_swap.h"
#include "fold/bcl_fold_mutate_protein_model_sse_swap_body.h"
#include "fold/bcl_fold_mutate_protein_model_sse_swap_with_pool.h"
#include "fold/bcl_fold_mutate_protein_model_sse_swap_with_pool_overlap.h"
#include "fold/bcl_fold_mutate_sse_bend_random.h"
#include "fold/bcl_fold_mutate_tree.h"
#include "fold/bcl_fold_placement_sse_into_body.h"
#include "fold/bcl_fold_score_weight_set.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_mutate_move_wrapper.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_printer_body_assignment.h"
#include "restraint/bcl_restraint_contains_body_origin.h"
#include "restraint/bcl_restraint_handler_body.h"
#include "restraint/bcl_restraint_mutate_transformation_matrix_3d_rotate.h"
#include "score/bcl_score_body_assignment.h"
#include "score/bcl_score_body_connectivity_density.h"
#include "score/bcl_score_body_extent_agreement.h"
#include "score/bcl_score_restraint_body_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProtocolEM::ProtocolEM()
    {
    }

    //! @brief Clone function
    //! @return pointer to new ProtocolEM
    ProtocolEM *ProtocolEM::Clone() const
    {
      return new ProtocolEM( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief access to only instance
    //! @return reference to only instance
    ProtocolEM &ProtocolEM::GetInstance()
    {
      static ProtocolEM s_protocol_instance;
      return s_protocol_instance;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProtocolEM::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ProtocolEM::GetAlias() const
    {
      static const std::string s_name( "ProtocolEM");
      return s_name;
    }

  ///////////
  // flags //
  ///////////

    //! @brief returns all flags that are specialized for this protocol
    //! @return all flags that are specialized for this protocol
    const util::ShPtrVector< command::FlagInterface> &ProtocolEM::GetAllFlags() const
    {
      // initialize static ShPtrVector of FlagInterfaces to form the comment line
      static util::ShPtrVector< command::FlagInterface> s_all_flags_vector;

      // if the flag vector is initialize for the first time
      if( s_all_flags_vector.IsEmpty())
      {
        // insert all the flags in the vector
        s_all_flags_vector.PushBack( GetFlagBodyRestraint());
        s_all_flags_vector.PushBack( GetFlagPrintBodyAssignment());
        s_all_flags_vector.PushBack( GetFlagScoreDensityAgreement());
        s_all_flags_vector.PushBack( GetFlagScoreDensityConnectivity());
        s_all_flags_vector.PushBack( GetFlagEMRefinement());
      }

      // end
      return s_all_flags_vector;
    }

    //! @brief return command line flag for using body restraints
    //! @return command line flag for using body restraints
    util::ShPtr< command::FlagInterface> &ProtocolEM::GetFlagBodyRestraint()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "body_restraint",
          "\tuse body restraint"
        )
      );

      // initialize parameters
      static util::ShPtr< command::ParameterInterface> s_filename
      (
        new command::Parameter( "body_restraint_filename", "\tfull path and name of restraint file", "body_restraint_file.pdb")
      );
      static util::ShPtr< command::ParameterInterface> s_upper_tolerance_width
      (
        new command::Parameter
        (
          "upper_tolerance_width",
          "bodies which have a length less than this tolerance longer than restraint body will still 100% agree with restraint\t",
          "0.0"
        )
      );
      static util::ShPtr< command::ParameterInterface> s_lower_tolerance_width
      (
        new command::Parameter
        (
          "lower_tolerance_width",
          "bodies which have a length less than this tolerance shorter than restraint body will still 100% agree with restraint\t",
          "0.0"
        )
      );
      static util::ShPtr< command::ParameterInterface> s_upper_transition_width
      (
        new command::Parameter
        (
          "upper_transtion_width",
          "agreement above restraint will transition from 100% agreement to 0% agreement according to cosine curve across this width\t",
          "0.0"
        )
      );
      static util::ShPtr< command::ParameterInterface> s_lower_transition_width
      (
        new command::Parameter
        (
          "lower_transtion_width",
          "agreement below restraint will transition from 100% agreement to 0% agreement according to cosine curve across this width\t",
          "0.0"
        )
      );
      static util::ShPtr< command::ParameterInterface> s_energy_well_depth
      (
        new command::Parameter
        (
          "energy_well_depth",
          "the amount of bonus that should be given to energy if restraint is completely fulfilled\t",
          "-1.0"
        )
      );

      // if the flag is initialized for the first time
      if( s_flag->GetParameterList().IsEmpty())
      {
        util::ShPtr< command::FlagStatic> flag( s_flag);
        // insert parameters
        flag->PushBack( s_filename);
        flag->PushBack( s_upper_tolerance_width);
        flag->PushBack( s_lower_tolerance_width);
        flag->PushBack( s_upper_transition_width);
        flag->PushBack( s_lower_transition_width);
        flag->PushBack( s_energy_well_depth);
      }
      // end
      return s_flag;
    }

    //! @brief return command line flag for writing body assignments when using EMFold
    //! @return command line flag for writing body assignments when using EMFold
    util::ShPtr< command::FlagInterface> &ProtocolEM::GetFlagPrintBodyAssignment()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "print_body_assignment", "\tuse enables writing of body assignment when used in EM-Fold protocol"
        )
      );

      // end
      return s_flag;
    }

    //! @brief return command line flag for using density agreement score
    //! @return command line flag for using density agreement score
    util::ShPtr< command::FlagInterface> &ProtocolEM::GetFlagScoreDensityAgreement()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStaticAndDynamic
        (
          "score_density_agreement",
          "\tflag to enable use of density agreement scores",
          command::Parameter
          (
            "agreement_score",
            "the agreement objective to be used",
            command::ParameterCheckEnumerate< density::ProteinAgreements>()
          ),
          0,
          density::GetProteinAgreements().GetEnumCount()
        )
      );

      // initialize parameters
      static util::ShPtr< command::ParameterInterface> s_filename
      (
        new command::Parameter
        (
          "density_map_filename",
          "\tfilename for density map (mrc/ccp4) file",
          command::ParameterCheckFileExistence(),
          ""
        )
      );
      static util::ShPtr< command::ParameterInterface> s_resolution
      (
        new command::Parameter
        (
          "resolution_density_map",
          "\tresolution of given electron density map [A] - will be used to simulate density and calculate correlation",
          command::ParameterCheckRanged< double>( 0.0, 100.0),
          "0.0"
        )
      );
      static util::ShPtr< command::ParameterInterface> s_simulator
      (
        new command::Parameter
        (
          "simulator",
          "the simulator to be used to generate density maps from a atom structure",
          command::ParameterCheckEnumerate< density::Simulators>(),
          density::GetSimulators().e_Gaussian.GetName()
        )
      );

      // if the flag is initialized for the first time
      if( s_flag->GetParameterList().IsEmpty())
      {
        util::ShPtr< command::FlagStaticAndDynamic> flag( s_flag);
        // insert parameters
        flag->PushBack( s_filename);
        flag->PushBack( s_resolution);
        flag->PushBack( s_simulator);
      }

      // end
      return s_flag;
    }

    //! @brief return command line flag for using density connectivity score
    //! @return command line flag for using density connectivity score
    util::ShPtr< command::FlagInterface> &ProtocolEM::GetFlagScoreDensityConnectivity()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "score_density_connectivity",
          "\tflag to enable use of density connectivity score that rewards models that connect SSEs in high density regions"
        )
      );

      // initialize parameters
      static util::ShPtr< command::ParameterInterface> s_filename
      (
        new command::Parameter
        (
          "density_map_filename", "\tfilename for density map (mrc) file", command::ParameterCheckExtension( ".mrc"), ""
        )
      );

      // if the flag is initialized for the first time
      if( s_flag->GetParameterList().IsEmpty())
      {
        util::ShPtr< command::FlagStatic> flag( s_flag);
        // insert parameters
        flag->PushBack( s_filename);
      }

      // end
      return s_flag;
    }

    //! @brief return command line flag for using em fold refinement protocol
    //! @return command line flag for using em fold refinement protocol
    util::ShPtr< command::FlagInterface> &ProtocolEM::GetFlagEMRefinement()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "em_refinement",
          "\tenable use of refinement protocol for EM fold"
        )
      );

      // end
      return s_flag;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief modifies the start model
    //! @param START_MODEL Start model to be modified
    void ProtocolEM::ModifyStartModel( assemble::ProteinModel &START_MODEL) const
    {
      // make sure that if there is one sse that it is properly placed in a density rod
      if( START_MODEL.GetSSEs().GetSize() == 1)
      {
        // set extent tolerances for the length extent picker
        // these tolerances are set to discriminate between helices and strands only, no length discrimination
        linal::Vector3D extent_tolerances( 0.1, 0.1, 50.0);

        // set the length extent picker
        const find::PickBodyExtent picker_extent( extent_tolerances, extent_tolerances);

        // create ShPtr to a PickCriteriaInterface "pick" which is used to determine how a restraint body is
        // selected
        util::ShPtr
        <
          find::PickCriteriaInterface
          <
            util::ShPtr< assemble::SSEGeometryInterface>,
            util::ShPtrVector< assemble::SSEGeometryInterface>, assemble::SSEGeometryInterface
          >
        > pick
        (
          // picker_random.Clone()
          picker_extent.Clone()
        );

        // random sse locator
        const util::ShPtr< assemble::LocatorSSERandom> sp_locator_sse_random( new assemble::LocatorSSERandom());

        const MutateProteinModelSSESwapBody mutate
        (
          GetBodyRestraint(),
          pick,
          sp_locator_sse_random,
          "sse_em_move_into_empty_body"
        );

        math::MutateResult< assemble::ProteinModel> result( mutate( START_MODEL));
        BCL_Assert
        (
          result.GetArgument().IsDefined(),
          "Unable to place initial SSE into body, make sure density map body SSE content is similar to SSE pool "
          "used (i.e. don't use strands in pool if none in body restraint)"
        );

        START_MODEL = ( *result.GetArgument());
      }
    }

    //! @brief initialize the scores and add them to Scores enumerator
    void ProtocolEM::InitializeScores()
    {
      // if the body restraint flag is given
      if( GetFlagBodyRestraint()->GetFlag())
      {
        // body restraints score
        e_BodyRestraintScore = GetScores().AddScore( GetScoreBodyRestraint());
      }

      // if connectivity score flag was provided
      if( GetFlagScoreDensityConnectivity()->GetFlag())
      {
        // body connectivity density (connectivity score)
        e_BodyConnectivityDensityScore = GetScores().AddScore( GetScoreBodyConnectivityDensity());
      }

      // density agreement score
      if( GetFlagScoreDensityAgreement()->GetFlag())
      {
        const util::ShPtrVector< score::ProteinModel> scores( GetScoresDensityAgreement());
        // density agreement score
        for
        (
          util::ShPtrVector< score::ProteinModel>::const_iterator itr( scores.Begin()), itr_end( scores.End());
          itr != itr_end;
          ++itr
        )
        {
          m_DensityAgreementScores.Insert( GetScores().AddScore( *itr));
          BCL_MessageDbg( "added score: " + ( --m_DensityAgreementScores.End())->GetName());
        }
      }
    }

    //! @brief modify the score weight set
    //! @param SCORE_WEIGHT_SET Score weight set
    void ProtocolEM::ModifyScoreWeightSet( ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
      // weights for EM specific score
      SCORE_WEIGHT_SET.SetWeight( e_BodyRestraintScore, 300.0);
      SCORE_WEIGHT_SET.SetWeight( e_BodyConnectivityDensityScore, 450.0);

      // iterate over density agreement scores
      for
      (
        storage::Set< Score>::const_iterator
          itr( m_DensityAgreementScores.Begin()), itr_end( m_DensityAgreementScores.End());
        itr != itr_end;
        ++itr
      )
      {
        SCORE_WEIGHT_SET.SetWeight( *itr, 1.0);
      }

      // if secondary structure predictions are provided
      for
      (
        storage::Map< sspred::Method, storage::VectorND< 2, Score> >::const_iterator
          sspred_itr( DefaultScores::GetInstance().m_SSPredScores.Begin()),
          sspred_itr_end( DefaultScores::GetInstance().m_SSPredScores.End());
        sspred_itr != sspred_itr_end;
        ++sspred_itr
      )
      {
        SCORE_WEIGHT_SET.SetWeight( sspred_itr->second.First(), 1.0);
        SCORE_WEIGHT_SET.SetWeight( sspred_itr->second.Second(), 1.0);
      }

      // density agreements 
      for
      (
        storage::Set< Score>::const_iterator
          itr( m_DensityAgreementScores.Begin()), itr_end( m_DensityAgreementScores.End());
        itr != itr_end;
        ++itr
      )
      {
        SCORE_WEIGHT_SET.SetWeight( *itr, 1.0);
      }
    }

    //! @brief modify the terminate object
    //! @param CRITERION which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolEM::ModifyCriterion
    (
      opti::CriterionCombine< assemble::ProteinModel, double> &CRITERION,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief modify the printer object
    //! @param PRINTER which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolEM::ModifyPrinter
    (
      mc::PrinterCombined< assemble::ProteinModel, double> &PRINTER,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief modify the pdb factory object
    //! @param FACTORY pdb factory to be modified
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolEM::ModifyFactory
    (
      util::ShPtr< pdb::Factory> &FACTORY,
      const mc::Stage &STAGE
    ) const
    {
      // write body sse assignment was selected
      if( GetFlagPrintBodyAssignment()->GetFlag())
      {
        FACTORY->AppendPrinter
        (
          util::ShPtr< util::FunctionInterface< assemble::ProteinModel, util::ShPtrList< pdb::Line> > >
          (
            new pdb::PrinterBodyAssignment( GetBodyRestraints())
          )
        );
      }
    }

    //! @brief initialize the mutates and add them to Mutates enumerator
    void ProtocolEM::InitializeMutates()
    {
      // if the body restraint flag is given
      if( GetFlagBodyRestraint()->GetFlag())
      {
        InitializeBodyAddMutates();
      }

      // initialize first collector for use in locator 1
      // this collector collects all the sses that agree in type with the criterion sse
      const util::ShPtr
      <
        find::CollectorCriteriaInterface
        <
          util::SiPtrList< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::SSE
        >
      > collector_same_sstype
      (
        // set up collector combined with single criterion (match type)
        new find::CollectorCriteriaCombined< assemble::SSE>
        (
          util::ShPtr< util::BinaryFunctionInterface< assemble::SSE, assemble::SSE, bool> >
          (
            new assemble::SSECompareType()
          )
        )
      );

      // initialize random picker (this picker will be used with both collectors)
      const util::ShPtr
      <
        find::PickCriteriaInterface
        <
          util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::SSE
        >
      > picker
      (
        new find::PickCriteriaWrapper
        <
          util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::SSE
        >
        (
          // set up picker with random picker
          assemble::PickSSERandom()
        )
      );

      // create locator to locate sses of same type
      const find::LocatorCriteria
      <
        util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::SSE,
        util::SiPtrList< const assemble::SSE>
      > locator_ss_type( collector_same_sstype, picker);

      // random sse locator
      const util::ShPtr< assemble::LocatorSSERandom> sp_locator_sse_random( new assemble::LocatorSSERandom());

      // remove
      GetMutates().AddMutate( MutateProteinModelSSERemove( sp_locator_sse_random, "remove_em_random"));

      // swaps
      GetMutates().AddMutate( MutateProteinModelSSESwap( locator_ss_type, false, "swap_em_sses"));
      GetMutates().AddMutate( MutateProteinModelSSESwapWithPool( locator_ss_type, false, "swap_em_sse_with_pool"));

      // swap move with single swap only
      GetMutates().AddMutate
      (
        MutateProteinModelSSESwapWithPoolOverlap
        (
          collector_same_sstype, true, false, "swap_em_sse_with_pool_overlap_single_swap_only"
        )
      );

      // set the random picker
      const find::PickCriteriaWrapper
      <
        util::ShPtr< assemble::SSEGeometryInterface>,
        util::ShPtrVector< assemble::SSEGeometryInterface>,
        assemble::SSEGeometryInterface
      >
      picker_random =
        find::PickCriteriaWrapper
        <
          util::ShPtr< assemble::SSEGeometryInterface>,
          util::ShPtrVector< assemble::SSEGeometryInterface>,
          assemble::SSEGeometryInterface
        >( find::PickBodyRandom());

      // set extent tolerances for the length extent picker
      // these tolerances are set to discriminate between helices and strands only, no length discrimination
      linal::Vector3D extent_tolerances( 0.1, 0.1, 50.0);

      // important this is only half the length!!
//      linal::Vector3D extent_tolerances( 0.5, 0.5, 5.0);

      // set the length extent picker
      const find::PickBodyExtent picker_extent( extent_tolerances, extent_tolerances);

      // create ShPtr to a PickCriteriaInterface "pick" which is used to determine how a restraint body is
      // selected
      util::ShPtr
      <
        find::PickCriteriaInterface
        <
          util::ShPtr< assemble::SSEGeometryInterface>,
          util::ShPtrVector< assemble::SSEGeometryInterface>, assemble::SSEGeometryInterface
        >
      > pick
      (
        // picker_random.Clone()
        picker_extent.Clone()
      );

      if( GetFlagBodyRestraint()->GetFlag())
      {
        GetMutates().AddMutate
        (
          MutateProteinModelSSESwapBody
          (
            GetBodyRestraint(),
            pick,
            sp_locator_sse_random,
            "sse_em_move_into_empty_body"
          )
        );
      }

      GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_sse_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveTranslateRandom( linal::Vector3D( 0.0, 0.0, 1.0), true)), false)), //for EM-Fold refinement step
          "sse_em_translate_random_z"
        )
      );

      GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_sse_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveRotateRandom( linal::Vector3D( math::g_Pi / 8.0, math::g_Pi / 8.0, 0), true)), false)), //for EM-Fold refinement step
          "sse_em_rotate_random_xy"
        )
      );

      GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_sse_random,
          util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( util::CloneToShPtr( coord::MoveRotateRandom( linal::Vector3D( 0.0, 0.0, math::g_Pi / 2.0), true)), false)), //for EM-Fold refinement step
          "sse_em_rotate_random_z"
        )
      );

      // get flips
      const util::ShPtr< coord::MoveInterface> x_flip( coord::MoveRotateDefined::GetFlipMove( coord::GetAxes().e_X).Clone());
      const util::ShPtr< coord::MoveInterface> y_flip( coord::MoveRotateDefined::GetFlipMove( coord::GetAxes().e_Y).Clone());
      const util::ShPtr< coord::MoveInterface> z_flip( coord::MoveRotateDefined::GetFlipMove( coord::GetAxes().e_Z).Clone());

      GetMutates().AddMutate( MutateProteinModelSSE( sp_locator_sse_random, util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( x_flip, false)), "sse_em_flip_x"));
      GetMutates().AddMutate( MutateProteinModelSSE( sp_locator_sse_random, util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( y_flip, false)), "sse_em_flip_y"));
      GetMutates().AddMutate( MutateProteinModelSSE( sp_locator_sse_random, util::CloneToShPtr( math::MutateMoveWrapper< assemble::SSE>( z_flip, false)), "sse_em_flip_z"));

      GetMutates().AddMutate
      (
        util::ShPtr< math::MutateInterface< assemble::ProteinModel> >
        (
          new MutateProteinModelSSE
          (
            sp_locator_sse_random,
            util::CloneToShPtr
            (
              MutateSSEBendRandom
              (
                math::Range< double>( math::Angle::Radian( -5.0), math::Angle::Radian( 5.0)),
                math::Range< double>( math::Angle::Radian( -5.0), math::Angle::Radian( 5.0))
              )
            ),
            "sse_em_bend_random_small"
          )
        )
      );

    }

    //! @brief modify the mutate tree used
    //! @param MUTATE_TREE MutateTree to be modified
    void ProtocolEM::ModifyMutateTree( MutateTree &MUTATE_TREE) const
    {
      // reset all the probabilities
      MUTATE_TREE.Reset();

      // if refinement flag given
      if( GetFlagEMRefinement()->GetFlag())
      {
        // first choose probabilities for add, remove, swap, move
        MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_Add,    0.00);
        MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_Remove, 0.00);
        MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_Swap,   0.00);
        MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_SSE,   1.00);

        // now assign probabilities for individual mutates
        MUTATE_TREE.SetMutateProbability( "sse_em_bend_random_small",  10.0);
        MUTATE_TREE.SetMutateProbability( "sse_em_translate_random_z",  7.0);
        MUTATE_TREE.SetMutateProbability( "sse_em_rotate_random_xy",    7.0);
        MUTATE_TREE.SetMutateProbability( "sse_em_rotate_random_z",     7.0);
        MUTATE_TREE.SetMutateProbability( "sse_em_flip_z",              1.0);
      }

      // else if regular em protocol
      else
      {
        // first choose probabilities for add, remove, swap, move
        MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_Add,    0.08);
        MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_Remove, 0.15);
        MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_Swap,   0.5);
        MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_SSE,   1.5);

        // now assign probabilities for individual mutates
        // adds
        // if the body restraint flag is given
        if( GetFlagBodyRestraint()->GetFlag())
        {
          MUTATE_TREE.SetMutateProbability( "add_em_sse_into_body_null",   1.0);
          MUTATE_TREE.SetMutateProbability( "add_em_sse_into_body_flip_x", 1.0);
          MUTATE_TREE.SetMutateProbability( "add_em_sse_into_body_flip_y", 1.0);
          MUTATE_TREE.SetMutateProbability( "add_em_sse_into_body_flip_z", 1.0);
        }
        // removes
        MUTATE_TREE.SetMutateProbability( "remove_em_random",   0.5);
        // swaps
        MUTATE_TREE.SetMutateProbability( "swap_em_sses",          0.7);
        MUTATE_TREE.SetMutateProbability( "swap_em_sse_with_pool", 0.3);
        MUTATE_TREE.SetMutateProbability( "swap_em_sse_with_pool_overlap_single_swap_only", 0.3);
        // moves
        MUTATE_TREE.SetMutateProbability( "sse_em_flip_x",          1.0);
        MUTATE_TREE.SetMutateProbability( "sse_em_flip_y",          1.0);
        MUTATE_TREE.SetMutateProbability( "sse_em_flip_z",          1.0);
        MUTATE_TREE.SetMutateProbability( "sse_em_rotate_random_z", 1.0);

        // if the body restraint flag is given
        if( GetFlagBodyRestraint()->GetFlag())
        {
          MUTATE_TREE.SetMutateProbability( "sse_em_move_into_empty_body", 3.0);
        }

        if( sspred::Methods::GetFlagReadSSPredictions()->GetFlag())
        {
          MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSEResizeCTerm, 10.0);
          MUTATE_TREE.SetMutateProbability( DefaultMutates::GetInstance().e_SSEResizeNTerm, 10.0);
        }
      }
    }

    //! @brief get the mutate tree associated with this protocol
    //! @return the mutate tree associated with this protocol
    util::ShPtr< MutateTree> ProtocolEM::GetMutateTree() const
    {
      util::ShPtr< MutateTree> sp_mutate_tree( new MutateTree());
      ModifyMutateTree( *sp_mutate_tree);
      return sp_mutate_tree;
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void ProtocolEM::MergeAndModifyMutateTree( MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProtocolEM::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Protocol for Using EM Density Maps");
      serializer.AddInitializer
        (
         "body restraint score",
         "score associated with the body restraint",
         io::Serialization::GetAgent( &e_BodyRestraintScore)
         );
      serializer.AddInitializer
        (
         "body connectivity score",
         "score associated with body connectivity",
         io::Serialization::GetAgent( &e_BodyConnectivityDensityScore)
         );
      serializer.AddInitializer
        (
         "density agreement score",
         "score associated with density agreement",
         io::Serialization::GetAgent( &m_DensityAgreementScores)
         );
      serializer.AddInitializer
        (
         "density map",
         "density map associated with the density agreement scores",
         io::Serialization::GetAgent( &m_AgreementScoreDensityMap)
         );

      return serializer;
    }

  ////////////
  // readme //
  ////////////

    //! @brief returns string containing short description of the protocol
    //! @return string containing short description of the protocol
    const std::string &ProtocolEM::GetDescription() const
    {
      // initialize description string
      static const std::string s_description
      (
        "EM protocol folds proteins into medium resolution density maps"
      );

      // end
      return s_description;
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ProtocolEM::GetReadMe() const
    {
      // create a static string to hold readme information
      static const std::string s_readme_text
      (
        "The EM protocol, or BCL::EMFold, folds proteins into medium resolution density maps. "
        "The Monte Carlo algorithm randomly samples moves that place predicted helices "
        "from the pool into these density rods. The results of the moves are scored by "
        "knowledge-based scores as well as scores that evaluate the agreement of the model "
        "with the experimental density map.\n\n"
        "When using BCL::EMFold in a publication, please cite the following "
        "publication describing the application's development:\n\n"
        "Lindert, S., Staritzbichler, R., Wotzel, N., Karakas, M., Stewart, P.L., "
        "and Meiler, J. (2009). EM-fold: De novo folding of alpha-helical proteins "
        "guided by intermediate-resolution electron microscopy density maps. "
        "Structure 17, 990-1003.\n\n"
        "Example command line flags for assembly:            \n"
        "-protocols EM                                       \n"
        "-nmodels 1                                          \n"
        "-fasta ????A.fasta                                  \n"
        "-pool ????.pool                                     \n"
        "-mc_number_iterations 2000 400                      \n"
        "-mc_temperature_fraction 0.25 0.05                  \n"
        "-protein_storage ./pdbs/                            \n"
        "-body_restraint ????.cst_body 2.5 2.5 5.0 5.0 -1.0  \n"
        "-print_body_assignment                              \n"
        "-message_level Critical                             \n"
        "-score_weightset_read ????.score                    \n"
        "-score_density_connectivity ????.mrc                \n"
        "-sspred JUFO PSIPRED                                \n"
        "-sspred_path_prefix ./sspred/ 1X91                  \n\n"
        "Example command line flags for refinement:          \n"
        "-protocols EM                                                            \n"
        "-nmodels 1                                                               \n"
        "-native ????.pdb                                                         \n"
        "-use_native_pool                                                         \n"
        "-quality RMSD                                                            \n"
        "-start_model ????_start.pdb                                              \n"
        "-score_density_agreement ????.mrc 4.6 TrilinearInterpolationGaussian CCC \n"
        "-fasta ????A.fasta                                                       \n"
        "-mc_number_iterations 2000 400                                           \n"
        "-mc_temperature_fraction 0.25 0.05                                       \n"
        "-protein_storage ./pdbs/                                                 \n"
        "-score_weightset_read ????.score                                         \n"
        "-em_refinement                                                           \n\n"
        "Speficific flags:\n"
        "-body_restraint ????.cst_body 500 1.0 1.0 2.0 2.0 -1.0 true \nThe file containing "
        "the restraints extracted from the density map, the weight for the occupancy score, "
        "and parameters describing the shape of the potential used in the occupancy score. "
        "-score_weightset weights_assembly.score The score file detailing which scores to use. "
        "-mc_number_iterations 2000 500 The number of total steps (2000) as well as rejected "
        "steps in a row (500) before the Monte Carlo search is stopped. \n\n"
        "The restraint file encodes the density rods as helices and looks "
        "like the following:\n\n"
        "SEQRES   1 A  152  PRO PRO LYS TRP LYS VAL LYS LYS GLN LYS LEU ALA GLU\n"
        "SEQRES   2 A  152  LYS ALA ALA ARG GLU ALA GLU LEU THR ALA LYS LYS ALA\n"
        "SEQRES   3 A  152  GLN ALA ARG GLN ALA LEU SER ILE TYR LEU ASN LEU PRO\n"
        "SEQRES   4 A  152  THR LEU ASP GLU ALA VAL ASN THR LEU LYS PRO TRP TRP\n"
        "SEQRES   5 A  152  PRO GLY LEU PHE ASP GLY ASP THR PRO ARG LEU LEU ALA\n"
        "SEQRES   6 A  152  CYS GLY ILE ARG ASP VAL LEU LEU GLU ASP VAL ALA GLN\n"
        "SEQRES   7 A  152  ARG ASN ILE PRO LEU SER HIS LYS LYS LEU ARG ARG ALA\n"
        "SEQRES   8 A  152  MET LYS ALA ILE THR ARG SER GLU SER TYR LEU CYS ALA\n"
        "SEQRES   9 A  152  MET LYS ALA GLY ALA CYS ARG TYR ASP THR GLU GLY TYR\n"
        "SEQRES  10 A  152  VAL THR GLU HIS ILE SER GLN GLU GLU GLU VAL TYR ALA\n"
        "SEQRES  11 A  152  ALA GLU ARG LEU ASP LYS ILE ARG ARG GLN ASN ARG ILE\n"
        "SEQRES  12 A  152  LYS ALA GLU LEU GLN ALA VAL LEU ASP\n"
        "HELIX    1   1 PRO A    2  SER A   33  1                                  32\n"
        "HELIX    5   5 GLY A   67  ARG A   79  1                                  13\n"
        "HELIX    6   6 SER A   84  ARG A   97  1                                  14\n"
        "HELIX    8   8 SER A  123  ASP A  152  1                                  30\n"
        "ATOM      8  N   PRO A   2     -34.913  10.141  21.737  1.00 45.99           N\n"
        "ATOM      9  CA  PRO A   2     -34.039  10.544  20.632  1.00 46.41           C\n"
        "ATOM     10  C   PRO A   2     -32.544  10.299  20.814  1.00 47.36           C\n"
        "ATOM     11  O   PRO A   2     -32.078   9.922  21.893  1.00 48.61           O\n"
        "ATOM     12  CB  PRO A   2     -34.341  12.030  20.469  1.00 45.93           C\n"
        "ATOM     13  CG  PRO A   2     -35.724  12.173  21.015  1.00 46.71           C\n"
        "ATOM     14  CD  PRO A   2     -35.693  11.293  22.226  1.00 46.13           C\n"
        "ATOM     15  N   LYS A   3     -31.807  10.533  19.730  1.00 47.24           N\n"
        "....\n\n"
        "The score file specifies used scores and should look like the following:\n\n"
        "bcl::storage::Table<double>         nr         co       rgyr     aadist   aasmooth "
        "aaneigh      aavmd     aanvec    annsasa      aaols    loop   ssepack_fr  strand_fr "
        "sum       rmsd    rmsd100\n"
        "weights                              0          0          0          0          0  "
        "0          0          0          0          0     265          0         0          "
        "0          0          0"
      );

      // return readme information
      return s_readme_text;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProtocolEM::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProtocolEM::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief creates the mutates for adding SSEs into the Protein model into body restraints
    void ProtocolEM::InitializeBodyAddMutates()
    {
      // initialize ssepool picker
      const util::ShPtr
      <
        find::PickCriteriaWrapper
        <
          util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::DomainInterface
        >
      > sp_picker_pool
      (
        new find::PickCriteriaWrapper
        <
          util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::DomainInterface
        >( assemble::PickSSERandom())
      );

      // set the random picker
      const find::PickCriteriaWrapper
      <
        util::ShPtr< assemble::SSEGeometryInterface>,
        util::ShPtrVector< assemble::SSEGeometryInterface>,
        assemble::SSEGeometryInterface
      >
      picker_random =
        find::PickCriteriaWrapper
        <
          util::ShPtr< assemble::SSEGeometryInterface>,
          util::ShPtrVector< assemble::SSEGeometryInterface>,
          assemble::SSEGeometryInterface
        >( find::PickBodyRandom());

      // set extent tolerances for the length extent picker
      // these tolerances are set to discriminate between helices and strands only, no length discrimination
      linal::Vector3D extent_tolerances( 0.1, 0.1, 50.0);

      // important this is only half the length!!
//      linal::Vector3D extent_tolerances( 0.5, 0.5, 5.0);

      // set the length extent picker
      const find::PickBodyExtent picker_extent( extent_tolerances, extent_tolerances);

      // create ShPtr to a PickCriteriaInterface "pick" which is used to determine how a restraint body is
      // selected
      util::ShPtr
      <
        find::PickCriteriaInterface
        <
          util::ShPtr< assemble::SSEGeometryInterface>,
          util::ShPtrVector< assemble::SSEGeometryInterface>, assemble::SSEGeometryInterface
        >
      > pick
      (
        // picker_random.Clone()
        picker_extent.Clone()
      );

      // create MutateTransformationMatrix3DNull for placing the sse into the body with the same orientation as the
      // body
      util::ShPtr< math::MutateInterface< math::TransformationMatrix3D> > orientation_null
      (
        new restraint::MutateTransformationMatrix3DNull()
      );

      // create MutateTransformationMatrix3DRotate for placing the sse into the body with the opposite orientation
      //  as the body according to the x-axis
      util::ShPtr< math::MutateInterface< math::TransformationMatrix3D> > orientation_flip_x_axis
      (
        new restraint::MutateTransformationMatrix3DRotate( coord::GetAxes().e_X, math::g_Pi, math::g_Pi)
      );

      // create MutateTransformationMatrix3DRotate for placing the sse into the body with the opposite orientation
      //  as the body according to the Y-axis
      util::ShPtr< math::MutateInterface< math::TransformationMatrix3D> > orientation_flip_y_axis
      (
        new restraint::MutateTransformationMatrix3DRotate( coord::GetAxes().e_Y, math::g_Pi, math::g_Pi)
      );

      // create MutateTransformationMatrix3DRotate for placing the sse into the body with orientation
      //  rotated with respect to the z-axis
      util::ShPtr< math::MutateInterface< math::TransformationMatrix3D> > orientation_rotate_z_axis
      (
        new restraint::MutateTransformationMatrix3DRotate( coord::GetAxes().e_Z, 2 * math::g_Pi, 0.0)
      );

    ////////////////
    // placements //
    ////////////////

      // create placement which just places an sse into a random body with matching orientations
      util::ShPtr< PlacementInterface< assemble::SSE, assemble::ProteinModel> > placement_null
      (
        new PlacementSSEIntoBody( GetBodyRestraints()->FirstElement(), pick, orientation_null)
      );

      // create placement which places an sse into a random body with orientation flipped by x-axis
      util::ShPtr< PlacementInterface< assemble::SSE, assemble::ProteinModel> > placement_flip_x
      (
        new PlacementSSEIntoBody( GetBodyRestraints()->FirstElement(), pick, orientation_flip_x_axis)
      );

      // create placement which places an sse into a random body with orientation flipped by y-axis
      util::ShPtr< PlacementInterface< assemble::SSE, assemble::ProteinModel> > placement_flip_y
      (
        new PlacementSSEIntoBody( GetBodyRestraints()->FirstElement(), pick, orientation_flip_y_axis)
      );

      // create placement which places an sse into a random body with orientation rotated about z-axis
      util::ShPtr< PlacementInterface< assemble::SSE, assemble::ProteinModel> > placement_rotate_z
      (
        new PlacementSSEIntoBody( GetBodyRestraints()->FirstElement(), pick, orientation_rotate_z_axis)
      );

      GetMutates().AddMutate
      (
        util::ShPtr< math::MutateInterface< assemble::ProteinModel> >
        (
          new MutateProteinModelSSEAdd( sp_picker_pool, placement_null, "add_em_sse_into_body_null")
        )
      );

      GetMutates().AddMutate
      (
        util::ShPtr< math::MutateInterface< assemble::ProteinModel> >
        (
          new MutateProteinModelSSEAdd( sp_picker_pool, placement_flip_x, "add_em_sse_into_body_flip_x")
        )
      );

      GetMutates().AddMutate
      (
        util::ShPtr< math::MutateInterface< assemble::ProteinModel> >
        (
          new MutateProteinModelSSEAdd( sp_picker_pool, placement_flip_y, "add_em_sse_into_body_flip_y")
        )
      );

      GetMutates().AddMutate
      (
        util::ShPtr< math::MutateInterface< assemble::ProteinModel> >
        (
          new MutateProteinModelSSEAdd( sp_picker_pool, placement_rotate_z, "add_em_sse_into_body_flip_z")
        )
      );
    }

    //! @brief return body restraint
    //! @return body restraint
    util::ShPtr< restraint::Body> ProtocolEM::GetBodyRestraint()
    {
    /////////////////////////////////////////////////////////////////////
    // move into empty density rod (specific to EM-Fold assembly step) //
    /////////////////////////////////////////////////////////////////////

      // create stream to restraint file
      io::IFStream read;
      io::File::MustOpenIFStream( read, GetFlagBodyRestraint()->GetFirstParameter()->GetValue());
      BCL_MessageStd
      (
        "read restraint file: " + GetFlagBodyRestraint()->GetFirstParameter()->GetValue()
      );

      // construct handler body in order to read in the restraints
      // (handler body is constructed with a method for determining if a restraint body is occupied)
      restraint::HandlerBody handler_body
      (
        // HandlerBody is constructed from ShPrt to FunctionInterface
        // (is method to determine if a restraint body is occupied, takes two coord::Bodies and returns bool)
        util::ShPtr< util::BinaryFunctionInterface< assemble::SSEGeometryInterface, assemble::SSE, bool> >
        (
          // this is derived from FunctionInterface and is method to determine if a restraint body is occupied
          new restraint::ContainsBodyOrigin()
        )
      );
      // returns the restraint::Body that HandlerBody created
      // although there is a ShPtrVector, there is only one restraint::Body per HandlerBody
      // this is the size where the fragment lengths are changed!!
      util::ShPtrVector< restraint::Body> density_map( handler_body.CreateRestraintsBody( read));
      io::File::CloseClearFStream( read);
      return density_map.FirstElement();
    }

    //! @brief GetScoresDensityAgreement creates the score objects for density agreement scores
    util::ShPtrVector< score::ProteinModel>
    ProtocolEM::GetScoresDensityAgreement()
    {
      const util::ShPtr< command::FlagStaticAndDynamic> agreement_scores_flag( GetFlagScoreDensityAgreement());

      // all parameters
      const util::ShPtrVector< command::ParameterInterface> &parameters
      (
        agreement_scores_flag->GetParameterList()
      );

      // mrc filename
      const std::string mrc_filename( parameters( 0)->GetValue());

      // resolution
      const double resolution_density_map( parameters( 1)->GetNumericalValue< double>());

      // density simulator
      const density::Simulator simulator( parameters( 2)->GetValue());

      // density agreement score
      storage::Set< density::ProteinAgreement> agreement_scores;
      for( util::ShPtrVector< command::ParameterInterface>::const_iterator itr( parameters.Begin() + 3), itr_end( parameters.End()); itr != itr_end; ++itr)
      {
        agreement_scores.Insert( density::ProteinAgreement( ( *itr)->GetValue()));
      }

      // protein model scores
      util::ShPtrVector< score::ProteinModel> protein_model_scores;

      if( agreement_scores.IsEmpty())
      {
        return protein_model_scores;
      }

      // instantiate Densitymap from mrc file
      io::IFStream read;
      m_AgreementScoreDensityMap = density::Map();
      BCL_MessageStd( "open density map: " + mrc_filename);
      BCL_Assert
      (
        io::File::TryOpenIFStream( read, mrc_filename, std::ios::binary),
        "unable to open file " + util::Format()( mrc_filename)
      );
      m_AgreementScoreDensityMap.ReadMRC( read, 0);
      io::File::CloseClearFStream( read);
      BCL_MessageStd( "density map has been read: " + mrc_filename);

      // SiPtr to m_AgreementScoreDensityMap
      util::SiPtr< density::Map> sp_map( new density::Map( m_AgreementScoreDensityMap));
      // iterate over all requested agreement scores and construct a protein model score
      for
      (
        storage::Set< density::ProteinAgreement>::const_iterator
          itr( agreement_scores.Begin()), itr_end( agreement_scores.End());
        itr != itr_end;
        ++itr
      )
      {
        protein_model_scores.PushBack
        (
          density::GetProteinAgreements().CreateProteinAgreement
          (
            *itr, simulator, sp_map, resolution_density_map
          )
        );
      }

      // return all scores
      return protein_model_scores;
    }

    //! @brief GetScoreBodyConnectivityDensity creates the score object for density connectivity score
    util::ShPtr< score::ProteinModel>
    ProtocolEM::GetScoreBodyConnectivityDensity()
    {
      // get restraints
      const util::ShPtr< restraint::Body> body_restraints( GetBodyRestraint());

      // instantiate Densitymap from mrc file
      io::IFStream read;
      util::ShPtr< density::Map> density_map( new density::Map);
      const std::string mrc_filename( GetFlagScoreDensityConnectivity()->GetFirstParameter()->GetValue());

      BCL_MessageStd( "open density map: " + mrc_filename);
      BCL_Assert
      (
        io::File::TryOpenIFStream( read, mrc_filename, std::ios::binary),
        "unable to open file " + util::Format()( mrc_filename)
      );
      density_map->ReadMRC( read, 0);
      io::File::CloseClearFStream( read);

      BCL_MessageStd( "density map has been read: " + mrc_filename);

      // construct density connectivity score
      util::ShPtr< score::ProteinModel> density_connectivity_score
      (
        new score::BodyConnectivityDensity( body_restraints, density_map)
      );

      //return density connectivity score
      return density_connectivity_score;
    }

    //! @brief GetBodyRestraintScore creates the scoring object for use with the body restraint
    //! @return the scoring object for use with the body restraint
    util::ShPtr< score::ProteinModel> ProtocolEM::GetScoreBodyRestraint()
    {
      const double upper_tolerance( GetFlagBodyRestraint()->GetParameterList()( 1)->GetNumericalValue< double>());
      const double lower_tolerance( GetFlagBodyRestraint()->GetParameterList()( 2)->GetNumericalValue< double>());
      const double upper_transition( GetFlagBodyRestraint()->GetParameterList()( 3)->GetNumericalValue< double>());
      const double lower_transition( GetFlagBodyRestraint()->GetParameterList()( 4)->GetNumericalValue< double>());
      const double energy_well_depth( GetFlagBodyRestraint()->GetParameterList()( 5)->GetNumericalValue< double>());

      // create score::BodyAssignment "assignment_score" as method for scoring agreement with a body restraint
      score::BodyAssignment assignment_score
      (
        util::ShPtr< math::BinaryFunctionInterface< assemble::SSEGeometryInterface, assemble::SSE, double> >
        (
          new score::BodyExtentAgreement
          (
            lower_tolerance,     //< lower tolerance
            lower_transition,    //< lower transition width
            upper_tolerance,     //< upper tolerance
            upper_transition,    //< upper transition width
            energy_well_depth,   //< energy well depth
            coord::GetAxes().e_Z //< axis of interest
          )
        )
      );

      // create score::Restraint "score_restraints" and return it
      return
      util::ShPtr< score::ProteinModel>
      (
        new score::RestraintBodyProteinModel
        (
          GetBodyRestraints(),
          // initialize with ShPtr to a function interface "assignment_score"
          assignment_score
        )
      );
    }

    //! @brief singleton function for creating static body restraints to be used during folding
    const util::ShPtr< util::ShPtrVector< restraint::Body> >
    &ProtocolEM::GetBodyRestraints()
    {
      // declare static variable restraints
      static
      util::ShPtr< util::ShPtrVector< restraint::Body> > restraints;

      // if restraints are not initialized
      if( !restraints.IsDefined())
      {
        // create HandlerAtomDistanceAssigned "handler" as method for determining if a restraint body is occupied
        restraint::HandlerBody handler
        (
          util::ShPtr< util::BinaryFunctionInterface< assemble::SSEGeometryInterface, assemble::SSE, bool> >
          (
            new restraint::ContainsBodyOrigin()
          )
        );

        // create string "restraint_filename" which has path for example restraint file
        const std::string restraint_filename( GetFlagBodyRestraint()->GetFirstParameter()->GetValue());

        // create stream to restraint file
        io::IFStream read;
        io::File::MustOpenIFStream( read, restraint_filename);
        BCL_MessageStd( "read restraint file: " + restraint_filename);

        // create restraints
        // initialize ShPtrVector with restraints from "restraint_filename"
        restraints = util::CloneToShPtr( handler.CreateRestraintsBody( read));

        io::File::CloseClearFStream( read);
      }

      // return
      return restraints;
    }

  } // namespace fold
} // namespace bcl
