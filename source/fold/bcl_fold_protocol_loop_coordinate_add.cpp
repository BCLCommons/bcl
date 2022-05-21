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
#include "fold/bcl_fold_protocol_loop_coordinate_add.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_collector_sse.h"
#include "assemble/bcl_assemble_pick_sse_random.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "find/bcl_find_locator.h"
#include "fold/bcl_fold_collector_unconnected_sses.h"
#include "fold/bcl_fold_default_scores.h"
#include "fold/bcl_fold_handler_locator_loop_domain.h"
#include "fold/bcl_fold_locator_loop_domain.h"
#include "fold/bcl_fold_mutate_protein_model_loop_resize.h"
#include "fold/bcl_fold_mutate_protein_model_move_aa.h"
#include "fold/bcl_fold_mutate_protein_model_sse.h"
#include "fold/bcl_fold_mutate_protein_model_sse_seed.h"
#include "fold/bcl_fold_mutate_sse_bend_random.h"
#include "fold/bcl_fold_mutate_tree.h"
#include "fold/bcl_fold_phi_psi_generator_ramachandran.h"
#include "fold/bcl_fold_score_weight_set.h"
#include "io/bcl_io_file.h"
#include "score/bcl_score_aa_pair_atom_clash.h"
#include "score/bcl_score_aa_sequence.h"
#include "score/bcl_score_aa_sequence_pair.h"
#include "score/bcl_score_loop_closure.h"
#include "score/bcl_score_phi_psi.h"
#include "score/bcl_score_protein_model_sse.h"
#include "score/bcl_score_protein_model_sse_neighbors.h"
#include "score/bcl_score_protein_model_sse_pairs.h"
#include "score/bcl_score_sse_pair_connectivity.h"
#include "score/bcl_score_sse_pair_gap.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProtocolLoopCoordinateAdd::ProtocolLoopCoordinateAdd()
    {
    }

    //! @brief Clone function
    //! @return pointer to new ProtocolLoopCoordinateAdd
    ProtocolLoopCoordinateAdd *ProtocolLoopCoordinateAdd::Clone() const
    {
      return new ProtocolLoopCoordinateAdd( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief access to only instance
    //! @return reference to only instance
    ProtocolLoopCoordinateAdd &ProtocolLoopCoordinateAdd::GetInstance()
    {
      static ProtocolLoopCoordinateAdd s_protocol_instance;
      return s_protocol_instance;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProtocolLoopCoordinateAdd::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ProtocolLoopCoordinateAdd::GetAlias() const
    {
      static const std::string s_name( "ProtocolLoopClose");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProtocolLoopCoordinateAdd::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "protocol for closing loops");
      serializer.AddInitializer
        (
         "mutate coil bend n term phi psi 45",
         "phi psi mutate for n terminus 45 degree",
         io::Serialization::GetAgent( &e_MutateCoilBendNTermPhiPsi45)
         );
      serializer.AddInitializer
        (
         "mutate coil bend n term phi 45",
         "phi mutate for n terminus 45 degree",
         io::Serialization::GetAgent( &e_MutateCoilBendNTermPhi45)
         );
      serializer.AddInitializer
        (
         "mutate coil bend n term psi 45",
         "psi mutate for n terminus 45 degree",
         io::Serialization::GetAgent( &e_MutateCoilBendNTermPsi45)
         );
      serializer.AddInitializer
        (
         "mutate coil bend c term phi psi 45",
         "phi psi mutate for c terminus 45 degree",
         io::Serialization::GetAgent( &e_MutateCoilBendCTermPhiPsi45)
         );
      serializer.AddInitializer
        (
         "mutate coil bend c term phi 45",
         "phi mutate for c terminus 45 degree",
         io::Serialization::GetAgent( &e_MutateCoilBendCTermPhi45)
         );
      serializer.AddInitializer
        (
         "mutate coil bend c term psi 45",
         "psi mutate for c terminus 45 degree",
         io::Serialization::GetAgent( &e_MutateCoilBendCTermPsi45)
         );
      serializer.AddInitializer
        (
         "mutate coil bend n term phi psi 5",
         "phi psi mutate for n terminus 5 degree",
         io::Serialization::GetAgent( &e_MutateCoilBendNTermPhiPsi5)
         );
      serializer.AddInitializer
        (
         "mutate coil bend c term phi psi 5",
         "phi psi mutate for c terminus 5 degree",
         io::Serialization::GetAgent( &e_MutateCoilBendCTermPhiPsi5)
         );
      serializer.AddInitializer
        (
         "mutate coil resize N term",
         "resize mutate for coil N terminus",
         io::Serialization::GetAgent( &e_MutateCoilResizeNTerm)
         );
      serializer.AddInitializer
        (
         "mutate coil resize C term",
         "resize mutate for coil C terminus",
         io::Serialization::GetAgent( &e_MutateCoilResizeCTerm)
         );
      serializer.AddInitializer
        (
         "mutate helix strand resize n term",
         "resize mutate for helix strand N terminus",
         io::Serialization::GetAgent( &e_MutateHelixStrandResizeNTerm)
         );
      serializer.AddInitializer
        (
         "mutate helix strand resize c term",
         "resize mutate for helix strand C terminus",
         io::Serialization::GetAgent( &e_MutateHelixStrandResizeCTerm)
         );
      serializer.AddInitializer
        (
         "mutate SSE seed N term",
         "add a coil as seed to an sse on left side",
         io::Serialization::GetAgent( &e_MutateSSESeedNTerm)
         );
      serializer.AddInitializer
      (
        "mutate SSE seed C term",
        "add a coil as seed to an sse on right side",
        io::Serialization::GetAgent( &e_MutateSSESeedCTerm)
       );
      serializer.AddInitializer
      (
        "mutate SSE seed cut N term",
        "add a coil as seed to an sse on left side cutting into the sse",
        io::Serialization::GetAgent( &e_MutateSSESeedCutNTerm)
       );
      serializer.AddInitializer
      (
        "mutate SSE seed cut C term",
        "add a coil as seed to an sse on right side cutting into the sse",
        io::Serialization::GetAgent( &e_MutateSSESeedCutCTerm)
       );
      serializer.AddInitializer
      (
        "mutate move AA N term",
        "move aas from one sse to the neighboring in nterm direction",
        io::Serialization::GetAgent( &e_MutateMoveAANTerm)
       );
      serializer.AddInitializer
      (
        "mutate move AA C term",
        "move aas from one sse to the neighboring in cterm direction",
        io::Serialization::GetAgent( &e_MutateMoveAACTerm)
       );
      serializer.AddInitializer
      (
        "score SSE interatom clash",
        "atom clash inter sse",
        io::Serialization::GetAgent( &e_ScoreSSEInterAtomClash)
       );
      serializer.AddInitializer
      (
        "score SSE intraatom clash",
        "atom clash intra sse",
        io::Serialization::GetAgent( &e_ScoreSSEIntraAtomClash)
       );
      serializer.AddInitializer
      (
        "score loop closure",
        "score the distance between two sse ends",
        io::Serialization::GetAgent( &e_ScoreSSEIntraAtomClash)
       );
      serializer.AddInitializer
      (
        "score phi psi coil",
        "phi psi",
        io::Serialization::GetAgent( &e_ScorePhiPsiCoil)
       );
      serializer.AddInitializer
      (
        "score sse pair gap",
        "scores the gap between two neighboring sses",
        io::Serialization::GetAgent( &e_ScoreSSEPairGap)
       );
      serializer.AddInitializer
      (
        "score sse pair connectivity",
        "scores the connectivity between two neighboring sses",
        io::Serialization::GetAgent( &e_ScoreSSEPairConnectivity)
       );

      return serializer;
    }

    //! @brief returns all flags that are specialized for this protocol
    //! @return all flags that are specialized for this protocol
    const util::ShPtrVector< command::FlagInterface> &ProtocolLoopCoordinateAdd::GetAllFlags() const
    {
      // initialize static ShPtrVector of FlagInterfaces to form the comment line
      static util::ShPtrVector< command::FlagInterface> s_all_flags_vector;

      // if the flag vector is initialize for the first time
      if( s_all_flags_vector.IsEmpty())
      {
        // insert all the flags in the vector
        s_all_flags_vector.PushBack( GetFlagLoopDomainFilename());
      }

      // end
      return s_all_flags_vector;
    }

    //! @brief return command line flag for specifying loop domains via a file
    //! @return command line flag for specifying loop domains via a file
    util::ShPtr< command::FlagInterface> &ProtocolLoopCoordinateAdd::GetFlagLoopDomainFilename()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "loop_coord_add",
          "\tSpecific loop segment(s) should be built. If this flag is not given, then all loops missing coordinates will be built."
        )
      );

      // initialize parameters
      static util::ShPtr< command::ParameterInterface> s_filename
      (
        new command::Parameter
        (
          "loop_domain_filename", "\tfilename for loop domain file", "loop_domain.txt"
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

  ////////////////
  // operations //
  ////////////////

    //! @brief modifies the start model
    //! @param START_MODEL Start model to be modified
    void ProtocolLoopCoordinateAdd::ModifyStartModel( assemble::ProteinModel &START_MODEL) const
    {
      if( !GetFlagLoopDomainFilename()->GetFlag())
      {
        return;
      }

      // read the given file
      io::IFStream read;
      io::File::MustOpenIFStream( read, GetFlagLoopDomainFilename()->GetFirstParameter()->GetValue());

      // get the sses locators read from file
      util::ShPtrList< LocatorLoopDomain> locator_loop_domain
      (
        HandlerLocatorLoopDomain().HandleReadMultiple( read, START_MODEL)
      );

      BCL_Assert( !locator_loop_domain.IsEmpty(), "locator_loop_domain is empty");
      // get all the sse locators from all of the domains
      for
      (
        util::ShPtrList< LocatorLoopDomain>::const_iterator
          domain_itr( locator_loop_domain.Begin()), domain_itr_end( locator_loop_domain.End());
        domain_itr != domain_itr_end;
        ++domain_itr
      )
      {
        // remove all sses from start model that are defined in the domain locator
        util::SiPtrVector< const assemble::LocatorSSE> sse_locators( ( *domain_itr)->GetLocatorSSEs());
        for( util::SiPtrVector< const assemble::LocatorSSE>::const_iterator itr( sse_locators.Begin()), itr_end( sse_locators.End()); itr != itr_end; ++itr)
        {
          // locate the sse and remove it if it is defined
          const util::SiPtr< const assemble::SSE> located_sse( ( *itr)->Locate( START_MODEL));
          if( located_sse.IsDefined())
          {
            START_MODEL.Remove( *located_sse);
          }
        }
      }
    }

    //! @brief initialize the scores and add them to Scores enumerator
    void ProtocolLoopCoordinateAdd::InitializeScores()
    {
      // initialize score only once
      if( e_ScoreSSEInterAtomClash.IsDefined())
      {
        return;
      }

      // between sses aa clash checking all atoms of each aa
      e_ScoreSSEInterAtomClash = GetScores().AddScore
      (
        Scores::WrapCacheProteinModelScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::ProteinModelSSEPairs
            (
              Scores::WrapCacheSSEPairScore
              (
                util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, assemble::SSE, double> >
                (
                  new score::AASequencePair
                  (
                    score::AAPairAtomClash( 0.2, 0, "sse_inter_atom_clash"), false
                  )
                ),
                false // symmetric
              ),
              false // no normalization
            )
          )
        )
      );

      // sse internal clashes between all atoms of the residues
      e_ScoreSSEIntraAtomClash = GetScores().AddScore
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
                  new score::AASequence
                  (
                    util::ShPtr< score::AAPairDistanceInterface>
                    (
                      new score::AAPairAtomClash( 0.2, 0, "sse_intra_atom_clash")
                    ),
                    false
                  )
                )
              ),
              false // no normalization
            )
          )
        )
      );

      // loop closure score
      e_ScoreLoopClosure = GetScores().AddScore
      (
        Scores::WrapCacheProteinModelScore
        (
          util::ShPtr< score::ProteinModelSSEPairs>
          (
            new score::ProteinModelSSEPairs
            (
              Scores::WrapCacheSSEPairScore
              (
                util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, assemble::SSE, double> >
                (
                  new score::LoopClosure( 0, 15.0, 0.8, true)
                ),
                false // non symmetric
              ),
              false // no normalization
            )
          )
        )
      );

      // phi psi for loop building
      e_ScorePhiPsiCoil = GetScores().AddScore
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
                  new score::PhiPsi
                  (
                    "phipsi_coil",
                    score::PhiPsi::GetDefaultHistogramFilename(),
                    storage::Set< biol::SSType>( biol::GetSSTypes().COIL)
                  )
                )
              ),
              false
            )
          )
        )
      );

      // between sses aa clash checking all atoms of each aa
      e_ScoreSSEPairGap = GetScores().AddScore
      (
        Scores::WrapCacheProteinModelScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::ProteinModelSSENeighbors
            (
              Scores::WrapCacheSSEPairScore
              (
                util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, assemble::SSE, double> >
                (
                  new score::SSEPairGap( 2.0, "sse_pair_gap")
                ),
                false // non symmetric
              ),
              false // no normalization
            )
          )
        )
      );

      // sse connectivity to help guide grown loops towards their correct neighboring sse
      e_ScoreSSEPairConnectivity = GetScores().AddScore
      (
        Scores::WrapCacheProteinModelScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::ProteinModelSSENeighbors
            (
              Scores::WrapCacheSSEPairScore
              (
                util::CloneToShPtr( score::SSEPairConnectivity( score::SSEPairConnectivity::GetDefaultScheme())),
                false // non symmetric
              ),
              false // normalization
            )
          )
        )
      );

      DefaultScores::GetInstance().InitializeScores();
    }

    //! @brief modify the score weight set
    //! @param SCORE_WEIGHT_SET Score weight set
    void ProtocolLoopCoordinateAdd::ModifyScoreWeightSet( ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
      SCORE_WEIGHT_SET.Reset();

      // default scores
      SCORE_WEIGHT_SET.SetWeight( DefaultScores::GetInstance().e_ScoreAAPairClash        ,  10.0);
      SCORE_WEIGHT_SET.SetWeight( DefaultScores::GetInstance().e_ScoreAAPairDistance     ,   0.5);
      SCORE_WEIGHT_SET.SetWeight( DefaultScores::GetInstance().e_ScoreAANeighborCount    ,  50.0);

      SCORE_WEIGHT_SET.SetWeight(                              e_ScoreSSEInterAtomClash  ,  10.0);
      SCORE_WEIGHT_SET.SetWeight(                              e_ScoreSSEIntraAtomClash  ,  10.0);
      SCORE_WEIGHT_SET.SetWeight(                              e_ScorePhiPsiCoil         , 100.0);
      SCORE_WEIGHT_SET.SetWeight(                              e_ScoreSSEPairGap         , 100.0);
      SCORE_WEIGHT_SET.SetWeight(                              e_ScoreSSEPairConnectivity,   2.5);
      SCORE_WEIGHT_SET.SetWeight(                              e_ScoreLoopClosure        , 500.0);
    }

    //! @brief modify the terminate object
    //! @param CRITERION which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolLoopCoordinateAdd::ModifyCriterion
    (
      opti::CriterionCombine< assemble::ProteinModel, double> &CRITERION,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief modify the printer object
    //! @param PRINTER which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolLoopCoordinateAdd::ModifyPrinter
    (
      mc::PrinterCombined< assemble::ProteinModel, double> &PRINTER,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief modify the pdb factory object
    //! @param FACTORY pdb factory to be modified
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolLoopCoordinateAdd::ModifyFactory
    (
      util::ShPtr< pdb::Factory> &FACTORY,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief initialize the mutates and add them to Mutates enumerator
    void ProtocolLoopCoordinateAdd::InitializeMutates()
    {
      if( e_MutateCoilResizeNTerm.IsDefined())
      {
        return;
      }

      typedef find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> LocatorUnconnectedSSE;
      const util::ShPtr< assemble::PickSSERandom> sp_pick_sse_random( new assemble::PickSSERandom());

      const util::ShPtr< CollectorUnconnectedSSE> sp_collector_nterminal_flexible
      (
        new CollectorUnconnectedSSE
        (
          biol::AASequenceFlexibility::e_NTerminal,
          true,
          storage::Set< biol::SSType>( biol::GetSSTypes().COIL),
          false
        )
      );
      const util::ShPtr< LocatorUnconnectedSSE> sp_locator_nterminal_flexible
      (
        new find::Locator< util::SiPtr< const assemble::SSE>, assemble::DomainInterface, util::SiPtrList< const assemble::SSE> >
        (
          sp_collector_nterminal_flexible,
          sp_pick_sse_random
        )
      );
      const util::ShPtr< CollectorUnconnectedSSE> sp_collector_cterminal_flexible
      (
        new CollectorUnconnectedSSE
        (
          biol::AASequenceFlexibility::e_CTerminal,
          true,
          storage::Set< biol::SSType>( biol::GetSSTypes().COIL),
          false
        )
      );
      const util::ShPtr< LocatorUnconnectedSSE> sp_locator_cterminal_flexible
      (
        new find::Locator< util::SiPtr< const assemble::SSE>, assemble::DomainInterface, util::SiPtrList< const assemble::SSE> >
        (
          sp_collector_cterminal_flexible,
          sp_pick_sse_random
        )
      );

      const util::ShPtr< CollectorUnconnectedSSE> sp_coil_collector_nterminal
      (
        new CollectorUnconnectedSSE
        (
          biol::AASequenceFlexibility::e_NTerminal,
          false,
          storage::Set< biol::SSType>( biol::GetSSTypes().COIL),
          false
        )
      );
      const util::ShPtr< LocatorUnconnectedSSE> sp_locator_nterminal
      (
        new find::Locator< util::SiPtr< const assemble::SSE>, assemble::DomainInterface, util::SiPtrList< const assemble::SSE> >
        (
          sp_coil_collector_nterminal,
          sp_pick_sse_random
        )
      );
      const util::ShPtr< CollectorUnconnectedSSE> sp_coil_collector_cterminal
      (
        new CollectorUnconnectedSSE
        (
          biol::AASequenceFlexibility::e_CTerminal,
          false,
          storage::Set< biol::SSType>( biol::GetSSTypes().COIL),
          false
        )
      );
      const util::ShPtr< LocatorUnconnectedSSE> sp_locator_cterminal
      (
        new find::Locator< util::SiPtr< const assemble::SSE>, assemble::DomainInterface, util::SiPtrList< const assemble::SSE> >
        (
          sp_coil_collector_cterminal,
          sp_pick_sse_random
        )
      );

      const util::ShPtr< CollectorUnconnectedSSE> sp_collector_nterminal_helix_strand
      (
        new CollectorUnconnectedSSE
        (
          biol::AASequenceFlexibility::e_NTerminal,
          false,
          storage::Set< biol::SSType>( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND),
          false
        )
      );
      const util::ShPtr< LocatorUnconnectedSSE> sp_locator_nterminal_helix_strand
      (
        new find::Locator< util::SiPtr< const assemble::SSE>, assemble::DomainInterface, util::SiPtrList< const assemble::SSE> >
        (
          sp_collector_nterminal_helix_strand,
          sp_pick_sse_random
        )
      );
      const util::ShPtr< CollectorUnconnectedSSE> sp_collector_cterminal_helix_strand
      (
        new CollectorUnconnectedSSE
        (
          biol::AASequenceFlexibility::e_CTerminal,
          false,
          storage::Set< biol::SSType>( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND),
          false
        )
      );
      const util::ShPtr< LocatorUnconnectedSSE> sp_locator_cterminal_helix_strand
      (
        new find::Locator< util::SiPtr< const assemble::SSE>, assemble::DomainInterface, util::SiPtrList< const assemble::SSE> >
        (
          sp_collector_cterminal_helix_strand,
          sp_pick_sse_random
        )
      );

      const util::ShPtr< math::FunctionInterfaceSerializable< MutationResidue, storage::VectorND< 2, double> > > phi_psi_generator
      (
        new PhiPsiGeneratorRamachandran( biol::Ramachandran::GetDefaultHistogramFilename())
      );

      const double extend_probability( 0.75);
      const math::Range< size_t> min_max_change_range( 1, 2);

      const storage::Map< biol::SSType, size_t> min_sse_sizes
      (
        storage::Map< biol::SSType, size_t>::Create( std::pair< biol::SSType, size_t>( biol::GetSSTypes().COIL, 1))
      );

      const storage::Map< biol::SSType, size_t> max_sse_sizes
      (
        storage::Map< biol::SSType, size_t>::Create( std::pair< biol::SSType, size_t>( biol::GetSSTypes().COIL, 999))
      );

      // grow shrink n terminal loops
      e_MutateCoilResizeNTerm = GetMutates().AddMutate
      (
        MutateProteinModelLoopResize
        (
          sp_locator_nterminal,
          extend_probability,
          min_max_change_range,
          min_sse_sizes,
          max_sse_sizes,
          biol::AASequenceFlexibility::e_NTerminal,
          phi_psi_generator,
          true,
          MutateTree::GetMutateTypeName( MutateTree::e_Add) + "_nterminal_resize"
        )
      );

      // grow shrink c terminal loops
      e_MutateCoilResizeCTerm = GetMutates().AddMutate
      (
        MutateProteinModelLoopResize
        (
          sp_locator_cterminal,
          extend_probability,
          min_max_change_range,
          min_sse_sizes,
          max_sse_sizes,
          biol::AASequenceFlexibility::e_CTerminal,
          phi_psi_generator,
          true,
          MutateTree::GetMutateTypeName( MutateTree::e_Add) + "_cterminal_resize"
        )
      );

      const storage::Map< biol::SSType, size_t> min_helix_strand_sizes
        (
          storage::Map< biol::SSType, size_t>::Create
          (
            std::pair< biol::SSType, size_t>( biol::GetSSTypes().HELIX, 5),
            std::pair< biol::SSType, size_t>( biol::GetSSTypes().STRAND, 3)
          )
        );

      // extend probability
      const double sse_extend_prob( 0.0);

      // grow shrink n terminal helix or strand
      e_MutateHelixStrandResizeNTerm = GetMutates().AddMutate
      (
        MutateProteinModelSSEResize
        (
          sp_locator_nterminal_helix_strand,
          sse_extend_prob,
          math::Range< size_t>( 1, 1),
          biol::AASequenceFlexibility::e_NTerminal,
          false,
          min_helix_strand_sizes,
          MutateTree::GetMutateTypeName( MutateTree::e_Add) + "_nterm_sse_resize"
        )
      );

      // grow shrink c terminal helix or strand
      e_MutateHelixStrandResizeCTerm = GetMutates().AddMutate
      (
        MutateProteinModelSSEResize
        (
          sp_locator_cterminal_helix_strand,
          sse_extend_prob,
          math::Range< size_t>( 1, 1),
          biol::AASequenceFlexibility::e_CTerminal,
          false,
          min_helix_strand_sizes,
          MutateTree::GetMutateTypeName( MutateTree::e_Add) + "_cterm_sse_resize"
        )
      );

      // permitted change ranges for phi and psi
      const math::Range< double> empty_change_range( 0, 0);
      const math::Range< double> phi_angle_change_range45( -0.25 * math::g_Pi, 0.25 * math::g_Pi);
      const math::Range< double> psi_angle_change_range45( -0.25 * math::g_Pi, 0.25 * math::g_Pi);
      const math::Range< double> phi_angle_change_range5( - 1.0 / 36.0 * math::g_Pi, 1.0 / 36.0 * math::g_Pi);
      const math::Range< double> psi_angle_change_range5( - 1.0 / 36.0 * math::g_Pi, 1.0 / 36.0 * math::g_Pi);

      // benders
      const util::ShPtr< MutateSSEBendRandom> sp_bend_phi_psi_nterm45
      (
        new MutateSSEBendRandom( phi_angle_change_range45, psi_angle_change_range45, biol::AASequenceFlexibility::e_NTerminal)
      );
      const util::ShPtr< MutateSSEBendRandom> sp_bend_phi_nterm45
      (
        new MutateSSEBendRandom( phi_angle_change_range45, empty_change_range, biol::AASequenceFlexibility::e_NTerminal)
      );
      const util::ShPtr< MutateSSEBendRandom> sp_bend_psi_nterm45
      (
        new MutateSSEBendRandom( empty_change_range, psi_angle_change_range45, biol::AASequenceFlexibility::e_NTerminal)
      );
      const util::ShPtr< MutateSSEBendRandom> sp_bend_phi_psi_cterm45
      (
        new MutateSSEBendRandom( phi_angle_change_range45, psi_angle_change_range45, biol::AASequenceFlexibility::e_CTerminal)
      );
      const util::ShPtr< MutateSSEBendRandom> sp_bend_phi_cterm45
      (
        new MutateSSEBendRandom( phi_angle_change_range45, empty_change_range, biol::AASequenceFlexibility::e_CTerminal)
      );
      const util::ShPtr< MutateSSEBendRandom> sp_bend_psi_cterm45
      (
        new MutateSSEBendRandom( empty_change_range, psi_angle_change_range45, biol::AASequenceFlexibility::e_CTerminal)
      );
      const util::ShPtr< MutateSSEBendRandom> sp_bend_phi_psi_nterm5
      (
        new MutateSSEBendRandom( phi_angle_change_range5, psi_angle_change_range5, biol::AASequenceFlexibility::e_NTerminal)
      );
      const util::ShPtr< MutateSSEBendRandom> sp_bend_phi_psi_cterm5
      (
        new MutateSSEBendRandom( phi_angle_change_range5, psi_angle_change_range5, biol::AASequenceFlexibility::e_CTerminal)
      );

      // add bending random n terminal loops
      e_MutateCoilBendNTermPhiPsi45 = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_nterminal_flexible,
          sp_bend_phi_psi_nterm45,
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_sse_bend_nterm_phi_psi_45"
        )
      );
      e_MutateCoilBendNTermPhiPsi5 = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_nterminal_flexible,
          sp_bend_phi_psi_nterm5,
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_sse_bend_nterm_phi_psi_5"
        )
      );
      e_MutateCoilBendNTermPhi45 = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_nterminal_flexible,
          sp_bend_phi_nterm45,
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_sse_bend_nterm_phi_45"
        )
      );
      e_MutateCoilBendNTermPsi45 = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_nterminal_flexible,
          sp_bend_psi_nterm45,
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_sse_bend_nterm_psi_45"
        )
      );

      // add bending random nonrigid
      e_MutateCoilBendCTermPhiPsi45 = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_cterminal_flexible,
          sp_bend_phi_psi_cterm45,
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_sse_bend_cterm_phi_psi_45"
        )
      );
      e_MutateCoilBendCTermPhiPsi5 = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_cterminal_flexible,
          sp_bend_phi_psi_cterm5,
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_sse_bend_cterm_phi_psi_5"
        )
      );
      e_MutateCoilBendCTermPhi45 = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_cterminal_flexible,
          sp_bend_phi_cterm45,
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_sse_bend_cterm_phi_45"
        )
      );
      e_MutateCoilBendCTermPsi45 = GetMutates().AddMutate
      (
        MutateProteinModelSSE
        (
          sp_locator_cterminal_flexible,
          sp_bend_psi_cterm45,
          MutateTree::GetMutateTypeName( MutateTree::e_SSE) + "_sse_bend_cterm_psi_45"
        )
      );

      // cut in range
      const math::Range< size_t> no_cut_range( 0, 0);
      const math::Range< size_t> cut_in_range( 1, 2);
      const math::Range< size_t> seed_length_range( 1, 2);
      const math::Range< size_t> seed_length_range_cut( 0, 1);
      e_MutateSSESeedNTerm = GetMutates().AddMutate
      (
        MutateProteinModelSSESeed
        (
          sp_locator_nterminal_helix_strand,
          seed_length_range,
          biol::AASequenceFlexibility::e_NTerminal,
          sp_bend_phi_psi_nterm45,
          no_cut_range,
          MutateTree::GetMutateTypeName( MutateTree::e_Add) + "_seed_coil_nterm"
        )
      );
      e_MutateSSESeedCTerm = GetMutates().AddMutate
      (
        MutateProteinModelSSESeed
        (
          sp_locator_cterminal_helix_strand,
          seed_length_range,
          biol::AASequenceFlexibility::e_CTerminal,
          sp_bend_phi_psi_cterm45,
          no_cut_range,
          MutateTree::GetMutateTypeName( MutateTree::e_Add) + "_seed_coil_cterm"
        )
      );
      e_MutateSSESeedCutNTerm = GetMutates().AddMutate
      (
        MutateProteinModelSSESeed
        (
          sp_locator_nterminal_helix_strand,
          seed_length_range_cut,
          biol::AASequenceFlexibility::e_NTerminal,
          sp_bend_phi_psi_nterm45,
          cut_in_range,
          MutateTree::GetMutateTypeName( MutateTree::e_Add) + "_seedcut_coil_nterm"
        )
      );
      e_MutateSSESeedCutCTerm = GetMutates().AddMutate
      (
        MutateProteinModelSSESeed
        (
          sp_locator_cterminal_helix_strand,
          seed_length_range_cut,
          biol::AASequenceFlexibility::e_CTerminal,
          sp_bend_phi_psi_cterm45,
          cut_in_range,
          MutateTree::GetMutateTypeName( MutateTree::e_Add) + "_seedcut_coil_cterm"
        )
      );

      const util::ShPtr< assemble::CollectorSSE> sp_coil_collector( new assemble::CollectorSSE( biol::GetSSTypes().COIL));
      const util::ShPtr< find::Locator< util::SiPtr< const assemble::SSE>, assemble::DomainInterface, util::SiPtrList< const assemble::SSE> > > sp_coil_locator
      (
        new find::Locator< util::SiPtr< const assemble::SSE>, assemble::DomainInterface, util::SiPtrList< const assemble::SSE> >
        (
          sp_coil_collector,
          sp_pick_sse_random
        )
      );
      const math::Range< size_t> nr_res_move( 1, 2);
      storage::Map< biol::SSType, size_t> min_sse_sizes_move;
      min_sse_sizes_move[ biol::GetSSTypes().COIL] = 0;
      min_sse_sizes_move[ biol::GetSSTypes().HELIX] = 5;
      min_sse_sizes_move[ biol::GetSSTypes().STRAND] = 3;

      e_MutateMoveAANTerm = GetMutates().AddMutate
      (
        MutateProteinModelMoveAA
        (
          sp_locator_cterminal_flexible,
          nr_res_move,
          biol::AASequenceFlexibility::e_NTerminal,
          min_sse_sizes_move,
          sp_bend_phi_psi_cterm45,
          MutateTree::GetMutateTypeName( MutateTree::e_Model) + "_move_aa_cterm"
        )
      );
      e_MutateMoveAACTerm = GetMutates().AddMutate
      (
        MutateProteinModelMoveAA
        (
          sp_locator_nterminal_flexible,
          nr_res_move,
          biol::AASequenceFlexibility::e_CTerminal,
          min_sse_sizes_move,
          sp_bend_phi_psi_nterm45,
          MutateTree::GetMutateTypeName( MutateTree::e_Model) + "_move_aa_nterm"
        )
      );
    }

    //! @brief modify the mutate tree used
    //! @param MUTATE_TREE MutateTree to be modified
    void ProtocolLoopCoordinateAdd::ModifyMutateTree( MutateTree &MUTATE_TREE) const
    {
      // reset all the probabilities
      MUTATE_TREE.Reset();
      MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_SSE     , 3.00);
      MUTATE_TREE.SetMutateProbability( e_MutateCoilBendNTermPhiPsi45, 0.33);
      MUTATE_TREE.SetMutateProbability( e_MutateCoilBendNTermPhi45   , 0.33);
      MUTATE_TREE.SetMutateProbability( e_MutateCoilBendNTermPsi45   , 0.33);
      MUTATE_TREE.SetMutateProbability( e_MutateCoilBendCTermPhiPsi45, 0.33);
      MUTATE_TREE.SetMutateProbability( e_MutateCoilBendCTermPhi45   , 0.33);
      MUTATE_TREE.SetMutateProbability( e_MutateCoilBendCTermPsi45   , 0.33);
      MUTATE_TREE.SetMutateProbability( e_MutateCoilBendNTermPhiPsi5 , 0.66);
      MUTATE_TREE.SetMutateProbability( e_MutateCoilBendCTermPhiPsi5 , 0.66);

      MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_Add         , 1.00);
      MUTATE_TREE.SetMutateProbability( e_MutateCoilResizeNTerm       , 1.00);
      MUTATE_TREE.SetMutateProbability( e_MutateCoilResizeCTerm       , 1.00);
      MUTATE_TREE.SetMutateProbability( e_MutateSSESeedCTerm          , 0.40);
      MUTATE_TREE.SetMutateProbability( e_MutateSSESeedNTerm          , 0.40);
      MUTATE_TREE.SetMutateProbability( e_MutateSSESeedCutCTerm       , 0.20);
      MUTATE_TREE.SetMutateProbability( e_MutateSSESeedCutNTerm       , 0.20);
//      MUTATE_TREE.SetMutateProbability( e_MutateHelixStrandResizeNTerm, 0.20);
//      MUTATE_TREE.SetMutateProbability( e_MutateHelixStrandResizeCTerm, 0.20);

      MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_Model   , 0.25);
      MUTATE_TREE.SetMutateProbability( e_MutateMoveAACTerm       , 0.125);
      MUTATE_TREE.SetMutateProbability( e_MutateMoveAANTerm       , 0.125);
    }

    //! @brief get the mutate tree associated with this protocol
    //! @return the mutate tree associated with this protocol
    util::ShPtr< MutateTree> ProtocolLoopCoordinateAdd::GetMutateTree() const
    {
      util::ShPtr< MutateTree> sp_mutate_tree( new MutateTree());
      ModifyMutateTree( *sp_mutate_tree);
      return sp_mutate_tree;
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void ProtocolLoopCoordinateAdd::MergeAndModifyMutateTree( MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }

  ////////////
  // readme //
  ////////////

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ProtocolLoopCoordinateAdd::GetDescription() const
    {
      static const std::string s_description( "Protocol Loop Coordinate Add");
      return s_description;
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ProtocolLoopCoordinateAdd::GetReadMe() const
    {
      static const std::string s_readme
      (
        "The loop coordinate add protocol adds in coordinates for missing residues in a protein sequence. These"
        " residues can be specified or, if not specified, all missing residues will be added. For each stretch of"
        " missing residues, the coordinates are dynamically added and conformations of the sequence are sampled as"
        " additional residues are added to complete each stretch."
      );

      return s_readme;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProtocolLoopCoordinateAdd::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProtocolLoopCoordinateAdd::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return the stream
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
