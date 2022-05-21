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
#include "fold/bcl_fold_protocol_loop_close.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_back_bone_completer.h"
#include "biol/bcl_biol_aa_complete.h"
#include "command/bcl_command_flag_static.h"
#include "fold/bcl_fold_collector_loop_domain.h"
#include "fold/bcl_fold_collector_loop_domain_random.h"
#include "fold/bcl_fold_criterion_loop_closure.h"
#include "fold/bcl_fold_handler_locator_loop_domain.h"
#include "fold/bcl_fold_loop_domain.h"
#include "fold/bcl_fold_mutate_loop_domain_dihedral.h"
#include "fold/bcl_fold_mutate_protein_model_grow_sse.h"
#include "fold/bcl_fold_mutate_protein_model_loop_domain.h"
#include "fold/bcl_fold_mutate_protein_model_loop_domain_ccd.h"
#include "fold/bcl_fold_phi_psi_generator_ramachandran.h"
#include "fold/bcl_fold_protocol_loop_coordinate_add.h"
#include "io/bcl_io_file.h"
#include "mc/bcl_mc_stage.h"
#include "pdb/bcl_pdb_printer_loop_closure.h"
#include "score/bcl_score_phi_psi_with_sspred.h"
#include "score/bcl_score_protein_model_loop_domain_closure.h"
#include "score/bcl_score_protein_model_sse.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProtocolLoopClose::ProtocolLoopClose()
    {
    }

    //! @brief Clone function
    //! @return pointer to new ProtocolLoopClose
    ProtocolLoopClose *ProtocolLoopClose::Clone() const
    {
      return new ProtocolLoopClose( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief access to only instance
    //! @return reference to only instance
    ProtocolLoopClose &ProtocolLoopClose::GetInstance()
    {
      static ProtocolLoopClose s_protocol_instance;
      return s_protocol_instance;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProtocolLoopClose::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ProtocolLoopClose::GetAlias() const
    {
      static const std::string s_name( "ProtocolLoopClose");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProtocolLoopClose::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "protocol for closing loops");
      serializer.AddInitializer
        (
         "mutate loop domain CCD",
         "ccd loop domain move to bring end closer to the other side",
         io::Serialization::GetAgent( &e_MutateLoopDomainCCD)
         );
      serializer.AddInitializer
        (
         "mutate loop domain ramachandran",
         "loop domain phi psi change base on ramachandan",
         io::Serialization::GetAgent( &e_MutateLoopDomainRamachandran)
         );
      serializer.AddInitializer
        (
         "score loop domain closure",
         "score associated with loop domain closure",
         io::Serialization::GetAgent( &e_ScoreLoopDomainClosure)
         );
      serializer.AddInitializer
        (
         "score phi psi SSPred close",
         "score associated with the phi psi SS prediction",
         io::Serialization::GetAgent( &e_ScorePhiPsiSSPredClose)
         );

      return serializer;
    }

  ///////////
  // flags //
  ///////////

    //! @brief returns all flags that are specialized for this protocol
    //! @return all flags that are specialized for this protocol
    const util::ShPtrVector< command::FlagInterface> &ProtocolLoopClose::GetAllFlags() const
    {
      // initialize static ShPtrVector of FlagInterfaces to form the comment line
      static util::ShPtrVector< command::FlagInterface> s_all_flags_vector;

      // if the flag vector is initialize for the first time
      if( s_all_flags_vector.IsEmpty())
      {
        // insert all the flags in the vector
        s_all_flags_vector.PushBack( LoopDomain::GetFlagLoopDomainFilename());
        s_all_flags_vector.PushBack( GetFlagRamachandranMutateProbability());
        s_all_flags_vector.PushBack( GetFlagUseFullFraction());
      }

      // end
      return s_all_flags_vector;
    }

    //! @brief return command line flag for specifying the probability of the mutate which selects a phi or psi change
    //!        from the ramachandran probability map (as opposed to the ccd mutate)
    //! @return command line flag for specifying the probability of the ramachandran phi psi mutate
    util::ShPtr< command::FlagInterface> &ProtocolLoopClose::GetFlagRamachandranMutateProbability()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "loop_rama_mutate_prob",
          "\tthe probability of the ramachandran phi psi mutate occurance"
        )
      );

      // initialize parameters
      static util::ShPtr< command::ParameterInterface> s_probability
      (
        new command::Parameter
        (
          "probability", "\tthe probability of the ramachandran phi psi mutate occuring", "0.1"
        )
      );

      // if the flag is initialized for the first time
      if( s_flag->GetParameterList().IsEmpty())
      {
        util::ShPtr< command::FlagStatic> flag( s_flag);
        // insert parameters
        flag->PushBack( s_probability);
      }

      // end
      return s_flag;
    }

    //! @brief return command line flag for specifying if random fraction is to be used
    //! @return command line flag for specifying if random fraction is to be used
    util::ShPtr< command::FlagInterface> &ProtocolLoopClose::GetFlagUseFullFraction()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "ccd_fraction",
          "\tuse a fraction of the optimal rotation in the cyclic coordinate descent",
          command::Parameter
          (
            "fraction",
            "multiply the angle with a random number in the given range",
            math::Range< double>( 1.0, 1.0).GetString()
          )
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
    void ProtocolLoopClose::ModifyStartModel( assemble::ProteinModel &START_MODEL) const
    {
      // complete the model with coordinates in missing places
      AddLoopCoordinates( START_MODEL);

      // split the coils, where they are not peptide bonded
      SplitCoilsAtNonPetideBond( START_MODEL);

      // add backbone hydrogens
      START_MODEL = *AddNitrogenHydrogens( START_MODEL);

      // initialize empty locator
      util::ShPtr< util::ShPtrList< LocatorLoopDomain> > sp_locator( new util::ShPtrList< LocatorLoopDomain>());

      // if a loop file was provided
      if( LoopDomain::GetFlagLoopDomainFilename()->GetFlag())
      {
        // read the given file
        io::IFStream read;
        io::File::MustOpenIFStream( read, LoopDomain::GetFlagLoopDomainFilename()->GetFirstParameter()->GetValue());

        // update sp_locator with the locators read from file
        *sp_locator = HandlerLocatorLoopDomain().HandleReadMultiple( read, START_MODEL);
      }
      else
      {
        // otherwise create locators using the SSE definitions from the model
        sp_locator = HandlerLocatorLoopDomain().CreateBidirectionalLocatorsForInteriorCoil( START_MODEL);
      }

      // print out the locator information
      for
      (
        util::ShPtrList< LocatorLoopDomain>::const_iterator
          domain_itr( sp_locator->Begin()), domain_itr_end( sp_locator->End());
        domain_itr != domain_itr_end;
        ++domain_itr
      )
      {
        BCL_MessageDbg( ( *domain_itr)->GetIdentification());
      }

      // get the protein model data
      util::ShPtr< assemble::ProteinModelData> sp_model_data( START_MODEL.GetProteinModelData());
      BCL_Assert( sp_model_data.IsDefined(), "sp_model_data is not defined");

      // try to see if the start model already has it
      util::ShPtr< util::ShPtrList< LocatorLoopDomain> > sp_loop_data
      (
        sp_model_data->GetData( assemble::ProteinModelData::e_LoopDomainLocators)
      );

      // if not defined insert it
      if( !sp_loop_data.IsDefined())
      {
        sp_model_data->Insert( assemble::ProteinModelData::e_LoopDomainLocators, sp_locator);
      }
      // otherwise replace it
      else
      {
        sp_model_data->Replace( assemble::ProteinModelData::e_LoopDomainLocators, sp_locator);
      }
    }

    //! @brief initialize the scores and add them to Scores enumerator
    void ProtocolLoopClose::InitializeScores()
    {
      // score to help align the atoms desired to be superimposed in order for the loop to be closed
      e_ScoreLoopDomainClosure = GetScores().AddScore
      (
        Scores::WrapCacheProteinModelScore
        (
          util::ShPtr< score::ProteinModel>
          (
            new score::ProteinModelLoopDomainClosure()
          )
        )
      );

      // phi psi for loop building with ss prediction
      e_ScorePhiPsiSSPredClose = GetScores().AddScore
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
                  new score::PhiPsiWithSSPred
                  (
                    sspred::Methods::GetCommandLineMethods(),
                    "phipsi_sspred_close"
                  )
                )
              ),
              false
            )
          )
        )
      );
      ProtocolLoopCoordinateAdd::GetInstance().InitializeScores();
    }

    //! @brief modify the score weight set
    //! @param SCORE_WEIGHT_SET Score weight set
    void ProtocolLoopClose::ModifyScoreWeightSet( ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
      SCORE_WEIGHT_SET.Reset();

//      SCORE_WEIGHT_SET.SetWeight( ProtocolLoopCoordinateAdd::GetInstance().e_ScoreSSEPairConnectivity, 200.0);
      SCORE_WEIGHT_SET.SetWeight( ProtocolLoopCoordinateAdd::GetInstance().e_ScorePhiPsiCoil         , 100.0);
      SCORE_WEIGHT_SET.SetWeight( e_ScoreLoopDomainClosure,                         5.0);

      if( sspred::Methods::GetFlagReadSSPredictions()->GetFlag())
      {
        SCORE_WEIGHT_SET.SetWeight( e_ScorePhiPsiSSPredClose, 5.0);
      }
    }

    //! @brief modify the terminate object
    //! @param CRITERION which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolLoopClose::ModifyCriterion
    (
      opti::CriterionCombine< assemble::ProteinModel, double> &CRITERION,
      const mc::Stage &STAGE
    ) const
    {
      // construct the loop closure criterion and insert it into the terminate combine
      CRITERION.InsertCriteria
      (
        CriterionLoopClosure< assemble::ProteinModel, double>
        (
          0.08
        )
      );
    }

    //! @brief modify the printer object
    //! @param PRINTER which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolLoopClose::ModifyPrinter
    (
      mc::PrinterCombined< assemble::ProteinModel, double> &PRINTER,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief modify the pdb factory object
    //! @param FACTORY pdb factory to be modified
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolLoopClose::ModifyFactory
    (
      util::ShPtr< pdb::Factory> &FACTORY,
      const mc::Stage &STAGE
    ) const
    {
      FACTORY->AppendPrinter
      (
        util::ShPtr< util::FunctionInterface< assemble::ProteinModel, util::ShPtrList< pdb::Line> > >
        (
          new pdb::PrinterLoopClosure
          (
            0.08
          )
        )
      );
    }

    //! @brief initialize the mutates and add them to Mutates enumerator
    void ProtocolLoopClose::InitializeMutates()
    {
      // was initialized?
      if( e_MutateLoopDomainCCD.IsDefined())
      {
        return;
      }
      // construct collector for unclosed loops
      util::ShPtr< find::CollectorInterface< util::ShPtrList< LoopDomain>, assemble::DomainInterface> >
      sp_unclosed_loop_collector
      (
        new CollectorLoopDomain
        (
          true,
          0.08
        )
      );

      math::Range< double> ccd_fraction;
      std::stringstream range_stream( GetFlagUseFullFraction()->GetFirstParameter()->GetValue());
      std::stringstream error_stream;
      BCL_Assert
      (
        ccd_fraction.FromStream( range_stream, error_stream),
        "Error while processing range from ccd flag: " + GetFlagUseFullFraction()->GetFirstParameter()->GetValue() +
        " " + error_stream.str()
      );

      // add Cyclic descent mutate
      e_MutateLoopDomainCCD = GetMutates().AddMutate
      (
        util::ShPtr< math::MutateInterface< assemble::ProteinModel> >
        (
          new MutateProteinModelLoopDomainCCD
          (
            sp_unclosed_loop_collector,
            random::GetGlobalRandom(),
            ccd_fraction
          )
        )
      );

      // add ramachandran mutate
      e_MutateLoopDomainRamachandran = GetMutates().AddMutate
      (
        util::ShPtr< math::MutateInterface< assemble::ProteinModel> >
        (
          new MutateProteinModelLoopDomain
          (
            sp_unclosed_loop_collector,
            util::ShPtr< math::MutateInterface< LoopDomain> >
            (
              new MutateLoopDomainDihedral
              (
                util::ShPtr< find::CollectorInterface< storage::List< MutationResidue>, LoopDomain> >
                (
                  new CollectorLoopDomainRandom( 1)
                ),
                util::ShPtr< math::FunctionInterfaceSerializable< MutationResidue, storage::VectorND< 2, double> > >
                (
                  PhiPsiGeneratorRamachandran::GetDefaultInstance().Clone()
                )
              )
            )
          )
        )
      );

      ProtocolLoopCoordinateAdd::GetInstance().InitializeMutates();
    }

    //! @brief modify the mutate tree used
    //! @param MUTATE_TREE MutateTree to be modified
    void ProtocolLoopClose::ModifyMutateTree( MutateTree &MUTATE_TREE) const
    {
      // reset all the probabilities
      MUTATE_TREE.Reset();

      MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_SSE, 1.0);
      MUTATE_TREE.SetMutateProbability( e_MutateLoopDomainCCD, 1.0);
      MUTATE_TREE.SetMutateProbability
      (
        e_MutateLoopDomainRamachandran,
        GetFlagRamachandranMutateProbability()->GetFirstParameter()->GetNumericalValue< double>()
      );

      MUTATE_TREE.SetMutateProbability
      (
        ProtocolLoopCoordinateAdd::GetInstance().e_MutateCoilBendCTermPhiPsi5,
        0.25
      );
      MUTATE_TREE.SetMutateProbability
      (
        ProtocolLoopCoordinateAdd::GetInstance().e_MutateCoilBendNTermPhiPsi5,
        0.25
      );
    }

    //! @brief get the mutate tree associated with this protocol
    //! @return the mutate tree associated with this protocol
    util::ShPtr< MutateTree> ProtocolLoopClose::GetMutateTree() const
    {
      util::ShPtr< MutateTree> sp_mutate_tree( new MutateTree());
      ModifyMutateTree( *sp_mutate_tree);
      return sp_mutate_tree;
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void ProtocolLoopClose::MergeAndModifyMutateTree( MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }

  ////////////
  // readme //
  ////////////

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ProtocolLoopClose::GetDescription() const
    {
      static const std::string s_description( "Protocol Loop Close");
      return s_description;
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ProtocolLoopClose::GetReadMe() const
    {
      static const std::string s_readme
      (
        "The loop close protocol is for tasks related to loops and closing loops using CCD. It can take a protein model "
        "with missing density and add the missing sequence and close the chain. It can do this for a single loop. It "
        "can also take a sequentially consecutive domain of sses and connect this loop domain. The purpose of this "
        "would be to explore other conformations that a rigid domain can adopt relative to flexible linking loops."
      );

      return s_readme;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProtocolLoopClose::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProtocolLoopClose::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief adds coordinates for residues that are in the protein model chain sequences, but are not within SSEs in
    //!        the protein model
    //! @param MODEL the protein model to be changed
    void ProtocolLoopClose::AddLoopCoordinates( assemble::ProteinModel &MODEL)
    {
      // create a phi psi generator
      static const util::ShPtr< math::FunctionInterfaceSerializable< MutationResidue, storage::VectorND< 2, double> > >
      phi_psi_generator
      (
        PhiPsiGeneratorRamachandran::GetDefaultInstance().Clone()
      );

      // grow all the missing coils to some initial conformation
      MODEL = *MutateProteinModelGrowSSE::GrowAllMissingCoilsBidirectional( MODEL, phi_psi_generator);

      // test for peptide bonding
      MODEL.Join( biol::GetSSTypes().COIL, true);
    }

    //! @brief adds coordinates for nitrogen hydrogens to residues by converting them to AAComplete and calcualting ideal N-H positions
    //! @param MODEL the protein model
    util::ShPtr< assemble::ProteinModel> ProtocolLoopClose::AddNitrogenHydrogens( const assemble::ProteinModel &MODEL)
    {
      // create a new model
      util::ShPtr< assemble::ProteinModel> new_model( MODEL.Empty());

      // set the data
      new_model->SetProteinModelData( MODEL.GetProteinModelData());

      // iterate over all chains
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
          chain_itr( MODEL.GetChains().Begin()), chain_itr_end( MODEL.GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        util::ShPtr< assemble::Chain> sp_chain( new assemble::Chain( ( *chain_itr)->GetSequence()));

        // pointer to last aa of previous sse
        util::SiPtr< const biol::AABase> last_aa_previous_sse;

        // iterate over SSEs
        for( assemble::Chain::const_ierator sse_itr( ( *chain_itr)->GetData().Begin()), sse_itr_end( ( *chain_itr)->GetData().End()); sse_itr != sse_itr_end; ++sse_itr)
        {
          biol::AASequence sequence( biol::GetAAClasses().e_AAComplete, 0, sp_chain->GetChainID());
          // iterate over sse
          biol::AASequence::const_iterator aa_itr( ( *sse_itr)->Begin()), aa_itr_end( ( *sse_itr)->End());

          // construct hydrogen for the first aa
          {
            util::ShPtr< biol::AAComplete> sp_aa_first( new biol::AAComplete( **aa_itr));

            // construct first hydrogen
            if( last_aa_previous_sse.IsDefined() && last_aa_previous_sse->DoesPrecede( *sp_aa_first) && biol::AABase::AreAminoAcidsPeptideBonded( *last_aa_previous_sse, *sp_aa_first, true))
            {
              const biol::Atom hydrogen( biol::AABackBoneCompleter::GenerateHydrogen( *sp_aa_first, last_aa_previous_sse->GetAtom( biol::GetAtomTypes().C).GetCoordinates()));
              sp_aa_first->SetAtom( hydrogen);
            }
            else
            {
              double phi( ( *sse_itr)->GetType()->GetIdealPhi());
              if( !util::IsDefined( phi))
              {
                phi = biol::GetSSTypes().HELIX->GetIdealPhi();
              }
              phi += math::g_Pi;

              const biol::Atom hydrogen
              (
                linal::CoordinatesDihedral
                (
                  sp_aa_first->GetAtom( biol::GetAtomTypes().N).GetCoordinates(),
                  sp_aa_first->GetAtom( biol::GetAtomTypes().CA).GetCoordinates(),
                  sp_aa_first->GetAtom( biol::GetAtomTypes().C).GetCoordinates(),
                  biol::GetAtomTypes().N->GetBondLength( biol::GetAtomTypes().H),
                  2 * math::g_Pi / 3,
                  phi
                ),
                biol::GetAtomTypes().H
              );
              sp_aa_first->SetAtom( hydrogen);
            }

            // insert the first amino acid
            sequence.PushBack( sp_aa_first);
          }

          // loop over all amino acids
          for( biol::AASequence::const_iterator aa_itr_next( aa_itr + 1); aa_itr_next != aa_itr_end; ++aa_itr, ++aa_itr_next)
          {
            util::ShPtr< biol::AAComplete> sp_aa( new biol::AAComplete( **aa_itr_next));

            // construct first hydrogen
            const biol::Atom hydrogen( biol::AABackBoneCompleter::GenerateHydrogen( *sp_aa, ( *aa_itr)->GetAtom( biol::GetAtomTypes().C).GetCoordinates()));
            sp_aa->SetAtom( hydrogen);
            sequence.PushBack( sp_aa);
          }
          last_aa_previous_sse = util::SiPtr< const biol::AABase>( *( *sse_itr)->GetLastAA());
          util::ShPtr< assemble::SSE> sp_sse( new assemble::SSE( sequence, ( *sse_itr)->GetType()));
          sp_chain->Insert( sp_sse);
        }

        new_model->Insert( sp_chain);
      }

      // end
      return new_model;
    }

    //! @brief split coils at non-petide connection
    //! @param MODEL the protein model to be changed
    void ProtocolLoopClose::SplitCoilsAtNonPetideBond( assemble::ProteinModel &MODEL)
    {
      // get all coils
      const util::SiPtrVector< const assemble::SSE> coils( MODEL.GetSSEs( biol::GetSSTypes().COIL));

      // iterate over all coils
      for( util::SiPtrVector< const assemble::SSE>::const_iterator itr( coils.Begin()), itr_end( coils.End()); itr != itr_end; ++itr)
      {
        const assemble::SSE &current( **itr);
        if( current.GetSize() < 2)
        {
          continue;
        }

        const util::ShPtrList< assemble::SSE> sses( SplitSSEAtNonPetideBond( current));
        MODEL.Remove( current);
        for( util::ShPtrList< assemble::SSE>::const_iterator coil_itr( sses.Begin()), coil_itr_end( sses.End()); coil_itr != coil_itr_end; ++coil_itr)
        {
          MODEL.Insert( *coil_itr);
        }
      }
    }

    //! @brief split sse into coils that are not connected
    //! @param SSE the sse to split
    //! @return list of coils
    util::ShPtrList< assemble::SSE> ProtocolLoopClose::SplitSSEAtNonPetideBond( const assemble::SSE &SS_ELEMENT)
    {
      util::ShPtrList< assemble::SSE> sses;

      assemble::SSE::const_iterator aa_itr( SS_ELEMENT.Begin()), aa_itr_end( SS_ELEMENT.End());

      while( aa_itr != aa_itr_end)
      {
        biol::AASequence new_seq( util::ShPtrVector< biol::AABase>(), SS_ELEMENT.GetChainID());

        // amino acids with undefined coordinates into new sse
        while( aa_itr != aa_itr_end && !( *aa_itr)->GetFirstSidechainAtom().GetCoordinates().IsDefined())
        {
          new_seq.PushBack( *aa_itr);
          ++aa_itr;
        }

        // found aa with undefined coordinates
        if( !new_seq.GetData().IsEmpty())
        {
          sses.PushBack( util::ShPtr< assemble::SSE>( new assemble::SSE( new_seq, SS_ELEMENT.GetType())));
          continue;
        }

        // peptide bond between amino acids
        new_seq.PushBack( *aa_itr);
        assemble::SSE::const_iterator aa_itr_next( aa_itr + 1);
        while( aa_itr_next != aa_itr_end && biol::AABase::AreAminoAcidsPeptideBonded( **aa_itr, **aa_itr_next, true))
        {
          ++aa_itr;
          ++aa_itr_next;
          new_seq.PushBack( *aa_itr);
        }

        // store coil
        sses.PushBack( util::ShPtr< assemble::SSE>( new assemble::SSE( new_seq, biol::GetSSTypes().COIL)));
        ++aa_itr;
      }

      // end
      return sses;
    }

  } // namespace fold
} // namespace bcl
