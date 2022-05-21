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
#include "assemble/bcl_assemble_collector_protein_model_conformation_by_score.h"
#include "assemble/bcl_assemble_printer_protein_model_ensemble.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "fold/bcl_fold_protocol_ensemble.h"
#include "fold/bcl_fold_setup.h"
#include "mc/bcl_mc_printer_combined.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProtocolEnsemble::s_Instance
    (
      GetObjectInstances().AddInstance( new ProtocolEnsemble())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProtocolEnsemble::ProtocolEnsemble()
    {
    }

    //! @brief Clone function
    //! @return pointer to new ProtocolEnsemble
    ProtocolEnsemble *ProtocolEnsemble::Clone() const
    {
      return new ProtocolEnsemble( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProtocolEnsemble::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ProtocolEnsemble::GetAlias() const
    {
      static const std::string s_name( "ProtocolEnsemble");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProtocolEnsemble::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "protocol to be applied to fold ensemble of models");

      return serializer;
    }

//    util::ShPtr< math::FunctionInterfaceSerializable< assemble::ProteinModel, double> > &ProtocolEnsemble::GetScore()
//    {
//      static util::ShPtr< math::FunctionInterfaceSerializable< assemble::ProteinModel, double> > s_score( new score::ProteinModelScoreSum());
//
//      return s_score;
//    }

//    util::ShPtr< math::MutateInterface< assemble::ProteinModel> >
//    ProtocolEnsemble::GetMutateProteinModelSwitchConformation
//    (
//      const util::ShPtr< math::FunctionInterfaceSerializable< assemble::ProteinModel, double> > &SCORE
//    )
//    {
//      find::Locator
//      <
//        util::SiPtr< const assemble::ProteinModel>,
//        assemble::ProteinModel,
//        util::SiPtrList< const assemble::ProteinModel>
//      >
//      worst_score_locator
//      (
//        util::ShPtr< assemble::CollectorProteinModelConformationByScore> //< collector
//        (
//          new assemble::CollectorProteinModelConformationByScore
//          (
//            SCORE,
//            1, //< collect single worst
//            util::ShPtr< util::BinaryFunctionInterface< double, double, bool> >
//            (
//              new util::BinaryFunctionSTLWrapper< std::greater< double> >() //< sort by worst score
//            ),
//            true //< do consider current conformation
//          )
//        ),
//        util::ShPtr< assemble::PickProteinModelConformationRandom>( new assemble::PickProteinModelConformationRandom())//< picker
//      );
//
//      util::ShPtr< math::MutateInterface< assemble::ProteinModel> > mutate
//      (
//        new MutateProteinModelSwitchConformation
//        (
//          util::ShPtr
//          <
//            find::LocatorInterface< util::SiPtr< const assemble::ProteinModel>, assemble::ProteinModel>
//          >
//          (
//            worst_score_locator.Clone()
//          ),
//          "ensemble_swtch_conf_wrst_scr"
//        )
//      );
//
//      return mutate;
//    }
//
//    util::ShPtr< math::MutateInterface< assemble::ProteinModel> >
//    ProtocolEnsemble::GetMutateCombineProteinModelSwitchConformation
//    (
//      const util::ShPtr< math::FunctionInterfaceSerializable< assemble::ProteinModel, double> > &SCORE,
//      const util::ShPtr< math::MutateInterface< assemble::ProteinModel> > &MUTATE
//    )
//    {
//      util::ShPtrList< math::MutateInterface< assemble::ProteinModel> > mutates;
//      mutates.PushBack( GetMutateProteinModelSwitchConformation( SCORE));
//      mutates.PushBack( MUTATE);
//
//      util::ShPtr< math::MutateInterface< assemble::ProteinModel> > combine
//      (
//        new math::MutateCombine< assemble::ProteinModel>( mutates, true, "ensemble_combine")
//      );
//
//      return combine;
//    }

    //! @brief access to only instance
    //! @return reference to only instance
    ProtocolEnsemble &ProtocolEnsemble::GetInstance()
    {
      static ProtocolEnsemble s_protocol_instance;
      return s_protocol_instance;
    }

    //! @brief return command line flag for specifying the domains
    //! @return command line flag for specifying domains
    util::ShPtr< command::FlagInterface> &ProtocolEnsemble::GetFlagEnsembleSize()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "ensemble_size",
          "\tThe size of the ensemble that should be used at the start of folding. "
          "i.e. how many models should be in the ensemble to start with.",
          command::Parameter( "number", "\tstart size of ensemble", "1")
        )
      );

      return s_flag;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns all flags that are specialized for this protocol
    //! @return all flags that are specialized for this protocol
    const util::ShPtrVector< command::FlagInterface> &ProtocolEnsemble::GetAllFlags() const
    {
      // initialize static ShPtrVector of FlagInterfaces to form the comment line
      static util::ShPtrVector< command::FlagInterface> s_all_flags_vector;

      // if the flag vector is initialize for the first time
      if( s_all_flags_vector.IsEmpty())
      {
        // insert all the flags in the vector
        s_all_flags_vector.PushBack( GetFlagEnsembleSize());
      }

      // end
      return s_all_flags_vector;
    }

    //! @brief modifies the start model
    //! @param START_MODEL Start model to be modified
    void ProtocolEnsemble::ModifyStartModel( assemble::ProteinModel &START_MODEL) const
    {
      const size_t start_ensemble_size( GetFlagEnsembleSize()->GetFirstParameter()->GetNumericalValue< size_t>() - 1);
      if( !START_MODEL.GetConformationalEnsemble().IsEmpty())
      {
        BCL_MessageCrt
        (
          "conformation ensemble is not empty so not modifying start model\nSize is "
          + util::Format()( START_MODEL.GetConformationalEnsemble().GetSize())
        );
        return;
      }

      util::ShPtr< assemble::ProteinModelData> pmd_here( START_MODEL.GetProteinModelData());

      BCL_Assert( pmd_here->GetData( assemble::ProteinModelData::e_Pool).IsDefined(), "pool is not defined");

      BCL_MessageCrt( "model type is " + START_MODEL.GetClassIdentifier());

      // ensemble for folding
      assemble::ProteinEnsemble ensemble;

      // insert as the start model as many times as desired into the ensemble
      for( size_t num( 0); num < start_ensemble_size; ++num)
      {
        util::ShPtr< assemble::ProteinModel> new_model( START_MODEL.HardCopy());

        util::ShPtr< assemble::ProteinModelData> pmd( new_model->GetProteinModelData());
        BCL_Assert( pmd_here->GetData( assemble::ProteinModelData::e_Pool).IsDefined(), "START_MODEL pool is not defined");
        BCL_Assert
        (
          pmd->GetData( assemble::ProteinModelData::e_Pool).IsDefined(),
          "pool is not defined during multiplying creation"
        );

        ensemble.InsertElement( new_model);
      }

      START_MODEL.SetConformationalEnsemble( ensemble);

      for
      (
        assemble::ProteinEnsemble::const_iterator
          itr( START_MODEL.GetConformationalEnsemble().Begin()), itr_end( START_MODEL.GetConformationalEnsemble().End());
        itr != itr_end;
        ++itr
      )
      {
        util::ShPtr< assemble::ProteinModelData> pmd( ( *itr)->GetProteinModelData());

        BCL_Assert( pmd->GetData( assemble::ProteinModelData::e_Pool).IsDefined(), "pool is not defined after ensemble creation");
      }

      BCL_MessageCrt( "ensemble size is " + util::Format()( START_MODEL.GetConformationalEnsemble().End() - START_MODEL.GetConformationalEnsemble().Begin()));
    }

    //! @brief initialize the scores and add them to Scores enumerator
    void ProtocolEnsemble::InitializeScores()
    {
    }

    //! @brief modify the score weight set
    //! @param SCORE_WEIGHT_SET Score weight set
    void ProtocolEnsemble::ModifyScoreWeightSet( ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
    }

    //! @brief modify the terminate object
    //! @param CRITERION which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolEnsemble::ModifyCriterion
    (
      opti::CriterionCombine< assemble::ProteinModel, double> &CRITERION,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief modify the printer object
    //! @param PRINTER which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolEnsemble::ModifyPrinter
    (
      mc::PrinterCombined< assemble::ProteinModel, double> &PRINTER,
      const mc::Stage &STAGE
    ) const
    {
      util::ShPtr< mc::PrintInterface< assemble::ProteinModel, double> > ensemble_printer
      (
        new assemble::PrinterProteinModelEnsemble
        (
          GetSetup().GetPrefix(),
          GetSetup().GetStorage(),
          GetSetup().GetSuperimposeMeasure()
        )
      );
      PRINTER.Insert( ensemble_printer);
    }

    //! @brief modify the pdb factory object
    //! @param FACTORY pdb factory to be modified
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolEnsemble::ModifyFactory
    (
      util::ShPtr< pdb::Factory> &FACTORY,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief initialize the mutates and add them to Mutates enumerator
    void ProtocolEnsemble::InitializeMutates()
    {

//
//      e_MutateEnsembleCombineSwitchConformationWorstScore = GetMutates().AddMutate
//      (
//        GetMutateCombineProteinModelSwitchConformation
//        (
//          util::ShPtr< math::FunctionInterfaceSerializable< assemble::ProteinModel, double> >( new score::ProteinModelScoreSum()),
//          util::ShPtr< math::MutateInterface< assemble::ProteinModel> >( new math::MutateDecisionNode< assemble::ProteinModel>())
//        )
//      );
    }

    //! @brief modify the mutate tree used
    //! @param MUTATE_TREE MutateTree to be modified
    void ProtocolEnsemble::ModifyMutateTree( MutateTree &MUTATE_TREE) const
    {
//      const storage::Vector< double> mutate_type_probs( MUTATE_TREE.GetMutateTypeProbabilities().GetMappedValues());
//
//      const double total_probabilities( math::Statistics::Sum( mutate_type_probs.Begin(), mutate_type_probs.End()));
//
//      MUTATE_TREE.SetMutateTypeProbability( MutateTree::e_Ensemble, total_probabilities);
//
//      MUTATE_TREE.SetMutateProbability( e_MutateEnsembleSwitchConformationWorstScore, 1);
    }

    //! @brief get the mutate tree associated with this protocol
    //! @return the mutate tree associated with this protocol
    util::ShPtr< MutateTree> ProtocolEnsemble::GetMutateTree() const
    {
      util::ShPtr< MutateTree> sp_mutate_tree( new MutateTree());
      ModifyMutateTree( *sp_mutate_tree);
      return sp_mutate_tree;
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void ProtocolEnsemble::MergeAndModifyMutateTree( MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }

  ////////////
  // readme //
  ////////////

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ProtocolEnsemble::GetDescription() const
    {
      // initialize string to store the description
      static const std::string s_description
      (
        "protocol to be applied to fold ensemble of models"
      );

      // end
      return s_description;
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ProtocolEnsemble::GetReadMe() const
    {
      // TODO: add readme for this protocol
      // initialize string to store the readme information
      static const std::string s_readme
      (
        "readme for ensemble folding"
      );

      // end
      return s_readme;
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProtocolEnsemble::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProtocolEnsemble::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////
  
  } // namespace fold
} // namespace bcl
