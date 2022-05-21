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
#include "fold/bcl_fold_protocol_ensemble_filter.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_collector_protein_model_conformation_by_score.h"
#include "assemble/bcl_assemble_printer_protein_model_ensemble.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_parameter.h"
#include "fold/bcl_fold_score_weight_set.h"
#include "fold/bcl_fold_setup.h"
#include "io/bcl_io_file.h"
#include "mc/bcl_mc_printer_combined.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProtocolEnsembleFilter::s_Instance
    (
      GetObjectInstances().AddInstance( new ProtocolEnsembleFilter())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProtocolEnsembleFilter::ProtocolEnsembleFilter()
    {
    }

    //! @brief Clone function
    //! @return pointer to new ProtocolEnsembleFilter
    ProtocolEnsembleFilter *ProtocolEnsembleFilter::Clone() const
    {
      return new ProtocolEnsembleFilter( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief access to only instance
    //! @return reference to only instance
    ProtocolEnsembleFilter &ProtocolEnsembleFilter::GetInstance()
    {
      static ProtocolEnsembleFilter s_protocol_instance;
      return s_protocol_instance;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProtocolEnsembleFilter::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ProtocolEnsembleFilter::GetAlias() const
    {
      static const std::string s_name( "ProtocolEnsembleFilter");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProtocolEnsembleFilter::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Protocol to be applied for filtering an ensemble.");

      return serializer;
    }

  ///////////
  // flags //
  ///////////

    //! @brief returns all flags that are specialized for this protocol
    //! @return all flags that are specialized for this protocol
    const util::ShPtrVector< command::FlagInterface> &ProtocolEnsembleFilter::GetAllFlags() const
    {
      // initialize static ShPtrVector of FlagInterfaces to form the comment line
      static util::ShPtrVector< command::FlagInterface> s_all_flags_vector;

      // if the flag vector is initialize for the first time
      if( s_all_flags_vector.IsEmpty())
      {
        // insert all the flags in the vector
        s_all_flags_vector.PushBack( GetFlagFilterScoreFunctions());
        s_all_flags_vector.PushBack( GetFlagFilterNumberToKeep());
      }

      // end
      return s_all_flags_vector;
    }

    //! @brief return command line flag for specifying the domains
    //! @return command line flag for specifying domains
    util::ShPtr< command::FlagInterface> &ProtocolEnsembleFilter::GetFlagFilterScoreFunctions()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagDynamic
        (
          "ensemble_filter_score_weights",
          "\tFilenames of score weight table formated files that will be used to construct score functions. The score"
          " functions will be used to filter the conformational ensemble. Multiple filenames can be provided. This"
          " flag should be used in conjunction with the \"ensemble_filter_num_to_keep\" flag. One weight file per"
          " number provided to that flag. According to each score function, the best N models will be kept.",
          command::Parameter( "score_weight_table_filename", "\tscore weight table filenames", "score_weights.tbl")
        )
      );

      return s_flag;
    }

    //! @brief return command line flag for specifying the domains
    //! @return command line flag for specifying domains
    util::ShPtr< command::FlagInterface> &ProtocolEnsembleFilter::GetFlagFilterNumberToKeep()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagDynamic
        (
          "ensemble_filter_number_to_keep",
          "\tThe number of models to keep based on the scoring weight table provided. The score"
          " functions will be used to filter the conformational ensemble. Multiple numbers can be provided. This"
          " flag should be used in conjunction with the \"ensemble_filter_score_weights\" flag. One weight file per"
          " number provided to that flag. According to each score function, the best N models will be kept.",
          command::Parameter( "number_to_keep", "\tsize_t", "1")
        )
      );

      return s_flag;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief modifies the start model
    //! @param START_MODEL Start model to be modified
    void ProtocolEnsembleFilter::ModifyStartModel( assemble::ProteinModel &START_MODEL) const
    {
      // get the vector that tells the number of models to keep for each score weight table
      const storage::Vector< int> numbers_to_keep( GetFlagFilterNumberToKeep()->GetNumericalList< int>());

      // get the vector of score weight table files
      const storage::Vector< std::string> score_weight_filenames( GetFlagFilterScoreFunctions()->GetStringList());

      BCL_Assert
      (
        numbers_to_keep.GetSize() == score_weight_filenames.GetSize(),
        "pass the same number of parameters to the flag "
        + GetFlagFilterNumberToKeep()->GetName() + " as passed to the flag "
        + GetFlagFilterScoreFunctions()->GetName()
        + "\nThe number of parameters passed were " + util::Format()( numbers_to_keep.GetSize()) + " and "
        + util::Format()( score_weight_filenames.GetSize()) + ", respectively"
      );

      // to hold all conformations that fulfill scoring filters
      util::SiPtrList< const assemble::ProteinModel> all_conformations;

      // iterate over the number of models to keep and the corresponding score weight table files
      for
      (
        storage::Vector< std::string>::const_iterator
          file_itr( score_weight_filenames.Begin()), file_itr_end( score_weight_filenames.End());
        file_itr != file_itr_end;
        ++file_itr
      )
      {
        // get the number of models to keep for the current score weight table
        const int num_to_keep( *( numbers_to_keep.Begin() + int( file_itr - score_weight_filenames.Begin())));

        BCL_MessageStd
        (
          "keeping " + util::Format()( num_to_keep) + " models for score weight table"
          " file " + *file_itr
        );

        // read in the score weight table
        io::IFStream read;
        io::File::MustOpenIFStream( read, *file_itr);
        storage::Table< double> score_weight_table;
        score_weight_table.ReadFormatted( read);

        // create a score weightset
        const ScoreWeightSet weightset( score_weight_table);

        // create score object to filter the ensemble
        const util::ShPtr< score::ProteinModelScoreSum> score_function( weightset.ConstructScoreSum());

        // create collector to get the models that fulfill the score function criteria
        const assemble::CollectorProteinModelConformationByScore collector
        (
          score_function,
          num_to_keep,
          util::ShPtr< util::BinaryFunctionInterface< double, double, bool> >
          (
            new util::BinaryFunctionSTLWrapper< std::less< double> >() //< collect best
          ),
          true
        );

        // get the conformations that fulfill the score criteria
        const util::SiPtrList< const assemble::ProteinModel> conformations( collector.Collect( START_MODEL));

        // append the conformations to the overall list that are desired
        all_conformations.Append( conformations);
      }

      BCL_MessageStd
      (
        "number conformations passing score filter "
        + util::Format()( all_conformations.GetSize())
      );

      // create new model using first model in the collected conformations
      util::ShPtr< assemble::ProteinModel> new_model( all_conformations.FirstElement()->Clone());

      // remove the first element from the collected conformations
      all_conformations.RemoveElement( all_conformations.Begin());

      // ensemble to hold the rest of the conformations
      assemble::ProteinEnsemble ensemble;

      BCL_MessageStd
      (
        "number conformations to add to ensemble "
        + util::Format()( all_conformations.GetSize())
      );

      // iterate through the rest of the conformations to put them into the ensemble
      for
      (
        util::SiPtrList< const assemble::ProteinModel>::const_iterator
          model_itr( all_conformations.Begin()), model_itr_end( all_conformations.End());
        model_itr != model_itr_end;
        ++model_itr
      )
      {
        // copy the model into a sh ptr
        util::ShPtr< assemble::ProteinModel> current_model( ( *model_itr)->Clone());

        // reset the conformational ensemble of the new model
        current_model->SetConformationalEnsemble( assemble::ProteinEnsemble());

        // insert the new model into the ensemble
        ensemble.InsertElement( current_model);
      }

      // set the conformational ensemble of the new model
      new_model->SetConformationalEnsemble( ensemble);

      START_MODEL = ( *new_model);
    }

    //! @brief initialize the scores and add them to Scores enumerator
    void ProtocolEnsembleFilter::InitializeScores()
    {
    }

    //! @brief modify the score weight set
    //! @param SCORE_WEIGHT_SET Score weight set
    void ProtocolEnsembleFilter::ModifyScoreWeightSet( ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
    }

    //! @brief modify the terminate object
    //! @param CRITERION which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolEnsembleFilter::ModifyCriterion
    (
      opti::CriterionCombine< assemble::ProteinModel, double> &CRITERION,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief modify the printer object
    //! @param PRINTER which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolEnsembleFilter::ModifyPrinter
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
    void ProtocolEnsembleFilter::ModifyFactory
    (
      util::ShPtr< pdb::Factory> &FACTORY,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief initialize the mutates and add them to Mutates enumerator
    void ProtocolEnsembleFilter::InitializeMutates()
    {
    }

    //! @brief modify the mutate tree used
    //! @param MUTATE_TREE MutateTree to be modified
    void ProtocolEnsembleFilter::ModifyMutateTree( MutateTree &MUTATE_TREE) const
    {
    }

    //! @brief get the mutate tree associated with this protocol
    //! @return the mutate tree associated with this protocol
    util::ShPtr< MutateTree> ProtocolEnsembleFilter::GetMutateTree() const
    {
      util::ShPtr< MutateTree> sp_mutate_tree( new MutateTree());
      ModifyMutateTree( *sp_mutate_tree);
      return sp_mutate_tree;
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void ProtocolEnsembleFilter::MergeAndModifyMutateTree( MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }

  ////////////
  // readme //
  ////////////

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ProtocolEnsembleFilter::GetDescription() const
    {
      // initialize string to store the description
      static const std::string s_description
      (
        "Protocol to be applied for filtering an ensemble."
      );

      // end
      return s_description;
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ProtocolEnsembleFilter::GetReadMe() const
    {
      // TODO: add readme for this protocol
      // initialize string to store the readme information
      static const std::string s_readme
      (
        "readme for protocol to be applied for filtering an ensemble."
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
    std::istream &ProtocolEnsembleFilter::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProtocolEnsembleFilter::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace fold
  
} // namespace bcl
