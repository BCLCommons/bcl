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
#include "fold/bcl_fold_protocol_default.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_printer_protein_model_movie.h"
#include "assemble/bcl_assemble_protein_model_inverter.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "fold/bcl_fold_default_flags.h"
#include "fold/bcl_fold_default_mutates.h"
#include "fold/bcl_fold_default_scores.h"
#include "mc/bcl_mc_movie_printers.h"
#include "mc/bcl_mc_printer_combined.h"
#include "mc/bcl_mc_stage.h"
#include "opti/bcl_opti_criterion_number_iterations.h"
#include "opti/bcl_opti_criterion_unimproved.h"
#include "pdb/bcl_pdb_printer_score.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProtocolDefault::ProtocolDefault()
    {
    }

    //! @brief Clone function
    //! @return pointer to new ProtocolDefault
    ProtocolDefault *ProtocolDefault::Clone() const
    {
      return new ProtocolDefault( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief access to only instance
    //! @return reference to only instance
    ProtocolDefault &ProtocolDefault::GetInstance()
    {
      static ProtocolDefault s_protocol_instance;
      return s_protocol_instance;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProtocolDefault::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ProtocolDefault::GetAlias() const
    {
      static const std::string s_name( "ProtocolDefault");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProtocolDefault::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "default folding protocol");

      return serializer;
    }

  ///////////
  // flags //
  ///////////

    //! @brief returns all flags that are specialized for this protocol
    //! @return all flags that are specialized for this protocol
    const util::ShPtrVector< command::FlagInterface> &ProtocolDefault::GetAllFlags() const
    {
      return DefaultFlags::GetAllFlags();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief modifies the start model
    //! @param START_MODEL Start model to be modified
    void ProtocolDefault::ModifyStartModel( assemble::ProteinModel &START_MODEL) const
    {
      // if a specified start model was not given and the model already does not have any SSEs
      if
      (
        !DefaultFlags::GetFlagStartModel()->GetFlag() &&
        START_MODEL.GetNumberSSEs() == 0
      )
      {
        AddRandomSSEToStartModel( START_MODEL);
      }
    }

    //! @brief initialize the scores and add them to Scores enumerator
    void ProtocolDefault::InitializeScores()
    {
      DefaultScores::GetInstance().InitializeScores();
    }

    //! @brief modify the score weight set
    //! @param SCORE_WEIGHT_SET Score weight set
    void ProtocolDefault::ModifyScoreWeightSet( ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
      DefaultScores::GetInstance().ModifyScoreWeightSet( SCORE_WEIGHT_SET);
    }

    //! @brief modify the terminate object
    //! @param CRITERION which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolDefault::ModifyCriterion
    (
      opti::CriterionCombine< assemble::ProteinModel, double> &CRITERION,
      const mc::Stage &STAGE
    ) const
    {
      CRITERION.InsertCriteria
      (
        opti::CriterionNumberIterations< assemble::ProteinModel, double>
        (
          STAGE.GetMaxNumberIterations()
        )
      );

      CRITERION.InsertCriteria
      (
        opti::CriterionUnimproved< assemble::ProteinModel, double>
        (
          STAGE.GetMaxNumberUnimprovedIterations()
        )
      );
    }

    //! @brief modify the printer object
    //! @param PRINTER which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolDefault::ModifyPrinter
    (
      mc::PrinterCombined< assemble::ProteinModel, double> &PRINTER,
      const mc::Stage &STAGE
    ) const
    {
      // construct model printer
      util::ShPtr< assemble::PrinterProteinModel> sp_model_printer
      (
        new assemble::PrinterProteinModel
        (
          GetSetup().GetPrefix(),
          GetSetup().GetStorage(),
          GetSetup().GetSuperimposeMeasure()
        )
      );
      PRINTER.Insert( sp_model_printer);

      // add movie printer if movie flag was given
      if( mc::MoviePrinterInterface::GetFlagMoviePrinter()->GetFlag())
      {
        // initialize mc movie printer
        util::ShPtr< mc::MoviePrinterInterface> sp_movie_printer
        (
          mc::MoviePrinter( mc::MoviePrinterInterface::GetParameterMoviePrinterType()->GetValue())
        );
        storage::Vector< std::string> row_names( STAGE.GetScoreFunction()->GetFunctionSchemes());
        row_names.PushBack( GetSetup().GetSuperimposeMeasure().GetName());
        sp_movie_printer->Initialize
        (
          mc::MoviePrinterInterface::GetParameterMoviePrefix()->GetValue(),
          storage::Vector< std::string>::Create
          (
            GetStaticClassName< storage::Table< double> >(),
            math::SumFunctionMixin< score::ProteinModel>::GetValueTableVerticalColumnNames()( 2)
          ),
          row_names,
          mc::MoviePrinterInterface::GetParameterMovieWidth()->GetNumericalValue< size_t>(),
          mc::MoviePrinterInterface::GetParameterMovieHeight()->GetNumericalValue< size_t>(),
          bool( mc::MoviePrinterInterface::GetParameterRayTrace()->GetNumericalValue< size_t>() > 0)
        );

        // construct protein model movie printer
        const util::ShPtr< assemble::PrinterProteinModelMovie> sp_model_movie_printer
        (
          new assemble::PrinterProteinModelMovie
          (
            mc::MoviePrinterInterface::GetParameterMoviePrefix()->GetValue(),
            sp_movie_printer,
            STAGE.GetScoreFunction(),
            GetSetup().GetStepStatuses(),
            GetSetup().GetSuperimposeMeasure(),
            GetSetup().GetQualityMeasures()
          )
        );

        // if pdb was given, use this as the starting structure
        if( DefaultFlags::GetFlagNativeModel()->GetFirstParameter()->GetWasSetInCommandLine())
        {
          sp_movie_printer->SetStartFrame( DefaultFlags::GetFlagNativeModel()->GetFirstParameter()->GetValue(), true);
        }

        // set the water mark
        sp_movie_printer->SetWaterMark();

        // add printer
        PRINTER.Insert( sp_model_movie_printer);
      }

    }

    //! @brief modify the pdb factory object
    //! @param FACTORY pdb factory to be modified
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolDefault::ModifyFactory
    (
      util::ShPtr< pdb::Factory> &FACTORY,
      const mc::Stage &STAGE
    ) const
    {
      FACTORY->AppendPrinter
      (
        util::ShPtr< util::FunctionInterface< assemble::ProteinModel, util::ShPtrList< pdb::Line> > >
        (
          new pdb::PrinterScore( STAGE.GetScoreFunction(), GetSetup().GetQualityMeasures())
        )
      );
    }

    //! @brief initialize the mutates and add them to Mutates enumerator
    void ProtocolDefault::InitializeMutates()
    {
      DefaultMutates::GetInstance().InitializeMutates();
    }

    //! @brief modify the mutate tree used
    //! @param MUTATE_TREE MutateTree to be modified
    void ProtocolDefault::ModifyMutateTree( MutateTree &MUTATE_TREE) const
    {
      DefaultMutates::GetInstance().ModifyMutateTree( MUTATE_TREE);
    }

    //! @brief reset the protocol (i.e. empty cache of scoring function)
    void ProtocolDefault::Reset()
    {
      // get the inverter
      util::ShPtr< assemble::ProteinModelInverter> sp_inverter( DefaultScores::GetInstance().GetProteinInverter());

      // reset it
      sp_inverter->Reset();
    }

    //! @brief get the mutate tree associated with this protocol
    //! @return the mutate tree associated with this protocol
    util::ShPtr< MutateTree> ProtocolDefault::GetMutateTree() const
    {
      util::ShPtr< MutateTree> sp_mutate_tree( new MutateTree());
      ModifyMutateTree( *sp_mutate_tree);
      return sp_mutate_tree;
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void ProtocolDefault::MergeAndModifyMutateTree( MutateTree &MUTATE_TREE) const
    {
      ModifyMutateTree( MUTATE_TREE);
    }

  ////////////
  // readme //
  ////////////

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ProtocolDefault::GetDescription() const
    {
      // initialize string to store the description
      static const std::string s_description
      (
        "default protocol for de novo folding"
      );

      // end
      return s_description;
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ProtocolDefault::GetReadMe() const
    {
      // initialize string to store the readme information
      static const std::string s_readme
      (
        "The Default protocol initializes all mutates and scores. These are then modified by other protocols to suit "
        "specific needs."
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
    std::istream &ProtocolDefault::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProtocolDefault::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief adds a random SSE from the pool to the start model
    //! @param START_MODEL starting model
    void ProtocolDefault::AddRandomSSEToStartModel( assemble::ProteinModel &START_MODEL)
    {
      // get the pool from start model
      util::ShPtr< assemble::SSEPool> sp_pool
      (
        START_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Pool)
      );
      BCL_Assert( sp_pool.IsDefined(), "SSE pool is not initialized!");

      // insert a random sse from the pool
      START_MODEL.Insert( random::GetGlobalRandom().Iterator( sp_pool->Begin(), sp_pool->End(), sp_pool->GetSize())->HardCopy());
    }

  } // namespace fold
} // namespace bcl
