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
#include "fold/bcl_fold_protocol_template.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_collector_sheet.h"
#include "fold/bcl_fold_mutate_protein_model_add_sheet_from_template.h"
#include "fold/bcl_fold_mutate_protein_model_domain.h"
#include "fold/bcl_fold_mutate_sheet_fit_to_template.h"
#include "fold/bcl_fold_mutate_tree.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProtocolTemplate::ProtocolTemplate()
    {
    }

    //! @brief returns a pointer to a new ProtocolTemplate
    //! @return pointer to a new ProtocolTemplate
    ProtocolTemplate *ProtocolTemplate::Clone() const
    {
      return new ProtocolTemplate( *this);
    }

  ///////////
  // flags //
  ///////////

    //! @brief returns all flags that are specialized for this protocol
    //! @return all flags that are specialized for this protocol
    const util::ShPtrVector< command::FlagInterface> &ProtocolTemplate::GetAllFlags() const
    {
      // initialize static ShPtrVector of FlagInterfaces to form the comment line
      static util::ShPtrVector< command::FlagInterface> s_all_flags_vector;

      // if the flag vector is initialize for the first time
      if( s_all_flags_vector.IsEmpty())
      {
        // insert all the flags in the vector
//        s_all_flags_vector.PushBack( assemble::FoldTemplateHandler::GetFlagFoldTemplates());
//        s_all_flags_vector.PushBack( assemble::SheetTemplateHandler::GetFlagSheetTemplates());
//        s_all_flags_vector.PushBack( GetFlagFoldSSEPairTemplates());
      }

      // end
      return s_all_flags_vector;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the only instance of this class
    //! @return the only instance of this class
    ProtocolTemplate &ProtocolTemplate::GetInstance()
    {
      static ProtocolTemplate s_protocol_instance;
      return s_protocol_instance;
    }

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &ProtocolTemplate::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ProtocolTemplate::GetAlias() const
    {
      static const std::string s_name( "ProtocolTemplate");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProtocolTemplate::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "protocol for de novo folding of proteins using fold templates");
      serializer.AddInitializer
        (
         "mutate sheet fit to template",
         "selects a sheet from the protein model and fits it to a sheet template from the database",
         io::Serialization::GetAgent( &e_MutateSheetFitToTemplate)
         );
      serializer.AddInitializer
        (
         "mutate add sheet from template",
         "selects strands from the SSE pool and fits them into sheet template from the database",
         io::Serialization::GetAgent( &e_MutateSheetFitToTemplate)
         );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief modifies the start model
    //! @param START_MODEL Start model to be modified
    void ProtocolTemplate::ModifyStartModel( assemble::ProteinModel &START_MODEL) const
    {
      // create the mutate that adds a sheet template to the model
      util::ShPtr< math::MutateInterface< assemble::ProteinModel> > sp_mutate
      (
        new MutateProteinModelAddSheetFromTemplate()
      );

      // apply the mutate on the model
      math::MutateResult< assemble::ProteinModel> result( sp_mutate->operator()( START_MODEL));
      if( result.GetArgument().IsDefined())
      {
        START_MODEL = *result.GetArgument();
      }
    }

    //! @brief initialize the scores and add them to Scores enumerator
    void ProtocolTemplate::InitializeScores()
    {
    }

    //! @brief modify the score weight set
    //! @param SCORE_WEIGHT_SET Score weight set
    void ProtocolTemplate::ModifyScoreWeightSet( ScoreWeightSet &SCORE_WEIGHT_SET) const
    {
    }

    //! @brief modify the terminate object
    //! @param CRITERION which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolTemplate::ModifyCriterion
    (
      opti::CriterionCombine< assemble::ProteinModel, double> &CRITERION,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief modify the printer object
    //! @param PRINTER which will be modified by protocols
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolTemplate::ModifyPrinter
    (
      mc::PrinterCombined< assemble::ProteinModel, double> &PRINTER,
      const mc::Stage &STAGE
    ) const
    {
    }

    //! @brief modify the pdb factory object
    //! @param FACTORY pdb factory to be modified
    //! @param STAGE Stage in which this terminate will be used
    void ProtocolTemplate::ModifyFactory( util::ShPtr< pdb::Factory> &FACTORY, const mc::Stage &STAGE) const
    {
    }

    //! @brief initialize the mutates and add them to Mutates enumerator
    void ProtocolTemplate::InitializeMutates()
    {
//      // mutate that bends a SSE according to a template
//      e_MutateSSEBendTemplate = GetMutates().AddMutate
//      (
//        MutateProteinModelSSE
//        (
//          util::CloneToShPtr( assemble::LocatorSSERandom()),
//          util::CloneToShPtr( MutateSSEBendTemplate( assemble::SSEGeometryWithinSizeTolerance( 0, 0))),
//          "sse_bend_template"
//        )
//      );

      // mutate that fits a sheet to template
      if( !GetMutates().HaveEnumWithName( "sheet_fit_to_template"))
      {
        e_MutateSheetFitToTemplate = GetMutates().AddMutate
        (
          MutateProteinModelDomain( assemble::CollectorSheet(), MutateSheetFitToTemplate(), "sheet_fit_to_template")
        );
      }

      if( !GetMutates().HaveEnumWithName( MutateProteinModelAddSheetFromTemplate().GetScheme()))
      {
        // mutate that forms a sheet from strands in the SSE pool and fits it to a template
        e_MutateAddSheetFromTemplate = GetMutates().AddMutate( MutateProteinModelAddSheetFromTemplate());
      }

    }

    //! @brief modify the mutate tree used
    //! @param MUTATE_TREE MutateTree to be modified
    void ProtocolTemplate::ModifyMutateTree( MutateTree &MUTATE_TREE) const
    {
      MUTATE_TREE.SetMutateProbability( e_MutateSheetFitToTemplate, 0.75);
      MUTATE_TREE.SetMutateProbability( e_MutateAddSheetFromTemplate, 0.25);
    }

    //! @brief get the mutate tree associated with this protocol
    //! @return the mutate tree associated with this protocol
    util::ShPtr< MutateTree> ProtocolTemplate::GetMutateTree() const
    {
      util::ShPtr< MutateTree> sp_mutate_tree( new MutateTree());
      ModifyMutateTree( *sp_mutate_tree);
      return sp_mutate_tree;
    }

    //! @brief merges this protocol's mutate tree into given mutate tree
    //! @param MUTATE_TREE tree into which to merge this protocol's tree
    void ProtocolTemplate::MergeAndModifyMutateTree( MutateTree &MUTATE_TREE) const
    {
      util::ShPtr< MutateTree> sp_mutate_tree( GetMutateTree());
      MUTATE_TREE.Merge( *sp_mutate_tree);
    }

  ////////////
  // readme //
  ////////////

    //! @brief returns readme information
    //! @return string containing information about this protocol
    const std::string &ProtocolTemplate::GetDescription() const
    {
      // initialize string to store the description
      static const std::string s_description
      (
        "protocol for using structural fold templates in the structure prediction process"
      );

      // end
      return s_description;
    }

    //! @brief returns readme information
    //! @return string containing information about this protocol
    const std::string &ProtocolTemplate::GetReadMe() const
    {
      // initialize string to store the readme information
      static const std::string s_readme
      (
        "readme for template protocol"
      );

      // end
      return s_readme;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read object from input stream
    //! @param ISTREAM input stream to read object from
    //! @return input stream which was read from
    std::istream &ProtocolTemplate::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write object into  output stream
    //! @param OSTREAM output stream to write object into
    //! @param INDENT number of indentations to separate members
    //! @return output stream object was written into
    std::ostream &ProtocolTemplate::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
