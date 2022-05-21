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
#include "bcl_app_template_model.h"

// includes from bcl - sorted alphabetically
#include "align/bcl_align_handler_classes.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_quality.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "fold/bcl_fold_default_flags.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_mutate_result.h"
#include "pdb/bcl_pdb_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {
  //////////
  // data //
  //////////

    const ApplicationType TemplateModel::TemplateModel_Instance
    (
      GetAppGroups().AddAppToGroup( new TemplateModel(), GetAppGroups().e_Protein)
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    TemplateModel::TemplateModel() :
      m_FlagPrintProblemAAs
      (
        new command::FlagStatic
        (
          "print_problem_aas",
          "\tif set, potentially structurally problematic residues will be printed to logger"
        )
      )
    {
    }

    //! @brief Clone function
    //! @return pointer to new TemplateModel
    TemplateModel *TemplateModel::Clone() const
    {
      return new TemplateModel( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &TemplateModel::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string TemplateModel::GetDescription() const
    {
      return "Uses an alignment to assign coordinates from a known structure onto a sequence of unknown structure";
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> TemplateModel::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // print zeros for undefined coordinates
      sp_cmd->AddFlag( pdb::Factory::GetFlagWriteZeroCoordinatesForUndefinedAminoAcids());

      // input alignments for chains
      sp_cmd->AddFlag( assemble::Quality::GetFlagAlignments());

      // flag for input format of the alignments
      sp_cmd->AddFlag( align::HandlerClasses< biol::AABase>::GetFlagInputFormat());

      // flag for specifying the pdb whose coordinates will be used
      sp_cmd->AddFlag( fold::DefaultFlags::GetFlagStartModel());

      // flag for specifying the prefix to prepend to the output filename
      sp_cmd->AddFlag( fold::DefaultFlags::GetFlagPrefix());

      // flag to output the list of residues which should be checked manually
      sp_cmd->AddFlag( m_FlagPrintProblemAAs);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &TemplateModel::GetReadMe() const
    {
      static const std::string readme
      (
        DefaultSectionSeparator() +
        "Uses an alignment to assign coordinates from a known structure onto a sequence of unknown structure. "
        "A sequence alignment between a sequence with known coordinates and a sequence with unknown coordinates is "
        "used to assign the coordinates from the given protein model to the second sequence in the alignment. The "
        "protein model must have the same sequence as the first sequence in the alignment. The provided start"
        "model is used as the template whose coordinates are used. The first sequence in the alignment must"
        "correspond to the sequence of the template start model. The output file generated is "
        "\"<prefix>threaded_<template_pdb_name_without_path>\"\n"
        + DefaultSectionSeparator()
      );
      return readme;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int TemplateModel::Main() const
    {
      // get the mutate for threading
      const fold::MutateProteinModelThreadSequence mutate( GetMutate());

      // get the template pdb
      const io::DirectoryEntry template_pdb
      (
        fold::DefaultFlags::GetFlagStartModel()->GetFirstParameter()->GetValue()
      );

      // get the template model
      const assemble::ProteinModel template_model
      (
        pdb::Factory().ProteinModelFromPDBFilename( template_pdb.GetFullName())
      );

      // mutate the protein model - i.e. get the threaded model
      const math::MutateResult< assemble::ProteinModel> mutate_result( mutate( template_model));
      BCL_Assert( mutate_result.GetArgument().IsDefined(), "result of mutate is not defined");

      // the name of the pdb file that will be created
      const std::string out_pdb_file
      (
        fold::DefaultFlags::GetFlagPrefix()->GetFirstParameter()->GetValue() + "threaded_" + template_pdb.GetName()
      );

      // write the pdb
      io::OFStream write;
      io::File::MustOpenOFStream( write, out_pdb_file);
      pdb::Factory().WriteModelToPDB( *mutate_result.GetArgument(), write);
      io::File::CloseClearFStream( write);

      //successful end
      return 0;
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
    std::istream &TemplateModel::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &TemplateModel::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief gives the mutate for threading sequence
    //! @return mutate for threading sequence onto template structure
    fold::MutateProteinModelThreadSequence TemplateModel::GetMutate() const
    {
      // the handler to read in the alignments
      const align::HandlerClasses< biol::AABase>::HandlerClass handler
      (
        align::GetHandlerClasses< biol::AABase>().GetEnumFromName
        (
          align::HandlerClasses< biol::AABase>::GetFlagInputFormat()->GetFirstParameter()->GetValue()
        )
      );

      // get the alignment files
      const storage::Vector< io::DirectoryEntry> align_files
      (
        assemble::Quality::GetFlagAlignments()->GetObjectList< io::DirectoryEntry>()
      );

      // get the map of chain ids to the corresponding alignment
      const storage::Map< char, util::ShPtr< align::AlignmentInterface< biol::AABase> > > chain_align
      (
        assemble::Quality::GetChainAlignments( **handler, align_files)
      );

      // create the mutate
      const fold::MutateProteinModelThreadSequence mutate( chain_align, "thread mutate", m_FlagPrintProblemAAs->GetFlag());

      // return the mutate
      return mutate;
    }

  } // namespace app
} // namespace bcl
