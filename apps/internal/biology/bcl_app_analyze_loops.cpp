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
#include "bcl_app_analyze_loops.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const ApplicationType AnalyzeLoops::AnalyzeLoopss_Instance
    (
      GetAppGroups().AddAppToGroup( new AnalyzeLoops(), GetAppGroups().e_InternalBiol)
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AnalyzeLoops::AnalyzeLoops() :
      m_InputPDB
      (
        new command::FlagStatic
        (
          "input_pdb",
          "\tpaths to the PDB for which to analyze the loop regions",
          command::Parameter( "input_pdb_filename", "\tfile name of the input PDB")
        )
      ),
      m_OutputPrefix
      (
        new command::FlagStatic
        (
          "output_prefix",
          "\tprefix for the output files",
          command::Parameter( "filename", "\tprefix for the output files", "")
        )
      ),
      m_IgnoreTermLoops
      (
        new command::FlagStatic
        (
          "ignore_terminal_loops",
          "\tignore terminal loops"
        )
      )
    {
    }

    //! @brief clone function
    //! @return pointer to a new AnalyzeLoops
    AnalyzeLoops *AnalyzeLoops::Clone() const
    {
      return new AnalyzeLoops( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &AnalyzeLoops::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns a shared pointer to the command object
    //! @return shared pointer to the command object
    util::ShPtr< command::Command> AnalyzeLoops::InitializeCommand() const
    {
      // add the BCL and the application specific flags to the command object
      util::ShPtr< command::Command> sp_cmd( new command::Command());
      sp_cmd->AddFlag( m_InputPDB);
      sp_cmd->AddFlag( m_OutputPrefix);
      sp_cmd->AddFlag( m_IgnoreTermLoops);
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      return sp_cmd;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief the main function of this application
    //! @return exit code - 0 for success
    int AnalyzeLoops::Main() const
    {
      // read in the protein model, which will be analyzed
      const std::string pdb_filename( m_InputPDB->GetFirstParameter()->GetValue());
      BCL_MessageTop( "analyzing " + pdb_filename);
      const pdb::Factory factory( biol::GetAAClasses().e_AABackBone);
      const assemble::ProteinModel model( factory.ProteinModelFromPDBFilename( pdb_filename));

      // // construct loop regions for the protein models
      // const fold::AnalyzeLoopsConstruction loop_builder( sp_library);
      // const size_t num_models( m_NumModels->GetFirstParameter()->GetNumericalValue< size_t>());
      // const std::string output_prefix( m_OutputPrefix->GetFirstParameter()->GetValue());
      // for( size_t index( 0); index < num_models; ++index)
      // {
      //  // construct loop regions for the input model
      //  util::ShPtr< assemble::ProteinModel> sp_model_with_loops( loop_builder.ConstructAnalyzeLoopss( model));

      //  // write out the protein model with the constructed loop regions
      //  const std::string new_file_name( output_prefix + util::Format()( index) + ".pdb");
      //  printer.Print( *sp_model_with_loops, new_file_name);
      // }

      return 0;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads members from a given input stream
    //! @param ISTREAM input stream to read members from
    //! @return input stream which members were read from
    std::istream &AnalyzeLoops::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief writes members into a given output stream
    //! @param OSTREAM output stream to write members into
    //! @param INDENT number of indentations
    //! @return output stream into which members were written
    std::ostream &AnalyzeLoops::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace app
} // namespace bcl
