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
#include "bcl_app_optimize.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_ensemble.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "fold/bcl_fold_protocols.h"
#include "io/bcl_io_file.h"
#include "mc/bcl_mc_optimization_mcm.h"
#include "opti/bcl_opti_ensemble_node.h"
#include "opti/bcl_opti_pipeline.h"
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
    const ApplicationType Optimize::Optimizes_Instance
    (
      GetAppGroups().AddAppToGroup( new Optimize(), GetAppGroups().e_InternalBiol)
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Optimize::Optimize() :
      m_PDBList
      (
        new command::FlagStatic
        (
          "pdb_list",
          "\tfile containing paths to the PDBs that should be optimized",
          command::Parameter( "pdb_list_filename", "\tfile name of the list of input PDBs")
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
      m_OptimizerFileName
      (
        new command::FlagStatic
        (
          "optimizer",
          "\tfile defining the optimizer used for the optimization or prediction",
          command::Parameter( "optimizer_filename", "\tfile name of the optimizer")
        )
      )
    {
    }

    //! @brief clone function
    //! @return pointer to a new Optimize
    Optimize *Optimize::Clone() const
    {
      return new Optimize( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &Optimize::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns a shared pointer to the command object
    //! @return shared pointer to the command object
    util::ShPtr< command::Command> Optimize::InitializeCommand() const
    {
      // add the BCL and the application specific flags to the command object
      util::ShPtr< command::Command> sp_cmd( new command::Command());
      sp_cmd->AddFlag( m_PDBList);
      sp_cmd->AddFlag( m_OutputPrefix);
      sp_cmd->AddFlag( m_OptimizerFileName);
      sp_cmd->AddFlag( sspred::GetMethods().GetFlagReadSSPredictions());
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      return sp_cmd;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief the main function of this application
    //! @return exit code - 0 for success
    int Optimize::Main() const
    {
      mc::OptimizationMCM::s_Instance.IsDefined();

      // initialize all scores from all protocols
      for( fold::Protocols::iterator itr( fold::GetProtocols().Begin()), itr_end( fold::GetProtocols().End()); itr != itr_end; ++itr)
      {
        ( **itr)->InitializeScores();
      }

      // read in the list of protein models to optimize
      io::IFStream pdb_list;
      io::File::MustOpenIFStream( pdb_list, m_PDBList->GetFirstParameter()->GetValue());
      storage::Vector< std::string> pdb_paths( util::StringLineListFromIStream( pdb_list));
      io::File::CloseClearFStream( pdb_list);

      // create the initial ensemble
      pdb::Factory pdb_factory;
      assemble::Ensemble< assemble::ProteinModel> ensemble;
      for( auto pdb_it( pdb_paths.Begin()), pdb_it_end( pdb_paths.End()); pdb_it != pdb_it_end; ++pdb_it)
      {
        // create a protein model from the current PDB file and add it to the ensemble
        BCL_MessageStd( "Reading PDB file " + *pdb_it);
        const std::string &pdb_path( *pdb_it);
        assemble::ProteinModel model( pdb_factory.ProteinModelFromPDBFilename( pdb_path));
        ensemble.AddElement( model);
      }

      // create the optimizer for the protein ensemble
      util::Implementation< opti::OptimizationInterface< assemble::Ensemble< assemble::ProteinModel> > > optimizer;
      io::IFStream optimizer_file;
      io::File::MustOpenIFStream( optimizer_file, m_OptimizerFileName->GetFirstParameter()->GetValue());
      optimizer_file >> optimizer;
      io::File::CloseClearFStream( optimizer_file);

      // optimize the protein ensemble
      ( *optimizer)( ensemble);

      // write out the results
      size_t model_index( 0);
      const std::string output_prefix( m_OutputPrefix->GetFirstParameter()->GetValue());
      for( auto result_it( ensemble.Begin()), result_it_end( ensemble.End()); result_it != result_it_end; ++result_it)
      {
        const std::string new_file_name( output_prefix + util::Format()( model_index) + ".pdb");
        io::OFStream write;
        io::File::MustOpenOFStream( write, new_file_name);
        pdb_factory.WriteModelToPDB( result_it->GetElement(), write);
        io::File::CloseClearFStream( write);
        ++model_index;
      }

      return 0;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads members from a given input stream
    //! @param ISTREAM input stream to read members from
    //! @return input stream which members were read from
    std::istream &Optimize::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief writes members into a given output stream
    //! @param OSTREAM output stream to write members into
    //! @param INDENT number of indentations
    //! @return output stream into which members were written
    std::ostream &Optimize::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace app
} // namespace bcl
