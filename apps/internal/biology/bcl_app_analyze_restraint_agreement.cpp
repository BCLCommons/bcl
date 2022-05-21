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
#include "app/bcl_app_apps.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_analyze_protein_ensemble_interface.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_protein_model_data.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "find/bcl_find_locator_coordinates_average.h"
#include "fold/bcl_fold_default_flags.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_head.h"
#include "sspred/bcl_sspred_method_handler.h"
#include "sspred/bcl_sspred_pdb.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{

  namespace find
  {
    // instantiate s_Instance
    template<>
    const util::SiPtr< const util::ObjectInterface> LocatorCoordinatesAverage< assemble::ProteinModel>::s_Instance
    (
      util::Enumerated< LocatorCoordinatesInterface< assemble::ProteinModel> >::AddInstance( new LocatorCoordinatesAverage< assemble::ProteinModel>())
    );
  } // namespace find
  namespace app
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinStatistics
    //! @brief analyzes a protein ensemble in the methods given by the command line
    //! @details
    //!
    //! @author alexanns
    //! @date 07/28/2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinStatistics :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      //! flag for specifying analysis types
      util::ShPtr< command::FlagInterface> m_AnalysisTypeFlag;

      //! flag for model identifiers within the ensemble
      util::ShPtr< command::FlagInterface> m_EnsembleModelIdentifiers;

      //! idealize flag
      util::ShPtr< command::FlagInterface> m_IdealizeFlag;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      ProteinStatistics();

    public:

      //! @brief Clone function
      //! @return pointer to new BenchmarkSSEPool
        ProteinStatistics *Clone() const
      {
        return new ProteinStatistics( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      // This function is used to retrieve prior names for applications that are still valid (for backwards compatibility)
      // but which should not be displayed, e.g. for help
      storage::Vector< std::string> GetDeprecatedAppNames() const
      {
        return storage::Vector< std::string>::Create( "AnalyzeRestraintAgreement");
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const
      {
        // initialize ShPtr to a new command object
        util::ShPtr< command::Command> sp_cmd( new command::Command());

        sp_cmd->AddFlag( restraint::GetFlagRestraintsTypes());

        sp_cmd->AddFlag( restraint::GetFlagRestraintsFilePrefix());

        sp_cmd->AddFlag( assemble::ProteinEnsemble::GetFlagPDBList());

        sp_cmd->AddFlag( assemble::AnalyzeProteinEnsembleInterface::GetFlagOutFilePrefix());

        sp_cmd->AddFlag( m_AnalysisTypeFlag);
        sp_cmd->AddFlag( m_EnsembleModelIdentifiers);
        sp_cmd->AddFlag( m_IdealizeFlag);

        sp_cmd->AddFlag( fold::DefaultFlags::GetFlagPDBIDNumbering());

        // Add All Flags from pdb::Factory
        for
        (
          util::ShPtrVector< command::FlagInterface>::const_iterator
           itr( pdb::Factory::GetAllFlags().Begin()), itr_end( pdb::Factory::GetAllFlags().End());
          itr != itr_end;
          ++itr
        )
        {
          sp_cmd->AddFlag( *itr);
        }

        // add default bcl parameters
        command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

        // return assembled Command object
        return sp_cmd;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      static const ApplicationType s_Instance;

    }; // class ProteinStatistics

    //! @brief the Main function
    //! @return error code - 0 for success
    int ProteinStatistics::Main() const
    {
      // get the analyses that are desired to be done
      storage::Vector< util::Implementation< assemble::AnalyzeProteinEnsembleInterface> > analyses
      (
        m_AnalysisTypeFlag->GetObjectList< util::Implementation< assemble::AnalyzeProteinEnsembleInterface> >()
      );

      // get the protein ensemble
      const std::string pdb_list_file( assemble::ProteinEnsemble::GetFlagPDBList()->GetFirstParameter()->GetValue());
      const size_t pdb_column
      (
        assemble::ProteinEnsemble::GetFlagPDBList()->GetParameterList()( 1)->GetNumericalValue< size_t>()
      );

      // read the flag from the command line and create the correct AAClass
      biol::AAClass aa_class( pdb::Factory::GetFlagAAClass()->GetFirstParameter()->GetValue());

      // Build the ensemble from the correct AAClass
      assemble::ProteinEnsemble ensemble( pdb_list_file, pdb_column, aa_class, std::string(), true);

      // iterate over the protein ensemble to set environment type for each protein model
      for
      (
        util::ShPtrVector< assemble::ProteinModel>::iterator protein_model_itr( ensemble.Begin()), protein_model_itr_end( ensemble.End());
          protein_model_itr != protein_model_itr_end;
        ++protein_model_itr
      )
      {

        // get the current protein model
        assemble::ProteinModel &current_protein_model( **protein_model_itr);

        if( m_IdealizeFlag->GetFlag())
        {
          current_protein_model.SetToIdealConformation();
        }
      }

      if( m_EnsembleModelIdentifiers->GetFlag())
      {
        BCL_Assert
        (
          ensemble.SetIdentifiers( m_EnsembleModelIdentifiers->GetStringList()),
          "unable to set identifiers list within ensemble"
        );
      }

      // iterate through the analysis types and write the analysis
      for
      (
        storage::Vector< util::Implementation< assemble::AnalyzeProteinEnsembleInterface> >::const_iterator
          analysis_itr( analyses.Begin()), analysis_itr_end( analyses.End());
        analysis_itr != analysis_itr_end;
        ++analysis_itr
      )
      {
        ( *analysis_itr)->WriteAnalysisFile
        (
          assemble::AnalyzeProteinEnsembleInterface::GetFlagOutFilePrefix()->GetFirstParameter()->GetValue(), ensemble
        );
      }

      return 0;
    }

    //! default constructor
    ProteinStatistics::ProteinStatistics() :
      m_AnalysisTypeFlag
      (
        new command::FlagDynamic
        (
          "analysis_type_enumerated",
          "method to analyze protein ensemble",
          command::Parameter
          (
            "analysis type",
            "method to analyze protein ensemble",
            command::ParameterCheckSerializable( util::Implementation< assemble::AnalyzeProteinEnsembleInterface>())
          ),
          0,
          5
        )
      ),
      m_EnsembleModelIdentifiers
      (
        new command::FlagDynamic
        (
          "ensemble_model_identifiers",
          "identifier for each model in the ensemble, for easier reading of the pymol object names",
          command::Parameter( "identifier", "short identifier for the model in the ensemble, e.g. PDB 4 letter code", pdb::Head::GetBCLPdbID()),
          0, util::GetUndefined< size_t>()
        )
      ),
      m_IdealizeFlag
      (
        new command::FlagStatic
        (
          "idealize",
          "set model to have idealized sses prior to computing statistics"
        )
      )
    {
    }

    const ApplicationType ProteinStatistics::s_Instance
    (
      GetAppGroups().AddAppToGroup( new ProteinStatistics(), GetAppGroups().e_Restraint)
    );

  } // namespace app
} // namespace bcl
