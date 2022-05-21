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
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container_generator_protein_model.h"
#include "assemble/bcl_assemble_aa_neighbor_vector.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_extension.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "fold/bcl_fold_default_flags.h"
#include "pdb/bcl_pdb_factory.h"
#include "restraint/bcl_restraint_epr_accessibility_data.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RestraintSimulateAccessibility
    //! @brief application for creating accessibility restraint for use in the bcl from a protein model
    //!
    //! @author alexanns
    //! @date Feb 1, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RestraintSimulateAccessibility :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      //! flag for pdb file from which the restraints will be simulated
      util::ShPtr< command::FlagInterface> m_PDBFilename;

      //! flag for specifying the output file path and name
      util::ShPtr< command::FlagInterface> m_OutputFile;

      //! flag for specifying the environment types the calculated exposures should be labeled as
      util::ShPtr< command::FlagInterface> m_AccessibilityEnvironments;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      RestraintSimulateAccessibility();

      //! @brief Clone function
      //! @return pointer to new RestraintSimulateAccessibility
      RestraintSimulateAccessibility *Clone() const
      {
        return new RestraintSimulateAccessibility( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      // This function is used to retrieve prior names for applications that are still valid (for backwards compatibility)
      // but which should not be displayed, e.g. for help
      storage::Vector< std::string> GetDeprecatedAppNames() const
      {
        return storage::Vector< std::string>( size_t( 1), "SimulateAccessibilityRestraints");
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

    ///////////////
    // operators //
    ///////////////

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
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      static const ApplicationType RestraintSimulateAccessibility_Instance;

    }; // class RestraintSimulateAccessibility

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RestraintSimulateAccessibility::RestraintSimulateAccessibility() :
      m_PDBFilename
      (
        new command::FlagStatic
        (
          pdb::GetDefaultFileExtension(),
          "\tpdb file for which restraints will be created",
          command::Parameter
          (
            "pdb_filename",
            "\tfilename for which restraints will be calculated",
              command::ParameterCheckExtension( "." + pdb::GetDefaultFileExtension()),
            ""
          )
        )
      ),
      m_OutputFile
      (
        new command::FlagStatic
        (
          "output_file",
          "Path and name of the output file which will hold the simulated restraints.",
          command::Parameter
          (
            "string",
            "\tpath and name of the outputfile",
            "restraints" + restraint::EPRAccessibilityData().GetDefaultExtension()
          )
        )
      ),
      m_AccessibilityEnvironments
      (
        new command::FlagDynamic
        (
          "accessibility_environments",
          "\tenvironments that exposures should be simulated as.",
          command::Parameter
          (
            "environment_type",
            "\tenvironment type the calculated exposure will be labeled as",
            command::ParameterCheckSerializable( restraint::AccessibilityAA::EnvironmentEnum())
          ),
          0,
          restraint::AccessibilityAA::s_NumberEnvironmentTypes
        )
      )
    {
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> RestraintSimulateAccessibility::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // add flag
      sp_cmd->AddFlag( m_PDBFilename);

      // add flag
      sp_cmd->AddFlag( m_OutputFile);

      // add flag
      sp_cmd->AddFlag( m_AccessibilityEnvironments);

      // add flag for min sse size of created model
      sp_cmd->AddFlag( pdb::Factory::GetFlagMinSSESize());

      sp_cmd->AddFlag( fold::DefaultFlags::GetFlagPDBIDNumbering());

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int RestraintSimulateAccessibility::Main() const
    {
      pdb::Factory factory;
      // create ProteinModel "protein_model" and initialize with protein of "pdb_filename"
      const assemble::ProteinModel
        protein_model( factory.ProteinModelFromPDBFilename( m_PDBFilename->GetFirstParameter()->GetValue()));

      const util::ShPtr< assemble::AAExposureInterface> exposure_calculator( new assemble::AANeighborVector());

      // generate the neighbor list from the protein
      const util::ShPtr< math::FunctionInterfaceSerializable< assemble::ProteinModel, assemble::AANeighborListContainer> > neighbor_list_generator
      (
        assemble::AANeighborListContainerGeneratorProteinModel::AANeighborListGenerator
        (
          exposure_calculator->GetDistanceCutoff(),
          exposure_calculator->GetMinimalSequenceSeparation(),
          true,
          false
        )
      );

      // make the neighbor lists
      const assemble::AANeighborListContainer neighbor_lists( neighbor_list_generator->operator()( protein_model));

      storage::List< restraint::AccessibilityAA> restraints;

      const storage::Vector< restraint::AccessibilityAA::EnvironmentEnum> environments
      (
        m_AccessibilityEnvironments->GetObjectList< restraint::AccessibilityAA::EnvironmentEnum>()
      );

      // iterate over the neighbor lists
      for
      (
        assemble::AANeighborListContainer::const_iterator
          neighbor_list_itr( neighbor_lists.Begin()), neighbor_list_itr_end( neighbor_lists.End());
        neighbor_list_itr != neighbor_list_itr_end;
        ++neighbor_list_itr
      )
      {
        // the current neighbor list
        const assemble::AANeighborList &current_neighbor_list( neighbor_list_itr->second);

        const util::SiPtr< const biol::AABase> &residue( current_neighbor_list.GetCenterAminoAcid());

        // true if the current neighbor list is not valid
        if( residue == util::SiPtr< const biol::AABase>( assemble::AANeighborList::GetDefaultCenterAA()))
        {
          BCL_MessageStd( "current neighbor list is being skipped");
          continue;
        }

        // the exposure of the current residue in the neighbor list
        const double exposure( exposure_calculator->operator()( current_neighbor_list));

        // locator for the current residue
        const util::ShPtr
        <
          assemble::LocatorAA
        > aa_locator( new assemble::LocatorAA( residue->GetChainID(), residue->GetSeqID()));

        storage::Map< restraint::AccessibilityAA::EnvironmentEnum, double> data;

        // add the accessibility data for each environment type to data
        for
        (
          storage::Vector< restraint::AccessibilityAA::EnvironmentEnum>::const_iterator
            environment_itr( environments.Begin()), environment_itr_end( environments.End());
          environment_itr != environment_itr_end;
          ++environment_itr
        )
        {
          data.Insert( std::make_pair( *environment_itr, exposure));
        }

        restraints.PushBack( restraint::AccessibilityAA( data, aa_locator, exposure_calculator));
      }

      const restraint::AccessibilityProfile profile( restraints);

      restraint::HandlerAccessibilityAA handler( exposure_calculator);

      io::OFStream ofstream;
      io::File::MustOpenOFStream( ofstream, m_OutputFile->GetFirstParameter()->GetValue());
      restraint::HandlerAccessibilityAA::WriteRestraints( ofstream, profile);

      return 0;
    }

    const ApplicationType RestraintSimulateAccessibility::RestraintSimulateAccessibility_Instance
    (
      GetAppGroups().AddAppToGroup( new RestraintSimulateAccessibility(), GetAppGroups().e_Restraint)
    );

  } // namespace app
} // namespace bcl
