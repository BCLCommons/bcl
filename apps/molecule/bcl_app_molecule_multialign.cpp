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
#include "bcl_app_molecule_multialign.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_aa_fragment_complete.h"
#include "chemistry/bcl_chemistry_conformation_comparison_interface.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_thunk_job.h"
#include "sched/bcl_sched_unary_function_job_with_data.h"
#include "storage/bcl_storage_template_instantiations.h"
#include "util/bcl_util_implementation.h"
namespace bcl
{
  namespace app
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief standard constructor
    MoleculeMultiAlign::MoleculeMultiAlign() :
      m_EnsembleASize( 0),
      m_MultiAlignScore( math::GetHighestBoundedValue< double>()),
      m_InputMolSDF
      (
        new command::Parameter
        (
          "fragment_filename",
          "filename for input sdf of fragment ensemble",
          command::ParameterCheckFileExistence()
        )
      ),
      m_InputPocketPDB
      (
        new command::Parameter
        (
          "input_filename",
          "filename for protein binding pocket(s) .sdf file. If not given, no bounding"
          "box will be placed around the multiple alignment",
          command::ParameterCheckFileExistence(),
          ""
        )
      ),
      m_ConformerComparerFlag
      (
        new command::FlagStatic
        (
          "method",
          "method to compare molecules with",
          command::Parameter
          (
            "",
            "",
            command::ParameterCheckSerializable( chemistry::ConformationComparisonMultiAlign())
          )
        )
      ),
      m_StartA
      (
        new command::FlagStatic
        (
          "ensemble_a_start",
          "flag for indicating which mol from ensemble a to load in first",
          command::Parameter
          (
            "index",
            "index of ensemble a molecules to start with",
            command::ParameterCheckRanged< size_t>( 0, math::GetHighestBoundedValue< size_t>()),
            "0"
          )
        )
      ),
      m_MaxMolsA
      (
        new command::FlagStatic
        (
          "ensemble_a_max",
          "flag for indicating maximum number of molecules to take from ensemble a",
          command::Parameter
          (
            "max",
            "max number of molecules from ensemble a",
            command::ParameterCheckRanged< size_t>( 0, math::GetHighestBoundedValue< size_t>()),
            util::Format()( math::GetHighestBoundedValue< size_t>())
          )
        )
      )
    {
    }

    //! @brief Clone function
    //! @return pointer to new FoldProtein
    MoleculeMultiAlign *MoleculeMultiAlign::Clone() const
    {
      return new MoleculeMultiAlign( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MoleculeMultiAlign::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string MoleculeMultiAlign::GetDescription() const
    {
      return "compare molecules by spatial, property, fingerprint, or substructural features";
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> MoleculeMultiAlign::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // ensemble containing the fragments
      sp_cmd->AddParameter( m_InputMolSDF);

      // ensemble containing the protein binding pocket residues
      sp_cmd->AddParameter( m_InputPocketPDB);

      // flag indication of first molecule to load from ensembles A and B
      sp_cmd->AddFlag( m_StartA);

      // flag indication of last molecule to load from ensembles A and B
      sp_cmd->AddFlag( m_MaxMolsA);

      // comparison
      sp_cmd->AddFlag( m_ConformerComparerFlag);

      // strict checking
      sp_cmd->AddFlag( chemistry::ConformationComparisonInterface::GetDisableStrictAtomBondTypeCheckingFlag());

      // hydrogen flags
      sdf::AddMoleculeIOPrefFlags( *sp_cmd);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &MoleculeMultiAlign::GetReadMe() const
    {
      util::Implementation< descriptor::Base< chemistry::AtomConformationalInterface, float> >::SetHaveDisplayedHelp();
      static io::FixedLineWidthWriter writer;
      static std::string s_read_me
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of how to run molecule:MultiAlign, terms of use, "
        "appropriate citation.\n\n"
        + DefaultSectionSeparator() +
        "II. WHAT IS BCL::MoleculeMultiAlign?\n"
        "molecule:MultiAlign is a C++ based application, created by Vanderbilt University's Meiler Laboratory, which is part "
        "of a larger library of applications called BCL::Commons.  molecule:MultiAlign allows simultaneous flexible "
        "alignment of small molecules given an input ensemble. \n"
        "\n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        + DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. APPROPRIATE CITATIONS FOR USING BCL::MoleculeMultiAlign.\n"
        "\n"

        "\n"
        + DefaultSectionSeparator() +
        "V. INSTALLATION PROCEDURES.\n"
        + DefaultInstallationProcedure() +
        "\n"
        + DefaultSectionSeparator() +
        "VI. RUNNING molecule:MultiAlign.\n"
        "\n"
        "Example command line:\n"
        "Information concerning syntax and flags can be obtained by typing bcl.exe molecule:MultiAlign -help\n"
        "\n"
        "For more general information about the product, type bcl.exe molecule:MultiAlign -readme\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        + DefaultTechnicalSupportString() +
        "\n"
        + DefaultSectionSeparator() +
        "VIII. FUTURE DEVELOPMENT OF BCL::MoleculeMultiAlign\n"
        "BCL::MoleculeMultiAlign is under ongoing further development.\n"
        "\n"
        + DefaultSectionSeparator()
      );
      util::Implementation< descriptor::Base< chemistry::AtomConformationalInterface, float> >::ResetHaveDisplayedHelp();
      return s_read_me;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int MoleculeMultiAlign::Main() const
    {
      // read in ensemble_a
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_InputMolSDF->GetValue());

      // remove any hydrogens because they slow down the scaffold search and are unnecessary
      math::Range< size_t> ens_a_load_rng( size_t( 0), math::GetHighestBoundedValue< size_t>());
      ens_a_load_rng.SetMin( util::ConvertStringToNumericalValue< size_t>( m_StartA->GetFirstParameter()->GetValue()));
      const size_t n_to_load_a( util::ConvertStringToNumericalValue< size_t>( m_MaxMolsA->GetFirstParameter()->GetValue()));
      ens_a_load_rng.SetMax( math::GetHighestBoundedValue< size_t>() - n_to_load_a > ens_a_load_rng.GetMin() ? ens_a_load_rng.GetMin() + n_to_load_a - 1 : math::GetHighestBoundedValue< size_t>());
      chemistry::FragmentEnsemble ensemble_a( input, sdf::GetCommandLineHydrogensPref(), ens_a_load_rng);
      io::File::CloseClearFStream( input);

      // check to make sure molecule ensemble_A is not empty
      if( !ensemble_a.GetSize())
      {
        BCL_MessageCrt( "The provided SDF file did not contain any molecules; no output will be given");
        return 0;
      }

      // check for second input filename
      chemistry::FragmentEnsemble ensemble_b;
      if( m_InputPocketPDB->GetWasSetInCommandLine())
      {
        // read in protein
        assemble::ProteinModel protein_model;
        const pdb::Factory factory;
        protein_model = factory.ProteinModelFromPDBFilename( m_InputPocketPDB->GetValue());
        chemistry::AAFragmentComplete aa_fragment( protein_model.GetAminoAcids(), true);
        aa_fragment.RemoveH();
        chemistry::FragmentEnsemble ensemble_b_temp( storage::List< chemistry::FragmentComplete>( 1, aa_fragment));

        // check to make sure molecule ensemble_B_temp is not empty
        if( !ensemble_b_temp.GetSize())
        {
          BCL_MessageCrt( "The provided binding pocket SDF file did not contain any molecules; no output will be given");
          return 0;
        }
        ensemble_b = ensemble_b_temp;
      }
      //otherwise just
      else
      {
        ensemble_b = ensemble_a;
      }

      m_Comparers.AssertRead( util::ObjectDataLabel( m_ConformerComparerFlag->GetFirstParameter()->GetValue()));
      BCL_MessageStd( "MultiAlignScore = " + util::Format()( m_Comparers( ensemble_a, ensemble_b)));
      return 0;
    }

//    //! @brief compare the molecules given by the indices in a vector
//    void MoleculeMultiAlign::RunThread( const size_t &THREAD_ID) const
//    {
//      m_MultiAlignScore = m_Comparers->operator ()( *m_EnsembleA);
//    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MoleculeMultiAlign::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &MoleculeMultiAlign::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    const ApplicationType MoleculeMultiAlign::MoleculeMultiAlign_Instance
    (
      GetAppGroups().AddAppToGroup( new MoleculeMultiAlign(), GetAppGroups().e_Molecule)
    );

  } // namespace app
} // namespace bcl

