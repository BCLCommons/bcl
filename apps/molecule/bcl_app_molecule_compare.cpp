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
#include "bcl_app_molecule_compare.h"

// includes from bcl - sorted alphabetically
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
    MoleculeCompare::MoleculeCompare() :
      m_InputFilenameA
      (
        new command::Parameter
        (
          "fragment_filename",
          "filename for input sdf of fragment.  All of them will be compared with one another, unless the second parameter"
          "(input_filename) is given, in which case molecules in fragment_filename will be compared with input_filename",
          command::ParameterCheckFileExistence()
        )
      ),
      m_InputFilenameB
      (
        new command::Parameter
        (
          "input_filename",
          "filename for input sdf of molecules.  If none is given, all molecules in the fragment file will be "
          "compared with all others in the same file",
          command::ParameterCheckFileExistence(),
          ""
        )
      ),
      m_OutputFilenameFlag
      (
        new command::FlagStatic
        (
          "output",
          "filename for similarity matrix",
          command::Parameter( "output_filename", "filename for similarity matrix")
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
            command::ParameterCheckSerializable( util::Implementation< chemistry::ConformationComparisonInterface>()),
            "RMSD"
          )
        )
      ),
      m_OutputTableFormat
      (
        new command::FlagStatic
        (
          "bcl_table_format",
          "flag for outputting the scores in bcl table format"
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
      m_StartB
      (
        new command::FlagStatic
        (
          "ensemble_b_start",
          "flag for indicating which mol from ensemble a to load in first",
          command::Parameter
          (
            "index",
            "index of ensemble b molecules to start with",
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
      ),
      m_MaxMolsB
      (
        new command::FlagStatic
        (
          "ensemble_b_max",
          "flag for indicating maximum number of molecules to take from ensemble b",
          command::Parameter
          (
            "max",
            "max number of molecules from ensemble b",
            command::ParameterCheckRanged< size_t>( 0, math::GetHighestBoundedValue< size_t>()),
            util::Format()( math::GetHighestBoundedValue< size_t>())
          )
        )
      )
    {
    }

    //! @brief Clone function
    //! @return pointer to new FoldProtein
    MoleculeCompare *MoleculeCompare::Clone() const
    {
      return new MoleculeCompare( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MoleculeCompare::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string MoleculeCompare::GetDescription() const
    {
      return "compare molecules by spatial, property, fingerprint, or substructural features";
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> MoleculeCompare::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // ensemble containing the fragments
      sp_cmd->AddParameter( m_InputFilenameA);

      // ensemble containing the molecules to be checked for the fragments
      sp_cmd->AddParameter( m_InputFilenameB);

      // flag indication of first molecule to load from ensembles A and B
      sp_cmd->AddFlag( m_StartA);
      sp_cmd->AddFlag( m_StartB);

      // flag indication of last molecule to load from ensembles A and B
      sp_cmd->AddFlag( m_MaxMolsA);
      sp_cmd->AddFlag( m_MaxMolsB);

      // matrix file that will be written out
      sp_cmd->AddFlag( m_OutputFilenameFlag);

      // comparison
      sp_cmd->AddFlag( m_ConformerComparerFlag);

      // format
      sp_cmd->AddFlag( m_OutputTableFormat);

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
    const std::string &MoleculeCompare::GetReadMe() const
    {
      util::Implementation< descriptor::Base< chemistry::AtomConformationalInterface, float> >::SetHaveDisplayedHelp();
      static io::FixedLineWidthWriter writer;
      static std::string s_read_me
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of how to run molecule:Compare, terms of use, "
        "appropriate citation.\n\n"
        + DefaultSectionSeparator() +
        "II. WHAT IS BCL::MoleculeCompare?\n"
        "molecule:Compare is a C++ based application, created by Vanderbilt University's Meiler Laboratory, which is part "
        "of a larger library of applications called BCL::Commons.  molecule:Compare allows molecules in an ensemble to "
        "be compared with one another (if only one file is given) or against another ensemble of molecules (if two are "
        "given). Comparison algorithms vary from tanimoto/substructure based, to fingerprint/hash-based, "
        "3D-coordinate-based, steric/electronic-based, and hybrid methods.\n"
        "\n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        + DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. APPROPRIATE CITATIONS FOR USING BCL::MoleculeCompare.\n"
        "\n"
        "Kothiwale S, Mendenhall JL, Meiler J. "
        "BCL::Conf: small molecule conformational sampling using a knowledge based rotamer library. "
        "J Cheminformatics, 2015, Sep 30\n"
        "Journal link: http://www.jcheminf.com/content/7/1/47\n"
        "\n"
        + DefaultSectionSeparator() +
        "V. INSTALLATION PROCEDURES.\n"
        + DefaultInstallationProcedure() +
        "\n"
        + DefaultSectionSeparator() +
        "VI. RUNNING molecule:Compare.\n"

        "Input: * 1 - 2 sdf files to compare molecules within. \n"
        "  * If only one file is given, all pairs of molecules in that file will be compared\n"
        "    The output will be a square matrix of distances/similarities. \n"
        "  * If two files (A,B) are given, all molecules in file A will be compared with all\n"
        "    molecules in B, yielding a rectangular matrix. \n"
        "\n"
        "Example command line:\n"
        " bcl.exe A.sdf -output A.rmsd -method RMSD -scheduler PThread 4 -remove_h\n"
        "  Compare all molecules in A.sdf by 3D-RMSD, allowing for super-imposition. Hydrogens will be ignored\n"
        " -add_h can be given instead of -remove_h to always consider hydrogens.\n"
        " The -method flag takes many options (shown in detail below). Be sure to . Complete help for these commands is "
        " given below: "
        + util::Implementation< chemistry::ConformationComparisonInterface>::WriteInstancesHelp( writer).String()
        + "\n"
        "\nADDITIONAL INFORMATION\n"
        "\n"
        "Runtime performance considerations: -remove_h should ordinarily be used when using substructure-based "
        " metrics like LargestCommonSubstructure because it decreases the runtime by a factor of 10-100 or more on large "
        "(>60 heavy atom) molecules. Only molecules of size < 20-30 are ordinarily comparable using disconnected "
        " substructure search with hydrogens present. Connected substructure search is often several times faster. "
        "\n"
        "An example property field correlation flag that we routinely use for comparing similarity of molecules is:\n"
        "Combine( Multiply(Atom_SigmaCharge, Atom_VDWVolume), "
        "Multiply(GreaterEqual(lhs=BondTypeCount(property=IsAromatic,value=1),rhs=2),0.694), "
        "Multiply(Atom_HbondAcceptors, Not(Atom_HbondDonors), 0.943), "
        "Multiply(Atom_HbondDonors, 1.98)) . The property lengths are often set of 0.75. This function emphasizes both "
        " steric and electronic considerations when comparing molecules that have been pre-aligned spatially.\n"
        "\n"
        "Information concerning syntax and flags can be obtained by typing bcl.exe molecule:Compare -help\n"
        "\n"
        "For more general information about the product, type bcl.exe molecule:Compare -readme\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        + DefaultTechnicalSupportString() +
        "\n"
        + DefaultSectionSeparator() +
        "VIII. FUTURE DEVELOPMENT OF BCL::MoleculeCompare\n"
        "BCL::MoleculeCompare is under ongoing further development.\n"
        "\n"
        + DefaultSectionSeparator()
      );
      util::Implementation< descriptor::Base< chemistry::AtomConformationalInterface, float> >::ResetHaveDisplayedHelp();
      return s_read_me;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int MoleculeCompare::Main() const
    {
      // read in ensemble_a
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_InputFilenameA->GetValue());
      // remove any hydrogens because they slow down the scaffold search and are unnecessary
      math::Range< size_t> ens_a_load_rng( size_t( 0), math::GetHighestBoundedValue< size_t>());
      ens_a_load_rng.SetMin( util::ConvertStringToNumericalValue< size_t>( m_StartA->GetFirstParameter()->GetValue()));
      const size_t n_to_load_a( util::ConvertStringToNumericalValue< size_t>( m_MaxMolsA->GetFirstParameter()->GetValue()));
      ens_a_load_rng.SetMax( math::GetHighestBoundedValue< size_t>() - n_to_load_a > ens_a_load_rng.GetMin() ? ens_a_load_rng.GetMin() + n_to_load_a - 1 : math::GetHighestBoundedValue< size_t>());
//      ens_a_load_rng.SetMax( math::GetHighestBoundedValue< size_t>() >= ens_a_load_rng.GetMin() + n_to_load_a ? ens_a_load_rng.GetMin() + n_to_load_a - 1 : math::GetHighestBoundedValue< size_t>());
      chemistry::FragmentEnsemble ensemble_a( input, sdf::GetCommandLineHydrogensPref(), ens_a_load_rng);
      io::File::CloseClearFStream( input);
      m_EnsembleA = util::SiPtrVector< const chemistry::ConformationInterface>( ensemble_a.Begin(), ensemble_a.End());
      m_EnsembleASize = m_EnsembleA.GetSize();

      // test whether only one ensemble was given or whether the same filename was given for both ensembles
      m_IdenticalEnsembles =
        !m_InputFilenameB->GetWasSetInCommandLine()
        ||
        (
          ( io::DirectoryEntry( m_InputFilenameB->GetValue()).GetFullName() == io::DirectoryEntry( m_InputFilenameA->GetValue()).GetFullName())
          && m_StartA->GetFirstParameter()->GetValue() == m_StartB->GetFirstParameter()->GetValue()
          && m_MaxMolsA->GetFirstParameter()->GetValue() == m_MaxMolsB->GetFirstParameter()->GetValue()
        );
      chemistry::FragmentEnsemble ensemble_b;

      if( !m_IdenticalEnsembles)
      {
        // read in ensemble_b
        io::File::MustOpenIFStream( input, m_InputFilenameB->GetValue());
        // remove any hydrogens because they slow down the scaffold search and are unnecessary
        math::Range< size_t> ens_b_load_rng( size_t( 0), math::GetHighestBoundedValue< size_t>());
        ens_b_load_rng.SetMin( util::ConvertStringToNumericalValue< size_t>( m_StartB->GetFirstParameter()->GetValue()));
        const size_t n_to_load_b( util::ConvertStringToNumericalValue< size_t>( m_MaxMolsB->GetFirstParameter()->GetValue()));
        ens_b_load_rng.SetMax( math::GetHighestBoundedValue< size_t>() - n_to_load_b > ens_b_load_rng.GetMin() ? ens_b_load_rng.GetMin() + n_to_load_b - 1 : math::GetHighestBoundedValue< size_t>());
//        ens_b_load_rng.SetMax( math::GetHighestBoundedValue< size_t>() >= ens_b_load_rng.GetMin() + n_to_load_b ? ens_b_load_rng.GetMin() + n_to_load_b - 1 : math::GetHighestBoundedValue< size_t>());
        ensemble_b.ReadMoreFromMdl( input, sdf::GetCommandLineHydrogensPref(), ens_b_load_rng);
        io::File::CloseClearFStream( input);
        m_EnsembleB = util::SiPtrVector< const chemistry::ConformationInterface>( ensemble_b.Begin(), ensemble_b.End());
        m_EnsembleBSize = m_EnsembleB.GetSize();
        // reset # of pairs to consider
        m_NumberPairsToConsider = m_EnsembleASize * m_EnsembleBSize;
      }
      else
      {
        m_EnsembleB = m_EnsembleA;
        m_EnsembleBSize = m_EnsembleASize;
        // reset # of pairs to consider
        m_NumberPairsToConsider = m_EnsembleASize * ( m_EnsembleASize + 1) / 2;
      }

      if( m_EnsembleASize == 0 || m_EnsembleBSize == 0)
      {
        BCL_MessageCrt( "One of the provided files did not contain any molecules; no output will be given");
        return 0;
      }

      m_NumberPairsConsidered = 0;
      m_LastAssignedMoleculeIds = storage::Pair< size_t, size_t>( 0, 0);
      m_ComparisonMatrix = linal::Matrix< double>( m_EnsembleASize, m_EnsembleBSize, 0.0);

      // get the number of processors
      const size_t n_processors( sched::GetNumberCPUs());

      util::Implementation< chemistry::ConformationComparisonInterface> comparer_templ
      (
        m_ConformerComparerFlag->GetFirstParameter()->GetValue()
      );
      comparer_templ->PrepareEnsemble( ensemble_a);
      comparer_templ->PrepareEnsemble( ensemble_b);
      m_Comparers.Resize( n_processors, comparer_templ);
      BCL_MessageStd( "check molecule compare");

      // create a vector to hold the jobs
      util::ShPtrVector< sched::JobInterface> jobs;

      // make the job vector big enough to hold all the jobs
      jobs.AllocateMemory( n_processors);
      const size_t group_id( 0);
      storage::Vector< size_t> index_vec( storage::CreateIndexVector( n_processors));
      for( size_t processor_number( 0); processor_number < n_processors; ++processor_number)
      {
        // create the job
        jobs.PushBack
        (
          util::ShPtr< sched::JobInterface>
          (
            new sched::UnaryFunctionJobWithData< const size_t, void, MoleculeCompare>
            (
              group_id,
              *this,
              &MoleculeCompare::RunThread,
              index_vec( processor_number),
              sched::JobInterface::e_READY,
              NULL
            )
          )
        );

        // submit it to the scheduler
        sched::GetScheduler().SubmitJob( jobs.LastElement());
      }

      // join all the jobs
      for( size_t processor_number( 0); processor_number < n_processors; ++processor_number)
      {
        sched::GetScheduler().Join( jobs( processor_number));
      }

      io::OFStream output;
      // open the output file so that we can write out scaffolds as we go
      io::File::MustOpenOFStream( output, m_OutputFilenameFlag->GetFirstParameter()->GetValue());
      // initialize bcl table format
      if( m_OutputTableFormat->GetFlag())
      {
        output << "bcl::storage::Table<double>";
        for( size_t column_id( 0), total( m_EnsembleBSize); column_id < total; ++column_id)
        {
          output << "\t" << util::Format()( column_id);
        }
        output << "\n";

        for( size_t row_id( 0); row_id < m_EnsembleASize; ++row_id)
        {
          output << util::Format()( row_id) << "\t";
          for( size_t column_id( 0); column_id < m_EnsembleBSize; ++column_id)
          {
            output << util::Format()( m_ComparisonMatrix( row_id, column_id)) << '\t';
          }
          output << '\n';
        }
      }
      else
      {
        output << m_ComparisonMatrix;
      }
      // close the file stream
      io::File::CloseClearFStream( output);

      return 0;
    }

    //! @brief get the next pair to compare
    //! @param IS_FIRST_PAIR whether this is the initial request from this thread for a pair
    storage::Pair< size_t, size_t> MoleculeCompare::GetNextPairToCompare( const bool &IS_FIRST_PAIR) const
    {
      m_GetNextPairMutex.Lock();

      // make a copy of the pair that will be returned
      storage::Pair< size_t, size_t> pair_to_compare( m_LastAssignedMoleculeIds);

      // test whether all pairs have been assigned to a thread to compare
      // if not, increment the second index
      if
      (
        m_LastAssignedMoleculeIds.First() < m_EnsembleASize
        && ++m_LastAssignedMoleculeIds.Second() >= m_EnsembleBSize
      )
      {
        if( !m_IdenticalEnsembles)
        {
          // second index reached end of file, increment first index and reset second
          m_LastAssignedMoleculeIds.Second() = 0;
          ++m_LastAssignedMoleculeIds.First();
        }
        else
        {
          // identical ensembles, only do upper triangle of comparisons
          m_LastAssignedMoleculeIds.Second() = ++m_LastAssignedMoleculeIds.First();
        }
      }

      if( !IS_FIRST_PAIR)
      {
        // log the status
        util::GetLogger().LogStatus
        (
          "% complete: "
          + util::Format().FFP( 2)( 100.0 * double( ++m_NumberPairsConsidered) / double( m_NumberPairsToConsider))
        );
        // molecule IDs in pair
        BCL_MessageStd( "First molecule in pair: " + util::Format()( m_LastAssignedMoleculeIds.First()));
        BCL_MessageStd( "Second molecule in pair: " + util::Format()( m_LastAssignedMoleculeIds.Second()));
      }

      // return a pair that is out of bounds, at which point the thread will stop
      m_GetNextPairMutex.Unlock();
      return pair_to_compare;
    }

    //! @brief compare the molecules given by the indices in a vector
    void MoleculeCompare::RunThread( const size_t &THREAD_ID) const
    {
      // keep track of what molecule # we are on
      for
      (
        storage::Pair< size_t, size_t> pair_to_compare( GetNextPairToCompare( true));
        pair_to_compare.First() < m_EnsembleASize;
        pair_to_compare = GetNextPairToCompare()
      )
      {
        const size_t mols_a_indice( pair_to_compare.First());
        const size_t mols_b_indice( pair_to_compare.Second());
        m_ComparisonMatrix( mols_a_indice, mols_b_indice)
          = m_Comparers( THREAD_ID)->operator ()( *m_EnsembleA( mols_a_indice), *m_EnsembleB( mols_b_indice));
        if( m_IdenticalEnsembles)
        {
          m_ComparisonMatrix( mols_b_indice, mols_a_indice) = m_ComparisonMatrix( mols_a_indice, mols_b_indice);
        }
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MoleculeCompare::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &MoleculeCompare::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    const ApplicationType MoleculeCompare::MoleculeCompare_Instance
    (
      GetAppGroups().AddAppToGroup( new MoleculeCompare(), GetAppGroups().e_Molecule)
    );

  } // namespace app
} // namespace bcl

