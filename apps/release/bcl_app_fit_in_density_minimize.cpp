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
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_quality.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "command/bcl_command_parameter_check_extension.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "density/bcl_density_fit_protein_minimizers.h"
#include "density/bcl_density_map.h"
#include "density/bcl_density_protein_agreements.h"
#include "density/bcl_density_simulators.h"
#include "io/bcl_io_file.h"
#include "mc/bcl_mc_temperature_accepted.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "storage/bcl_storage_table.h"
#include "util/bcl_util_stopwatch.h"

#undef SetPrinter // visual studio defines that macro

namespace bcl
{
  namespace app
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FitInDensityMinimize
    //! @brief Class is an application for fitting a protein structure into a density map using an initial fit
    //! @details Using a monte carlo minimization, a protein model is fitted into the density map using an initial
    //! position, that should be close. The correlation between the given experimental density and the simulated
    //! density for the protein model is optimized.
    //!
    //! @see @link example_app_fit_in_density_minimize.cpp @endlink
    //! @author woetzen
    //! @date 07/18/2009
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FitInDensityMinimize :
      public InterfaceRelease
    {
    private:

    //////////
    // data //
    //////////

      //!input pdb and mrc file parameter
      util::ShPtr< command::ParameterInterface> m_PdbFilenameParam;

      //!input density filename
      util::ShPtr< command::ParameterInterface> m_MrcFilenameParam;

      //!flag to write the complete minimization
      util::ShPtr< command::FlagInterface> m_WriteMinimizationFlag;

      //! flag for choice of density agreement function
      util::ShPtr< command::FlagInterface> m_DensityAgreementFlag;

      //! flag or choice of density simulator
      util::ShPtr< command::FlagInterface> m_SimulatorFlag;

      //! flag and parameter for adjusting the resolution of the given and simulated density
      util::ShPtr< command::FlagInterface> m_MrcResolutionFlag;

      //! flag and parameter for the output prefix of files that are written
      util::ShPtr< command::FlagInterface> m_OutputPrefix;

      //! flags for monte carlo minimization
      //! max number of rejected steps and max iterations
      util::ShPtr< command::FlagStatic> m_MaxIterationsFlag;

      //! mutate parameters
      util::ShPtr< command::FlagStatic> m_MaxTransRotFlag;
      util::ShPtr< command::ParameterInterface> m_MaxTransParam;
      util::ShPtr< command::ParameterInterface> m_MaxRotParam;

      //! approximator flag - which to use
      util::ShPtr< command::FlagStatic> m_ApproximatorFlag;

      //! vector containing all column names for result table
      static const storage::TableHeader s_TableHeader;

      //! density map
      mutable util::ShPtr< density::Map> m_SpMap;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      FitInDensityMinimize();

      //! @brief Clone function
      //! @return pointer to new FitInDensityMinimize
      FitInDensityMinimize *Clone() const
      {
        return new FitInDensityMinimize( *this);
      }

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief return the bcl::commons name
      //! @return string for the bcl::commons name of that application
      std::string GetBCLScopedName() const
      {
        static const std::string s_bcl_commons_name( "BCL::EMFitMinimize");
        return s_bcl_commons_name;
      }

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const
      {
        return "Protein model refinement using a density map. "
               "Optimizes the correlation between the given experimental density and the simulated density";
      }

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief main
      int Main() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

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

    private:

      //! @brief the stats table
      //! @param START_MODEL model before minimzation
      //! @param BEST_MODEL model after minimization
      storage::Table< double> CalculateStatsTable( const assemble::ProteinModel &START_MODEL, const assemble::ProteinModel &BEST_MODEL) const;

      static const ApplicationType FitInDensityMinimize_Instance;

    }; // class FitInDensityMinimize

    //! vector containing all column names for result table
    const storage::TableHeader FitInDensityMinimize::s_TableHeader
    (
      storage::Vector< std::string>::Create( "ccc_start", "ccc_min", "RMSD")
    );

    //! @brief main
    int FitInDensityMinimize::Main() const
    {
      //initialize read stream object
      io::IFStream read;

      //instantiate DensityMap from mrc file
      BCL_MessageStd( "read DensityMap from mrc file");
      io::File::MustOpenIFStream( read, m_MrcFilenameParam->GetValue(), std::ios::binary);
      m_SpMap = util::ShPtr< density::Map>( new density::Map());
      m_SpMap->ReadMRC( read);
      io::File::CloseClearFStream( read);
      BCL_MessageStd( "read Density map from mrc file done");

      //construct density map
      const std::string mrc_name( io::File::RemovePath( m_MrcFilenameParam->GetValue()));
      BCL_MessageStd( "name of mrc: " + mrc_name);

      //open pdb file
      io::File::MustOpenIFStream( read, m_PdbFilenameParam->GetValue());
      pdb::Handler pdb( read);
      io::File::CloseClearFStream( read);

      BCL_MessageStd( "fitting: " + m_PdbFilenameParam->GetValue());

      //create protein model from pdb
      pdb::Factory factory;
      const assemble::ProteinModel model( factory.ProteinModelFromPDB( pdb));

      BCL_MessageStd( "minimization started");

      util::Stopwatch minimize_time( "minimization " + m_ApproximatorFlag->GetFirstParameter()->GetValue(), util::Message::e_Standard, true);
      const density::FitProteinMinimizer minimizer_enum( m_ApproximatorFlag->GetFirstParameter()->GetValue());
      if( !minimizer_enum.IsDefined() && !minimizer_enum->IsDefined())
      {
        return -1;
      }

      util::ShPtr< density::FitProteinMinimizerInterface> sp_minimizer( minimizer_enum->HardCopy());
      sp_minimizer->SetResolution( m_MrcResolutionFlag->GetFirstParameter()->GetNumericalValue< double>());
      sp_minimizer->SetProteinAgreement( density::ProteinAgreement( m_DensityAgreementFlag->GetFirstParameter()->GetValue()));
      sp_minimizer->SetSimulator( density::Simulator( m_SimulatorFlag->GetFirstParameter()->GetValue()));
      sp_minimizer->SetMaxIterations( m_MaxIterationsFlag->GetFirstParameter()->GetNumericalValue< size_t>());
      sp_minimizer->SetMaxTranslationAndRotation
      (
        m_MaxTransParam->GetNumericalValue< double>(),
        m_MaxRotParam->GetNumericalValue< double>()
      );

      const assemble::ProteinModel best_fit( sp_minimizer->operator ()( model, *m_SpMap));

      // write minimized pdb to file
      {
        io::OFStream write;
        const std::string filename
        (
          m_OutputPrefix->GetFirstParameter()->GetValue() + "transformed_min.pdb"
        );
        BCL_MessageStd( "write transformed protein model in file " + filename);

        io::File::MustOpenOFStream( write, filename);
        factory.WriteModelToPDB( best_fit, write);
        io::File::CloseClearFStream( write);
      }

      BCL_MessageStd( "minimization finished");
      // print results
      //counts, calculated correlation factors, and rmsd
      const storage::Table< double> corr_corrmin( CalculateStatsTable( model, best_fit));
      corr_corrmin.WriteFormatted( util::GetLogger());
      BCL_MessageStd( "fit in density finished");

      //end
      return 0;
    } // end Main

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &FitInDensityMinimize::GetReadMe() const
    {
      // create a static string to hold readme information
      static const std::string s_readme_text
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of BCL::EMFitMinimize, terms of use, appropriate citation, installation "
        "procedures, BCL::EMFitMinimize execution, technical support, and future research directions.\n"
        "\n"
        + DefaultSectionSeparator() +
        "II. WHAT IS BCL::EMFitMinimize?\n"
        "BCL::EMFitMinimize is a C++ based application, created by Vanderbilt University's Meiler Laboratory, which is part of "
        "a larger library of applications called BCL::Commons.\n"
        "BCL::EMFitMinimize fits a given atomic protein structure in a medium resolution electron density map given an "
        "initial placement. It utilizes either a Monte Carlo/Metropolis simulated annealing optimization protocol or "
        "a Powell gradient based minimizer to minimze the cross correlation coefficient between the density map and "
        "the simulated density map for the protein. The key advantage is, that with a graphics card with OpenCL "
        "double precision support, significant speed ups can be achieved.\n"
        "\n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        + DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. APPROPRIATE CITATIONS FOR USING BCL::EMFitMinimize.\n"
        "When using BCL::EMFitMinimize in a publication, please cite the publications describing the application's "
        "development, of which one is currently under review. Any news will be published at http://www.meilerlab.org\n"
        "The publication for the OpenCL implementation of the Powell optimizer is:\n"
        "\n"
        "N. Woetzel, E. W. Lowe and J. Meiler, Poster: GPU-accelerated rigid body fitting of atomic structures into "
        "electron density maps, in 2011 IEEE 1st International Conference on Computational Advances in Bio and Medical "
        "Sciences (ICCABS), 2011, pp. 265–265."
        "\n"
        + DefaultSectionSeparator() +
        "V. INSTALLATION PROCEDURES.\n"
        + DefaultInstallationProcedure() +
        "Additionally, an Opencl driver is required which can be obtained from AMD (ATI_SDK), NVIDIA (CUDA_SDK) or IBM; "
        "they need to be able to support double precision.\n"
        "\n"
        + DefaultSectionSeparator() +
        "VI. RUNNING BCL::EMFitMinimize.\n"
        "Running BCL::EMFitMinimize consists of three main steps.\n"
        "\n"
        "1) Identify the density map of interest in which a protein will be fitted.\n"
        "The only format currently supported is the mrc/ccp4 format.  If sections of the density maps are to be "
        "considered, they will have to be cut manually. Any format conversions will have to be performed with other "
        "software.\n"
        "\n"
        "2) Create the protein structure of interest using the pdb format and generate an initial placement.\n"
        "If biomolecules (multimers) have to be fitted, BCL::PDBConvert or other tools to generate the appropriate pdb "
        " structure and store it in the pdb. An initial fit can be obtained using BCL::EMFit or any other tool.\n"
        "\n"
        "3) Run BCL::EMFitMinimize:\n"
        "Please read the publication on the different parameters and there meaning; also call BCL::EMFitMinimize on the "
        "command line with \"-help\" which documents the meaning of each flag. A sample command line could be:\n"
        "\n"
        "bcl_em_fit_minimize.exe initial.pdb groel.mrc -mrc_resolution 5.4 -write_minimization -density_agreement "
        "ProteinAgreementCCCOpencl -density_simulator OpenclGaussianSphere -approximator powellopencl "
        "-max_iterations_max_unimproved 120 10 -max_trans_rot 1.0 0.034 -opencl NVIDIA_CUDA TYPE_GPU 0\n"
        "\n"
        "PARAMETERS:\n"
        "initial.pdb groel.mrc are the protein and the density maps of interest\n"
        "\n"
        "FLAGS:\n"
        "\n"
        "-mrc_resolution 5.4 -> is the experimental resolution of the density map\n"
        "-write_minimization -> output a pdb for each minimization iteration, which can be used for movies\n"
        "-density_agreement ProteinAgreementCCCOpencl -> the objective function to minimize on; CCC on the opencl Device\n"
        "-density_simulator OpenclGaussianSphere -> the simulator to generate a density map from the atomic structure\n"
        "-approximator PowellPpencl -> to optimizer to use\n"
        "-max_iterations 120 -> the maximal number of iterations without improvement\n"
        "-max_trans_rot 1.0 0.034 -> the maximal step size for one iteration applied to the protein position in the "
        "density map, used for monte carlo and powell optimization (one line search)\n"
        "-opencl NVIDIA_CUDA TYPE_GPU 0 -> run opencl kernels on the nvidia platform, on GPU devices, device 0\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        + DefaultTechnicalSupportString() +
        "\n"
        + DefaultSectionSeparator() +
        "VIII. FUTURE DEVELOPMENT OF BCL::EMFitMinimize.\n"
        "BCL::EMFitMinimize is under ongoing further development. For current research please refer to www.meilerlab.org\n"
        + DefaultSectionSeparator()
      );

      // return readme information
      return s_readme_text;
    }

    util::ShPtr< command::Command> FitInDensityMinimize::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // pdb reader min sse size
      pdb::Factory::GetFlagMinSSESize()->GetParameterList()( biol::GetSSTypes().COIL)->SetDefaultParameter( "0");
      sp_cmd->AddFlag( pdb::Factory::GetFlagMinSSESize());

      // pdb reader default amino acid
      pdb::Factory::GetFlagAAClass()->GetParameterList()( 0)->SetDefaultParameter( biol::GetAAClasses().e_AAComplete);
      sp_cmd->AddFlag( pdb::Factory::GetFlagAAClass());

      //input pdb and mrc file parameter
      sp_cmd->AddParameter( m_PdbFilenameParam);

      //input density filename
      sp_cmd->AddParameter( m_MrcFilenameParam);

      // density map resolution parameters
      sp_cmd->AddFlag( m_MrcResolutionFlag);

      // flag for output prefix
      sp_cmd->AddFlag( m_OutputPrefix);

      sp_cmd->AddFlag( m_WriteMinimizationFlag);

      // flag for adjusting the agreement and density simulation to be used
      sp_cmd->AddFlag( m_DensityAgreementFlag);
      sp_cmd->AddFlag( m_SimulatorFlag);

      // flags for monte carlo minimization
      // adjust start and end condition
      sp_cmd->AddFlag( mc::TemperatureAccepted::GetFlagTemperature());

      //max number of rejected steps and max iterations
      sp_cmd->AddFlag( m_MaxIterationsFlag);

      //mutate parameters
      sp_cmd->AddFlag( m_MaxTransRotFlag);

      // approximator
      sp_cmd->AddFlag( m_ApproximatorFlag);

      // sses from backbone if no sse is given in pdb
      sp_cmd->AddFlag( pdb::Factory::GetFlagSSEsFromBackBone());

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief the stats table
    //! @param START_MODEL model before minimzation
    //! @param BEST_MODEL model after minimization
    storage::Table< double> FitInDensityMinimize::CalculateStatsTable( const assemble::ProteinModel &START_MODEL, const assemble::ProteinModel &BEST_MODEL) const
    {
      // construct a function that calculates the deviation between simulated density from atoms and the given density map
      util::ShPtr< density::ProteinAgreementInterface> density_agreement
      (
        density::GetProteinAgreements().CreateProteinAgreement
        (
          density::ProteinAgreement( m_DensityAgreementFlag->GetFirstParameter()->GetValue()),
          density::Simulator( m_SimulatorFlag->GetFirstParameter()->GetValue()),
          m_SpMap,
          m_MrcResolutionFlag->GetFirstParameter()->GetNumericalValue< double>()
        )
      );

      // counts, calculated correlation factors and rmsd
      // transform all coordinates and correlate them for search hashmap
      storage::Table< double> count_corr_rmsd_corrmin_rmsdmin( s_TableHeader);

      storage::Row< double> &current_row( count_corr_rmsd_corrmin_rmsdmin.InsertRow( "protein"));
      current_row[ "ccc_start"] = -density_agreement->operator()( START_MODEL);
      current_row[ "ccc_min"]   = -density_agreement->operator()( BEST_MODEL);
      current_row[ "RMSD"]      = assemble::Quality::Calculate( quality::GetMeasures().e_RMSD_NoSuperimposition, BEST_MODEL, START_MODEL, storage::Set< biol::AtomType>( biol::GetAtomTypes().CA));

      return count_corr_rmsd_corrmin_rmsdmin;
    }

    FitInDensityMinimize::FitInDensityMinimize() :
      m_PdbFilenameParam
      (
        new command::Parameter
        (
          "pdb_filename",
          "\tfilename for input pdb to be fitted in mrc electron density map",
          command::ParameterCheckExtension( ".pdb")
        )
      ),
      m_MrcFilenameParam
      (
        new command::Parameter
        (
          "mrc_filename",
          "\tfilename for input mrc to be used for fitting",
          command::ParameterCheckExtension( ".mrc")
        )
      ),
      m_WriteMinimizationFlag
      (
        new command::FlagStatic( "write_minimization", "\twrite pdbs of entire minimization")
      ),
      m_DensityAgreementFlag
      (
        new command::FlagStatic
        (
          "density_agreement",
          "choice of density agreement objectives",
          command::Parameter
          (
            "agreement_score",
            "the agreement to be used for quality of fit",
            command::ParameterCheckEnumerate< density::ProteinAgreements>(),
            density::GetProteinAgreements().e_CCC.GetName()
          )
        )
      ),
      m_SimulatorFlag
      (
        new command::FlagStatic
        (
          "density_simulator",
          "choice of the way the density is simulated",
          command::Parameter
          (
            "simulator",
            "the simulator to be used to generate density maps from a atom sturcture",
            command::ParameterCheckEnumerate< density::Simulators>(),
            density::GetSimulators().e_Gaussian.GetName()
          )
        )
      ),
      m_MrcResolutionFlag
      (
        new command::FlagStatic
        (
          "mrc_resolution",
          "\tresoltuion of given experimental density",
          command::Parameter
          (
            "density map resolution",
            "\tresolution of given electron density [A] - will be used to simulate density and calculate correlation"
          )
        )
      ),
      m_OutputPrefix
      (
        new command::FlagStatic
        (
          "prefix",
          "\toutput prefix for files written",
          command::Parameter
          (
            "output_prefix",
            "\tprefix that is prepended to each file that is written - can also contain directories",
            "./"
          )
        )
      ),
      m_MaxIterationsFlag
      (
        new command::FlagStatic
        (
          "max_iterations",
          "\tmodify the max number of iterations",
          command::Parameter
          (
            "max_iterations",
            "\tmaximal number of minimization iterations",
            command::ParameterCheckRanged< size_t>( 0, 1000),
            "250"
          )
        )
      ),
      m_MaxTransRotFlag
      (
        new command::FlagStatic
        (
          "max_trans_rot",
          "\tmaximal translation and rotation per iteration"
        )
      ),
      m_MaxTransParam
      (
        new command::Parameter
        (
          "translation",
          "\tmaximal translation [A]",
          command::ParameterCheckRanged< double>( 0.0, 10.0), // restriction
          "1.0" // default
        )
      ),
      m_MaxRotParam
      (
        new command::Parameter
        (
          "rotation",
          "\t\tmaximal rotation [rad]",
          command::ParameterCheckRanged< double>( 0, math::g_Pi), // restriction
          util::Format()( 2.0 / 180 * math::g_Pi) // default 2 degree
        )
      ),
      m_ApproximatorFlag
      (
        new command::FlagStatic
        (
          "approximator",
          "choosing approximator",
          command::Parameter
          (
            "approximator",
            "the approximator to use for fit minimization",
            command::ParameterCheckEnumerate< density::FitProteinMinimizers>(),
            density::GetFitProteinMinimizers().e_MC.GetName()
          )
        )
      )
    {
      // attach parameters to flags
      m_MaxTransRotFlag->PushBack( m_MaxTransParam);
      m_MaxTransRotFlag->PushBack( m_MaxRotParam);
    }

    const ApplicationType FitInDensityMinimize::FitInDensityMinimize_Instance
    (
      GetAppGroups().AddAppToGroup( new FitInDensityMinimize(), GetAppGroups().e_Density)
    );

  } // namespace app
} // namespace bcl
