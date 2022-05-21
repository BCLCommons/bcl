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
#include "biol/bcl_biol_atom.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "command/bcl_command_parameter_check_extension.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "density/bcl_density_map.h"
#include "density/bcl_density_simulators.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_template_instantiations.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "util/bcl_util_stopwatch.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DensityFromPDB
    //! @brief Class generates a density map from a pdb
    //! @details Class reads the given pdb to and generates a density map for the given protein model
    //! has options to print out the fasta and bcl style pdb output for the given pdb
    //!
    //! @author woetzen, karakam
    //! @date 03/24/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DensityFromPDB :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      //! input pdb file parameter
      util::ShPtr< command::ParameterInterface> m_PdbFilenameParam;

      //! flag for changing the resolution requested for the produced density map
      util::ShPtr< command::FlagInterface> m_ResolutionFlag;

      //! flag for changing the voxel size requested for the produced density map
      util::ShPtr< command::FlagInterface> m_VoxelSizeFlag;

      //! flag for changing the smoothing kernel used in the process
      util::ShPtr< command::FlagInterface> m_SmoothingKernelFlag;

      //! flag for adding noise
      util::ShPtr< command::FlagInterface> m_NoiseFlag;

      //! flag for changing the output path and prefix
      util::ShPtr< command::FlagInterface> m_OutputPrefixFlag;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      DensityFromPDB();

    public:

      //! @brief Clone function
      //! @return pointer to new DensityFromPDB
      DensityFromPDB *Clone() const
      {
        return new DensityFromPDB( *this);
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
        return storage::Vector< std::string>( size_t( 1), "PDBToDensity");
      }

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

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      static const ApplicationType DensityFromPDB_Instance;

    }; // class DensityFromPDB

    //! @brief the Main function
    //! @return error code - 0 for success
    int DensityFromPDB::Main() const
    {
      const double resolution( m_ResolutionFlag->GetFirstParameter()->GetNumericalValue< double>());
      const double voxel_size( m_VoxelSizeFlag->GetFirstParameter()->GetNumericalValue< double>());
      const util::ShPtr< density::SimulateInterface> simulator
      (
        density::GetSimulators().CreateSimulator
        (
          density::Simulators::EnumType( m_SmoothingKernelFlag->GetFirstParameter()->GetValue()),
          linal::Vector3D( voxel_size), resolution
        )
      );

      //instantiate ofstream for writing mrc and pdb to file
      io::OFStream write;

      //read pdbfile
      io::IFStream read;
      BCL_MessageStd( "read pdbfile: " + util::Format()( m_PdbFilenameParam->GetValue()));
      io::File::MustOpenIFStream( read, m_PdbFilenameParam->GetValue());
      pdb::Handler pdb( read);
      io::File::CloseClearFStream( read);

      pdb::Factory factory;

      //build models from pdb sequences
      BCL_MessageStd( "building models from pdb chains and sse information");
      assemble::ProteinModel model( factory.ProteinModelFromPDB( pdb));

    //////////////////////////
    // simulate density map //
    //////////////////////////

      // density map calculated from SimpleAtoms, with target resolution, voxelsize( usually sigma = resolution/3) and
      // smoothingkernel of choice
      util::Stopwatch simulate_watch( m_SmoothingKernelFlag->GetFirstParameter()->GetValue(), util::Message::e_Standard, false);
      const density::Map density_map( simulator->operator ()( model.GetAtoms()));
      simulate_watch.WriteMessage();

      // pdb input path and filename
      const storage::VectorND< 2, std::string> pdb_path_filename
      (
        io::File::SplitToPathAndFileName( m_PdbFilenameParam->GetValue())
      );

      //write simulated density map to file
      const std::string mrc_filename
      (
        m_OutputPrefixFlag->GetFirstParameter()->GetValue() +
        io::File::RemoveLastExtension( pdb_path_filename.Second()) + "_res_" + util::Format().W( 3).FFP( 1)( resolution) +
        "voxelsize_" + util::Format().W( 5).FFP( 3)( voxel_size) + m_SmoothingKernelFlag->GetFirstParameter()->GetValue() +
        ".mrc"
      );

      // if no noise is to be added, write plain mrc map
      if( !m_NoiseFlag->GetFlag())
      {
        io::File::MustOpenOFStream( write, mrc_filename, std::ios::binary);
        density_map.WriteMRC( write);
        io::File::CloseClearFStream( write);
      }
      // add noise
      else
      {
        const double target_ccc( m_NoiseFlag->GetFirstParameter()->GetNumericalValue< double>());
        double noise_ccc( density_map.CrossCorrelationCoefficient( density_map, density_map.GetMinimum()));
        density::Map noised_density;
        for
        (
          size_t i( 0);
          i < 100 &&
            ( noise_ccc > target_ccc || util::Format().W( 4).FFP( 3)( noise_ccc) == util::Format().W( 4).FFP( 3)( target_ccc));
          ++i
        )
        {
          // reset noised density
          noised_density = density_map;
          noised_density.AddNoise( random::GetGlobalRandom(), double( 0.0), double( i) * density_map.GetRmsd() / double( 20));
          noise_ccc = noised_density.CrossCorrelationCoefficient( density_map, density_map.GetMinimum());
        }

        if( noise_ccc > target_ccc)
        {
          BCL_MessageCrt
          (
            "unable to generate noised map with CCC to original map of: " + util::Format()( target_ccc) +
            " after 100 iterations"
          );

          return 1;
        }

        // write noised density map to file
        const std::string mrc_noised_filename
        (
          m_OutputPrefixFlag->GetFirstParameter()->GetValue() +
          io::File::RemoveLastExtension( pdb_path_filename.Second()) + "_res_" +
          util::Format().W( 3).FFP( 1)( resolution) + "voxelsize_" + util::Format().W( 5).FFP( 3)( voxel_size) +
          m_SmoothingKernelFlag->GetFirstParameter()->GetValue() + "_noiseccc_" + util::Format().W( 5).FFP( 3)( noise_ccc) + ".mrc"
        );

        io::File::MustOpenOFStream( write, mrc_noised_filename, std::ios::binary);
        noised_density.WriteMRC( write);
        io::File::CloseClearFStream( write);
      }

      return 0;
    } // end Main

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> DensityFromPDB::InitializeCommand() const
    {
      // initialize ShPtr to a new command object
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // input pdb file parameter
      sp_cmd->AddParameter( m_PdbFilenameParam);

      // flag for changing the resolution requested for the produced density map
      sp_cmd->AddFlag( m_ResolutionFlag);

      // flag for changing the voxel size requested for the produced density map
      sp_cmd->AddFlag( m_VoxelSizeFlag);

      // flag for changing the smoothing kernel used in the process
      sp_cmd->AddFlag( m_SmoothingKernelFlag);

      // flag for adding noise
      sp_cmd->AddFlag( m_NoiseFlag);

      // flag for prefix
      sp_cmd->AddFlag( m_OutputPrefixFlag);

      // flag to change the aa class to be used
      pdb::Factory::GetFlagAAClass()->GetParameterList()( 0)->SetDefaultParameter
      (
        biol::GetAAClasses().e_AAComplete.GetName()
      );
      sp_cmd->AddFlag( pdb::Factory::GetFlagAAClass());
      sp_cmd->AddFlag( pdb::Factory::GetFlagMinSSESize());

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! default constructor
    DensityFromPDB::DensityFromPDB() :
      m_PdbFilenameParam
      (
        new command::Parameter
        (
          "pdb_filename",
          "\tfilename for input pdb to be fitted in mrc electron density map",
          command::ParameterCheckExtension( ".pdb")
        )
      ),
      m_ResolutionFlag
      (
        new command::FlagStatic
        (
          "resolution",
          "\t\tthis is the resolution of the simulated density map",
          command::Parameter
          (
            "resolution_value",
            "\tchoose the resolution for the density",
            "6.9"
          )
        )
      ),
      m_VoxelSizeFlag
      (
        new command::FlagStatic
        (
          "voxel_size",
          "\t\tthis is the voxel size of the simulated density map",
          command::Parameter
          (
            "voxel_size_value",
            "\tchoose the voxel size for the density",
            "1.5"
          )
        )
      ),
      m_SmoothingKernelFlag
      (
        new command::FlagStatic
        (
          "kernel",
          "\tchange the smoothing kernel",
          command::Parameter
          (
            "kernel",
            "\tchoice of possible smoothing kernel",
            command::ParameterCheckEnumerate< density::Simulators>(),
            density::GetSimulators().e_Gaussian
          )
        )
      ),
      m_NoiseFlag
      (
        new command::FlagStatic
        (
          "noise",
          "add gaussian noise to the simulated map",
          command::Parameter
          (
            "noise_ccc",
            "overall ccc of the simulated map to the noised map",
            command::ParameterCheckRanged< double>( double( 0.0), double( 1.0)),
            "0.5"
          )
        )
      ),
      m_OutputPrefixFlag
      (
        new command::FlagStatic
        (
          "prefix",
          "output prefix",
          command::Parameter
          (
            "output_prefix",
            "path with and prefix for written file like /home/user/run1_",
            ""
          )
        )
      )
    {
    }

    const ApplicationType DensityFromPDB::DensityFromPDB_Instance
    (
      GetAppGroups().AddAppToGroup( new DensityFromPDB(), GetAppGroups().e_Density)
    );

  } // namespace app
} // namespace bcl
