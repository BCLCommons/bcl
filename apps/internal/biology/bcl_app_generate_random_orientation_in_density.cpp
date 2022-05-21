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
#include "biol/bcl_biol_aa_back_bone.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "density/bcl_density_map.h"
#include "density/bcl_density_mask_3d.h"
#include "density/bcl_density_simulators.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_running_average_sd.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "storage/bcl_storage_template_instantiations.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class GenerateRandomOrientation
    //! @brief Class randomly orients protein structures into a simulated density maps
    //! @details Class takes protein structures to randomly orient them into a simulated density maps for a given set
    //! of protein.
    //!
    //! @author woetzen, bitterd
    //! @date 03/24/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API GenerateRandomOrientation :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      //! input pdbs file flag
      util::ShPtr< command::FlagInterface> m_PdbFilenamesFlag;

      //! input pdbs to synthesize density for
      util::ShPtr< command::FlagInterface> m_PdbFilenamesForDensityFlag;

      //! input move object for generating random orientation
      util::ShPtr< command::FlagInterface> m_MutateFlag;

      //! flag for changing the resolution requested for the produced density map
      util::ShPtr< command::FlagInterface> m_ResolutionFlag;

      //! flag for adding noise
      util::ShPtr< command::FlagInterface> m_NoiseFlag;

      //! flag for the number of random orientations to be probed
      util::ShPtr< command::FlagInterface> m_NumberOrientationsFlag;

      //! flag for high res or low resolution mask mode
      util::ShPtr< command::FlagInterface> m_MaskModeFlag;

      //! flag or choice of density simulator
      util::ShPtr< command::FlagInterface> m_SimulatorFlag;

      //! flag for a dynamic list of atoms to be fitted
      util::ShPtr< command::FlagInterface> m_AtomListFlag;

      //! mutate to orient protein structure
      mutable util::ShPtr< math::MutateInterface< assemble::ProteinModel> > m_Mutate;

      //! density simlulator object
      mutable util::ShPtr< density::SimulateInterface> m_SpSimulator;

      //! atom types to be considered
      mutable storage::Set< biol::AtomType> m_AtomTypes;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      GenerateRandomOrientation();

    public:

      //! @brief Clone function
      //! @return pointer to new GenerateRandomOrientation
      GenerateRandomOrientation *Clone() const
      {
        return new GenerateRandomOrientation( *this);
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
        return storage::Vector< std::string>( size_t( 1), "GenerateRandomOrientationInDensity");
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

      //! @brief collect orientation ccc statistics for all models within a simulated map of a given pdb file
      //! @param MODELS the models to randomly orient
      //! @param PDB_MRC_FILENAME filename of pdb to simulate density from
      //! @param STATISTICS DataSetStatisticInterface to add statistic to
      void CorrStatistics( const storage::List< assemble::ProteinModel> &MODELS, const std::string &PDB_MRC_FILENAME, math::RunningAverageSD< double> &STATISTICS) const
      {
        // demsity
        util::ShPtr< density::Map> sp_map_sim( PrepareDensity( PDB_MRC_FILENAME));

        // get the mean intensity and standard deviation of simulated map
        const double map_sim_mean( sp_map_sim->GetMean());
        const double map_sim_sd( sp_map_sim->GetRmsd());

        BCL_MessageStd
        (
          "mean and sd of intensities in map: " + util::Format()( map_sim_mean) + " " + util::Format()( map_sim_sd)
        );

        // statistics object
        math::RunningAverageSD< double> mean_sd_corr_models_all;
        math::RunningAverageSD< double> mean_sd_corr_models_mask;

        // high res mask modus?
        const bool high_res_modus( m_MaskModeFlag->GetFlag());
        const double masking_distance( high_res_modus ? 5.0 : 8.0);
        BCL_MessageStd( "starting to generate random orientations");

        // iterate over all models
        for( storage::List< assemble::ProteinModel>::const_iterator model_itr( MODELS.Begin()), model_itr_end( MODELS.End()); model_itr != model_itr_end; ++model_itr)
        {
          // statistics object
          math::RunningAverageSD< double> mean_sd_corr_current_all;
          math::RunningAverageSD< double> mean_sd_corr_current_mask;

          // for original orientation
          {
            const assemble::ProteinModel &new_orientation( *model_itr);

            // generate density for new orientation
            const density::Map map_orient_all( m_SpSimulator->operator()( high_res_modus ? new_orientation.GetAtoms() : new_orientation.GetAtoms( storage::Set< biol::AtomType>( biol::GetAtomTypes().CA))));

            // calculate ccc for entire map
            const double ccc_all( sp_map_sim->CrossCorrelationCoefficient( map_orient_all, 0.0));
            if( util::IsDefined( ccc_all))
            {
              mean_sd_corr_current_all +=  ccc_all;
              mean_sd_corr_models_all +=  ccc_all;
            }

            // get residues of model
            const util::SiPtrVector< const biol::AABase> residues( new_orientation.GetAminoAcids());
            const size_t nr_residues( residues.GetSize());
            for( size_t res( 0); res < nr_residues; ++res)
            {
              util::SiPtrVector< const linal::Vector3D> coords;
              // get all and neighboring atom coordinates in high res mode
              if( high_res_modus)
              {
                // skip first and last residue
                if( res == 0 || res == ( nr_residues - 1))
                {
                  continue;
                }
                coords.Append( residues( res - 1)->GetAtomCoordinates());
                coords.Append( residues( res    )->GetAtomCoordinates());
                coords.Append( residues( res + 1)->GetAtomCoordinates());
              }
              // get only ca atom in low res mode
              else
              {
                coords.Append( residues( res    )->GetAtomCoordinates( storage::Set< biol::AtomType>( biol::GetAtomTypes().CA)));
              }

              // current mask and ccc over that mask
              const density::Mask3d current_mask( coords, masking_distance, sp_map_sim->GetCellWidth(), sp_map_sim->GetOrigin());
              const double mask_ccc( current_mask.CrossCorrelationCoefficient( *sp_map_sim, map_orient_all));
              if( util::IsDefined( mask_ccc))
              {
                STATISTICS += mask_ccc;
                mean_sd_corr_current_mask +=  mask_ccc;
                mean_sd_corr_models_mask +=  mask_ccc;
              }
            }
          }
          BCL_MessageStd( "unchanged density and model over all  ccc: " + util::Format()( mean_sd_corr_current_all.GetAverage()));
          BCL_MessageStd( "unchanged density and model over mask ccc: " + util::Format()( mean_sd_corr_current_mask.GetAverage()));

          // randomly orient protein in density
          for( size_t round( 0), max_nr( m_NumberOrientationsFlag->GetFirstParameter()->GetNumericalValue< size_t>()); round < max_nr; ++round)
          {
            // print message every 10%
            if( round % int( ceil( max_nr / 10.0)) == 0)
            {
              BCL_MessageStd( "finished " + util::Format().W( 6).FFP( 2)( 100.0 * double( round) / max_nr) + "%");
            }

            // generate new protein model
            const math::MutateResult< assemble::ProteinModel> new_orientation_result( m_Mutate->operator()( *model_itr));
            if( !new_orientation_result.GetArgument().IsDefined())
            {
              BCL_MessageCrt( "unsuccessful generation of new orientation in round: " + util::Format()( round));
              continue;
            }
            const assemble::ProteinModel &new_orientation( *new_orientation_result.GetArgument());

            // generate density for new orientation
            const density::Map map_orient_all( m_SpSimulator->operator()( high_res_modus ? new_orientation.GetAtoms() : new_orientation.GetAtoms( storage::Set< biol::AtomType>( biol::GetAtomTypes().CA))));

            // calculate ccc for entire map
            const double ccc_all( sp_map_sim->CrossCorrelationCoefficient( map_orient_all, -std::numeric_limits< double>::max()));
            if( util::IsDefined( ccc_all))
            {
              mean_sd_corr_models_all +=  ccc_all;
              mean_sd_corr_current_all +=  ccc_all;
            }

            // get residues of model
            const util::SiPtrVector< const biol::AABase> residues( new_orientation.GetAminoAcids());
            const size_t nr_residues( residues.GetSize());
            for( size_t res( 0); res < nr_residues; ++res)
            {
              util::SiPtrVector< const linal::Vector3D> coords;
              // get all and neighboring atom coordinates in high res mode
              if( high_res_modus)
              {
                // skip first and last residue
                if( res == 0 || res == ( nr_residues - 1))
                {
                  continue;
                }
                coords.Append( residues( res - 1)->GetAtomCoordinates());
                coords.Append( residues( res    )->GetAtomCoordinates());
                coords.Append( residues( res + 1)->GetAtomCoordinates());
              }
              // get only ca atom in low res mode
              else
              {
                coords.Append( residues( res    )->GetAtomCoordinates( m_AtomTypes));
              }

              // current mask and ccc over that mask
              const density::Mask3d current_mask( coords, masking_distance, sp_map_sim->GetCellWidth(), sp_map_sim->GetOrigin());
              const double mask_ccc( current_mask.CrossCorrelationCoefficient( *sp_map_sim, map_orient_all));
              if( util::IsDefined( mask_ccc))
              {
                STATISTICS += mask_ccc;
                mean_sd_corr_current_mask +=  mask_ccc;
                mean_sd_corr_models_mask +=  mask_ccc;
              }
            }
            BCL_MessageVrb
            (
              util::Format().W( 5).Fill( '0').R()( round) + " density and model over all ccc mean sd: " +
              util::Format()( mean_sd_corr_current_all.GetAverage()) + " " + util::Format()( mean_sd_corr_current_all.GetStandardDeviation())
            );
            BCL_MessageVrb
            (
              util::Format().W( 5).Fill( '0').R()( round) + " density and model over mask ccc mean sd: " +
              util::Format()( mean_sd_corr_current_mask.GetAverage()) + " " + util::Format()( mean_sd_corr_current_mask.GetStandardDeviation())
            );
          }

          BCL_MessageTop
          (
            "current mean and sd for over all  ccc: " + util::Format()( mean_sd_corr_current_all.GetAverage()) + " " +
            util::Format()( mean_sd_corr_current_all.GetStandardDeviation())
          );
        }
        BCL_MessageTop
        (
          "models mean and sd for over all  ccc: " + util::Format()( mean_sd_corr_models_all.GetAverage()) + " " +
          util::Format()( mean_sd_corr_models_all.GetStandardDeviation())
        );
      }

      //! @brief prepare a density map for a given pdb filename
      //! @param PDB_MRC_FILENAME filename of pdb to simulate density from
      //! @return ShPtr to simulated density map
      util::ShPtr< density::Map> PrepareDensity( const std::string &PDB_MRC_FILENAME) const
      {
        io::IFStream read;
        BCL_MessageStd( "read pdbfile for density: " + PDB_MRC_FILENAME);
        io::File::MustOpenIFStream( read, PDB_MRC_FILENAME);
        pdb::Handler pdb( read);
        io::File::CloseClearFStream( read);

        // factory
        const pdb::Factory factory;

        // build models from pdb and translate center to origin
        assemble::ProteinModel density_model( factory.ProteinModelFromPDB( pdb));
        density_model.Translate( -density_model.GetCenterOfMass());

        // simulate density
        util::ShPtr< density::Map> sp_map_sim( new density::Map( m_SpSimulator->operator()( density_model.GetAtoms())));

        // add noise
        if( m_NoiseFlag->GetFlag())
        {
          BCL_MessageStd( "adding noise to start density");
          const double target_ccc( m_NoiseFlag->GetFirstParameter()->GetNumericalValue< double>());
          double noise_ccc( sp_map_sim->CrossCorrelationCoefficient( *sp_map_sim, 0.0));
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
            noised_density = *sp_map_sim;
            noised_density.AddNoise( random::GetGlobalRandom(), double( 0.0), double( i) * sp_map_sim->GetRmsd() / double( 20));
            noise_ccc = noised_density.CrossCorrelationCoefficient( *sp_map_sim, sp_map_sim->GetMinimum());
          }

          if( noise_ccc > target_ccc)
          {
            BCL_Exit
            (
              "unable to generate noised map with CCC to original map of: " + util::Format()( target_ccc) +
              " after 100 iterations",
              -1
            );
          }
          else
          {
            *sp_map_sim = noised_density;
          }
        }

        // end
        return sp_map_sim;
      }

    private:

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      static const ApplicationType GenerateRandomOrientation_Instance;

    }; // class GenerateRandomOrientation

    //! @brief the Main function
    //! @return error code - 0 for success
    int GenerateRandomOrientation::Main() const
    {
      // create storage::Set of BBAtomTypes to be matched - need to extract them from protein model
      m_AtomTypes = m_AtomListFlag->GetObjectSet< biol::AtomType>();

      const double resolution( m_ResolutionFlag->GetFirstParameter()->GetNumericalValue< double>());
      const linal::Vector3D grid_spacing( resolution / 3);

      // factory
      const pdb::Factory factory;

      // read pdbfiles
      io::IFStream read;
      storage::List< assemble::ProteinModel> models;
      {
        // list of pdb files
        const storage::Vector< std::string> pdb_filenames( m_PdbFilenamesFlag->GetObjectList< std::string>());
        for( storage::Vector< std::string>::const_iterator pdb_itr( pdb_filenames.Begin()), pdb_itr_end( pdb_filenames.End()); pdb_itr != pdb_itr_end; ++pdb_itr)
        {
          BCL_MessageStd( "read pdbfile to orient: " + util::Format()( *pdb_itr));
          io::File::MustOpenIFStream( read, *pdb_itr);
          pdb::Handler pdb( read);
          io::File::CloseClearFStream( read);
          // translate center of mass to origin
          assemble::ProteinModel model( factory.ProteinModelFromPDB( pdb));
          model.Translate( -model.GetCenterOfMass());
          models.PushBack( model);
        }
      }

      // read mutate object
      io::File::MustOpenIFStream( read, m_MutateFlag->GetFirstParameter()->GetValue());
      read >> m_Mutate;
      io::File::CloseClearFStream( read);

      // density simulator
      m_SpSimulator =
        density::GetSimulators().CreateSimulator
        (
          density::Simulator( m_SimulatorFlag->GetFirstParameter()->GetValue()),
          grid_spacing,
          resolution
        );

      // list of pdb mrc files
      const storage::Vector< std::string> pdb_mrc_filenames( m_PdbFilenamesForDensityFlag->GetObjectList< std::string>());

      // cross correlation mean and sd
      math::RunningAverageSD< double> mean_sd_corr_mask;

      // iterate over mrc files
      for( storage::Vector< std::string>::const_iterator mrc_itr( pdb_mrc_filenames.Begin()), mrc_itr_end( pdb_mrc_filenames.End()); mrc_itr != mrc_itr_end; ++mrc_itr)
      {
        CorrStatistics( models, *mrc_itr, mean_sd_corr_mask);
      }

      // write mean and sd
      BCL_MessageTop
      (
        "final mean and sd over mask ccc: " + util::Format()( mean_sd_corr_mask.GetAverage()) + " " +
        util::Format()( mean_sd_corr_mask.GetStandardDeviation())
      );

      return 0;
    } // end Main

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> GenerateRandomOrientation::InitializeCommand() const
    {
      // initialize ShPtr to a new command object
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // input pdb file parameter
      sp_cmd->AddFlag( m_PdbFilenamesFlag);

      // input pdb for synthesized density
      sp_cmd->AddFlag( m_PdbFilenamesForDensityFlag);

      // input move object
      sp_cmd->AddFlag( m_MutateFlag);

      // flag for changing the resolution requested for the produced density map
      sp_cmd->AddFlag( m_ResolutionFlag);

      // flag for adding noise
      sp_cmd->AddFlag( m_NoiseFlag);

      // flag for the number of orientations
      sp_cmd->AddFlag( m_NumberOrientationsFlag);

      // flag for switching between high or low resolution mask modus
      sp_cmd->AddFlag( m_MaskModeFlag);

      // flag for adjusting the used density simulator
      sp_cmd->AddFlag( m_SimulatorFlag);

      // flag for a dynamic list of atoms to be used
      sp_cmd->AddFlag( m_AtomListFlag);

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
    GenerateRandomOrientation::GenerateRandomOrientation() :
      m_PdbFilenamesFlag
      (
        new command::FlagDynamic
        (
          "pdb_filenames_orient",
          "\tfilenames for input pdbs to be oriented in mrc electron density maps",
          command::Parameter
          (
            "pdb_filename",
            "file name of pdb",
            command::ParameterCheckFileExistence()
          ),
          1
        )
      ),
      m_PdbFilenamesForDensityFlag
      (
        new command::FlagDynamic
        (
          "pdb_filenames_densities",
          "\tfilenames for pdbs which are used for synthesizing density maps",
          command::Parameter
          (
            "pdb_filename",
            "file name for pdb to synthesize map for",
            command::ParameterCheckFileExistence()
          ),
          1
        )
      ),
      m_MutateFlag
      (
        new command::FlagStatic
        (
          "mutate",
          "\t\tuse this flag to pass a custom mutate object",
          command::Parameter
          (
            "mutate_filename",
            "\tname of file, containing the mutate bcl object"
          )
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
      m_NumberOrientationsFlag
      (
        new command::FlagStatic
        (
          "orientations",
          "number_orientations",
          command::Parameter
          (
            "nr_orientations",
            "that many orientations will be generated",
            command::ParameterCheckRanged< size_t>( 1, 10000),
            "1000"
          )
        )
      ),
      m_MaskModeFlag
      (
        new command::FlagStatic
        (
          "mask_high_res",
          "use high resolution modus with 5A cutoff distance using all side chain atoms, vs the default low res modus"
          " with 8A cutoff using only CA atoms"
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
      m_AtomListFlag
      (
        new command::FlagDynamic
          (
          "atoms",
          "\t\tthis is list of atoms used for the fitting",
          command::Parameter
          (
            "backboneatom",
            "any backbone atom from the list",
            command::ParameterCheckAllowed
            (
              storage::Vector< std::string>( biol::AABackBone().GetTypesOfAtoms().Begin(), biol::AABackBone().GetTypesOfAtoms().End())
            )
          ),
          1,
          biol::AABackBone().GetTypesOfAtoms().GetSize()
        )
      )
    {
    }

    const ApplicationType GenerateRandomOrientation::GenerateRandomOrientation_Instance
    (
      GetAppGroups().AddAppToGroup( new GenerateRandomOrientation(), GetAppGroups().e_Density)
    );

  } // namespace app
} // namespace bcl
