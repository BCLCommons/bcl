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

// App header
#include "app/bcl_app_apps.h"

// include header for this application
#include "bcl_app_pharm_map.h"

// include headers from the bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "chemistry/bcl_chemistry_fragment_split_interface.h"
#include "chemistry/bcl_chemistry_molecule_feature_mapper.h"
#include "chemistry/bcl_chemistry_pharmacophore_mapper.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "descriptor/bcl_descriptor_prediction.h"
#include "descriptor/bcl_descriptor_structure_search.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_running_average_sd.h"
#include "math/bcl_math_template_instantiations.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_thunk_job.h"
#include "sdf/bcl_sdf_fragment_factory.h"
#include "sdf/bcl_sdf_mdl_handler.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_stopwatch.h"
#include "util/bcl_util_string_numeric_conversion.h"

using bcl::command::ParameterCheckFileExistence;

namespace bcl
{
  namespace app
  {

    // Strings for specifying methods
    static const storage::Vector< std::string> g_PharmMapMethodFlags = storage::Vector< std::string>::Create( "Analyze", "Score", "FeatureMap");

    namespace
    {
      struct Worker
      {
        util::SiPtr< const storage::Vector< chemistry::FragmentComplete> > m_Ensemble;
        size_t m_NumberMols;
        size_t m_ThreadNumber;
        size_t m_Begin;
        size_t m_End;
        chemistry::FragmentEnsemble m_LocalEnsemble;
        chemistry::ConformationGraphConverter m_GraphMaker;
        void RunThread();
      };
    }

    //! @brief initializes the command object for this application
    //! @return a ShPtr to a Command containing all of this applications parameters
    util::ShPtr< command::Command> PharmMap::InitializeCommand() const
    {
      // initialize command to be returned
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // insert all the flags and params

      // flags for input/output
      sp_cmd->AddFlag( m_MethodFlag);
      sp_cmd->AddFlag( m_ScaffoldFlag);
      sp_cmd->AddFlag( m_FragmentsFlag);
      sp_cmd->AddFlag( m_SampleSizeFlag);
      sp_cmd->AddFlag( m_DerivativesFlag);
      sp_cmd->AddFlag( m_PymolOutputFilenameFlag);
      sp_cmd->AddFlag( m_PropertiesToMapFlag);
      sp_cmd->AddFlag( m_ScorerFlag);
      sp_cmd->AddFlag( m_MoleculeOutputFlag);
      sp_cmd->AddFlag( m_RemoveZerodPropertiesFlag);
      sp_cmd->AddFlag( m_DetailsFilenameFlag);
      sp_cmd->AddFlag( m_PharmMapFileFlag);
      sp_cmd->AddFlag( m_BinaryDifferenceFlag);
      sp_cmd->AddFlag( m_ScoreToleranceFlag);
      sp_cmd->AddFlag( m_PercentCoverageForScaffoldFlag);

      // default flags
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    } // InitializeCommand

    //! @brief the Main function
    //! @return 0 for success
    int PharmMap::Main() const
    {
      // Run options
      bool analyze( true);
      bool pharmmap( true);
      bool feature_map( false);

      const std::string &method( m_MethodFlag->GetFirstParameter()->GetValue());
      if( method == g_PharmMapMethodFlags( 1))
      {
        analyze = false;
      }
      else if( method == g_PharmMapMethodFlags( 2))
      {
        pharmmap = false;
        feature_map = true;
      }

      // Get the scoring descriptor - needed no matter what we do
      util::ObjectDataLabel scorer_label( m_ScorerFlag->GetFirstParameter()->GetValue());

      // TODO just wrap things in a combine if they aren't already
      if( scorer_label.GetValue() != "Combine")
      {
        BCL_MessageStd( "May be unable to read activity descriptor.  If you experience problems, wrap it in a Combine() descriptor, e.g. Combine(PredictionMean(...))");
      }

      if( pharmmap)
      {
        //
        // User input
        //

        // Retrieve the scaffold
        io::IFStream input;
        io::File::MustOpenIFStream( input, m_ScaffoldFlag->GetFirstParameter()->GetValue());
        chemistry::FragmentComplete scaffold( sdf::FragmentFactory::MakeFragment( input, sdf::e_Maintain));
        io::File::CloseClearFStream( input);
        BCL_Assert( scaffold.GetNumberAtoms(), "Scaffold must have at least one atom in it");

        // Get the properties to map
        io::File::MustOpenIFStream( input, m_PropertiesToMapFlag->GetFirstParameter()->GetValue());
        util::ObjectDataLabel properties_label( input);
        io::File::CloseClearFStream( input);

        BCL_Assert
        (
          properties_label.GetValue() == "Combine",
          "Unable to read the desired properties.  Wrap them in a Combine() descriptor, e.g. Combine(Atom_Vcharge,Atom_PiEN)"
        );

        //
        // Set up the pharmacophore mapper
        //

        chemistry::PharmacophoreMapper mapper
        (
          scorer_label,
          properties_label
        );

        //
        // Pharmacophore generation/derivative analysis
        //

        // Make sure a molecules file was given
        BCL_Assert
        (
          m_DerivativesFlag->GetFlag() && !m_DerivativesFlag->GetParameterList().IsEmpty(),
          "You must give a file containing some molecules"
        );

        // Input derivatives
        chemistry::FragmentEnsemble derivs;
        BCL_MessageStd( "Reading molecules...");
        chemistry::FragmentFeed derivs_feed( m_DerivativesFlag->GetStringList(), sdf::e_Saturate);
        size_t n_mols( 0);
        for( ; derivs_feed.NotAtEnd(); ++derivs_feed)
        {
          derivs.PushBack( *derivs_feed);
          ++n_mols;

          // Write out a log status every 1000 molecules
          if( !( n_mols % 1000))
          {
            util::GetLogger().LogStatus( "Read " + util::Format()( n_mols) + " molecules");
          }
        }

        // Generate a pharmmap from the inputs
        if( false)
        {
        }
        else if( pharmmap)
        {

          // Calculate a pharm map
          /*
          mapper.CalculateMaps
          (
            scaffold,
            derivs,
            m_RemoveZerodPropertiesFlag->GetFlag(),
            m_ScoreToleranceFlag->GetFlag() ? m_ScoreToleranceFlag->GetFirstParameter()->GetNumericalValue< float>() : 0.0,
            m_BinaryDifferenceFlag->GetFlag()
          );
          */
          size_t tmp( 0);
          mapper.CalculateMaps
          (
            scaffold,
            derivs,
            tmp
          );

          std::string details_file;
          if( m_DetailsFilenameFlag->GetFlag())
          {
            details_file = m_DetailsFilenameFlag->GetFirstParameter()->GetValue();
          }

          // Write details to a file
          if( details_file.length())
          {
            mapper.WriteDetailsCSV( details_file);
            mapper.WriteDetailsMols( details_file);
          }

          // Write pymol script to display results
          mapper.WriteMapsPML
          (
            m_ScaffoldFlag->GetFirstParameter()->GetValue(),
            m_PymolOutputFilenameFlag->GetFirstParameter()->GetValue()
          );
        }
        else
        {
          //
          // Score molecules using the pharmmap
          //

          BCL_Assert
          (
            m_PharmMapFileFlag->GetFlag() && !m_PharmMapFileFlag->GetParameterList().IsEmpty(),
            "Must provide a pharmacophore map file (csv file that matches the given properties file)"
          );

          // Read in the coefficient map
          storage::Map< size_t, linal::Vector< float> > coeffs_map
          (
            ParsePharmMapCSV( m_PharmMapFileFlag->GetFirstParameter()->GetValue())
          );

          // Call the mapper
          storage::Map< size_t, storage::Vector< storage::Pair< size_t, linal::Vector< float> > > > results
          (
            mapper.ScoreMolecules( scaffold, derivs, coeffs_map)
          );

          BCL_Assert
          (
            m_DetailsFilenameFlag->GetFlag() && !m_DetailsFilenameFlag->GetParameterList().IsEmpty(),
            "Must specify a details output"
          );

          std::string details_file( m_DetailsFilenameFlag->GetFirstParameter()->GetValue());

          io::OFStream out;
          io::File::MustOpenOFStream( out, details_file + "_molecule_comparisons.csv");

          const storage::Vector< std::string> &prop_names( mapper.GetPropertyStrings());

          out << "growpt,mol_1,mol_2,";
          for( size_t i( 0); i < prop_names.GetSize(); ++i)
          {
            out << prop_names( i) << ",";
          }
          out << "pharmmap_score_diff,scorer_score_diff" << std::endl;

          for
          (
            storage::Map< size_t, storage::Vector< storage::Pair< size_t, linal::Vector< float> > > >::const_iterator
              itr_map( results.Begin()), itr_map_end( results.End());
            itr_map != itr_map_end;
            ++itr_map
          )
          {
            const size_t &first_mol_no( itr_map->first);
            const storage::Vector< storage::Pair< size_t, linal::Vector< float> > > &comp_vector( itr_map->second);
            for( size_t i( 0); i < comp_vector.GetSize(); ++i)
            {
              // Number of coefficients
              const size_t n_cols( comp_vector( i).Second().GetSize());

              // there should be exactly this many columns; one for grow point, one for pharmmap score diff, one for scorer score diff
              if( n_cols < prop_names.GetSize() + 3)
              {
                continue;
              }

              // first column is the grow point
              const size_t &growpt( comp_vector( i).Second()( 0));
              const size_t &second_mol_no( comp_vector( i).First());

              // intermediate columns are property differences
              out << growpt << "," << first_mol_no << "," << second_mol_no << ",";
              for( size_t j( 0); j < prop_names.GetSize(); ++j)
              {
                out << comp_vector( i).Second()( j + 1) << ",";
              }
              // last two columns are score differences
              out << comp_vector( i).Second()( n_cols - 2) << "," << comp_vector( i).Second()( n_cols - 1) << std::endl;
            }
          }

          io::File::CloseClearFStream( out);
          if( !analyze)
          {
          }
        }
      }
      else if( feature_map)
      {
        //
        // Get molecules
        //
        BCL_Assert
        (
          m_DerivativesFlag->GetFlag() && !m_DerivativesFlag->GetParameterList().IsEmpty(),
          "You must provide a molecules file to -" + m_DerivativesFlag->GetName()
        );

        chemistry::MoleculeFeatureMapper mapper( scorer_label);

        math::RunningAverageSD< float> natom_avg;
        chemistry::FragmentEnsemble mols;
        chemistry::FragmentFeed feed( m_DerivativesFlag->GetStringList(), sdf::e_Remove);
        for( ; feed.NotAtEnd(); ++feed)
        {
          if( feed->GetNumberAtoms() > 0)
          {
            natom_avg += feed->GetNumberAtoms();
            mols.PushBack( *feed);
          }
        }

        size_t n_mols( mols.GetSize());
/*
        util::ObjectDataLabel label( "Constant(0)");
        chemistry::PharmacophoreMapper pm( label, label);

        size_t min_scaff_size
        (
          std::max
          (
            size_t( natom_avg.GetAverage() * m_PercentCoverageForScaffoldFlag->GetFirstParameter()->GetNumericalValue< float>()),
            size_t( 1)
          )
        );
        BCL_MessageStd( "Finding scaffolds that have at least " + util::Format()( min_scaff_size) + " atoms");

        storage::Vector< storage::Pair< graph::ConstGraph< size_t, size_t>, storage::Vector< size_t> > > scaffolds
        (
          // on average find scaffolds that are at least 35% of the molecule size
          pm.FindScaffolds( mols, 0.025, min_scaff_size)
        );

        for( size_t i( 0); i < scaffolds.GetSize(); ++i)
        {
          BCL_MessageStd( "Scaffold " + util::Format()( i) + ": " + util::Format()( scaffolds( i).First().GetSize()) + " atoms, " + util::Format()( scaffolds( i).Second().GetSize()) + " matches");
        }
        */
        // Must saturate here instead of on read to make sure all hydrogens are at the end
        // this is important during the CSI searching below
        mols.SaturateWithH();

        io::OFStream out;

        //
        // individual molecules
        //

        std::vector< std::vector< std::map< size_t, float> > > maps;
        maps.reserve( n_mols);
        size_t cur_mol( 1);
        for
        (
          chemistry::FragmentEnsemble::const_iterator itr_mol( mols.Begin()), itr_mol_end( mols.End());
          itr_mol != itr_mol_end;
          ++itr_mol, ++cur_mol
        )
        {
          util::GetLogger().LogStatus
          (
            "Feature mapping molecule " + util::Format()( cur_mol) + " / " + util::Format()( n_mols) + " (" + util::Format().FFP( 2)( float( cur_mol) / float( n_mols) * 100) + " %)"
          );
          maps.push_back( mapper.Perturb( *itr_mol, false, false, storage::Vector< size_t>()));
        }

        std::string molecules_full_filename
        (
          io::File::MakeAbsolutePath( m_PymolOutputFilenameFlag->GetFirstParameter()->GetValue() + "_molecules.sdf")
        );
        io::File::MustOpenOFStream( out, molecules_full_filename);
        mols.WriteMDL( out);
        io::File::CloseClearFStream( out);

        std::string mol_pml_filename( m_PymolOutputFilenameFlag->GetFirstParameter()->GetValue() + ".pml");
        io::File::MustOpenOFStream( out, mol_pml_filename);
        out << "load " << molecules_full_filename << ", molecule" << std::endl;
        out << "alter molecule and state 1, b=0" << std::endl;
        out << "alter molecule and state 1, vdw=0" << std::endl;
        out << "set label_size, -0.4" << std::endl;
        out << "hide lines " << std::endl;
        out << "show sticks" << std::endl;
        out << "set stick_quality, 20" << std::endl;
        out << "set valence, 1" << std::endl;

        float min( 0), max( 0);
        for( size_t i( 0); i < n_mols; ++i)
        {
          if( maps[ i].empty())
          {
            continue;
          }
          std::map< size_t, float> &map( maps[ i][ 0]);
          size_t mol_no( i + 1);
          for
          (
            std::map< size_t, float>::const_iterator itr_map( map.begin()), itr_map_end( map.end());
            itr_map != itr_map_end;
            ++itr_map
          )
          {
            if( itr_map->second < min)
            {
              min = itr_map->second;
            }
            if( itr_map->second > max)
            {
              max = itr_map->second;
            }
            out << "alter \"molecule\" and state " << mol_no << " and id " << itr_map->first + 1;
            out << ", b=" << itr_map->second << std::endl;
          }
        }
        float range( std::max( std::fabs( min), std::fabs( max)));
        out << "spectrum b, blue_white_red, selection=molecule, ";
        out << "minimum=" << ( -1 * range) << ", maximum=" << range << std::endl;

        io::File::CloseClearFStream( out);

        mols.RemoveH();

        storage::Vector< chemistry::FragmentComplete> scaffold_mols;
        if( !m_ScaffoldFlag->GetFlag())
        {
          storage::Vector< chemistry::FragmentComplete> scaffold_mols = storage::Vector< chemistry::FragmentComplete>( mols.Begin(), mols.End());
          size_t n_scaffold_mols( scaffold_mols.GetSize());

          //
          // Find common scaffold
          //
          size_t n_threads( std::min( sched::GetNumberCPUs(), n_scaffold_mols / 2));
          std::vector< Worker> threads( n_threads);
          for( size_t i( 0); i < n_threads; ++i)
          {
            threads[ i].m_GraphMaker = chemistry::ConformationGraphConverter
                                       (
                                         //chemistry::ConformationGraphConverter::e_ElementType,
                                         chemistry::ConformationGraphConverter::e_ElementType,
                                         chemistry::ConfigurationalBondTypeData::e_BondOrderOrAromatic
                                       );
          }

          BCL_MessageStd( "Finding common substructures between molecules");
          while( ( n_scaffold_mols = scaffold_mols.GetSize()) > 1)
          {

            // Find the common scaffold amongst all molecules
            // TODO
            // set up threads

            n_threads = std::max< size_t>( std::min< size_t>( n_threads, n_scaffold_mols / 2), 1);

            BCL_MessageStd( "  Comparing " + util::Format()( n_scaffold_mols) + " molecules");

            size_t n_per_thread( n_scaffold_mols / n_threads);
            size_t overflow( n_scaffold_mols % n_threads);

            for( size_t i( 0); i < n_threads; ++i)
            {
              threads[ i].m_Ensemble = util::SiPtr< const storage::Vector< chemistry::FragmentComplete> >( &scaffold_mols);
              threads[ i].m_ThreadNumber = i;
              threads[ i].m_NumberMols = n_scaffold_mols;
              threads[ i].m_LocalEnsemble = chemistry::FragmentEnsemble();
              threads[ i].m_Begin = std::min( n_per_thread * i, n_scaffold_mols);
              threads[ i].m_End = std::min( threads[ i].m_Begin + n_per_thread, n_scaffold_mols);
            }
            threads[ n_threads - 1].m_End += overflow; // last thread gets the extras

            util::ShPtrVector< sched::JobInterface> jobs;
            jobs.AllocateMemory( n_threads);
            const size_t group( 1);

            for( size_t p( 0); p < n_threads; ++p)
            {
              jobs.PushBack
              (
                util::ShPtr< sched::JobInterface>
                (
                  new sched::ThunkJob< Worker, void>
                  (
                    group,
                    threads[ p],
                    &Worker::RunThread,
                    sched::JobInterface::e_READY,
                    NULL
                  )
                )
              );

              // submit the job
              sched::GetScheduler().RunJob( jobs.LastElement());
            }

            // wait for jobs to finish
            for( size_t p( 0); p < n_threads; ++p)
            {
              sched::GetScheduler().Join( jobs( p));
            }

            scaffold_mols = storage::Vector< chemistry::FragmentComplete>();
            for( size_t p( 0); p < n_threads; ++p)
            {
              scaffold_mols.InsertElements
              (
                scaffold_mols.End(), threads[ p].m_LocalEnsemble.Begin(), threads[ p].m_LocalEnsemble.End()
              );
            }
          }
        }
        else
        {
          io::IFStream in;
          io::File::MustOpenIFStream( in, m_ScaffoldFlag->GetFirstParameter()->GetValue());
          chemistry::FragmentEnsemble scaffs( in);
          scaffold_mols = storage::Vector< chemistry::FragmentComplete>( scaffs.Begin(), scaffs.End());
        }

        if( !scaffold_mols.IsEmpty())
        {
          for( size_t scaff_no( 0), end_no( scaffold_mols.GetSize()); scaff_no < end_no; ++scaff_no)
          {
            BCL_MessageStd( "Found a common scaffold between  molecules, mapping average features over this scaffold");
            chemistry::FragmentComplete &scaffold( scaffold_mols( scaff_no));

            // compare each molecule to each other molecule
            storage::Vector< storage::Triplet< size_t, size_t, float> > rmsd_vector( mapper.GetFeatureRMSDOnScaffold( scaffold, mols, maps));

            io::OFStream out( std::string( m_PymolOutputFilenameFlag->GetFirstParameter()->GetValue() + "_feature_rmsd.csv").c_str(), std::ios::out);
            for( size_t i( 0); i < rmsd_vector.GetSize(); ++i)
            {
              out << rmsd_vector( i).First() << "," << rmsd_vector( i).Second() << "," << rmsd_vector( i).Third() << std::endl;
            }
            io::File::CloseClearFStream( out);

            // Need to add hydrogens for the feature mapping
            std::vector< std::map< size_t, float> > scaff_features( mapper.AverageFeatures( scaffold, mols, maps));

            if( scaffold.GetNumberAtoms() > 0 && !scaff_features.empty())
            {
              std::string scaff_filename( m_PymolOutputFilenameFlag->GetFirstParameter()->GetValue() + "_scaffold" + util::Format()( scaff_no) + ".sdf");
              io::File::MustOpenOFStream( out, scaff_filename);
              scaffold.WriteMDL( out);
              io::File::CloseClearFStream( out);

              std::string scaff_full_filename
              (
                io::File::MakeAbsolutePath( scaff_filename)
              );

              std::map< size_t, float> &map( scaff_features[ 0]);

              std::string scaff_pymol_filename( m_PymolOutputFilenameFlag->GetFirstParameter()->GetValue() + "_scaffold" + util::Format()( scaff_no) + ".pml");
              io::File::MustOpenOFStream( out, scaff_pymol_filename);
              out << "load " << scaff_full_filename << ", molecule" << std::endl;
              out << "alter molecule and state 1, b=0" << std::endl;
              out << "alter molecule and state 1, vdw=0" << std::endl;
              out << "set label_size, -0.4" << std::endl;
              out << "hide lines " << std::endl;
              out << "show sticks" << std::endl;
              out << "set stick_quality, 20" << std::endl;
              out << "set valence, 1" << std::endl;

              for
              (
                std::map< size_t, float>::const_iterator itr_map( map.begin()), itr_map_end( map.end());
                itr_map != itr_map_end;
                ++itr_map
              )
              {
                out << "alter \"molecule\" and state 1 and id " << itr_map->first + 1;
                out << ", b=" << itr_map->second << std::endl;
              }
              out << "spectrum b, blue_white_red, selection=molecule and state 1, ";
              out << "minimum=" << ( -1 * range) << ", maximum=" << range << std::endl;
              io::File::CloseClearFStream( out);
            }
          }
        }
        else
        {
          BCL_MessageStd( "Could not find a common scaffold between molecules");
        }

      }
      else
      {
        BCL_MessageStd( "Nothing to do...");
      }

      return 0;

    } // Main

    //! @brief helper function to parse a pharmmap coefficient file
    //! @details reads each line, skipping comments and string lines, and extracts coefficients by commas
    //!          the first field of each line must be the grow point that the coefficients correspond to
    //! @return a map from grow point (key) to a vector of property coefficients( values)
    storage::Map< size_t, linal::Vector< float> > PharmMap::ParsePharmMapCSV( const std::string &FILENAME) const
    {
      storage::Map< size_t, linal::Vector< float> > coeff_map;

      if( !FILENAME.length())
      {
        return coeff_map;
      }

      io::OFStream nullout;

      io::IFStream input;
      io::File::MustOpenIFStream( input, FILENAME);

      // Read each line of the file
      std::string line;
      while( std::getline( input, line))
      {
        // Skip comments and blank lines
        if( !line.length() || line[ 0] == '#')
        {
          continue;
        }

        storage::Vector< float> coeffs;

        // Split the line on commas
        std::stringstream l( line);
        std::string field;

        // Read fields as floats
        bool skip_line( false);
        while( std::getline( l, field, ','))
        {
          float val( 0);

          // Can't convert something on this line to a float, skip the line
          if( !util::TryConvertFromString( val, field, nullout))
          {
            skip_line = true;
            break;
          }
          coeffs.PushBack( val);
        }

        if( skip_line || coeffs.IsEmpty())
        {
          continue;
        }

        // Copy relevant data to the coefficients map
        size_t growpt( coeffs( 0));
        linal::Vector< float> &coeff( coeff_map[ growpt]);
        coeff = linal::Vector< float>( coeffs.GetSize() - 1);
        std::copy( coeffs.Begin() + 1, coeffs.End(), coeff.Begin());
      }

      io::File::CloseClearFStream( input);

      return coeff_map;
    }

    //! @brief standard constructor
    PharmMap::PharmMap() :
      m_MethodFlag
      (
        new command::FlagStatic
        (
          "method",
          "specifies whether pharmacophore map building or scoring is being done",
          command::Parameter
          (
            "method_type",
            "the method to use",
            command::ParameterCheckAllowed
            (
              g_PharmMapMethodFlags
            )
          )
        )
      ),
      m_ScaffoldFlag
      (
        new command::FlagDynamic
        (
          "scaffold",
          "filename containing the scaffold/base fragment",
          command::Parameter
          (
            "filename",
            "the name of an SDF file containing the scaffold",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_FragmentsFlag
      (
        new command::FlagDynamic
        (
          "fragments", "filename for (preferably small) fragments to add to the scaffold",
          command::Parameter
          (
            "filename", "file containing the fragments that will be added to grow points of the scaffold",
            command::ParameterCheckFileExistence()
          ),
          0,
          21
        )
      ),
      m_SampleSizeFlag
      (
        new command::FlagStatic
        (
          "sample_size", "the number of molecules to generate through Monte-Carlo fragment addition",
          command::Parameter
          (
            "number",
            "the number of molecules to generate",
            "50"
          )
        )
      ),
      m_DerivativesFlag
      (
        new command::FlagDynamic
        (
          "molecules", "filename for molecules which are derived from the scaffold that should be analyzed",
          command::Parameter
          (
            "filename", "file containing derivative molecules",
            command::ParameterCheckFileExistence()
          ),
          0,
          1
        )
      ),
      m_PymolOutputFilenameFlag
      (
        new command::FlagStatic
        (
          "output_prefix", "where to output the pymol script",
          command::Parameter
          (
            "output", "filename for where to write a pymol script illustrating the pharmacophore map"
          )
        )
      ),
      m_PropertiesToMapFlag
      (
        new command::FlagDynamic
        (
          "properties", "a file containing the properties to correlate with activity",
          command::Parameter
          (
            "properties",
            "filename containing a descriptor label (e.g. Combine(Atom_VCharge,...) )",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_ScorerFlag
      (
        new command::FlagStatic
        (
          "score",
          "how to compute the activity or score of a molecule",
          command::Parameter
          (
            "descriptor",
            "the descriptor to use to calculate activities (probably should be PredictionMean())",
            command::ParameterCheckSerializable
            (
              descriptor::Combine< chemistry::AtomConformationalInterface, float>()
            )
          )
        )
      ),
      m_MoleculeOutputFlag
      (
        new command::FlagDynamic
        (
          "output_molecules", "filename for molecules that are generated from the MC addition (if -fragments has been specified)",
          command::Parameter
          (
            "filename", "file to output generated molecules to"
          )
        )
      ),
      m_RemoveZerodPropertiesFlag
      (
        new command::FlagStatic
        (
          "skip_unchanged_props", "whether to exclude (set coefficients to zero) for all-zero parameters"
        )
      ),
      m_DetailsFilenameFlag
      (
        new command::FlagDynamic
        (
          "details", "where to write detailed information to",
          command::Parameter
          (
            "filename",
            "the filename to write to"
          )
        )
      ),
      m_PharmMapFileFlag
      (
        new command::FlagDynamic
        (
          "pharmmap",
          "file containing pharmmap coefficients",
          command::Parameter
          (
            "filename",
            "name of the file that contains coefficient data",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_BinaryDifferenceFlag
      (
        new command::FlagDynamic
        (
          "binary_score",
          "whether scores between pairs of molecules should be evaluated as a binary greater/less",
          command::Parameter
          (
            "binary_score",
            "if scores are within this value they are considered to be exactly the same",
            command::ParameterCheckRanged< size_t>( 0, 1),
            "0"
          )
        )
      ),
      m_ScoreToleranceFlag
      (
        new command::FlagDynamic
        (
          "tolerance",
          "score tolerance for whether two scores are considered equivalent",
          command::Parameter
          (
            "tolerance",
            "if scores are within this value they are considered to be exactly the same",
            command::ParameterCheckRanged< float>( 0, std::numeric_limits< float>::max()),
            "0.0"
          )
        )
      ),
      m_PercentCoverageForScaffoldFlag
      (
        new command::FlagDynamic
        (
          "coverage",
          "this percentage of the average number of atoms in provided molecules is how big scaffolds must be",
          command::Parameter
          (
            "coverage",
            "if scores are within this value they are considered to be exactly the same",
            command::ParameterCheckRanged< float>( 0, 1),
            "0.5"
          )
        )
      )
    {
    }

    void Worker::RunThread()
    {
      graph::CommonSubgraphIsomorphism< size_t, size_t> csi( graph::CommonSubgraphIsomorphismBase::e_Connected);

      if( !m_Ensemble.IsDefined() || !( m_End <= m_NumberMols && m_Begin < m_End))
      {
        BCL_MessageStd( "Bad parameters given to thread");
        return;
      }

      size_t first( m_Begin);
      size_t second( m_Begin + 1);

      size_t pairs_compared( 0);
      for( ; first < m_End; first += 2, second += 2, ++pairs_compared)
      {
        if( second >= m_End)
        {
          second = m_Begin;
        }
        const chemistry::FragmentComplete &first_mol( ( *m_Ensemble)( first));
        const chemistry::FragmentComplete &second_mol( ( *m_Ensemble)( second));

        graph::ConstGraph< size_t, size_t> first_graph( m_GraphMaker( first_mol));
        graph::ConstGraph< size_t, size_t> second_graph( m_GraphMaker( second_mol));

        csi.SetGraphA
        (
          util::OwnPtr< graph::ConstGraph< size_t, size_t> >( &first_graph, false)
        );
        csi.SetGraphB
        (
          util::OwnPtr< graph::ConstGraph< size_t, size_t> >( &second_graph, false)
        );

        csi.FindIsomorphism( csi.EstimateUpperBounds(), size_t( 2));

        storage::Vector< size_t> common_atoms
        (
          csi.GetIsomorphism().GetKeysAsVector()
        );
        if( !common_atoms.IsEmpty())
        {
          chemistry::AtomVector< chemistry::AtomComplete> new_vector( first_mol.GetAtomVector());
          new_vector.Reorder
          (
            common_atoms
          );
          m_LocalEnsemble.PushBack( chemistry::FragmentComplete( new_vector, ""));
        }
/*        else
        {
          io::OFStream out;
          io::File::MustOpenOFStream( out, "/tmp/fail.sdf", std::ios::app);
          first_mol.WriteMDL( out);
          second_mol.WriteMDL( out);
          chemistry::FragmentComplete newmol;
          newmol.WriteMDL( out);
        } */
      }
    }

    // Construct the static instance of this application, and add it to the ChemInfo group
    const ApplicationType PharmMap::PharmMap_Instance
    (
      GetAppGroups().AddAppToGroup( new PharmMap(), GetAppGroups().e_ChemInfo)
    );

  } // namespace app
} // namespace bcl
