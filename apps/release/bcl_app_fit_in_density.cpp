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
#include "assemble/bcl_assemble_printer_protein_model_movie.h"
#include "assemble/bcl_assemble_protein_storage_file.h"
#include "assemble/bcl_assemble_quality.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "command/bcl_command_parameter_check_extension.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "coord/bcl_coord_geometric_hash_storage_classes.h"
#include "coord/bcl_coord_geometric_hashing.h"
#include "coord/bcl_coord_move_rotate_random.h"
#include "coord/bcl_coord_move_transform_random.h"
#include "coord/bcl_coord_move_translate_random.h"
#include "coord/bcl_coord_point_cloud.h"
#include "coord/bcl_coord_point_to_key_classes.h"
#include "density/bcl_density_map.h"
#include "density/bcl_density_protein_agreements.h"
#include "density/bcl_density_simulators.h"
#include "fold/bcl_fold_mutate_protein_model.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_histogram.h"
#include "math/bcl_math_mutate_decision_node.h"
#include "mc/bcl_mc_approximator.h"
#include "mc/bcl_mc_movie_printers.h"
#include "mc/bcl_mc_printer_combined.h"
#include "mc/bcl_mc_temperature_accepted.h"
#include "opti/bcl_opti_criterion_combine.h"
#include "opti/bcl_opti_criterion_number_iterations.h"
#include "opti/bcl_opti_criterion_unimproved.h"
#include "pdb/bcl_pdb_handler.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_unary_function_job_with_data.h"
#include "score/bcl_score_protein_atom_density.h"
#include "util/bcl_util_stopwatch.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FitInDensity
    //! @brief Class to fit structures in density maps
    //! @details Application is for fitting atomic protein structure into cryoEM or other medium resolution (5-20 A) density
    //! maps.
    //!
    //! @see @link example_app_fit_in_density.cpp @endlink
    //! @author woetzen
    //! @date 02/04/2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FitInDensity :
      public InterfaceRelease
    {
    private:

    //////////
    // data //
    //////////

      //!input pdb and mrc file parameter
      util::ShPtr< command::ParameterInterface> m_PDBFilenameParam;

      //!input density filename
      util::ShPtr< command::ParameterInterface> m_MRCFilenameParam;

      //! flag to determine the prefix for output files
      util::ShPtr< command::FlagInterface> m_PrefixFlag;

      //! flag and parameter for resolution
      util::ShPtr< command::FlagStatic> m_QuantizationResolutionFlag;
      util::ShPtr< command::ParameterInterface> m_QuantizationResolutionAngularParam;
      util::ShPtr< command::ParameterInterface> m_QuantizationResolutionDistanceParam;

      //! flag and parameter for adjusting the resolution of the given and simulated density
      util::ShPtr< command::FlagInterface> m_MrcResolutionFlag;

      //! flag and parameter for the number of features used for representing the electron density map
      util::ShPtr< command::FlagInterface> m_NumberFeaturesFlag;

      //! flag and parameter for the multiple of the number of points derived from all the parameters and the protein to be fitted
      util::ShPtr< command::FlagInterface> m_MultipleOfExpectedPointsFlag;

      //! flag and parameter for the minimal distance between two features
      util::ShPtr< command::FlagInterface> m_FeatureDistanceFlag;

      //! flag and parameters for generating features from intensity and gradient in density map
      util::ShPtr< command::FlagInterface> m_RatioIntensityGradientFlag;

      //! flag for writing feature clouds to file
      util::ShPtr< command::FlagInterface> m_WriteFeatureCloudFlag;

      //! flag and two parameters for the threshold for the length of the sides of the triangular base in the geometric hashing algorithm
      util::ShPtr< command::FlagDynamic> m_ThresholdFlag;

      //! flag with parameter for the radius around the triangular base where coordinates are considered features
      util::ShPtr< command::FlagInterface> m_FeatureRadiusFlag;

      //! flag and parameter for the coordinate system to be used for the geometric hashing
      util::ShPtr< command::FlagInterface> m_CoordinateSystemFlag;

      //! flag and parameter for number of trials and number of best results
      util::ShPtr< command::FlagStatic> m_TrialsSavebestFlag;
      util::ShPtr< command::ParameterInterface> m_NumberTrialsParam;
      util::ShPtr< command::ParameterInterface> m_NumberSavebestParam;

      //! flag and two parameters for the allowed difference in rotation and translation between two hits
      util::ShPtr< command::FlagStatic> m_DiffRotTransFlag;
      util::ShPtr< command::ParameterInterface> m_DiffRotParam;
      util::ShPtr< command::ParameterInterface> m_DiffTransParam;

      //! flag for a dynamic list of atoms to be fitted
      util::ShPtr< command::FlagInterface> m_AtomListFlag;

      //! flags for Monte Carlo minimization
      //! flag to write the complete minimization
      util::ShPtr< command::FlagInterface> m_WriteMinimizationFlag;

      //! flag for writing to shared memory
      util::ShPtr< command::FlagInterface> m_WriteSharedMemoryFlag;

      //! flag for choice of density agreement function
      util::ShPtr< command::FlagInterface> m_DensityAgreementFlag;

      //! flag or choice of density simulator
      util::ShPtr< command::FlagInterface> m_SimulatorFlag;

      //! max number of rejected steps and max iterations
      util::ShPtr< command::FlagStatic> m_McMaxIterationsUnimprovedFlag;
      util::ShPtr< command::ParameterInterface> m_McMaxIterationsParam;
      util::ShPtr< command::ParameterInterface> m_McMaxStepsUnimprovedParam;

      //! mutate parameters
      util::ShPtr< command::FlagStatic> m_McMutateTransRotFlag;
      util::ShPtr< command::ParameterInterface> m_McMutateTransParam;
      util::ShPtr< command::ParameterInterface> m_McMutateRotParam;

      //! quantize flag
      util::ShPtr< command::FlagStatic> m_QuantizeFlag;
      util::ShPtr< command::ParameterInterface> m_CenterXParam;
      util::ShPtr< command::ParameterInterface> m_CenterYParam;
      util::ShPtr< command::ParameterInterface> m_CenterZParam;

      //! flag for removing features in low occupied regions
      util::ShPtr< command::FlagStatic> m_RemoveNoiseFlag;
      util::ShPtr< command::ParameterInterface> m_NumberNeighbors;
      util::ShPtr< command::ParameterInterface> m_MinNeighborDistance;

      //! protein storage
      mutable util::ShPtr< assemble::ProteinStorageFile> m_Storage;

      //! density map
      mutable util::ShPtr< density::Map> m_SpMap;

      //! vector containing all column names for result table
      static const storage::TableHeader s_TableHeader;

      //! @brief source string for storing initial fits
      static const std::string &GetSourceInitial()
      {
        static const std::string s_source_inital( "transformed");
        return s_source_inital;
      }

      //! @brief source string for storing minimized fits
      static const std::string &GetSourceMinimized()
      {
        static const std::string s_source_minimized( "transformed_min");
        return s_source_minimized;
      }

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      FitInDensity();

    public:

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      FitInDensity *Clone() const
      {
        return new FitInDensity( *this);
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

      //! @brief return the bcl::commons name
      //! @return string for the bcl::commons name of that application
      std::string GetBCLScopedName() const
      {
        return "BCL::EMFit";
      }

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const
      {
        return "fits atomic protein structure into cryoEM or other medium resolution (5-20 A) density maps";
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief initializes the command object for that application
      util::ShPtr< command::Command> InitializeCommand() const;

      //! main function
      int Main() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief create a minimizer from a density map
      //! @param NUMBER the number of this minimizer used to identify eventual written minimizations
      //! @return the minimizer
      util::ShPtr< mc::Approximator< assemble::ProteinModel, double> >
      CreateMinimizer( const size_t NUMBER) const;

      //! @brief minimizes the initial fits
      //! @param TRANSFORMS
      //! @param PROTEIN_MODEL
      //! @param MAX_NUMBER the maximal number of fits to be minimized, after sorting the initial fits by cross correlation
      //! @return a storage::Table
      storage::Table< double>
      MinimizeInitialFits
      (
        const storage::List< storage::Pair< math::TransformationMatrix3D, size_t> > &TRANSFORMS,
        const assemble::ProteinModel &PROTEIN_MODEL,
        const size_t MAX_NUMBER
      ) const;

      //! @brief initialize PointCloud - try to read from file, otherwise call the appropriate function
      //! @param FILENAME_PREFIX prefix for the file read from or written to
      //! @param FEATURE_DISTANCE the minimal distance between two features
      //! @param NUMBER_FEATURES the number of features that should be in the density map
      //! @return feature cloud either read from file if existent or generated from static parameters from command line
      coord::PointCloud InitializePointCloud
      (
        const std::string &FILENAME_PREFIX,
        const double FEATURE_DISTANCE,
        const size_t NUMBER_FEATURES
      ) const;

      //! @brief quantizes a feature cloud
      //! @param POINTCLOUD unquantized point cloud from density or pdb structure
      //! @param POINT_TO_KEY function to convert into a key
      //! @return feature cloud, containing the quantized feature cloud depending on the chosen parameteres
      coord::PointCloud QuantizeFeatureCloud
      (
        const coord::PointCloud &POINTCLOUD,
        const coord::PointToKeyInterface &POINT_TO_KEY
      ) const;

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

      //! instance of the application
      static const ApplicationType FitInDensity_Instance;

    }; // class FitInDensity

    //! vector containing all column names for result table
    const storage::TableHeader FitInDensity::s_TableHeader
    (
      storage::Vector< std::string>::Create( "hashscorerank", "counts", "corr", "RMSD", "min_corr", "min_RMSD")
    );

    //! main function
    int FitInDensity::Main() const
    {
      // initialize write and read stream objects
      io::OFStream write;
      io::IFStream read;

      // initialize protein storage
      m_Storage = assemble::ProteinStorageFile::GetDefaultStorage();

      // check if storage could be initialized
      BCL_Assert( m_Storage.IsDefined(), "unable to initialize the protein storage");

      // create a source string for the storage
      const std::string source( io::File::SplitToPathAndFileName( m_PrefixFlag->GetFirstParameter()->GetValue()).Second());

      // instantiate DensityMap from mrc file
      BCL_MessageStd( "read density map from mrc file");
      io::File::MustOpenIFStream( read, m_MRCFilenameParam->GetValue(), std::ios::binary);
      m_SpMap = util::ShPtr< density::Map>( new density::Map());
      m_SpMap->ReadMRC( read);
      io::File::CloseClearFStream( read);
      BCL_MessageStd( "read density map from mrc file done");

      // create simulator
      const util::ShPtr< density::SimulateInterface> simulator
      (
        density::GetSimulators().CreateSimulator
        (
          density::Simulator( m_SimulatorFlag->GetFirstParameter()->GetValue()),
          m_SpMap->GetCellWidth(), m_MrcResolutionFlag->GetFirstParameter()->GetNumericalValue< double>()
        )
      );

      // construct Hashmap
      const std::string mrc_name( io::File::RemovePath( m_MRCFilenameParam->GetValue()));
      BCL_MessageStd( "name of mrc: " + mrc_name);

      // open pdb file
      io::File::MustOpenIFStream( read, m_PDBFilenameParam->GetValue());
      pdb::Handler pdb( read);
      io::File::CloseClearFStream( read);

      BCL_MessageStd( "matching: " + m_PDBFilenameParam->GetValue());

      pdb::Factory factory;

      // create protein model from pdb
      assemble::ProteinModel model( factory.ProteinModelFromPDB( pdb));
      storage::Map< biol::SSType, size_t> sse_min_size;
      sse_min_size[ biol::GetSSTypes().HELIX] = 0;
      sse_min_size[ biol::GetSSTypes().STRAND] = 0;
      const assemble::ProteinModel all_sse_model( factory.ProteinModelFromPDB( pdb, sse_min_size));
      sse_min_size[ biol::GetSSTypes().COIL] = 0;
      const assemble::ProteinModel complete_model( factory.ProteinModelFromPDB( pdb, sse_min_size));

      // create storage::Vector of BBAtomTypes to be matched - need to extract them from protein model
      const storage::Set< biol::AtomType> bb_atoms( m_AtomListFlag->GetObjectSet< biol::AtomType>());

      // collect all atom coordinates
      const util::SiPtrVector< const linal::Vector3D> atomcoordinates( model.GetAtomCoordinates( bb_atoms));
      const util::SiPtrVector< const linal::Vector3D> atomcoordinates_helix_strand( all_sse_model.GetAtomCoordinates( bb_atoms));

      // Calculate the radius of gyration of the atoms to be fitted in the protein
      const double fit_radius_gyration( coord::RadiusOfGyration( atomcoordinates));
      const double helix_sheet_radius_of_gyration( coord::RadiusOfGyration( atomcoordinates_helix_strand));
      const double all_radius_of_gyration( coord::RadiusOfGyration( complete_model.GetAtomCoordinates( biol::GetAtomTypes().GetBackBoneAtomTypes())));

      double feature_radius( m_FeatureRadiusFlag->GetFirstParameter()->GetNumericalValue< double>());

      BCL_MessageStd( "number of atoms to be matched: " + util::Format()( atomcoordinates.GetSize()));
      BCL_MessageStd( "radius of gyration of coordinates to be matched: " + util::Format()( fit_radius_gyration));
      BCL_MessageStd( "radius of gyration of coordinates of helix-sheet model: " + util::Format()( helix_sheet_radius_of_gyration));
      BCL_MessageStd( "radius of gyration of coordinates of complete model: " + util::Format()( all_radius_of_gyration));

      // feature distance (minimal distance between two points in feature cloud) from command line
      double feature_distance( std::max( m_SpMap->GetCellWidth().X(), m_FeatureDistanceFlag->GetFirstParameter()->GetNumericalValue< double>()));

      // derive threshold for triangular base and the feature radius from the radius of gyration of the atoms to be matched
      if( !m_FeatureDistanceFlag->GetFirstParameter()->GetWasSetInCommandLine())
      {
        feature_distance = std::max( m_SpMap->GetCellWidth().X(), 0.15 * all_radius_of_gyration);
      }
      if( !m_FeatureRadiusFlag->GetFirstParameter()->GetWasSetInCommandLine())
      {
        feature_radius = 1.25 * all_radius_of_gyration;
      }

      // message to user
      BCL_MessageStd
      (
        "threshold and feature radius was derived from the radius of gyration - in cases where it was not set in the commandline\n"
        "feature_distance: " + util::Format()( feature_distance) + "\n"
        "feature_radius: " + util::Format()( feature_radius)
      );

      if( fit_radius_gyration > feature_radius)
      {
        BCL_MessageCrt
        (
          "feature radius should be at least as big as the radius of gyration of the atoms of the protein to be matched: "
          + util::Format()( fit_radius_gyration) + " > " + util::Format()( feature_radius)
        );
      }

      // a score object for calculating the atom density of backbone atoms of SSEs within a protein model
      const score::ProteinAtomDensity score_protein_atom_density
      (
        storage::Set< biol::SSType>( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND),
        biol::GetAtomTypes().GetBackBoneAtomTypes(),
        m_SpMap->GetCellWidth()
      );

      // number of filled voxels considering all atoms that are expected by sses and using the voxel size
      const size_t number_filled_boxes_voxel( coord::QuantizePoints( atomcoordinates, m_SpMap->GetCellWidth()).GetSize());
      BCL_MessageStd( "number of expectedly filled voxels: " + util::Format()( number_filled_boxes_voxel));

      // number of filled voxels considering all atoms that are expected by sses and using the voxel size
      // TODO: round up to a multiple of the cell width!!
      const linal::Vector3D cell_width_feature_distance
      (
        std::ceil( feature_distance / m_SpMap->GetCellWidth()( 0)) * m_SpMap->GetCellWidth()( 0),
        std::ceil( feature_distance / m_SpMap->GetCellWidth()( 1)) * m_SpMap->GetCellWidth()( 1),
        std::ceil( feature_distance / m_SpMap->GetCellWidth()( 2)) * m_SpMap->GetCellWidth()( 2)
      );
      const size_t number_filled_boxes_feature_distance( coord::QuantizePoints( atomcoordinates, cell_width_feature_distance).GetSize());
      BCL_MessageStd( "number of expected points with given feature distance: " + util::Format()( number_filled_boxes_feature_distance));

      // number of points that are expected per protein model
      const size_t number_points_in_model
      (
        std::min( number_filled_boxes_voxel, number_filled_boxes_feature_distance)
        //        number_filled_boxes_voxel * m_SpMap->GetVoxelVolume() * actual_target_point_density // protein_volume * density
      );

      // extract all atoms for later transformation and correlation calculation
      coord::PointCloud match_pdb_positions
      (
        coord::PointCloud( util::ConvertToStorageVector< linal::Vector3D, const linal::Vector3D>( atomcoordinates))
      );

      BCL_MessageStd( "number of features to match: " + util::Format()( match_pdb_positions.GetSize()));

      // derive the number of features in the feature cloud
      size_t number_features( 0);
      if( m_NumberFeaturesFlag->GetFirstParameter()->GetWasSetInCommandLine())
      {
        number_features = m_NumberFeaturesFlag->GetFirstParameter()->GetNumericalValue< size_t>();
      }
      else
      {
        // multiply the number of features per model with the multiplier from the commandline, for densities that contain a multiple of the protein to be matched
        number_features = size_t
        (
          m_MultipleOfExpectedPointsFlag->GetFirstParameter()->GetNumericalValue< double>() * number_points_in_model
        );
      }

      // create quantization function
      const util::ShPtr< coord::PointToKeyInterface> point_to_key
      (
        coord::GetPointToKeyClasses().GetPointToKeyFunction
        (
          m_CoordinateSystemFlag->GetFirstParameter()->GetValue(),
          m_QuantizationResolutionAngularParam->GetNumericalValue< double>(),
          m_QuantizationResolutionDistanceParam->GetWasSetInCommandLine() ?
            m_QuantizationResolutionDistanceParam->GetNumericalValue< double>() : feature_distance
        )
      );

      // initialize feature cloud
      const coord::PointCloud feature_cloud( InitializePointCloud( mrc_name, feature_distance, number_features));

      storage::VectorND< 4, double> thresholds;

      if( m_ThresholdFlag->GetParameterList().GetSize() < 4)
      {
        // determine thresholds, so that they are equally occupied with distance
        const storage::VectorND< 4, double> thresholds_equal_occupied
        (
          coord::GeometricHashing::CreateEqualOccupiedIntervals
          (
            coord::GeometricHashing::CalculateDistanceDistribution
            (
              util::ConvertToConstSiPtrVector< linal::Vector3D>( match_pdb_positions),
              feature_distance / 10.0,
              math::Range< double>( 0.0, 1.25 * all_radius_of_gyration)
            )
          )
        );
        thresholds = thresholds_equal_occupied;
        BCL_MessageStd
        (
          "threshold derived from points to be fitted: " + util::Format()( thresholds)
        );
      }
      else
      {
        for( size_t i( 0); i < 4; ++i)
        {
          thresholds( i) = m_ThresholdFlag->GetParameterList()( i)->GetNumericalValue< double>();
        }
        BCL_MessageStd
        (
          "threshold from command line: " + util::Format()( thresholds)
        );
      }

      // initialize hash storage and hashmap
      util::ShPtr< coord::GeometricHashStorageInterface> hash_storage
      (
        coord::GetGeometricHashStorageClasses().ConstructFromCommandline
        (
          mrc_name,
          m_MrcResolutionFlag->GetFirstParameter()->GetNumericalValue< double>(),
          m_SpMap->GetCellWidth().X(),
          storage::VectorND< 3, size_t>( m_SpMap->GetDimensions()),
          number_features,
          feature_distance,
          m_RatioIntensityGradientFlag->GetFirstParameter()->GetNumericalValue< double>(),
          thresholds,
          feature_radius,
          *point_to_key,
          m_NumberNeighbors->GetNumericalValue< size_t>(),
          m_MinNeighborDistance->GetNumericalValue< double>()
        )
      );
      // pass the storage to the geometric hashing object
      coord::GeometricHashing hashmap( hash_storage);

      // write the points that should be matched to file
      if( m_WriteFeatureCloudFlag->GetFlag())
      {
        // write the pointcloud pdb
        const std::string pdb_pointcloud_filename( m_PDBFilenameParam->GetValue() + "_points.pdb");
        match_pdb_positions.WriteToPDB( pdb_pointcloud_filename);
      }

      BCL_MessageStd( "target number of features in cloud: " + util::Format()( number_features));

      // quantize
      if( m_QuantizeFlag->GetFlag())
      {
        const std::string filename_appendix( "_" + m_QuantizationResolutionFlag->GetFirstParameter()->GetValue() + "_" + m_CoordinateSystemFlag->GetFirstParameter()->GetValue());
        // translate the protein's center to the center of the coordinate frame, since the spherical quantization would
        // not give neat picture
        BCL_MessageStd( "center of mass protein model: " + util::Format()( model.GetCenterOfMass()));
        const linal::Vector3D center
        (
          m_CenterXParam->GetNumericalValue< double>(),
          m_CenterYParam->GetNumericalValue< double>(),
          m_CenterZParam->GetNumericalValue< double>()
        );
        model.Translate( -center);

        const std::string centered_pdb_file_name
        (
          m_PrefixFlag->GetFirstParameter()->GetValue() + m_PDBFilenameParam->GetValue() + "_center.pdb"
        );

        // write centered protein
        BCL_MessageStd( "write centered protein to: " + centered_pdb_file_name);
        io::File::MustOpenOFStream( write, centered_pdb_file_name);
        factory.WriteModelToPDB( model, write);
        io::File::CloseClearFStream( write);

        // write selected atoms
        const std::string selected_pdb_file_name
        (
          m_PrefixFlag->GetFirstParameter()->GetValue() + m_PDBFilenameParam->GetValue() + "_selected.pdb"
        );
        match_pdb_positions.Translate( -center);
        BCL_MessageStd( "write selected atoms to: " + selected_pdb_file_name);
        match_pdb_positions.WriteToPDB( selected_pdb_file_name);

        // write quantized atoms
        const std::string quantized_pdb_file_name
        (
          m_PrefixFlag->GetFirstParameter()->GetValue() + "_protein_" +
          filename_appendix + "_quantized.pdb"
        );

        const coord::PointCloud quantized_features( QuantizeFeatureCloud( match_pdb_positions, *point_to_key));
        BCL_MessageStd( "write quantized atoms to: " + quantized_pdb_file_name);
        quantized_features.WriteToPDB( quantized_pdb_file_name);

        // center feature cloud
        coord::PointCloud feature_cloud( InitializePointCloud( mrc_name, feature_distance, number_features));
        BCL_MessageStd( "center of mass feature cloud: " + util::Format()( feature_cloud.GetCenter()));
        feature_cloud.Translate( -center);
        std::string feature_cloud_filename;
        feature_cloud_filename =
          m_PrefixFlag->GetFirstParameter()->GetValue() + "features_both11_" + util::Format()( number_features) +
          "_featuredist_" + util::Format()( feature_distance) + "A_feature_cloud_center.pdb";
        BCL_MessageStd( "write centered feature cloud to: " + feature_cloud_filename);
        feature_cloud.WriteToPDB( feature_cloud_filename);

        // quantize feature cloud
        const coord::PointCloud quantized_feature_cloud( QuantizeFeatureCloud( feature_cloud, *point_to_key));
        feature_cloud_filename =
          m_PrefixFlag->GetFirstParameter()->GetValue() + "features_both11_" + util::Format()( number_features) +
          "_featuredist_" + util::Format()( feature_distance) + "A" + filename_appendix +
          "_feature_cloud_quantized.pdb";
        BCL_MessageStd( "write quantized feature_cloud to: " + feature_cloud_filename);
        quantized_feature_cloud.WriteToPDB( feature_cloud_filename);

        return 0;
      }

      // start with the algorithm
      // if there was not such geometric hashmap before, calculate the feature cloud and build hashmap
      if( !hash_storage->IsReadOnly())
      {
        BCL_MessageStd
        (
          "actual number of features in cloud: " + util::Format()( feature_cloud.GetSize())
        );

        // build hash with normal feature cloud
        BCL_MessageStd
        (
          "building hash with thresholds: " + util::Format()( thresholds) + ", feature radius of: " +
          util::Format()( feature_radius) + " and quantization resolution: " +
          m_QuantizationResolutionFlag->GetFirstParameter()->GetValue() + " - can take a little time"
        );

        // building hashmap
        BCL_MessageStd( "calculate geometric hash started");
        util::Stopwatch stopwatchhash( "build hash", util::Message::e_Standard, true);
        hashmap.BuildHash( feature_cloud);
        BCL_MessageStd( "calculate geometric hash finished");
      }
      else
      {
        BCL_MessageStd
        (
          "no necessity to calculate feature cloud and building hash since map was already calculated"
        );
      }

      // search hash and return best pairs( with highest count in Hashmap) of matching transformation matrices
      storage::List< storage::Pair< math::TransformationMatrix3D, size_t> > transforms;

      // this is the geometric hash fitting
      BCL_MessageStd( "fitting started");
      {
        util::Stopwatch stopwatchfitting( "hash fitting", util::Message::e_Standard, true);
        transforms = hashmap.SearchTarget
                     (
                       match_pdb_positions,
                       m_NumberSavebestParam->GetNumericalValue< size_t>(),
                       m_NumberTrialsParam->GetNumericalValue< size_t>(),
                       storage::VectorND< 2, double>
                       (
                         m_DiffRotParam->GetNumericalValue< double>(),
                         m_DiffTransParam->GetNumericalValue< double>()
                       )
                     );
        // after this scope util::Stopwatch will output process duration
      }
      BCL_MessageStd( "minimization started");

      util::Stopwatch minimize_watch( "minimization", util::Message::e_Standard, true);
      // counts, calculated correlation factors and rmsd
      storage::Table< double>
        count_corr_rmsd_corrmin_rmsdmin
        (
          MinimizeInitialFits( transforms, complete_model, m_NumberSavebestParam->GetNumericalValue< size_t>())
        );

      BCL_MessageStd( "minimization finished");

      // print results
      count_corr_rmsd_corrmin_rmsdmin.SortByColumn( "min_corr", ( **math::Comparisons< double>::GetEnums().e_Greater));
      const std::string result_filename( m_PrefixFlag->GetFirstParameter()->GetValue() + "result.table");
      io::File::MustOpenOFStream( write, result_filename);
      count_corr_rmsd_corrmin_rmsdmin.WriteFormatted( write);
      io::File::CloseClearFStream( write);

      // print also in log file
      count_corr_rmsd_corrmin_rmsdmin.WriteFormatted( util::GetLogger());

      BCL_MessageStd( "fit in density finished");

      // end
      return 0;
    } // end Main

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &FitInDensity::GetReadMe() const
    {
      // create a static string to hold readme information
      static const std::string s_readme_text
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of BCL::EMFit, terms of use, appropriate citation, installation "
        "procedures, BCL::EMFit execution, technical support, and future research directions.\n"
        "\n"
        + DefaultSectionSeparator() +
        "II. WHAT IS BCL::EMFit?\n"
        "BCL::EMFit is a C++ based application, created by Vanderbilt University's Meiler Laboratory, which is part of "
        "a larger library of applications called BCL::Commons.\n"
        "BCL::EMFit fits a given atomic protein structure in a medium resolution electron density map. It utilizes an "
        "a geometric hashing algorithm on extracted features from the given density map and the protein to be fitted, "
        "to generate a list of initial fits. This initial list is minimized employing a Monte Carlo/Metropolis "
        "simulated annealing optimization protocol on the CCC between the experimental density map and the simulated "
        "map from the protein to fit.\n"
        "\n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        + DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. APPROPRIATE CITATIONS FOR USING BCL::EMFit.\n"
        "When using BCL::EMFit in a publication, please cite the following publications describing the application's "
        "development:\n"
        "\n"
        "N. Woetzel, S. Lindert, P. L. Stewart, and J. Meiler, BCL::EM-Fit: Rigid body fitting of atomic structures "
        "into density maps using geometric hashing and real space refinement. J. Struct. Biol., May 2011."
        "\n"
        + DefaultSectionSeparator() +
        "V. INSTALLATION PROCEDURES.\n"
        + DefaultInstallationProcedure() +
        "\n"
        + DefaultSectionSeparator() +
        "VI. RUNNING BCL::EMFit.\n"
        "Running BCL::EMFit consists of three main steps.\n"
        "\n"
        "1) Identify the density map of interest in which a protein will be fitted.\n"
        "The only format currently supported is the mrc/ccp4 format.  If sections of the density maps are to be "
        "considered, they will have to be cut manually. Any format conversions will have to be performed with other "
        "software.\n"
        "\n"
        "2) Create the protein structure of interest using the pdb format.\n"
        "If biomolecules (multimers) have to be fitted, BCL::PDBConvert or other tools to generate the appropriate pdb "
        " structure and store it in the pdb.\n"
        "\n"
        "3) Run BCL::EMFit:\n"
        "Please read the publication on the different parameters and there meaning; also call BCL::EMFit on the "
        "command line with \"-help\" which documents the meaning of each flag. A sample command line could be:\n"
        "\n"
        "bcl_em_fit.exe 1OELG.pdb groel.mrc -mrc_resolution 5.4 -multiple_features 14 -atoms CA -min_sse_size 0 999 "
        "999 -number_bases_refine 1000 25 -diff_rot_trans 3 5 -protein_storage File fit4/ Create -prefix fit4 "
        "-scheduler PThread 8\n"
        "\n"
        "PARAMETERS:\n"
        "1OELG.pdb groel.mrc are the protein and the density maps of interest\n"
        "\n"
        "FLAGS:\n"
        "\n"
        "-mrc_resolution 5.4 -> is the experimental resolution of the density map\n"
        "-multiple_features 14 -> is the approximate number of copies of the given protein expected in the "
        "experimental density map\n"
        "-atoms CA -> a space separated list are of the backbone atoms of the protein that are used as features to "
        "hash in the map of the features for the experimental density map\n"
        "-min_sse_size 0 999 999 are the minimal sizes of the SSEs that are selected for matching, e.g. all helices, "
        "and virtually no strands or loops (is they are shorter than 999)\n"
        "-number_bases_refine <bases> <number initial fits to report and refine> bases is the number of bases that are "
        "picked from the protein features and hashed in the hash map; number initial fits to report and minimize is "
        "the list of best hash scoring fits, that are refined in the MCM protocol\n"
        "-protein_storage File fit4/ -> store all proteins in the folder \"fit4/\" as files\n"
        "-prefix fit4 -> all other files will have a prefix \"fit4\"\n"
        "-scheduler PThread 8 -> hash and minimize using 8 threads (processors)\n"
        "\n"
        "INPUT AND OUTPUT.\n"
        "\n"
        "BCL::EMFit requires two inputs, a pdb file and a mrc file\n"
        "the output are a list of initial fits and refined fits with the prefix \"transformed\"; additionally the "
        "hash score and CCC is written to the std output and to a <prefix>result.table file\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        + DefaultTechnicalSupportString() +
        "\n"
        + DefaultSectionSeparator() +
        "VIII. FUTURE DEVELOPMENT OF BCL::EMFit.\n"
        "BCL::EMFit is under ongoing further development.  For current research please refer to www.meilerlab.org\n"
        + DefaultSectionSeparator()
      );

      // return readme information
      return s_readme_text;
    }

    //! @brief create a minimizer from a density map
    //! @param NUMBER the number of this minimizer used to identify eventual written minimizations
    //! @return the minimizer
    util::ShPtr< mc::Approximator< assemble::ProteinModel, double> >
    FitInDensity::CreateMinimizer( const size_t NUMBER) const
    {
      util::ShPtr< score::ProteinModelScoreSum> sp_score( new score::ProteinModelScoreSum());
      sp_score->NewOperand
      (
        *density::GetProteinAgreements().CreateProteinAgreement
        (
          density::ProteinAgreement( m_DensityAgreementFlag->GetFirstParameter()->GetValue()),
          density::Simulator( m_SimulatorFlag->GetFirstParameter()->GetValue()),
          m_SpMap,
          m_MrcResolutionFlag->GetFirstParameter()->GetNumericalValue< double>()
        )
      );

      // mutate rotate
      util::ShPtr< math::MutateInterface< assemble::ProteinModel> > sp_mutate_rotate
      (
        new fold::MutateProteinModel
        (
          coord::MoveRotateRandom
          (
            m_McMutateRotParam->GetNumericalValue< double>()
          )
        )
      );

      // mutate translate
      util::ShPtr< math::MutateInterface< assemble::ProteinModel> > sp_mutate_translate
      (
        new fold::MutateProteinModel
        (
          coord::MoveTranslateRandom
          (
            m_McMutateTransParam->GetNumericalValue< double>()
          )
        )
      );

      // mutate function rot trans
      util::ShPtr< math::MutateInterface< assemble::ProteinModel> > sp_mutate_rot_trans
      (
        new fold::MutateProteinModel
        (
          coord::MoveTransformRandom
          (
            m_McMutateTransParam->GetNumericalValue< double>(),
            m_McMutateRotParam->GetNumericalValue< double>()
          )
        )
      );

      math::ObjectProbabilityDistribution< math::MutateInterface< assemble::ProteinModel> > mutates;
      mutates.PushBack
      (
          0.25,
          *sp_mutate_translate
      );
      mutates.PushBack
      (
          0.25,
          *sp_mutate_rotate
      );
      mutates.PushBack
      (
          0.25,
          *sp_mutate_rot_trans
      );

      // mutate function
      util::ShPtr< math::MutateDecisionNode< assemble::ProteinModel> > sp_mutate
      (
        new math::MutateDecisionNode< assemble::ProteinModel>( mutates)
      );

      // create shared pointer to PrinterCombined
      util::ShPtr< mc::PrinterCombined< assemble::ProteinModel, double> >
        sp_printer( new mc::PrinterCombined< assemble::ProteinModel, double>());

      // if user wished to print all intermediate minimization steps
      if( m_WriteMinimizationFlag->GetFlag())
      {
        // initialize storage::Set of step status to be used for printer
        storage::Set< opti::StepStatusEnum> step_status_set;

        // iterate over step statuses
        for
        (
          util::ShPtrVector< command::ParameterInterface>::const_iterator
            step_status_itr( m_WriteMinimizationFlag->GetParameterList().Begin()),
            step_status_itr_end( m_WriteMinimizationFlag->GetParameterList().End());
          step_status_itr != step_status_itr_end;
          ++step_status_itr
        )
        {
          // insert into the set the correct enumerators (e.g. opti::e_Improved) corresponding to the strings
          // (e.g. improved)
          step_status_set.Insert( opti::StepStatusEnum( ( *step_status_itr)->GetValue()));
        }

        // mc movie printer
        util::ShPtr< mc::MoviePrinterInterface> sp_movie_printer
        (
          mc::MoviePrinter( mc::MoviePrinterInterface::GetParameterMoviePrinterType()->GetValue())
        );

        const std::string &prefix( m_PrefixFlag->GetFirstParameter()->GetValue() + util::Format()( NUMBER));
        sp_movie_printer->Initialize
        (
          prefix,
          storage::Vector< std::string>::Create
          (
            GetStaticClassName< storage::Table< double> >(),
            math::SumFunctionMixin< score::ProteinModel>::GetValueTableVerticalColumnNames()( 2)
          ),
          sp_score->GetFunctionSchemes(),
          mc::MoviePrinterInterface::GetParameterMovieWidth()->GetNumericalValue< size_t>(),
          mc::MoviePrinterInterface::GetParameterMovieHeight()->GetNumericalValue< size_t>(),
          bool( mc::MoviePrinterInterface::GetParameterRayTrace()->GetNumericalValue< size_t>() > 0)
        );
        util::ShPtr< mc::PrintInterface< assemble::ProteinModel, double> > printer_movie
        (
          new assemble::PrinterProteinModelMovie
          (
            prefix,
            sp_movie_printer,
            sp_score,
            step_status_set,
            quality::GetSuperimposeMeasures().e_NoSuperimpose,
            storage::Set< quality::Measure>()
          )
        );

        // set the water mark
        sp_movie_printer->SetWaterMark();

        sp_printer->Insert( printer_movie);
      }

      // set round number
      sp_printer->Initialize( NUMBER);

      // create the temperature control
      util::ShPtr< mc::TemperatureInterface> sp_temperature
      (
        new mc::TemperatureAccepted
        (
          mc::TemperatureAccepted::GetParameterStartFraction()->GetNumericalValue< double>(),
          mc::TemperatureAccepted::GetParameterEndFraction()->GetNumericalValue< double>(),
          m_McMaxIterationsParam->GetNumericalValue< size_t>(),
          mc::TemperatureAccepted::GetParameterStartTemperature()->GetNumericalValue< double>(),
          mc::TemperatureAccepted::GetParameterUpdateInterval()->GetNumericalValue< size_t>()
        )
      );

      // create the metropolis
      util::ShPtr< mc::Metropolis< double> > sp_metropolis
      (
        new mc::Metropolis< double>( sp_temperature, true)
      );

      // create the termination criterion
      opti::CriterionCombine< assemble::ProteinModel, double> criterion_combine;

      // insert termination criterion that depends on the total number of mc iterations
      criterion_combine.InsertCriteria
      (
        opti::CriterionNumberIterations< assemble::ProteinModel, double>
        (
          m_McMaxIterationsParam->GetNumericalValue< size_t>()
        )
      );

      // insert termination criterion that depends on the number of unimproved steps in a row
      criterion_combine.InsertCriteria
      (
        opti::CriterionUnimproved< assemble::ProteinModel, double>
        (
          m_McMaxStepsUnimprovedParam->GetNumericalValue< size_t>()
        )
      );

      // create approximator
      util::ShPtr< mc::Approximator< assemble::ProteinModel, double> > sp_approximator
      (
        new mc::Approximator< assemble::ProteinModel, double>
        (
          *sp_score,
          *sp_mutate,
          *sp_metropolis,
          criterion_combine
        )
      );

      // end
      return sp_approximator;
    }

    //! minimizes the initial fits
    storage::Table< double>
    FitInDensity::MinimizeInitialFits
    (
      const storage::List< storage::Pair< math::TransformationMatrix3D, size_t> > &TRANSFORMS,
      const assemble::ProteinModel &PROTEIN_MODEL,
      const size_t MAX_NUMBER
    ) const
    {
      // make pointers of transforms so that we can have random access later
      const util::SiPtrVector< const storage::Pair< math::TransformationMatrix3D, size_t> >
        transforms( TRANSFORMS.Begin(), TRANSFORMS.End());
      const storage::Set< biol::AtomType> rmsd_atom_types( biol::GetAtomTypes().CA);

      // construct a function that calculates the deviation between simulated density from atoms and the given density map
      util::ShPtr< density::ProteinAgreementInterface> density_ccc
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

      size_t i( 0);
      for
      (
        storage::List< storage::Pair< math::TransformationMatrix3D, size_t> >::const_iterator
          itr( TRANSFORMS.Begin()), itr_end( TRANSFORMS.End());
        itr != itr_end;
        ++itr, ++i
      )
      {
        util::ShPtr< assemble::ProteinModel> transformed_model_ptr( PROTEIN_MODEL.HardCopy());
        assemble::ProteinModel &transformed_model( *transformed_model_ptr);
        transformed_model.Transform( itr->First());

        storage::Row< double> &current_row( count_corr_rmsd_corrmin_rmsdmin.InsertRow( GetSourceInitial() + util::Format()( i)));
        current_row[ "hashscorerank"]   = double( i);
        current_row[ "counts"] = double( itr->Second());
        current_row[ "corr"]   = -density_ccc->operator()( transformed_model);
        current_row[ "RMSD"]   =
          assemble::Quality::Calculate
          (
            quality::GetMeasures().e_RMSD_NoSuperimposition, transformed_model, PROTEIN_MODEL, rmsd_atom_types
          );
        // Initialize minimized values to the initial values
        current_row[ "min_corr"] = current_row[ "corr"];
        current_row[ "min_RMSD"] = current_row[ "RMSD"];
      }

      // sort by initial CCC
      count_corr_rmsd_corrmin_rmsdmin.SortByColumn( "corr", ( **math::Comparisons< double>::GetEnums().e_Greater));

      // check that there are enough initial fits to minimize
      if( MAX_NUMBER > count_corr_rmsd_corrmin_rmsdmin.GetSize())
      {
        BCL_MessageStd
        (
          "cannot minimize " + util::Format()( MAX_NUMBER) + " initial fits, since that count is: " +
          util::Format()( count_corr_rmsd_corrmin_rmsdmin.GetSize())
        );
      }
      else
      {
        storage::Table< double>::iterator itr( count_corr_rmsd_corrmin_rmsdmin.Begin());
        storage::AdvanceIterator( itr, count_corr_rmsd_corrmin_rmsdmin.End(), MAX_NUMBER);
        // cutoff at the maximal number that should be minimzed
        count_corr_rmsd_corrmin_rmsdmin.Remove( itr, count_corr_rmsd_corrmin_rmsdmin.End());
      }

      // list to store jobs and their result
      storage::List
      <
        storage::Pair
        <
          util::ShPtr< sched::JobInterface>,
          storage::Triplet< assemble::ProteinModel, util::ShPtr< storage::Pair< assemble::ProteinModel, double> >, util::SiPtr< storage::Pair< std::string, storage::Row< double> > > >
        >
      >
      jobs_result;

      // list for the minimzers used
      util::ShPtrList< mc::Approximator< assemble::ProteinModel, double> > minimzers;

      i = 0;
      for
      (
        storage::Table< double>::iterator table_itr( count_corr_rmsd_corrmin_rmsdmin.Begin()), table_itr_end( count_corr_rmsd_corrmin_rmsdmin.End());
        // only minimize to the max number of models that are to be minimized
        table_itr != table_itr_end && i < m_NumberSavebestParam->GetNumericalValue< size_t>();
        ++table_itr, ++i
      )
      {
        std::string &row_name( table_itr->First());
        const size_t number( size_t( table_itr->Second()[ "hashscorerank"]));
        BCL_MessageStd( "minimizing model with hashscorerank " + util::Format()( number));

        util::ShPtr< assemble::ProteinModel> transformed_model_ptr( PROTEIN_MODEL.HardCopy());
        assemble::ProteinModel &transformed_model( *transformed_model_ptr);
        transformed_model.Transform( transforms( number)->First());

        // create pdb file with initial fit
        const std::string key( assemble::ProteinStorageFile::KeyToString( i));
        const bool store_success( m_Storage->Store( transformed_model, GetSourceInitial(), key));
        BCL_Assert( store_success, "unable to store initial model to " + GetSourceInitial() + " with key: " + key);
        BCL_MessageStd( "wrote initial fit " + util::Format()( number) + " to " + GetSourceInitial() + " with key: " + key);
        row_name = GetSourceInitial() + key;

        // if no minimization is required got to next model
        if( m_McMaxIterationsParam->GetNumericalValue< size_t>() == 0)
        {
          continue;
        }

        // create new job trans_center result triplet
        jobs_result.PushBack
        (
          storage::Pair
          <
            util::ShPtr< sched::JobInterface>,
            storage::Triplet
            <
              assemble::ProteinModel,
              util::ShPtr< storage::Pair< assemble::ProteinModel, double> >,
              util::SiPtr< storage::Pair< std::string, storage::Row< double> > >
            >
          >()
        );

        // create reference to the just create triplet
        storage::Pair
        <
          util::ShPtr< sched::JobInterface>,
          storage::Triplet< assemble::ProteinModel, util::ShPtr< storage::Pair< assemble::ProteinModel, double> >, util::SiPtr< storage::Pair< std::string, storage::Row< double> > > >
        > &current_job_result( jobs_result.LastElement());

        // transform the model into the initial fit position
        // this is a very complicated way of hard copying the transformed model for threaded use
        // one could just retrieve it from m_Storage, but the min sse size is not necessarily 0 0 0
        {
          const pdb::Factory tmp_factory( biol::GetAAClasses().e_AAComplete);
          storage::Map< biol::SSType, size_t> min_sse_size;
          min_sse_size[ biol::GetSSTypes().HELIX]  = 0;
          min_sse_size[ biol::GetSSTypes().STRAND] = 0;
          min_sse_size[ biol::GetSSTypes().COIL]   = 0;
          std::stringstream pdb_stream;
          tmp_factory.WriteModelToPDB( transformed_model, pdb_stream);
          pdb::Handler handler( pdb_stream);
          current_job_result.Second().First() = tmp_factory.ProteinModelFromPDB( handler, min_sse_size);
        }

        // initialize result to start protein model
        current_job_result.Second().Second() = util::ShPtr< storage::Pair< assemble::ProteinModel, double> >( new storage::Pair< assemble::ProteinModel, double>);
        current_job_result.Second().Second()->First() = current_job_result.Second().First();

        // point to row
        current_job_result.Second().Third() = util::ToSiPtrNonConst( ( *table_itr));

        // instantiate the approximator
        util::ShPtr< mc::Approximator< assemble::ProteinModel, double> > approximator
        (
          CreateMinimizer( i)
        );
        minimzers.PushBack( approximator);

        // initialize the actual job object
        current_job_result.First() =
        util::ShPtr< sched::JobInterface>
        (
          new sched::UnaryFunctionJobWithData
          <
            const assemble::ProteinModel,
            util::ShPtr< storage::Pair< assemble::ProteinModel, double> >,
            mc::Approximator< assemble::ProteinModel, double>
          >
          (
            0,
            *approximator,
            &mc::Approximator< assemble::ProteinModel, double>::ApproximateReturnCopy,
            current_job_result.Second().First(),
            sched::JobInterface::e_READY,
            &current_job_result.Second().Second()
          )
        );

        // submit job to the scheduler
        sched::GetScheduler().SubmitJob( current_job_result.First());
      }

      i = 0;
      // iterate over jobs and join them, writing the results
      for
      (
        storage::List
        <
          storage::Pair
          <
            util::ShPtr< sched::JobInterface>,
            storage::Triplet< assemble::ProteinModel, util::ShPtr< storage::Pair< assemble::ProteinModel, double> >, util::SiPtr< storage::Pair< std::string, storage::Row< double> > > >
          >
        >::iterator itr( jobs_result.Begin()), itr_end( jobs_result.End()); itr != itr_end;
        ++itr, ++i
      )
      {
        // join
        sched::GetScheduler().Join( itr->First());

        std::string &row_name( itr->Second().Third()->First());
        storage::Row< double> &current_row( itr->Second().Third()->Second());
        const assemble::ProteinModel &min_model( itr->Second().Second()->First());

        // apply transformations
        current_row[ "min_corr"] = -density_ccc->operator()( min_model);
        current_row[ "min_RMSD"] =
          assemble::Quality::Calculate
          (
            quality::GetMeasures().e_RMSD_NoSuperimposition, min_model, PROTEIN_MODEL, rmsd_atom_types
          );

        // write minimized pdb to file
        const std::string key( assemble::ProteinStorageFile::KeyToString( i));
        const bool store_success( m_Storage->Store( itr->Second().Second()->First(), GetSourceMinimized(), key));
        BCL_Assert( store_success, "unable to store minimized model to " + GetSourceMinimized() + " with key: " + key);
        BCL_MessageStd( "wrote minimized fit " + util::Format()( row_name) + " to " + GetSourceMinimized() + " with key: " + key);
        row_name += "_" + GetSourceMinimized() + key;
      }

      // end
      return count_corr_rmsd_corrmin_rmsdmin;
    }

    //! @brief initializes the command object for that application
    util::ShPtr< command::Command> FitInDensity::InitializeCommand() const
    {
      // new command line object
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      //input pdb and mrc file parameter
      sp_cmd->AddParameter( m_PDBFilenameParam);

      //input density filename
      sp_cmd->AddParameter( m_MRCFilenameParam);

      // density map resolution parameters
      sp_cmd->AddFlag( m_MrcResolutionFlag);

      // flag to choose the geometric hash storage
      sp_cmd->AddFlag( coord::GetGeometricHashStorageClasses().GetGeometricHashStorageClassFlag());

      // flag for protein storage
      sp_cmd->AddFlag( assemble::ProteinStorageFile::GetDefaultStorageFlag());

      // flag for prefix - corresponds to source in protein storage
      sp_cmd->AddFlag( m_PrefixFlag);

      // flag and parameter for the number of features used for representing the electron density map
      sp_cmd->AddFlag( m_NumberFeaturesFlag);

      // flag to give a multiple of the number of features determined by a formula using all the other parameters
      sp_cmd->AddFlag( m_MultipleOfExpectedPointsFlag);

      // flag and parameters for generating features from intensity and gradient in density map
      sp_cmd->AddFlag( m_RatioIntensityGradientFlag);

      // flag and parameter for the min distance between two feature
      sp_cmd->AddFlag( m_FeatureDistanceFlag);

      // noise
      sp_cmd->AddFlag( m_RemoveNoiseFlag);

      // write feature cloud
      sp_cmd->AddFlag( m_WriteFeatureCloudFlag);

      // flag for the threshold for the length of the sides of the triangular base in the geometric hashing algorithm
      sp_cmd->AddFlag( m_ThresholdFlag);

      // quantize
      sp_cmd->AddFlag( m_QuantizeFlag);

      // flag with parameter for the radius around the triangular base where coordinates are considered features
      sp_cmd->AddFlag( m_FeatureRadiusFlag);

      // flag and parameter for the coordinatesystem to be used for the geometric hashing
      sp_cmd->AddFlag( m_CoordinateSystemFlag);

      // flag and parameter for resolution
      sp_cmd->AddFlag( m_QuantizationResolutionFlag);

      // sses from backbone if no sse is given in pdb
      sp_cmd->AddFlag( pdb::Factory::GetFlagSSEsFromBackBone());

      // pdb reader min sse size
      pdb::Factory::GetFlagMinSSESize()->GetParameterList()( biol::GetSSTypes().COIL)->SetDefaultParameter( "999");
      sp_cmd->AddFlag( pdb::Factory::GetFlagMinSSESize());

      // pdb reader default amino acid
      pdb::Factory::GetFlagAAClass()->GetParameterList()( 0)->SetDefaultParameter( biol::GetAAClasses().e_AAComplete.GetName());
      sp_cmd->AddFlag( pdb::Factory::GetFlagAAClass());

      // flag for a dynamic list of atoms to be fitted
      sp_cmd->AddFlag( m_AtomListFlag);

      // flag to write minimization
      sp_cmd->AddFlag( m_WriteMinimizationFlag);
      sp_cmd->AddFlag( m_WriteSharedMemoryFlag);

      // flag for choice of density agreement
      sp_cmd->AddFlag( m_DensityAgreementFlag);

      // flag for simulator
      sp_cmd->AddFlag( m_SimulatorFlag);

      // flag and parameter for number of trials and number of best results
      sp_cmd->AddFlag( m_TrialsSavebestFlag);

      // flag and two parameters for the allowed difference in rotation and translation between two hits
      sp_cmd->AddFlag( m_DiffRotTransFlag);

      // flags for monte carlo minimization
      sp_cmd->AddFlag( m_McMaxIterationsUnimprovedFlag);

      // adjust start and end condition
      sp_cmd->AddFlag( mc::TemperatureAccepted::GetFlagTemperature());

      sp_cmd->AddFlag( m_McMutateTransRotFlag);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief initialize PointCloud - try to read from file, otherwise call the appropriate function
    //! @param FILENAME_PREFIX prefix for the file read from or written to
    //! @param NUMBER_FEATURES the number of features that should be in the density map
    //! @return feature cloud either read from file if existent or generated from static parameters from command line
    coord::PointCloud FitInDensity::InitializePointCloud
    (
      const std::string &FILENAME_PREFIX,
      const double FEATURE_DISTANCE,
      const size_t NUMBER_FEATURES
    ) const
    {
      io::IFStream read;

      // try to read feature cloud from file
      std::string feature_cloud_filename;
      feature_cloud_filename = m_PrefixFlag->GetFirstParameter()->GetValue() + FILENAME_PREFIX +
        "_features_" + util::Format()( NUMBER_FEATURES) + "_mindist_" +
        util::Format()( FEATURE_DISTANCE) + "A" + "_ratio_" +
        m_RatioIntensityGradientFlag->GetFirstParameter()->GetValue() + "_neighbors_" +
        m_NumberNeighbors->GetValue() + "_dist_" + m_MinNeighborDistance->GetValue() + "A" + ".feature_cloud";
      coord::PointCloud feature_cloud;

      if( io::File::TryOpenIFStream( read, feature_cloud_filename))
      {
        BCL_MessageStd( "read point cloud coordinates from file " + feature_cloud_filename);
        read >> feature_cloud;
        io::File::CloseClearFStream( read);
      }
      else
      {
        // construct feature cloud from density map containing all features over a certain density value
        {
          BCL_MessageStd( "calculate point cloud started");
          util::Stopwatch stopwatch;
          feature_cloud = m_SpMap->CalculatePointCloud
                          (
                            NUMBER_FEATURES,
                            FEATURE_DISTANCE,
                            m_RatioIntensityGradientFlag->GetFirstParameter()->GetNumericalValue< double>()
                          );
          BCL_MessageStd
          (
            "calculate point cloud finished - extracted " + util::Format()( feature_cloud.GetSize()) + " features!"
          );
        }

        // removing noise by removing features form feature cloud that do not have specified amount of neighbors within
        // a certain distance
        if( m_RemoveNoiseFlag->GetFlag())
        {
          const size_t number_removed_features
          (
            feature_cloud.RemoveSingles
            (
              m_NumberNeighbors->GetNumericalValue< size_t>(),
              m_MinNeighborDistance->GetNumericalValue< double>()
            )
          );
          BCL_MessageStd
          (
            "removed " + util::Format()( number_removed_features) + " since they did not have " +
            m_NumberNeighbors->GetValue() + " neighbors within a distance of " +
            m_MinNeighborDistance->GetValue() + "A"
          );
        }

        // write the feature cloud to file
        if( m_WriteFeatureCloudFlag->GetFlag())
        {
          BCL_MessageStd( "write feature cloud coordinates in file " + feature_cloud_filename);
          io::OFStream write;
          io::File::MustOpenOFStream( write, feature_cloud_filename);
          write << feature_cloud;
          io::File::CloseClearFStream( write);
          //write feature cloud to pdb file
          feature_cloud_filename += ".pdb";
          feature_cloud.WriteToPDB( feature_cloud_filename);
        }
      }

      //end
      return feature_cloud;
    } // end InitializePointCloud

    //! @brief quantizes a feature cloud
    //! @param POINTCLOUD unquantized point cloud from density or pdb structure
    //! @param POINT_TO_KEY function to convert into a key
    //! @return feature cloud, containing the quantized feature cloud depending on the chosen parameteres
    coord::PointCloud FitInDensity::QuantizeFeatureCloud
    (
      const coord::PointCloud &POINTCLOUD,
      const coord::PointToKeyInterface &POINT_TO_KEY
    ) const
    {
      // feature cloude for qunatized features
      coord::PointCloud quantized_features;
      quantized_features.GetData().AllocateMemory( POINTCLOUD.GetSize());

      // iterate over argument point ccoud, qunatize poitns and insert them in the quantized_features
      for
      (
        coord::PointCloud::const_iterator point_itr( POINTCLOUD.Begin()), point_itr_end( POINTCLOUD.End());
          point_itr != point_itr_end;
        ++point_itr
      )
      {
        // convert to key and convert back to cartesian coordinate and insert
        quantized_features.PushBack( POINT_TO_KEY( POINT_TO_KEY( *point_itr)));
      }

      // end
      return quantized_features;
    }

    //! default constructor
    FitInDensity::FitInDensity() :
      m_PDBFilenameParam
      (
        new command::Parameter
            (
              "pdb_filename",
              "\tfilename for input pdb to be fitted in mrc electron density map",
              command::ParameterCheckExtension( ".pdb")
            )
      ),
      m_MRCFilenameParam
      (
        new command::Parameter
            (
              "mrc_filename",
              "\tfilename for input mrc to be used for fitting",
              command::ParameterCheckExtension( ".mrc")
            )
      ),
      m_PrefixFlag
      (
        new command::FlagStatic
            (
              "prefix",
              "prefix for written files",
              command::Parameter
              (
                "output_prefix",
                "prefix added in front of filenames for output",
                ""
              )
            )
      ),
      m_QuantizationResolutionFlag
      (
        new command::FlagStatic
            (
              "quantization_resolution",
              "\t\tresolution in which features are quantized for the geometric hash"
            )
      ),
      m_QuantizationResolutionAngularParam
      (
        new command::Parameter
        (
          "angular_resolution_value",
          "\tchoose the angular resolution for the features encoded in hash",
          command::ParameterCheckRanged< double>( 0, 720),
          "12.0"
        )
      ),
      m_QuantizationResolutionDistanceParam
      (
        new command::Parameter
        (
          "distance_resolution_value",
          "\tchoose the distance resolution for the features encoded in hash",
          command::ParameterCheckRanged< double>( 0, 100),
          "2.0"
        )
      ),
      m_MrcResolutionFlag
      (
        new command::FlagStatic
            (
              "mrc_resolution",
              "\tresolution of given density map",
              command::Parameter
              (
                "density map resolution",
                "\tresolution of given density map [A] - will be used to simulate density and calculate correlation",
                command::ParameterCheckRanged< double>( 0.0, 100.0)
              )
            )
      ),
      m_NumberFeaturesFlag
      (
        new command::FlagStatic
            (
              "number_features",
              "\t\tchoose the number of features to be extracted from the density map",
              command::Parameter
              (
                "number_features",
                "\tchoose the number of features",
                command::ParameterCheckRanged< size_t>( 0, 10000),
                "0"
              )
            )
      ),
      m_MultipleOfExpectedPointsFlag
      (
        new command::FlagStatic
            (
              "multiple_features",
              "\t\tchoose a multiple of the number of features to be extracted from the density map as determined by "
                "the number of points to be fitted and all the other parameters",
              command::Parameter
              (
                "mulitple_features",
                "\tchoose the a multiple of features",
                command::ParameterCheckRanged< double>( 0, 100),
                "1.0"
              )
            )
      ),
      m_FeatureDistanceFlag
      (
        new command::FlagStatic
        (
          "feature_distance",
          "\t\tthis is the minimal distance of two features representing the density map",
          command::Parameter
          (
            "feature_distance",
            "\tminimal distance [A] between two features in cloud",
            command::ParameterCheckRanged< double>( 1, 100),
            "3.8"
          )
        )
      ),
      m_RatioIntensityGradientFlag
      (
        new command::FlagStatic
            (
              "ratio_int_grad",
              "\t\tthis is the ratio between the intensity vs. gradient to assign a feature within the density map",
              command::Parameter( "ratio", "\tthe ratio intensity/gradient of a voxel", "1.0")
            )
      ),
      m_WriteFeatureCloudFlag
      (
        new command::FlagStatic
            (
              "write_feature_cloud",
              "\t\tuset this to write a file for the feature cloud - coordinate file and pdb file"
            )
      ),
      m_ThresholdFlag
      (
        new command::FlagDynamic
        (
          "threshold",
          "\t\t4 threshold [A] for the length of the sides of the triangular base in the geometric hashing algorithm, will de derived if not given",
          command::Parameter
          (
            "threshold [A]",
            "lower middle1, middle2, highest threshold",
            command::ParameterCheckRanged< double>( 0.0, 100.0),
            "0"
          ),
          0,
          4
        )
      ),
      m_FeatureRadiusFlag
      (
        new command::FlagStatic
        (
          "feature_radius",
          "\tradius around the base in which features are quantized and stored in hash",
          command::Parameter
          (
            "feature_radius", //flag name
            "\tradius around the base to include features", //description
            command::ParameterCheckRanged< double> //parameter restriction
            (
              0.0, 1000.0
            ),
            "57.0" //default
          )
        )
      ),
      m_CoordinateSystemFlag
      (
        new command::FlagStatic
        (
          "coordinatesystem",
          "\tcoordinatesystem for the quantization for the geometric hash",
          command::Parameter
          (
            "coordinatesystem", //parameter name
            "\ta choice of possible coordinate systems for quantizing features in the geometric hash", //description
            command::ParameterCheckEnumerate< coord::PointToKeyClasses>(),
            "Spherical"
//            coord::GetPointToKeyClasses().e_SphericalRadius //default
          )
        )
      ),
      m_TrialsSavebestFlag
      (
        new command::FlagStatic
            (
              "number_bases_refine",
              "\tnumber of bases considered and initial fits refined"
            )
      ),
      m_NumberTrialsParam
      (
        new command::Parameter
            (
              "number_bases",
              "\tnumber of bases in the atomic structure that are used for the geometric hash fit",
              command::ParameterCheckRanged< size_t>( 1, 2500),
              "200"
            )
      ),
      m_NumberSavebestParam
      (
         new command::Parameter
             (
               "number_refine",
               "\tnumber of best initial fits by hash score to be refined",
               command::ParameterCheckRanged< size_t>( 1, 100),
               "10"
             )
      ),
      m_DiffRotTransFlag
      (
        new command::FlagStatic
            (
              "diff_rot_trans",
              "\tallowed difference in rotation and translation between two fittings, if they are within those ranges, "
                "keep the one with better hash score"
            )
      ),
      m_DiffRotParam
      (
        new command::Parameter
            (
              "diff_rot",
              "\t\teffective rotation angle between two orientations [rad]",
              command::ParameterCheckRanged< double>( 0.0, math::g_Pi),
              "0.0"
            )
      ),
      m_DiffTransParam
      (
        new command::Parameter
            (
              "diff_trans",
              "\tdifference in translation in Angstroem",
              command::ParameterCheckRanged< double>( 0.0, 1000.0),
              "0.0"
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
                command::ParameterCheckAllowed( biol::GetAtomTypes().GetBackBoneAtomNames())
              ),
              0,
              biol::GetAtomTypes().GetBackBoneAtomNames().GetSize()
            )
      ),
      m_WriteMinimizationFlag
      (
        new command::FlagDynamic
        (
          "write_minimization", "write the minimization to files - only of given step statuses",
          command::Parameter
          (
            "stepstatuses",
            "any step status from the list",
            command::ParameterCheckSerializable( opti::StepStatusEnum()),
            opti::GetStepStatusName( opti::e_Improved)
          ),
          0,
          opti::s_NumberStepStatus
        )
      ),
      m_WriteSharedMemoryFlag
      (
        new command::FlagDynamic
        (
          "shared_memory", "write the minimization to shared memory - only of given step statuses",
          command::Parameter
          (
            "stepstatuses",
            "any step status from the list",
            command::ParameterCheckAllowed( opti::StepStatusEnum::GetStringVector())
          ),
          0,
          opti::s_NumberStepStatus
        )
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
      m_McMaxIterationsUnimprovedFlag
      (
        new command::FlagStatic
            (
              "mc_max_unimproved_steps",
              "\tmodify the number of total steps and max number of steps in a row without improvement"
            )
      ),
      m_McMaxIterationsParam
      (
        new command::Parameter
            (
              "max_iterations",
              "\tmaximal number of steps with improvement",
              command::ParameterCheckRanged< size_t>( 0, 1000),
              "250"
            )
      ),
      m_McMaxStepsUnimprovedParam
      (
        new command::Parameter
            (
               "steps_without_improvement",
               "\tmax number of steps without improvement before terminating",
               command::ParameterCheckRanged< size_t>( 0, 1000),
               "50"
            )
      ),
      m_McMutateTransRotFlag
      (
        new command::FlagStatic
        (
          "mc_mutate_trans_rot",
          "\tmaximal mutation step sizes for monte carlo refinement"
        )
      ),
      m_McMutateTransParam
      (
        new command::Parameter
        (
          "translation",
          "\tmaximal translational mutation [A]",
          command::ParameterCheckRanged< double>( 0.0, 10.0), // restriction
          "1.0" // default
        )
      ),
      m_McMutateRotParam
      (
        new command::Parameter
        (
          "rotation",
          "\t\tmaximal rotational mutation [rad]",
          command::ParameterCheckRanged< double>( 0, math::g_Pi), // restriction
          util::Format()( 2.0 / 180 * math::g_Pi) // default 2 degree
        )
      ),
      m_QuantizeFlag
      (
        new command::FlagStatic
        (
          "quantize",
          "\tquantize the centered pdb and write to file"
        )
      ),
      m_CenterXParam
      (
        new command::Parameter
        (
          "center x",
          "x component of translation vector to center",
          "0.0" //default
        )
      ),
      m_CenterYParam
      (
        new command::Parameter
        (
          "center y",
          "y component of translation vector to center",
          "0.0" //default
        )
      ),
      m_CenterZParam
      (
        new command::Parameter
        (
          "center z",
          "z component of translation vector to center",
          "0.0" //default
        )
      ),
      m_RemoveNoiseFlag
      (
        new command::FlagStatic
        (
          "remove_noise",
          "remove features that do not have certain amount of neighbors within a certain distance"
        )
      ),
      m_NumberNeighbors
      (
        new command::Parameter
        (
          "number neighbors",
          "minimal number of neighbors within given distance - otherwise points gets removed",
          "0"
        )
      ),
      m_MinNeighborDistance
      (
        new command::Parameter
        (
          "neighbor distance [A]",
          "distance in which a certain amount of neighbors should be in Angstroem",
          "0.0" // default
        )
      )
    {
      // attach parameters to flags
      m_QuantizeFlag->PushBack( m_CenterXParam);
      m_QuantizeFlag->PushBack( m_CenterYParam);
      m_QuantizeFlag->PushBack( m_CenterZParam);
      m_QuantizationResolutionFlag->PushBack( m_QuantizationResolutionAngularParam);
      m_QuantizationResolutionFlag->PushBack( m_QuantizationResolutionDistanceParam);
      m_TrialsSavebestFlag->PushBack( m_NumberTrialsParam);
      m_TrialsSavebestFlag->PushBack( m_NumberSavebestParam);
      m_DiffRotTransFlag->PushBack( m_DiffRotParam);
      m_DiffRotTransFlag->PushBack( m_DiffTransParam);

      // max number of rejected steps and max iterations
      m_McMaxIterationsUnimprovedFlag->PushBack( m_McMaxIterationsParam);
      m_McMaxIterationsUnimprovedFlag->PushBack( m_McMaxStepsUnimprovedParam);

      // mutate parameters
      m_McMutateTransRotFlag->PushBack( m_McMutateTransParam);
      m_McMutateTransRotFlag->PushBack( m_McMutateRotParam);

      // noise parameters
      m_RemoveNoiseFlag->PushBack( m_NumberNeighbors);
      m_RemoveNoiseFlag->PushBack( m_MinNeighborDistance);
    }

    //! instance of the application
    const ApplicationType FitInDensity::FitInDensity_Instance
    (
      GetAppGroups().AddAppToGroup( new FitInDensity(), GetAppGroups().e_Density)
    );

  } // namespace app
} // namespace bcl
