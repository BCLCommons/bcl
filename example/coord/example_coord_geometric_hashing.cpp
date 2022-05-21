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

// include example header
#include "example.h"
#include "example_proteins.h"
// include the header of the class which this example is for
#include "coord/bcl_coord_geometric_hashing.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_quality.h"
#include "biol/bcl_biol_atom.h"
#include "coord/bcl_coord_geometric_hash_storage_hash_map.h"
#include "coord/bcl_coord_point_cloud.h"
#include "coord/bcl_coord_point_to_key_classes.h"
#include "density/bcl_density_map.h"
#include "density/bcl_density_simulators.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_histogram.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_coord_geometric_hashing.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCoordGeometricHashing :
    public ExampleInterface
  {
  public:

    ExampleCoordGeometricHashing *Clone() const
    { return new ExampleCoordGeometricHashing( *this);}

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      //initialize filenames
      const std::string mrc_filename( AddExampleInputPathToFilename( e_Biology, "1ubi_res_6.6voxelsize_2.200Gaussian.mrc"));
      // create string "pdb_filename" which has path for example pdb file
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      std::string hashmap_filename;

      const size_t resolution( 12);
      const size_t number_of_points( 47);
      const double mindist( 3.8);
      const storage::VectorND< 2, double> threshold( 7.6, 22.9495);
      const double feature_radius( 17.2122);
      const size_t savebest( 5);
      const size_t trials( 200);
      const storage::VectorND< 2, double> difference( 0, 0);
      const storage::Set< biol::AtomType> atom_types( biol::GetAtomTypes().CA);

      BCL_MessageStd( "processing screening with following arguments");
      BCL_MessageStd
      (
        "-pdb " + pdb_filename
        + "\n-mrc "                      + util::Format()( mrc_filename)
        + "\n-resolution "               + util::Format()( resolution)
        + "\n-points "                   + util::Format()( number_of_points)
        + "\n-mindistance "              + util::Format()( mindist)
        + "\n-minthreshold "             + util::Format()( threshold.First())
        + "\n-maxthreshold "             + util::Format()( threshold.Second())
        + "\n-feature_radius "           + util::Format()( feature_radius)
        + "\n-savebest "                 + util::Format()( savebest)
        + "\n-trials "                   + util::Format()( trials)
        + "\n-diffrot "                  + util::Format()( difference.First())
        + "\n-difftrans "                + util::Format()( difference.Second())
        + "\n-matching following atoms " + util::Format()( atom_types)
      );

      //instantiate DensityMap from mrc file
      io::IFStream read;
      BCL_MessageStd( "read DensityMap from mrc file");
      BCL_ExampleMustOpenBinaryInputFile( read, mrc_filename);
      density::Map mrc;
      mrc.ReadMRC( read);
      io::File::CloseClearFStream( read);
      BCL_MessageStd( "read DensityMap from mrc file done");

      // create density simulator
      const util::ShPtr< density::SimulateInterface> simulator
      (
        density::GetSimulators().CreateSimulator
        (
          density::GetSimulators().e_Gaussian,
          mrc.GetCellWidth(),
          3 * mrc.GetCellWidth().X()
        )
      );

      // construct point cloud from DensityMap containing all points over a certain density value
      coord::PointCloud pointcloud( mrc.CalculatePointCloud( number_of_points, mindist, 1.0));

      // calculate the distance distribution
      const math::Histogram distance_distribution
      (
        coord::GeometricHashing::CalculateDistanceDistribution
        (
          util::ConvertToConstSiPtrVector< linal::Vector3D>( pointcloud), threshold, 20
        )
      );
      BCL_MessageStd
      (
        "this is the distribution of distances within the pointcloud: " + util::Format()( distance_distribution)
      );

      // determine thresholds, so that they are equally occupied with distance
      const storage::VectorND< 4, double> thresholds_equal_occupied( coord::GeometricHashing::CreateEqualOccupiedIntervals( distance_distribution));
      BCL_MessageStd( "thresholds that contain equal amounts of distances: " + util::Format()( thresholds_equal_occupied));

      // write the point cloud to file
      io::OFStream write;
      const std::string pointcloud_filename( AddExampleOutputPathToFilename( pointcloud, "points_both11_" + util::Format()( number_of_points) + "_threshold_" + util::Format().W( 3).FFP( 1)( mindist) + "A.pointcloud"));
      BCL_MessageStd( "write point cloud coordinates in file " + pointcloud_filename);
      BCL_ExampleMustOpenOutputFile( write, pointcloud_filename);

      write << pointcloud;
      io::File::CloseClearFStream( write);

      // build hash with normal point cloud
      BCL_MessageStd
      (
        "building hash with threshold range: " + util::Format()( threshold.First()) + " ... " +
        util::Format()( threshold.Second()) + " and resolution factor: " + util::Format()( resolution) +
        " - can take a little time"
      );

      if( !coord::GetPointToKeyClasses().HaveEnumWithName( "Spherical"))
      {
        BCL_MessageCrt
        (
          GetClassIdentifier() + " not tested completely because Spherical was not available"
        );
        return 0;
      }

      //construct hashmap and write to file
      coord::GeometricHashing hashmap
      (
        util::ShPtr< coord::GeometricHashStorageInterface>
        (
          new coord::GeometricHashStorageHashMap
          (
            mrc_filename,
            6.9,
            mindist,
            mrc.GetDimensions(),
            number_of_points,
            mindist,
            1.0,
            thresholds_equal_occupied,
            feature_radius,
            *coord::GetPointToKeyClasses().GetPointToKeyFunction
            (
              "Spherical", resolution, mindist
            ),
            0,
            0.0
          )
        )
      );
      hashmap.BuildHash( pointcloud);

      storage::Map< biol::SSType, size_t> sse_min_size;
      sse_min_size[ biol::GetSSTypes().HELIX] = 0;
      sse_min_size[ biol::GetSSTypes().STRAND] = 0;
      const assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AAComplete, sse_min_size)
      );

      BCL_MessageStd( "matching: " + pdb_filename);

      //extract all atoms for later transformation and correlation calculation
      coord::PointCloud match_pdb_positions( util::ConvertToStorageVector< linal::Vector3D, const linal::Vector3D>( protein_model.GetAtomCoordinates( atom_types)));

      //all atom coordinates

      BCL_MessageStd( "number of atoms to match: " + util::Format()( match_pdb_positions.GetSize()));

      //search hash and return best pairs( with highest count in Hashmap) of matching Transformationmatrices
      const storage::List< storage::Pair< math::TransformationMatrix3D, size_t> >
        transforms( hashmap.SearchTarget( match_pdb_positions, savebest, trials, difference));

      //calculate correlation factors
      linal::Vector< double> corrfactors( transforms.GetSize());
      //calculate correlation factors
      linal::Vector< double> rmsd( transforms.GetSize());

      sse_min_size[ biol::GetSSTypes().COIL] = 0;
      size_t i( 0);
      // transform all coordinates and correlate them for search hashmap1
      for
      (
        storage::List< storage::Pair< math::TransformationMatrix3D, size_t> >::const_iterator
          itr( transforms.Begin()), itr_end( transforms.End());
        itr != itr_end;
        ++itr, ++i
      )
      {
        assemble::ProteinModel current_model
        (
          Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AAComplete, sse_min_size)
        );

        // apply transformations
        current_model.Transform( itr->First());
        corrfactors( i) = mrc.Correlation( simulator->operator ()( current_model.GetAtoms()));
        rmsd( i) =
            assemble::Quality::Calculate
            (
              quality::GetMeasures().e_RMSD_NoSuperimposition,
              current_model, biol::GetAtomTypes().GetBackBoneAtomTypes()
            );

        const std::string filename( AddExampleOutputPathToFilename( hashmap, "transformed" + util::Format()( i) + ".pdb"));
        BCL_MessageStd( "hash score: " + util::Format()( itr->Second()));
        BCL_MessageStd( "write transformed coordinates in file " + filename);
        Proteins::WriteModelToPDB( current_model, filename);
      }

      if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Standard))
      {
        BCL_MessageStd( "CorrelationFactors");
        io::Serialize::Write( corrfactors, util::GetLogger(), 0, util::Format().W( 12).FFP( 9)) << '\n';
        BCL_MessageStd( "RMSD");
        io::Serialize::Write( rmsd, util::GetLogger(), 0, util::Format().W( 7).FFP( 3)) << '\n';
      }

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // end ExampleCoordGeometricHashing

  const ExampleClass::EnumType ExampleCoordGeometricHashing::s_Instance
  (
    GetExamples().AddEnum( ExampleCoordGeometricHashing())
  );

} // namespace bcl
