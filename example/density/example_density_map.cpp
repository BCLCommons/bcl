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
// include the header of the class which this example is for
#include "density/bcl_density_map.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_atom.h"
#include "coord/bcl_coord_point_cloud.h"
#include "density/bcl_density_simulators.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_density_map.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDensityMap :
    public ExampleInterface
  {
  public:

    ExampleDensityMap *Clone() const
    { return new ExampleDensityMap( *this);}

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
      const std::string example_pdb( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      // preparation
      // create pdb model
      storage::Map< biol::SSType, size_t> sse_min_size;
      sse_min_size[ biol::GetSSTypes().HELIX] = 0;
      sse_min_size[ biol::GetSSTypes().STRAND] = 0;
      sse_min_size[ biol::GetSSTypes().COIL] = 0;
      assemble::ProteinModel protein( Proteins::GetModel( example_pdb, biol::GetAAClasses().e_AAComplete, sse_min_size));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      density::Map density;

    /////////////////
    // data access //
    /////////////////

      //set parameters
      const double resolution( 6.6), voxelsize( 2.2);
      const util::ShPtr< density::SimulateInterface> simulator
      (
        density::GetSimulators().CreateSimulator
        (
          density::GetSimulators().e_Gaussian, linal::Vector3D( voxelsize), resolution
        )
      );

      // density map calculated from SimpleAtoms, with target resolution, voxelsize( usually sigma = resolution/3) and smoothingkernel of choice
      density = simulator->operator ()( protein.GetAtoms());

      // write simulated density map to file
      io::OFStream write;
      std::string mrc_filename( AddExampleOutputPathToFilename( density, "1ubi.pdb"));
      mrc_filename = io::File::RemoveLastExtension( mrc_filename) + "_res_" + util::Format().W( 3).FFP( 1)( resolution) + "voxelsize_" + util::Format().W( 5).FFP( 3)( voxelsize) + density::GetSimulators().e_Gaussian.GetName() + ".mrc";
      BCL_ExampleMustOpenBinaryOutputFile( write, mrc_filename);
      density.WriteMRC( write);
      io::File::CloseClearFStream( write);

      // instantiate Densitymap from mrc file
      density::Map mrc_read;
      BCL_MessageStd( "open density map: " + mrc_filename);
      io::IFStream read;
      BCL_ExampleMustOpenBinaryInputFile( read, mrc_filename);
      mrc_read.ReadMRC( read, 0);
      io::File::CloseClearFStream( read);

      // compare the written and read information
      BCL_Example_Check
      (
        density.GetAngles() == mrc_read.GetAngles()
        && density.GetAxis()      == mrc_read.GetAxis()
//        && density.GetCellWidth() == mrc_read.GetCellWidth()
//        && density.GetData() == mrc_read.GetData()
        && density.GetDimensions() == mrc_read.GetDimensions()
        && density.GetIndex() == mrc_read.GetIndex()
        && density.GetIntervals() == mrc_read.GetIntervals()
//        && density.GetMaximum() == mrc_read.GetMaximum()
//        && density.GetMean() == mrc_read.GetMean()
        && density.GetMinimum() == mrc_read.GetMinimum()
        && density.GetOrigin() == mrc_read.GetOrigin()
//        && density.GetRmsd() == mrc_read.GetRmsd()
        && density.GetSize() == mrc_read.GetSize()
//        && density.GetUnitCellLength() == mrc_read.GetUnitCellLength()
        ,
        "written and read map do not agree"
      );

      //this is a way to get a subdensity map from your density map by passing start indices and extension in each direction (in voxel)
      density::Map mrc1( density.SubMap( 5, 5, 5, 10, 10, 10));

      //write all header information to util::GetLogger()
      mrc1.WriteHeader( util::GetLogger());

      //write subdensity map to file
      BCL_MessageStd( "write subdensity map to file subdensity_555_101010.mrc");
      const std::string out_filename( AddExampleOutputPathToFilename( mrc1, "subdensity_555_101010.mrc"));
      BCL_ExampleMustOpenBinaryOutputFile( write, out_filename);
      mrc1.WriteMRC( write);
      io::File::CloseClearFStream( write);

      //construct pointcloud from density map containing all points over a certain density value
      coord::PointCloud pointcloud( density.CalculatePointCloud( 150, 1.0, 1.0));
      BCL_MessageStd( "Number of points in point cloud " + util::Format()( pointcloud.GetSize()));

      //calculate correlation between density map and simulated density map
      density = simulator->operator ()( protein.GetAtoms());
      const double ccc( density.Correlation( simulator->operator ()( protein.GetAtoms())));
      BCL_MessageStd( "correlation: " + util::Format()( ccc));

      // normalize
      BCL_MessageStd
      (
        " density. min=" + util::Format()( density.GetMinimum()) +
        " max=" + util::Format()( density.GetMaximum()) +
        " mean=" + util::Format()( density.GetMean()) +
        " rmsd=" + util::Format()( density.GetRmsd())
      );
      density.Normalize();
      BCL_MessageStd
      (
        " normalized density. min=" + util::Format()( density.GetMinimum()) +
        " max=" + util::Format()( density.GetMaximum()) +
        " mean=" + util::Format()( density.GetMean()) +
        " rmsd=" + util::Format()( density.GetRmsd())
      );

      // add noise
      const double cc_noise( density.AddNoise( random::GetGlobalRandom(), double( 0.0), 1.5 * density.GetRmsd()));
      density.Normalize();
      BCL_MessageStd
      (
        " added noise and normalized. min=" + util::Format()( density.GetMinimum()) +
        " max=" + util::Format()( density.GetMaximum()) +
        " mean=" + util::Format()( density.GetMean()) +
        " rmsd=" + util::Format()( density.GetRmsd())
      );
      BCL_MessageStd( "noise added to simulated density map. CCC=" + util::Format()( cc_noise));

      // reset filename
      mrc_filename = AddExampleOutputPathToFilename( density, "1ubi.pdb");
      mrc_filename = io::File::RemoveLastExtension( mrc_filename) + "_res_" + util::Format().W( 3).FFP( 1)( resolution) + "voxelsize_" + util::Format().W( 5).FFP( 3)( voxelsize) + density::GetSimulators().e_Gaussian.GetName() + "_noise.mrc";
      BCL_ExampleMustOpenBinaryOutputFile( write, mrc_filename);
      density.WriteMRC( write);
      io::File::CloseClearFStream( write);

    ////////////////
    // operations //
    ////////////////

      density = simulator->operator ()( protein.GetAtoms());
      density.EdgeDetectionSobel();
//      mrc_filename = "1ubi.pdb";
//      mrc_filename = mrc_filename.substr( 0, mrc_filename.find( ".")) + "_res_" + util::Format().W( 3, 1)( resolution) + "voxelsize_" + util::Format().W( 5).FFP( 3)( voxelsize) + density::g_SmoothingKernelsDescription[ kernel] + "_edge.mrc";
//      BCL_ExampleMustOpenBinaryOutputFile( write, AddExampleOutputPathToFilename( density, mrc_filename));
//      density.WriteMRC( write);
//      io::File::CloseClearFStream( write);

      // test the OrthogonalizeMap() function
      // instantiate filename of non-orthogonalized map
      std::string non_orthogonal_mrc_filename( AddExampleInputPathToFilename( e_Biology, "3kgv-sf_pdb_small.mrc"));

      // read in an experimental map that has angles 90, 105, 90
      density::Map mrc_non_orthogonal_read;
      BCL_MessageStd( "open density map: " + non_orthogonal_mrc_filename);
//      io::IFStream read;
      BCL_ExampleMustOpenBinaryInputFile( read, non_orthogonal_mrc_filename);
      mrc_non_orthogonal_read.ReadMRC( read, 0);
      io::File::CloseClearFStream( read);

      // orthogonalize map
      density::Map orthogonalized_map( mrc_non_orthogonal_read.OrthogonalizeMap());

      // write out the orthogonalized map
      std::string orthogonal_mrc_filename( AddExampleOutputPathToFilename( density, "orthogonalized_map.mrc"));
      BCL_ExampleMustOpenBinaryOutputFile( write, orthogonal_mrc_filename);
      orthogonalized_map.WriteMRC( write);
      io::File::CloseClearFStream( write);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // end ExampleDensityMap

  const ExampleClass::EnumType ExampleDensityMap::s_Instance
  (
    GetExamples().AddEnum( ExampleDensityMap())
  );
} // namespace bcl
