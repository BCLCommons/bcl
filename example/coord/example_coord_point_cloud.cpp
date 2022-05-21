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
#include "coord/bcl_coord_point_cloud.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom.h"
#include "density/bcl_density_map.h"
#include "density/bcl_density_simulators.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_coord_point_cloud.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCoordPointCloud :
    public ExampleInterface
  {
  public:

    ExampleCoordPointCloud *Clone() const
    { return new ExampleCoordPointCloud( *this);}

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
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      BCL_MessageStd( "reading pdb: " + pdb_filename);

      //set parameters
      const double resolution( 6.6), voxelsize( 2.2);
      const util::ShPtr< density::SimulateInterface> simulator( density::GetSimulators().CreateSimulator( density::GetSimulators().e_Gaussian, linal::Vector3D( voxelsize), resolution));

      //read atoms from pdb
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, pdb_filename);

      const pdb::Handler pdb( read);
      io::File::CloseClearFStream( read);

      //extract all atoms for conversion to densitymap
      const util::ShPtrList< pdb::Line> atom_lines( pdb.GetLines( pdb::GetLineTypes().ATOM));
      util::ShPtrVector< biol::Atom> atoms;
      for
      (
        util::ShPtrList< pdb::Line>::const_iterator
          line_itr( atom_lines.Begin()), line_itr_end( atom_lines.End());
        line_itr != line_itr_end; ++line_itr
      )
      {
        atoms.PushBack( pdb::Factory( biol::GetAAClasses().e_AABackBone).AtomFromLine( **line_itr));
      }

      BCL_MessageStd( "atoms in pdb: " + util::Format()( atoms.GetSize()));

      //density map calculated from SimpleAtoms, with target resolution, voxelsize( usually sigma = resolution/3) and smoothingkernel of choice
      const density::Map density_map( simulator->operator ()( atoms));

      // point cloud from densitymap
      coord::PointCloud point_cloud( density_map.CalculatePointCloud( 200, 2.0, 1.0));

      // write point_cloud to file
      point_cloud.WriteToPDB( AddExampleOutputPathToFilename( point_cloud, "pointcloud_200.pdb"));

      // remove all points, that do not have at least 2 neighbors within a range of 3
      BCL_MessageStd( "removed " + util::Format()( point_cloud.RemoveSingles( 2, 3.0)) +
        " points that did have less than 2 neighbors within a range of 3.0 A");

      // write thinner point_cloud to file
      point_cloud.WriteToPDB( AddExampleOutputPathToFilename( point_cloud, "pointcloud_200_without_singles.pdb"));

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // end ExampleCoordPointCloud

  const ExampleClass::EnumType ExampleCoordPointCloud::s_Instance
  (
    GetExamples().AddEnum( ExampleCoordPointCloud())
  );

} // namespace bcl
