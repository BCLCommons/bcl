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
#include "density/bcl_density_simulate_gaussian_sphere.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_atom.h"
#include "density/bcl_density_map.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_density_simulate_gaussian_sphere.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDensitySimulateGaussianSphere :
    public ExampleInterface
  {
  public:

    ExampleDensitySimulateGaussianSphere *Clone() const
    {
      return new ExampleDensitySimulateGaussianSphere( *this);
    }

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
      const double resolution( 6.9);
      const linal::Vector3D grid_spacing( resolution / 3);

      const std::string example_pdb( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      // preparation
      // create pdb model
      storage::Map< biol::SSType, size_t> sse_min_size;
      sse_min_size[ biol::GetSSTypes().HELIX] = 0;
      sse_min_size[ biol::GetSSTypes().STRAND] = 0;
      sse_min_size[ biol::GetSSTypes().COIL] = 0;
      assemble::ProteinModel protein( Proteins::GetModel( example_pdb, biol::GetAAClasses().e_AAComplete, sse_min_size));

      // all atoms
      const util::SiPtrVector< const biol::Atom> atoms( protein.GetAtoms());

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct from grid spacing and resolution
      density::SimulateGaussianSphere simulate_gs( grid_spacing, resolution);

      // copy constructor
      density::SimulateGaussianSphere simulate_gs_copy( simulate_gs);

      // clone
      util::ShPtr< util::ObjectInterface> ptr( simulate_gs.Clone());

    /////////////////
    // data access //
    /////////////////

      // check the class identifier
      BCL_Example_Check
      (
        GetStaticClassName( simulate_gs) == "bcl::density::SimulateGaussianSphere",
        "unexpected static class name: " + GetStaticClassName( simulate_gs) + " should be: bcl::density::SimulateGaussianSphere"
      );

      // class identifier
      BCL_Example_Check
      (
        ptr->GetClassIdentifier() == GetStaticClassName< density::SimulateGaussianSphere>(),
        "unexpected class identifier class name: " + ptr->GetClassIdentifier() + " should be: " + GetStaticClassName< density::SimulateGaussianSphere>()
      );

    ///////////////
    // operators //
    ///////////////

      // simulate density map
      const density::Map map_sim( simulate_gs( atoms));

      // write simulated density map to file
      io::OFStream write;
      std::string mrc_filename( AddExampleOutputPathToFilename( map_sim, "1ubi.pdb"));
      mrc_filename = mrc_filename.substr( 0, mrc_filename.find( ".")) + "_gaussian_sphere_" + util::Format().W( 3).FFP( 1)( resolution) + ".mrc";
      BCL_ExampleMustOpenBinaryOutputFile( write, mrc_filename);
      map_sim.WriteMRC( write);
      io::File::CloseClearFStream( write);

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities for density::SimulateGaussianSphere");
      // open output file
      BCL_MessageVrb( "write object");
      WriteBCLObject( simulate_gs);
      BCL_MessageVrb( "read object");
      density::SimulateGaussianSphere simulate_gs_read( grid_spacing / 2, resolution / 2);
      ReadBCLObject( simulate_gs_read);

      // compare the maps, by comparing their mean and sd and cross correlation
      const density::Map map_from_read( simulate_gs_read( atoms));

      // compare the objects
      BCL_Example_Check
      (
        math::EqualWithinTolerance( map_sim.GetMean(), map_from_read.GetMean()) &&
        math::EqualWithinTolerance( map_sim.GetRmsd(), map_from_read.GetRmsd()) &&
        math::EqualWithinTolerance( 1.0, map_sim.Correlation( map_from_read)),
        "read simulated map different than original simulator map\n" +
        util::Format()( map_sim.GetMean()) + " != " + util::Format()( map_from_read.GetMean()) + " or " +
        util::Format()( map_sim.GetRmsd()) + " != " + util::Format()( map_from_read.GetRmsd()) + " or " +
        "ccc " + util::Format()( map_sim.Correlation( map_from_read)) + " != 1"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDensitySimulateGaussianSphere

  const ExampleClass::EnumType ExampleDensitySimulateGaussianSphere::s_Instance
  (
    GetExamples().AddEnum( ExampleDensitySimulateGaussianSphere())
  );

} // namespace bcl
