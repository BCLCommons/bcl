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
#include "density/bcl_density_simulators.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_atom.h"
#include "density/bcl_density_map.h"
#include "io/bcl_io_file.h"
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_density_simulators.cpp
  //! @details this example demonstrates how different simulators will simulate a density map from a list of atoms
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDensitySimulators :
    public ExampleInterface
  {
  public:

    ExampleDensitySimulators *Clone() const
    {
      return new ExampleDensitySimulators( *this);
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

      // no constructors available, since only singleton class

    /////////////////
    // data access //
    /////////////////

      // class identifier
      BCL_ExampleCheck( density::GetSimulators().GetClassIdentifier(), GetStaticClassName< density::Simulators>());

    ////////////////
    // operations //
    ////////////////

      // iterator over simulators and generate density maps
      for( density::Simulators::const_iterator itr( density::GetSimulators().Begin()), itr_end( density::GetSimulators().End()); itr != itr_end; ++itr)
      {
        if( !( *itr)->IsDefined())
        {
          continue;
        }
        util::Stopwatch stopwatch( "simulate " + itr->GetName(), util::Message::e_Standard, true);
        // copy the current simulator
        util::ShPtr< density::SimulateInterface> current_simulator( density::GetSimulators().CreateSimulator( *itr, grid_spacing, resolution));

        // simulate density map
        const density::Map map_sim( current_simulator->operator()( atoms));

        // write simulated density map to file
        io::OFStream write;
        std::string mrc_filename( AddExampleOutputPathToFilename( map_sim, "1ubi.pdb"));
        mrc_filename = mrc_filename.substr( 0, mrc_filename.find( ".")) + itr->GetName() + ".mrc";
        BCL_ExampleMustOpenBinaryOutputFile( write, mrc_filename);
        map_sim.WriteMRC( write);
        io::File::CloseClearFStream( write);

        // report min max and sd
        BCL_MessageStd
        (
          itr->GetName() + " min, max, mean, rmsd: " +
          util::Format()( map_sim.GetMinimum()) + " " +
          util::Format()( map_sim.GetMaximum()) + " " +
          util::Format()( map_sim.GetMean()) + " " +
          util::Format()( map_sim.GetRmsd())
        );
      }

    //////////////////////
    // input and output //
    //////////////////////

      // simulators cannot be written or read - except via the enum read and write mechanism

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // class ExampleDensitySimulators

  const ExampleClass::EnumType ExampleDensitySimulators::s_Instance
  (
    GetExamples().AddEnum( ExampleDensitySimulators())
  );

} // namespace bcl
