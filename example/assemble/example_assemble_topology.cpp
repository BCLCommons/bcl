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
#include "assemble/bcl_assemble_topology.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_collector_topology_combined.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_topology.cpp
  //!
  //! @author nobody
  //! @date
  //! @remarks status empty
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleTopology :
    public ExampleInterface
  {
  public:

    ExampleAssembleTopology *Clone() const
    {
      return new ExampleAssembleTopology( *this);
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
      BCL_MessageStd( "Reading native structure and perturbed models");

      // initialize stream and pdb factory and minsse_sizes
      storage::Map< biol::SSType, size_t> sse_min_sizes;
      sse_min_sizes[ biol::GetSSTypes().HELIX] = 5;
      sse_min_sizes[ biol::GetSSTypes().STRAND] = 3;
      sse_min_sizes[ biol::GetSSTypes().COIL] = 999;

      // initialize native structure pdb filename
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
      //build models from pdb
      BCL_MessageStd( "building model from pdb file");
      assemble::ProteinModel native_structure
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, sse_min_sizes)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      assemble::Topology protein_topology
      (
        assemble::CollectorTopologyCombined().CalculateTopology( native_structure.GetSSEs())
      );

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, AddExampleOutputPathToFilename( protein_topology, "1ubi_topology.viz"));
      //protein_topology.WriteGraphVizScript( write, util::SiPtr< const math::FunctionInterfaceSerializable< assemble::SSEGeometryPacking, double> >());
      io::File::CloseClearFStream( write);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleTopology

  const ExampleClass::EnumType ExampleAssembleTopology::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleTopology())
  );

} // namespace bcl
