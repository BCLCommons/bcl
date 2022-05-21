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
#include "assemble/bcl_assemble_collector_topology_sheet.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_collector_topology_sheet.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleCollectorTopologySheet :
    public ExampleInterface
  {
  public:

    ExampleAssembleCollectorTopologySheet *Clone() const
    {
      return new ExampleAssembleCollectorTopologySheet( *this);
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
      // initialize pdb filename
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "2CRT.pdb"));
      // get the protein model
      assemble::ProteinModel model( Proteins::GetModel( pdb_filename));
      // join the adjacent strand if any exists
      model.Join( biol::GetSSTypes().STRAND, false);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct a collector
      assemble::CollectorTopologySheet collector;

      // test copy constructor
      BCL_MessageStd( "Testing copy constructor");
      assemble::CollectorTopologySheet collector_copy( collector);

      // test clone
      BCL_MessageStd( "Testing clone function");
      util::ShPtr< assemble::CollectorTopologySheet> sp_collector( collector.Clone());

    /////////////////
    // data access //
    /////////////////

      // example check
      BCL_ExampleCheck( collector.GetClassIdentifier(), GetStaticClassName( collector));

      // check the static variable s_MaximumDistance
      BCL_ExampleCheck( assemble::CollectorTopologySheet::s_MaximumDistance, 7.0);

      // check the static variable s_MinimumStrandStrandPairingWeight
      BCL_ExampleCheck( assemble::CollectorTopologySheet::s_MinimumStrandStrandPairingWeight, 0.5);

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleCollectorTopologySheet

  const ExampleClass::EnumType ExampleAssembleCollectorTopologySheet::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleCollectorTopologySheet())
  );

} // namespace bcl
