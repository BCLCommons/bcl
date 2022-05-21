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
#include "chemistry/bcl_chemistry_fragment_graph_marker.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_split_rings.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically
#include <stdio.h>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_fragment_graph_marker.cpp
  //!
  //! @author kothiwsk
  //! @date Oct 23, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryFragmentGraphMarker :
    public ExampleInterface
  {
  public:

    ExampleChemistryFragmentGraphMarker *Clone() const
    {
      return new ExampleChemistryFragmentGraphMarker( *this);
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

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      chemistry::FragmentSplitRings split_rings( true);
      chemistry::ConformationGraphConverter graph_maker;

      chemistry::FragmentGraphMarker map_rings( graph_maker, split_rings);

    /////////////////
    // data access //
    /////////////////

      // load an ensemble directly from an input stream
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "mGluR5_five_actives.sdf"));
      chemistry::FragmentEnsemble ensemble_a( input);
      io::File::CloseClearFStream( input);

      map_rings.Insert( ensemble_a);

      BCL_ExampleCheck( map_rings.GetConstitution().GetSize(), size_t( 4));

      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "taxol.sdf"));
      chemistry::FragmentEnsemble ensemble( input);
      io::File::CloseClearFStream( input);

      map_rings.Insert( ensemble);

      BCL_ExampleCheck( map_rings.GetConstitution().GetSize(), size_t( 5));

      graph::ConstGraph< size_t, size_t> marked_graph( map_rings( ensemble.GetMolecules().FirstElement()));

      linal::Matrix< size_t> edge_data_matrix( marked_graph.GetEdgeDataMatrix());

      std::string edge_data( AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), "fragment_graph_marker.mat"));

      io::OFStream output;
      BCL_ExampleMustOpenOutputFile( output, edge_data);
      edge_data_matrix.WriteObject( output, 1);
      io::File::CloseClearFStream( output);

      std::string correct_edge_data( edge_data + ".correct");
      if
      (
        BCL_ExampleIndirectCheck
        (
          io::File::FilesMatch( edge_data, correct_edge_data),
          true,
          "graph was highlighted correctly"
        )
      )
      {
        remove( edge_data.c_str());
      }

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
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryFragmentGraphMarker

  const ExampleClass::EnumType ExampleChemistryFragmentGraphMarker::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryFragmentGraphMarker())
  );

} // namespace bcl
