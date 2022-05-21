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
#include "chemistry/bcl_chemistry_fragment_connector.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically
#include <stdio.h>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_fragment_connector.cpp
  //!
  //! @author kothiwsk
  //! @date Oct 10, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryFragmentConnector :
    public ExampleInterface
  {
  public:

    ExampleChemistryFragmentConnector *Clone() const
    {
      return new ExampleChemistryFragmentConnector( *this);
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

      // create input stream for reading a fragment ensemble
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "1_3_pentadiene_E.sdf"));
      // read in ensemble
      util::ShPtr< chemistry::FragmentEnsemble> fragments
      (
        new chemistry::FragmentEnsemble( input, sdf::e_Remove)
      );
      io::File::CloseClearFStream( input);

      // read another set of molecules
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "benzene_prepared.sdf"));
      fragments->Append( chemistry::FragmentEnsemble( input, sdf::e_Remove));
      // close stream
      io::File::CloseClearFStream( input);

      // get reference to first and second molecules in the ensemble
      chemistry::FragmentComplete &first_mol( fragments->GetMolecules().FirstElement());
      chemistry::FragmentComplete &second_mol( fragments->GetMolecules().LastElement());

      // create a set of vertices to represent what atoms the second molecule represents in the final assembled molecule
      const linal::Vector< size_t> fill_vec
      (
        linal::FillVector< size_t>( second_mol.GetNumberAtoms(), size_t( 0), size_t( 1))
      );

      // get storage vector from linal vector
      storage::Vector< size_t> second_mol_vertices( fill_vec.Begin(), fill_vec.End());

      // construct FragmentConnector
      chemistry::FragmentConnector fragment_connector( second_mol, second_mol_vertices);

      // join first molecule to fragment held in FragmentConnector
      fragment_connector
      (
        first_mol,
        chemistry::ConfigurationalBondType( 1),
        storage::Pair< size_t, size_t>( 1, 0),
        storage::Vector< size_t>::Create( 6, 7, 8, 9, 10)
      );

      // update second_mol_vertices
      second_mol_vertices.Append( storage::Vector< size_t>::Create( 6, 7, 8, 9, 10));

      // store the connected molecule in an ensmble
      chemistry::FragmentEnsemble merged_molecule_ensemble;
      merged_molecule_ensemble.PushBack( fragment_connector.GetFragment());

      // read in vinyl to be added to moleclue held in fragment_connector object
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "corina_vinyl.sdf"));
      fragments->Append( chemistry::FragmentEnsemble( input, sdf::e_Remove));
      // close stream
      io::File::CloseClearFStream( input);

      chemistry::FragmentComplete &vinyl( fragments->GetMolecules().LastElement());

      // create a map to tell the fragment connector object the overlapping atoms between vinyl and fragment held by fragment connector object
      storage::Map< size_t, size_t> common_indices;
      common_indices.Insert( storage::Pair< size_t, size_t>( 6, 3));
      common_indices.Insert( storage::Pair< size_t, size_t>( 8, 1));

      // add vinyl to fragment held by fragment connector
      fragment_connector
      (
        vinyl, common_indices, storage::Vector< size_t>::Create( 11, 8, 12, 6)
      );

      second_mol_vertices.Append( storage::Vector< size_t>::Create( 11, 12));

      merged_molecule_ensemble.PushBack( fragment_connector.GetFragment());

      // check the vertices
      BCL_ExampleCheck( fragment_connector.GetVertices() == second_mol_vertices, true);

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // check if generated molecules are the same
      const std::string merge_fragments( AddExampleOutputPathToFilename( merged_molecule_ensemble, "connected_molecules.sdf"));
      io::OFStream output;
      BCL_ExampleMustOpenOutputFile( output, merge_fragments);
      merged_molecule_ensemble.WriteMDL( output);
      io::File::CloseClearFStream( output);

      std::string correct_filename( merge_fragments + ".correct");

      // check if the generated file is right
      if
      (
        BCL_ExampleCheck( io::File::FilesMatch( merge_fragments, correct_filename), true)
      )
      {
        remove( merge_fragments.c_str());
      }

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
      } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryFragmentConnector

  const ExampleClass::EnumType ExampleChemistryFragmentConnector::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryFragmentConnector())
  );

} // namespace bcl
