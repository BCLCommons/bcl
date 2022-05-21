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
#include "chemistry/bcl_chemistry_search_fragment_library_from_tree.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_directory.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically
#include <stdio.h>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_search_fragment_library_from_tree.cpp
  //!
  //! @author kothiwsk
  //! @date Jul 09, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistrySearchFragmentLibraryFromTree :
    public ExampleInterface
  {
  public:

    ExampleChemistrySearchFragmentLibraryFromTree *Clone() const
    {
      return new ExampleChemistrySearchFragmentLibraryFromTree( *this);
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

      // create input stream for reading an ensemble of fragments
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "mGluR5_rotamers.sdf"));

      // read in ensemble of fragments
      chemistry::FragmentEnsemble fragment_ensemble( input, sdf::e_Remove);

      // close stream
      io::File::CloseClearFStream( input);

      std::string output_path_file( AddExampleOutputPathToFilename( fragment_ensemble, "confgenerator_tree"));
      std::string output_path_file_a( AddExampleOutputPathToFilename( fragment_ensemble, "confgenerator_tree_search"));

      util::Implementation< chemistry::RotamerLibraryInterface> rotamer_lib
      (
        "File(prefix=" + output_path_file_a + ",number_of_files=10, compression=Uncompressed)"
      );

      rotamer_lib->Create( fragment_ensemble);

      // constructor
      chemistry::SearchFragmentLibraryFromTree search_lib( *rotamer_lib);

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // create input stream for reading molecule of interest
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "mGluR5_five_actives.sdf"));

      // read in the first molecule molecules
      chemistry::FragmentEnsemble molecules
      (
        storage::List< chemistry::FragmentComplete>
        (
          1,
          chemistry::FragmentEnsemble( input, sdf::e_Remove).GetMolecules().FirstElement()
        )
      );
      // close stream
      io::File::CloseClearFStream( input);

      // search for fragments of the molecule
      const util::ShPtrVector< chemistry::SmallMoleculeFragmentIsomorphism> molecule_isomorphism
      (
        search_lib.FindFragmentsOfMolecule( molecules.GetMolecules().FirstElement())
      );

      // print out fragments of the molecule to a file
      const std::string fragment_search( AddExampleOutputPathToFilename( molecules, "fragments_tree_contained.sdf"));
      io::OFStream output;
      BCL_ExampleMustOpenOutputFile( output, fragment_search);
      for
      (
        util::ShPtrVector< chemistry::SmallMoleculeFragmentIsomorphism>::const_iterator
          itr_ensemble( molecule_isomorphism.Begin()), itr_ensemble_end( molecule_isomorphism.End());
        itr_ensemble != itr_ensemble_end;
        ++itr_ensemble
      )
      {
        // write the fragment out to a file
        ( *itr_ensemble)->GetFragment().WriteMDL( output);
      }
      io::File::CloseClearFStream( output);
      std::string correct_filename( fragment_search + ".correct");

      std::string substructure( output_path_file_a + ".substructure.txt");
      std::string substructure_correct( output_path_file + ".substructure.txt.correct");

      std::string constitutions( output_path_file_a + ".constitutions.txt");
      std::string constitutions_correct( output_path_file + ".constitutions.txt.correct");

      std::string mapping( output_path_file_a + ".configuration_mapping.txt");
      std::string mapping_correct( output_path_file + ".configuration_mapping.txt.correct");

      // check if the generated file is right
      if
      (
        BCL_ExampleCheck( io::File::FilesMatch( fragment_search, correct_filename), true)
      )
      {
        remove( fragment_search.c_str());

        if
        (
          BCL_ExampleCheck( io::File::FilesMatch( substructure, substructure_correct), true) &&
          BCL_ExampleCheck( io::File::FilesMatch( constitutions, constitutions_correct), true) &&
          BCL_ExampleCheck( io::File::FilesMatch( mapping, mapping_correct), true)
        )
        {
          remove( std::string( output_path_file_a+".substructure.txt").c_str());
          remove( std::string( output_path_file_a+".constitutions.txt").c_str());
          remove( std::string( output_path_file_a+".configuration_mapping.txt").c_str());
          remove( std::string( output_path_file_a+".requirements.txt").c_str());
          io::Directory( output_path_file_a + "_conformations").Remove( true);
        }
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

  }; //end ExampleChemistrySearchFragmentLibraryFromTree

  const ExampleClass::EnumType ExampleChemistrySearchFragmentLibraryFromTree::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistrySearchFragmentLibraryFromTree())
  );

} // namespace bcl
