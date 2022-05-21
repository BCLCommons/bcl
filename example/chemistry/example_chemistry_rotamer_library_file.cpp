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
#include "chemistry/bcl_chemistry_rotamer_library_file.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_directory.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_factory.h"

// external includes - sorted alphabetically
#include <stdio.h>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_rotamer_library_file.cpp
  //!
  //! @author kothiwsk
  //! @date Jul 01, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryRotamerLibraryFile :
    public ExampleInterface
  {
  public:

    ExampleChemistryRotamerLibraryFile *Clone() const
    {
      return new ExampleChemistryRotamerLibraryFile( *this);
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

      // create input stream for reading an ensemble of fragments
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "mGluR5_rotamers.sdf"));

      // read in ensemble of fragments
      chemistry::FragmentEnsemble fragment_ensemble( input, sdf::e_Remove);

      // close stream
      io::File::CloseClearFStream( input);

    //////////////////
    // construction //
    //////////////////

      std::string output_path_file( AddExampleOutputPathToFilename( fragment_ensemble, "confgenerator_tree"));
      util::Implementation< chemistry::RotamerLibraryInterface> rotamer_lib
      (
        "File(prefix=" + output_path_file + ",number_of_files=10,compression=Uncompressed)"
      );

      rotamer_lib->Create( fragment_ensemble);

    /////////////////
    // data access //
    /////////////////

    ////////////////
    // operations //
    ////////////////

      std::string substructure( output_path_file + ".substructure.txt");
      std::string substructure_correct( output_path_file + ".substructure.txt.correct");

      std::string constitutions( output_path_file + ".constitutions.txt");
      std::string constitutions_correct( output_path_file + ".constitutions.txt.correct");

      std::string mapping( output_path_file + ".configuration_mapping.txt");
      std::string mapping_correct( output_path_file + ".configuration_mapping.txt.correct");

      if
      (
        BCL_ExampleCheck( io::File::FilesMatch( substructure, substructure_correct), true) &&
        BCL_ExampleCheck( io::File::FilesMatch( constitutions, constitutions_correct), true) &&
        BCL_ExampleCheck( io::File::FilesMatch( mapping, mapping_correct), true)
      )
      {
        remove( std::string( output_path_file+".substructure.txt").c_str());
        remove( std::string( output_path_file+".constitutions.txt").c_str());
        remove( std::string( output_path_file+".configuration_mapping.txt").c_str());
        remove( std::string( output_path_file+".requirements.txt").c_str());
        io::Directory( output_path_file + "_conformations").Remove( true);
      }

    //////////////////////
    // input and output //
    //////////////////////

    /////////////
    // cleanup //
    /////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryRotamerLibraryFile

  const ExampleClass::EnumType ExampleChemistryRotamerLibraryFile::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryRotamerLibraryFile())
  );

} // namespace bcl
