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
#include "descriptor/bcl_descriptor_molecule_complexity.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_molecule_complexity.cpp
  //!
  //! @author geanesar, mendenjl
  //! @date Dec 29, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMoleculeComplexity :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMoleculeComplexity *Clone() const
    {
      return new ExampleDescriptorMoleculeComplexity( *this);
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

      // default constructor
      descriptor::MoleculeComplexity complexity;

      // copy constructor
      descriptor::MoleculeComplexity complexity_copy( complexity);

      // Clone/destruction
      util::ShPtr< descriptor::MoleculeComplexity> copy( util::CloneToShPtr( complexity));

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck( complexity.GetAlias(), "MoleculeComplexity");

      BCL_ExampleCheck( complexity.GetString(), "MoleculeComplexity");

    ///////////////
    // operators //
    ///////////////

      storage::Map< std::string, float> filenames_to_expected_complexity;

      // files names with expected maximum distances calculated using pymol
      filenames_to_expected_complexity[ "diazepam.sdf"] = 0.77;
      filenames_to_expected_complexity[ "taxol.sdf"] = 4.68;
      filenames_to_expected_complexity[ "hexane.sdf"] = 0.05;

      for
      (
        storage::Map< std::string, float>::const_iterator
        itr_files( filenames_to_expected_complexity.Begin()),
        itr_files_end( filenames_to_expected_complexity.End());
        itr_files != itr_files_end;
        ++itr_files
      )
      {
        // create input stream for reading a smallmolecule ensemble
        io::IFStream input;
        BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, itr_files->first));
        // read in ensemble
        chemistry::FragmentEnsemble ensemble( input);
        // close stream
        io::File::CloseClearFStream( input);

        const chemistry::FragmentComplete &mol( *ensemble.Begin());

        BCL_ExampleIndirectCheckWithinAbsTolerance
        (
          complexity( mol)( 0),
          float( itr_files->second),
          0.01,
          "Calculating the complexity of molecules in " + itr_files->first
        );
      }

    //////////////////////
    // input and output //
    //////////////////////

      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( complexity, complexity_copy),
        true,
        "MoleculeComplexity I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMoleculeComplexity

  const ExampleClass::EnumType ExampleDescriptorMoleculeComplexity::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMoleculeComplexity())
  );

} // namespace bcl
