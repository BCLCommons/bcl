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
#include "descriptor/bcl_descriptor_molecule_entropy_qha.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "chemistry/bcl_chemistry_sample_conformations.h"
#include "io/bcl_io_file.h"

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_molecule_entropy_qha.cpp
  //!
  //! @author brownbp1, mendenjl
  //! @date Aug 22, 2019
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMoleculeEntropyQHA :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMoleculeEntropyQHA *Clone() const
    {
      return new ExampleDescriptorMoleculeEntropyQHA( *this);
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

      // Get a sample conformations object
      chemistry::RotamerLibraryFile rotamer_library_file;
      chemistry::SampleConformations sampler
      (
        rotamer_library_file,
        "",
        0.25,   // tolerance
        100,    // number of conformations
        1000,    // number of iterations
        false,  // change chirality
        0.0,    // random dihedral change weight
        false,  // generate 3d
        0.05    // clash tolerance
      );

      // Construction
      descriptor::MoleculeEntropyQHA orig( sampler);

      // Copying and destruction
      descriptor::MoleculeEntropyQHA copy( orig);
      util::ShPtr< descriptor::MoleculeEntropyQHA> sp_copy( orig.Clone());
      sp_copy = util::ShPtr< descriptor::MoleculeEntropyQHA>();

      // Data access
      BCL_ExampleCheck( orig.GetAlias(), "EntropyQHA");

      // operations
      chemistry::FragmentEnsemble ensemble;
      io::IFStream read;

      // read in diazepam
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Chemistry, "diazepam.sdf"));
      ensemble.ReadMoreFromMdl( read, sdf::e_Saturate);
      io::File::CloseClearFStream( read);

    ///////////////
    // operators //
    ///////////////

      storage::Vector< chemistry::FragmentComplete> mols( ensemble.Begin(), ensemble.End());

      // Check the calculated values
      linal::Vector< float> diffs( size_t( 1));
      diffs( 0) = 1290.68;
      for( size_t mol_no( 0); mol_no < mols.GetSize(); ++mol_no)
      {
        // generate entropy descriptor output
        linal::Vector< float> entropy( orig( mols( mol_no)));

        // check to make sure each value in entropy vector matches expected
        BCL_ExampleIndirectCheckWithinAbsTolerance
        (
          entropy( 2), diffs( mol_no), 100.0, "QHA Entropy difference (global - local)"
        );
      }

    //////////////////////
    // input and output //
    //////////////////////

      // Check for symmetry
      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( orig, copy),
        true,
        "MoleculeEntropyQHA I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMoleculeEntropyQHA

  const ExampleClass::EnumType ExampleDescriptorMoleculeEntropyQHA::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMoleculeEntropyQHA())
  );

} // namespace bcl
