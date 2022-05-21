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
#include "descriptor/bcl_descriptor_molecule_lipinski_violations.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_molecule_lipinski_violations.cpp
  //!
  //! @author geanesar
  //! @date Dec 29, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMoleculeLipinskiViolations :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMoleculeLipinskiViolations *Clone() const
    {
      return new ExampleDescriptorMoleculeLipinskiViolations( *this);
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

      // Construction
      descriptor::MoleculeLipinskiViolations orig;
      descriptor::MoleculeLipinskiViolations veber( descriptor::MoleculeLipinskiViolations::e_Veber);

      // Copying and destruction
      descriptor::MoleculeLipinskiViolations copy( orig);
      util::ShPtr< descriptor::MoleculeLipinskiViolations> sp_copy( orig.Clone());
      sp_copy = util::ShPtr< descriptor::MoleculeLipinskiViolations>();

      // Data access
      BCL_ExampleCheck( orig.GetAlias(), "LipinskiViolations");
      BCL_ExampleCheck( veber.GetAlias(), "LipinskiViolationsVeber");

      // operations
      chemistry::FragmentEnsemble ensemble;
      io::IFStream read;

      // read taxol - should have 2 original Ro5 violations
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Chemistry, "taxol.sdf"));
      ensemble.ReadMoreFromMdl( read, sdf::e_Saturate);
      io::File::CloseClearFStream( read);

      // read in diazepam - good drug molecule
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Chemistry, "diazepam.sdf"));
      ensemble.ReadMoreFromMdl( read, sdf::e_Saturate);
      io::File::CloseClearFStream( read);

    ///////////////
    // operators //
    ///////////////
      
      storage::Vector< chemistry::FragmentComplete> mols( ensemble.Begin(), ensemble.End());

      // Expected values (logP values are calculated, not experimental, so these may not match up with other sources)
      float expect_orig[] = { 3.0, 1.0};
      float expect_veber[] = { 2.0, 0.0};

      // Check the calculated values
      for( size_t mol_no( 0); mol_no < mols.GetSize(); ++mol_no)
      {
        BCL_ExampleIndirectCheck
        (
          orig( mols( mol_no))( 0),
          expect_orig[ mol_no],
          "Molecule " + util::Format()( mol_no) + " original lipinski check"
        );
        BCL_ExampleIndirectCheck
        (
          veber( mols( mol_no))( 0),
          expect_veber[ mol_no],
          "Molecule " + util::Format()( mol_no) + " Veber variant check"
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
        "MoleculeLipinskiViolations I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMoleculeLipinskiViolations

  const ExampleClass::EnumType ExampleDescriptorMoleculeLipinskiViolations::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMoleculeLipinskiViolations())
  );

} // namespace bcl
