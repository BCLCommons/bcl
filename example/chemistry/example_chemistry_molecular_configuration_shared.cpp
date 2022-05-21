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
#include "chemistry/bcl_chemistry_molecular_configuration_shared.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_molecular_configuration_shared.cpp
  //! @details Tests ChemistryMolecularConfigurationShared class which contains small molecule configuration data
  //!
  //! @author kothiwsk
  //! @date
  //! @remarks status complete
  //! @remarks
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryMolecularConfigurationShared :
    public ExampleInterface
  {
  public:

    ExampleChemistryMolecularConfigurationShared *Clone() const
    {
      return new ExampleChemistryMolecularConfigurationShared( *this);
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

      // read in molecule and construct configuration from conformation interface
      io::IFStream input_sdf;
      const std::string hexane_filename( AddExampleInputPathToFilename( e_Chemistry, "cyclohexane.sdf"));
      BCL_ExampleMustOpenInputFile( input_sdf, hexane_filename);

      // load information into small_mol_conformation
      chemistry::FragmentComplete small_mol_conformation;
      small_mol_conformation
        = sdf::FragmentFactory::MakeFragment( input_sdf, sdf::e_Saturate);
      // close the input stream
      io::File::CloseClearFStream( input_sdf);

      chemistry::MolecularConfigurationShared configuration_a( small_mol_conformation);

    /////////////////
    // data access //
    /////////////////

      // test GetBonds
      BCL_ExampleCheck( configuration_a.GetNumberBonds(), 18);

      // test GetNumberofAtoms
      BCL_ExampleCheck( configuration_a.GetNumberAtoms(), 18);

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // test the stream operators
      BCL_MessageStd
      (
        "Outputting configuration: " + util::Format()( configuration_a)
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSmallMolecule

  const ExampleClass::EnumType ExampleChemistryMolecularConfigurationShared::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryMolecularConfigurationShared())
  );

} // namespace bcl
