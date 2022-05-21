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
#include "sdf/bcl_sdf_fragment_factory.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_fragment_configuration_shared.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_sdf_fragment_factory.cpp
  //! @details Demonstrates how the FragmentFactory class is used to read and write SmallMolecule
  //! and SmallMoleculeEnsembles from files
  //!
  //! @author mendenjl
  //! @date Mar 12, 2012
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleSdfFragmentFactory :
    public ExampleInterface
  {
  public:

    ExampleSdfFragmentFactory *Clone() const
    {
      return new ExampleSdfFragmentFactory( *this);
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

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      //  read small_mol_conformation
      // let's get information for small_mol_conformation
      io::IFStream input_sdf;
      const std::string hexane_filename( AddExampleInputPathToFilename( e_Chemistry, "cyclohexane.sdf"));
      BCL_ExampleMustOpenInputFile( input_sdf, hexane_filename);

      // load information into small_mol_conformation
      chemistry::FragmentComplete small_mol_conformation( sdf::FragmentFactory::MakeFragment( input_sdf));

      BCL_ExampleCheck( small_mol_conformation.GetNumberAtoms(), 6);

      // close the input stream
      io::File::CloseClearFStream( input_sdf);
      BCL_ExampleMustOpenInputFile( input_sdf, hexane_filename);

      // load information into small_mol_conformation
      chemistry::FragmentConfigurationShared
        configuration( sdf::FragmentFactory::MakeConfiguration( input_sdf, sdf::e_Saturate));

      BCL_ExampleCheck( configuration.GetNumberAtoms(), 18);
      BCL_ExampleCheck
      (
        configuration.GetAtomTypesString(),
        "C_TeTeTeTe C_TeTeTeTe C_TeTeTeTe C_TeTeTeTe C_TeTeTeTe C_TeTeTeTe H_S H_S H_S H_S H_S H_S H_S H_S H_S H_S H_S H_S "
      );

      // close the input stream
      io::File::CloseClearFStream( input_sdf);

    //////////////////////
    // input and output //
    //////////////////////

//      // test the stream operators
//      BCL_MessageStd( "Outputting factory_a: " + util::Format()( factory_a));

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;
  }; //end ExampleSdfFragmentFactory

  const ExampleClass::EnumType ExampleSdfFragmentFactory::s_Instance
  (
    GetExamples().AddEnum( ExampleSdfFragmentFactory())
  );
} // namespace bcl
