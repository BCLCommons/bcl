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
#include "chemistry/bcl_chemistry_bond_configurational.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_configuration_shared.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_bond_configurational.cpp
  //!
  //! @author kothiwsk
  //! @date Dec 02, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryBondConfigurational :
    public ExampleInterface
  {
  public:

    ExampleChemistryBondConfigurational *Clone() const
    {
      return new ExampleChemistryBondConfigurational( *this);
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

      // construct from conformation interface
      io::IFStream input_sdf;
      const std::string hexane_filename( AddExampleInputPathToFilename( e_Chemistry, "cyclohexane.sdf"));
      BCL_ExampleMustOpenInputFile( input_sdf, hexane_filename);

      // load information into small_mol_conformation

      // create configuration to test bond_configuration
      chemistry::FragmentConfigurationShared configuration;
      configuration = sdf::FragmentFactory::MakeConfiguration( input_sdf);
      // close the input stream
      io::File::CloseClearFStream( input_sdf);

    /////////////////
    // data access //
    /////////////////

      // test GetBondType
      BCL_ExampleCheck
      (
        configuration.GetAtomsIterator()->GetBonds().FirstElement().GetBondType(),
        chemistry::GetConfigurationalBondTypes().e_NonConjugatedSingleBondInRing
      );

      // test GetTargetAtom
      iterate::Generic< const chemistry::AtomConfigurationalInterface> itr( configuration.GetAtomsIterator());
      BCL_ExampleCheck
      (
        &configuration.GetAtomsIterator()->GetBonds().FirstElement().GetTargetAtom()
        == &*++itr,
        true
      );

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

  }; //end ExampleChemistryBondConfigurational

  const ExampleClass::EnumType ExampleChemistryBondConfigurational::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryBondConfigurational())
  );

} // namespace bcl
