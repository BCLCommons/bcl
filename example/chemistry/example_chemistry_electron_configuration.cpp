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
#include "chemistry/bcl_chemistry_electron_configuration.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_element_types.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_electron_configuration.cpp
  //!
  //! @author mendenjl
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryElectronConfiguration :
    public ExampleInterface
  {
  public:

    ExampleChemistryElectronConfiguration *Clone() const
    {
      return new ExampleChemistryElectronConfiguration( *this);
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

      // get electron configuration for Br
      chemistry::ElectronConfiguration br_ec( chemistry::GetElementTypes().e_Bromine->GetElectronConfiguration());

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // make sure the 4p orbital of bromine has 5 electrons
      BCL_ExampleCheck
      (
        br_ec( chemistry::ElectronConfiguration::e_4, chemistry::ElectronConfiguration::e_P),
        5
      );

      // ensure bromine has 7 valence electrons when considering sp orbitals
      BCL_ExampleCheck
      (
        br_ec.ValenceElectronsSP(),
        7
      );

      // ensure bromine has 17 valence electrons when considering spd orbitals
      BCL_ExampleCheck
      (
        br_ec.ValenceElectronsSPD(),
        17
      );

    //////////////////////
    // input and output //
    //////////////////////

      // make sure that everything that an electron configuration writes out is read back in
      BCL_ExampleCheck( ( TestBCLObjectIOForSymmetry( br_ec, chemistry::ElectronConfiguration())), true);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryElectronConfiguration

  const ExampleClass::EnumType ExampleChemistryElectronConfiguration::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryElectronConfiguration())
  );

} // namespace bcl
