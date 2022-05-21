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
#include "chemistry/bcl_chemistry_element_types.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_element_types.cpp
  //!
  //! @author mendenjl
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryElementTypes :
    public ExampleInterface
  {
  public:

    ExampleChemistryElementTypes *Clone() const
    {
      return new ExampleChemistryElementTypes( *this);
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

      const chemistry::ElementType carbon( chemistry::GetElementTypes().ElementTypeLookup( "C"));
      const chemistry::ElementType calcium( chemistry::GetElementTypes().ElementTypeLookup( "Ca"));
      const chemistry::ElementType deuterium( chemistry::GetElementTypes().ElementTypeLookup( "2H"));
      const chemistry::ElementType bromine( chemistry::GetElementTypes().e_Bromine);

      // assert that each of the element types is defined
      BCL_ExampleAssert( carbon.IsDefined(), true);
      BCL_ExampleAssert( calcium.IsDefined(), true);
      BCL_ExampleAssert( deuterium.IsDefined(), true);
      BCL_ExampleAssert( bromine.IsDefined(), true);

    /////////////////
    // data access //
    /////////////////

      // check that each of the names is correct
      BCL_ExampleCheck( carbon->GetChemicalName(), "Carbon");
      BCL_ExampleCheck( calcium->GetChemicalName(), "Calcium");
      BCL_ExampleCheck( deuterium->GetChemicalName(), "Hydrogen");
      BCL_ExampleCheck( bromine->GetChemicalName(), "Bromine");

      // confirm that properties have been set properly
      BCL_ExampleCheck
      (
        bromine->GetElectronConfiguration().ValenceElectronsSP(),
        7
      );
      BCL_ExampleCheck
      (
        bromine->GetProperty( chemistry::ElementTypeData::e_ElectroNegativity),
        3.0
      );
      BCL_ExampleCheck
      (
        bromine->GetProperty( chemistry::ElementTypeData::e_BoilingPoint),
        332.25
      );
      BCL_ExampleCheck
      (
        bromine->GetProperty( chemistry::ElementTypeData::e_IonizationPotential),
        11.81
      );
      BCL_ExampleCheck
      (
        bromine->GetProperty( chemistry::ElementTypeData::e_CovalentRadius),
        1.14
      );
      BCL_ExampleCheck
      (
        bromine->GetProperty( chemistry::ElementTypeData::e_VDWaalsRadius),
        1.85
      );
      BCL_ExampleCheck
      (
        calcium->GetElectronConfiguration()
        (
          chemistry::ElectronConfiguration::e_1,
          chemistry::ElectronConfiguration::e_S
        ),
        2
      );
      BCL_ExampleCheck
      (
        bromine->GetElectronConfiguration()
        (
          chemistry::ElectronConfiguration::e_4,
          chemistry::ElectronConfiguration::e_P
        ),
        5
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

      BCL_MessageStd
      (
        " here is bromine's electron configuration: "
        + util::Format()( bromine->GetElectronConfiguration())
      )

      // test element type i/o
      BCL_ExampleCheck( ( TestBCLObjectIOForSymmetry( carbon, deuterium)), true);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;
  }; //end ExampleChemistryElementTypes

  const ExampleClass::EnumType ExampleChemistryElementTypes::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryElementTypes())
  );
} // namespace bcl
