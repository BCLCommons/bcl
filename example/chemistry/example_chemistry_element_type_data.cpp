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
#include "chemistry/bcl_chemistry_element_type_data.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_element_type_data.cpp
  //!
  //! @author mendenjl
  //! @date August 01, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryElementTypeData :
    public ExampleInterface
  {
  public:

    ExampleChemistryElementTypeData *Clone() const
    {
      return new ExampleChemistryElementTypeData( *this);
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

      #define NONSTANDARD_FORMAT
                                      //  atomic number  period  Grp Symbol    Element Name       mass             gyromagnetic ratio  covalent radiu  vdw radius  melting poin  boiling poi             electro negativity  ionization poten                                  sp valence electrons  spd valence electro   1s   1p   1d   1f   2s   2p   2d   2f   3s   3p   3d   3f   4s   4p   4d   4f   5s   5p   5d   5f   6s   6p   6d   6f   7s   7p   7d   7f     pymol rgb r   pymol rgb g  pymol rgb b
      const chemistry::ElementTypeData h(             1,      1,   1, "H"     , "Hydrogen"     ,    1.01,                     267510000,           0.32,       0.79,        14.03,       20.27,                           2.2,             13.6, 1.2, 0.66, chemistry::ElectronConfiguration(                             1,                   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ), 0.900000000f, 0.900000000f, 0.900000000f, chemistry::ElementStructureFactor());
      const chemistry::ElementTypeData e2(            2,      1,   2, "He"    , "Helium"       ,       4,                     203780000,           0.93,       0.49,         0.95,        4.22, util::GetUndefined< double>(),            24.59, 1.4, 0.66, chemistry::ElectronConfiguration(                             2,                   2,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ), 0.850980392f, 1.000000000f, 1.000000000f, chemistry::ElementStructureFactor());
      #undef NONSTANDARD_FORMAT

    /////////////////
    // data access //
    /////////////////

      // chemical name
      BCL_ExampleCheck( h.GetChemicalName(), "Hydrogen");
      BCL_ExampleCheck( h.GetChemicalSymbol(), "H");
      BCL_ExampleCheck( h.GetPeriod(), 1);
      BCL_ExampleCheck( h.GetAtomicNumber(), size_t( 1));

      for( size_t i( 0); i < size_t( chemistry::ElementTypeData::s_NumberOfProperties); i++)
      {
        BCL_MessageStd
        (
          chemistry::ElementTypeData::GetPropertyName( chemistry::ElementTypeData::Properties( i)) + ": "
          + util::Format()( h.GetProperty( chemistry::ElementTypeData::Properties( i)))
        );
      }

      BCL_MessageStd
      (
        "Electron configuration for H: " + util::Format()( h.GetElectronConfiguration())
      );
      BCL_MessageStd
      (
        "Pymol RGB code for H: " + util::Format()( h.GetPymolColorRGB()[ 0]) + ','
        + util::Format()( h.GetPymolColorRGB()[ 1]) + ',' + util::Format()( h.GetPymolColorRGB()[ 2])
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

      // Checks Read and Write
      BCL_ExampleIndirectCheck
      (
        ExampleInterface::TestBCLObjectIOForSymmetry( h, chemistry::ElementTypeData()),
        true,
        "I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;
  }; //end ExampleChemistryElementTypeData

  const ExampleClass::EnumType ExampleChemistryElementTypeData::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryElementTypeData())
  );

} // namespace bcl
