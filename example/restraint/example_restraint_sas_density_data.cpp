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
#include "restraint/bcl_restraint_sas_density_data.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_sas_density_data.cpp
  //! @details This example tests the reading and storing of SAXS Density Distribution Data
  //!
  //! @author putnamdk
  //! @date July 3, 2014
  //! @remarks status complete
  //! @remarks reviewed by putnamdk on July 3, 2014
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRestraintSasDensityData :
    public ExampleInterface
  {
  public:

    ExampleRestraintSasDensityData *Clone() const
    {
      return new ExampleRestraintSasDensityData( *this);
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

      // test default constructor
      restraint::SasDensityData sas_data;

      BCL_ExampleIndirectCheck( sas_data.GetDensityData().IsEmpty(), true, "default constructor");

    //////////////////////
    // input and output //
    //////////////////////

      // test reading in a data file
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "A10_2.mccd.gnom"));

      sas_data.ReadFromDataFile( read);

      BCL_ExampleIndirectCheck( sas_data.GetDensityData().GetSize(), 101, "ReadFromDataFile");
      io::File::CloseClearFStream( read);

      // test read and write
      WriteBCLObject( sas_data);

      // Create new object
      restraint::SasDensityData sas_data_read;

      // Read in new object
      ReadBCLObject( sas_data_read);

      // Compare original and read in object

      // Test the size of vector
      BCL_ExampleIndirectAssert
      (
        sas_data.GetDensityData().GetSize(), sas_data_read.GetDensityData().GetSize(),
        "read and write"
      );

      // Test the first element in the vector (r-value)
      BCL_ExampleIndirectCheckWithinTolerance
      (
        sas_data.GetDensityData()( 2).GetRvalue(),
        sas_data_read.GetDensityData()( 2).GetRvalue(),
        0.001,
        "read and write"
      );

      // Test the second element in the vector (density)
      BCL_ExampleIndirectCheckWithinTolerance
      (
        sas_data.GetDensityData()( 2).GetDensity(),
        sas_data_read.GetDensityData()( 2).GetDensity(),
        0.001,
        "read and write"
      );

      // Test the third element in the vector (error)
      BCL_ExampleIndirectCheckWithinTolerance
      (
        sas_data.GetDensityData()( 2).GetError(),
        sas_data_read.GetDensityData()( 2).GetError(),
        0.001,
        "read and write"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRestraintSasDensityData

  const ExampleClass::EnumType ExampleRestraintSasDensityData::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintSasDensityData())
  );

} // namespace bcl
