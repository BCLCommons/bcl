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
#include "restraint/bcl_restraint_sas_experimental_and_calculated_data.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_saxs_experimental_and_calculated_data.cpp
  //! @details This example tests the reading and storing of SAXS Experimental Data
  //!
  //! @author putnamdk, mendenjl
  //! @date Sep 04, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRestraintSasExperimentalAndCalculatedData :
    public ExampleInterface
  {
  public:

    ExampleRestraintSasExperimentalAndCalculatedData *Clone() const
    {
      return new ExampleRestraintSasExperimentalAndCalculatedData( *this);
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
      restraint::SasExperimentalAndCalculatedData saxs_data;

      BCL_ExampleIndirectCheck( saxs_data.GetExperimentalData().GetScatteringData().IsEmpty(), true, "default constructor");

    //////////////////////
    // input and output //
    //////////////////////

      // test reading in a data file
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "A10_2.mccd.gnom"));
      restraint::SasScatteringData saxs_experimental_data;

      saxs_experimental_data.ReadFromDataFile( read);
      saxs_data = restraint::SasExperimentalAndCalculatedData( saxs_experimental_data, saxs_experimental_data);

      BCL_ExampleIndirectCheck( saxs_data.GetExperimentalData().GetScatteringData().GetSize(), 508, "ReadFromDataFile");
      io::File::CloseClearFStream( read);

      // test read and write
      WriteBCLObject( saxs_data);

      // Create new object
      restraint::SasExperimentalAndCalculatedData saxs_data_read;

      // Read in new object
      ReadBCLObject( saxs_data_read);

      // Compare original and read in object

      // Test the size of vector
      BCL_ExampleIndirectAssert
      (
        saxs_data.GetExperimentalData().GetScatteringData().GetSize(),
        saxs_data_read.GetCalculatedData().GetScatteringData().GetSize(),
        "read and write"
      );

      // Test the first element in the vector (q-value)
      BCL_ExampleIndirectCheckWithinTolerance
      (
        saxs_data.GetExperimentalData().GetScatteringData().FirstElement().GetQvalue(),
        saxs_data_read.GetExperimentalData().GetScatteringData().FirstElement().GetQvalue(),
        0.001,
        "read and write"
      );

      // Test the second element in the vector (intensity)
      BCL_ExampleIndirectCheckWithinTolerance
      (
        saxs_data.GetExperimentalData().GetScatteringData().FirstElement().GetIntensity(),
        saxs_data_read.GetExperimentalData().GetScatteringData().FirstElement().GetIntensity(),
        0.001,
        "read and write"
      );

      // Test the third element in the vector (error)
      BCL_ExampleIndirectCheckWithinTolerance
      (
        saxs_data.GetExperimentalData().GetScatteringData().FirstElement().GetError(),
        saxs_data_read.GetExperimentalData().GetScatteringData().FirstElement().GetError(),
        0.001,
        "read and write"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRestraintSasExperimentalAndCalculatedData

  const ExampleClass::EnumType ExampleRestraintSasExperimentalAndCalculatedData::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintSasExperimentalAndCalculatedData())
  );

} // namespace bcl
