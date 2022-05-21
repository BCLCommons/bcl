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
#include "io/bcl_io_binary_serialize.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_io_binary_serialize.cpp
  //!
  //! @author mendenjl
  //! @date Feb 02, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleIoBinarySerialize :
    public ExampleInterface
  {
  public:

    ExampleIoBinarySerialize *Clone() const
    { return new ExampleIoBinarySerialize( *this);}

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
      // the actual code that shows and test the functions is implemented in the TestClass' Read and Write function
      const std::string filename( AddExampleOutputPathToFilename( io::GetNamespaceIdentifier(), "test.bin"));
      io::OFStream write;
      BCL_ExampleMustOpenBinaryOutputFile( write, filename);
      std::string test;
      for( size_t i( 0); i < 256; ++i)
      {
        // add each character
        test += char( i);
      }
      double dbl_test( 1.7e65);
      float flt_test( -234.45);
      io::BinarySerialize::Write( test, write);
      io::BinarySerialize::Write( true, write);
      io::BinarySerialize::Write( false, write);
      io::BinarySerialize::Write( dbl_test, write);
      io::BinarySerialize::Write( flt_test, write);
      io::File::CloseClearFStream( write);

      io::IFStream read;
      BCL_ExampleMustOpenBinaryInputFile( read, filename);
      std::string test_read;
      bool true_read, false_read;
      double dbl_read;
      float flt_read;

      io::BinarySerialize::Read( test_read, read);
      io::BinarySerialize::Read( true_read, read);
      io::BinarySerialize::Read( false_read, read);
      io::BinarySerialize::Read( dbl_read, read);
      io::BinarySerialize::Read( flt_read, read);
      io::File::CloseClearFStream( read);

      BCL_ExampleIndirectCheck( test_read.size(), 256, "I/O of binary string");
      BCL_ExampleIndirectCheck( test_read == test, true, "I/O of binary string; contains all characters (0-255)");
      BCL_ExampleIndirectCheck( true_read, true, "I/O of binary bool");
      BCL_ExampleIndirectCheck( false_read, false, "I/O of binary bool");
      BCL_ExampleIndirectCheck( dbl_read, dbl_test, "I/O of binary double");
      BCL_ExampleIndirectCheck( flt_read, flt_test, "I/O of binary float");

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleIoBinarySerialize

  const ExampleClass::EnumType ExampleIoBinarySerialize::s_Instance
  (
    GetExamples().AddEnum( ExampleIoBinarySerialize())
  );

} // namespace bcl
