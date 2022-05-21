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
#include "io/bcl_io_serialize.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_io_serialize.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleIoSerialize :
    public ExampleInterface
  {
  public:

    ExampleIoSerialize *Clone() const
    { return new ExampleIoSerialize( *this);}

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! example class that shows how to implement the write function for an object properly
    class TestClass :
      public util::ObjectInterface
    {
    public:
    //////////
    // data //
    //////////

      double m_Double;
      char m_Char;
      std::string m_String;
      storage::Pair< std::string, storage::Pair< size_t, double> > m_Data;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      TestClass() :
        m_Double( 0),
        m_Char( ' ')
      {
      }

      //! construct TestClass initializing all members
      TestClass
      (
        const double DOUBLE,
        const char CHAR,
        const std::string &STRING,
        const storage::Pair< std::string, storage::Pair< size_t, double> > &DATA
      ) :
        m_Double( DOUBLE),
        m_Char( CHAR),
        m_String( STRING),
        m_Data( DATA)
      {
      }

      //! @brief Clone function
      TestClass *Clone() const
      {
        return new TestClass( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // read member
        io::Serialize::Read( m_Double, ISTREAM);
        io::Serialize::Read( m_Char, ISTREAM);
        io::Serialize::Read( m_String, ISTREAM);
        io::Serialize::Read( m_Data, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write member
        io::Serialize::Write( m_Double, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Char, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_String, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Data, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // class TestClass

    int Run() const
    {
      // the actual code that shows and test the functions is implemented in the TestClass' Read and Write function

      // construct a TestClass object
      TestClass test1
      (
        5.5,
        ' ',
        "hello \"BCL friend\"",
        storage::Pair< std::string, storage::Pair< size_t, double> >
        (
          "nifty",
          storage::Pair< size_t, double>( 5, util::GetUndefined< double>())
        )
      );

      BCL_MessageStd( "this is a formatted written bcl object:\n" + util::Format()( test1));
      const std::string filename( AddExampleOutputPathToFilename( io::GetNamespaceIdentifier(), "test_class1.ser"));
      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, filename);
      write << test1;
      io::File::CloseClearFStream( write);

      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, filename);
      TestClass test_read;
      read >> test_read;
      io::File::CloseClearFStream( read);

      BCL_MessageStd( "test class after writing and reading:\n" + util::Format()( test_read));

      // check that the undefined daouble was read correctly
      BCL_Example_Check
      (
        !util::IsDefined( test_read.m_Data.Second().Second()),
        "read double should be undefined but is " + util::Format()( test_read.m_Data.Second().Second())
      );
      // check that every member was read correctly
      BCL_Example_Check
      (
        test1.m_Char == test_read.m_Char
        && test1.m_Double == test_read.m_Double
        && test1.m_String == test_read.m_String,
        "writing and reading was unsuccessful: " + util::Format()( test1) + "!=" + util::Format()( test_read)
      );

      std::stringstream ss;
      // write bool
      io::Serialize::Write( bool( true), ss, 0);
      io::Serialize::Write( bool( false), ss, 1);

      // write char
      io::Serialize::Write( char( ' '), ss, 2);

      // write double
      io::Serialize::Write( double( 0.00123), ss, 3);

      // write float
      io::Serialize::Write( float( 0.00321), ss, 4);

      // write int
      io::Serialize::Write( int( -5), ss, 5);

      // write long int
      io::Serialize::Write( ( long int)( -6), ss, 6);

      // write short int
      io::Serialize::Write( ( short int)( -7), ss, 7);

      // write size_t
      io::Serialize::Write( size_t( 8), ss, 8);

      // write std::pair
      io::Serialize::Write( std::pair< size_t, size_t>( 9, 9), ss, 9);

      // write std::string
      io::Serialize::Write( std::string( "10 10"), ss, 10);

      // write unsigned long long int
      io::Serialize::Write( ( unsigned long long int)( 11), ss, 11);

      BCL_MessageStd( "writing nearly every possible basic type " + ss.str());

      const std::string expected_string
      (
        "1"
        "  0"
        "    ' '"
        "      0.00123"
        "        0.00321"
        "          -5"
        "            -6"
        "              -7"
        "                8"
        "                  std::pair<size_t,size_t>\n"
        "                    9\n"
        "                    9"
        "                    \"10 10\""
        "                      11"
      );
      BCL_ExampleIndirectCheck( ss.str(), expected_string, "writing basic types");

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleIoSerialize

  const ExampleClass::EnumType ExampleIoSerialize::s_Instance
  (
    GetExamples().AddEnum( ExampleIoSerialize())
  );

} // namespace bcl
