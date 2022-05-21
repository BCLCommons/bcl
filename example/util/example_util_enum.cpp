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
#include "util/bcl_util_enum.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //! @example example_util_enum.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class SomeData :
    public util::ObjectInterface
  {
  private:

    double m_Value;

  public:

    //! @brief default constructor
    SomeData() :
      m_Value( util::GetUndefined< double>())
    {
    }

    SomeData( const double VALUE) :
      m_Value( VALUE)
    {}

    SomeData *Clone() const
    {
      return new SomeData( *this);
    }

    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    const double GetValue() const
    {
      return m_Value;
    }

    std::istream &Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_Value, ISTREAM);

      return ISTREAM;
    }

    std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_Value, OSTREAM, INDENT);

      return OSTREAM;
    }
  };

  class SomeDataTypes :
    public util::Enumerate< SomeData, SomeDataTypes>
  {
    friend class util::Enumerate< SomeData, SomeDataTypes>;

  public:
    const EnumType e_XXX;
    const EnumType e_ALA;
    const EnumType e_ARG;
    const EnumType e_ASN;
    const EnumType e_ASP;
    const EnumType e_CYS;
    const EnumType e_GLN;
    const EnumType e_GLU;

  private:

    SomeDataTypes() :
      e_XXX( AddEnum( "XXX", SomeData( 1.1))),
      e_ALA( AddEnum( "ALA", SomeData( 2.2))),
      e_ARG( AddEnum( "ARG", SomeData( 3.3))),
      e_ASN( AddEnum( "ASN", SomeData( 4.4))),
      e_ASP( AddEnum( "ASP", SomeData( 5.5))),
      e_CYS( AddEnum( "CYS", SomeData( 6.6))),
      e_GLN( AddEnum( "GLN", SomeData( 7.7))),
      e_GLU( AddEnum( "GLU", SomeData( 8.8)))
    {}

  public:

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  };

  SomeDataTypes &GetSomeDataTypes()
  {
    return SomeDataTypes::GetEnums();
  }

  class ExampleUtilEnum :
    public ExampleInterface
  {
  public:

    ExampleUtilEnum *Clone() const
    { return new ExampleUtilEnum( *this);}

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

      // Construct all data types by calling that function with a static variable in it
      GetSomeDataTypes();

      // construct undefined enums
      SomeDataTypes::EnumType test_undef1( util::GetUndefined< SomeDataTypes::EnumType>());
      SomeDataTypes::EnumType test_undef2( GetSomeDataTypes().e_Undefined);
      BCL_Example_Check
      (
        !test_undef1.IsDefined(), ""
      );
      BCL_Example_Check
      (
        !test_undef2.IsDefined(), ""
      );

      //construct Enum< biol::AATypeData> from name
      SomeDataTypes::EnumType test_name( "ALA");
      BCL_Example_Check
      (
        test_name.GetIndex() == 1, "enum constructed from string did not result in expected index of 1"
      );

      //construct Enum< biol::AATypeData> from index
      SomeDataTypes::EnumType test_index( 4);
      BCL_Example_Check
      (
        test_index.GetName() == "ASP", "enum constructed from index did not result in expected name of \"ASP\""
      );

    /////////////////
    // data access //
    /////////////////

      // GetEnumCount
      BCL_MessageDbg( "GetEnumCount=" + util::Format()( GetSomeDataTypes().GetEnumCount()));
      BCL_Example_Check
      (
        GetSomeDataTypes().GetEnumCount() == 8,
        "Enum< biol::AATypeData> does not return the expected number of enum values of 8!"
      );

      // get static class name
      BCL_MessageStd( "name of Enums: " + GetStaticClassName< SomeDataTypes>());

      // defined enumerator values - check if it is defined and ask for index and name
      BCL_MessageDbg
      (
        util::Format()( GetSomeDataTypes().e_ALA)
        + util::Format()( " IsDefined=") + util::Format()( GetSomeDataTypes().e_ALA.IsDefined())
        + util::Format()( " GetIndex=") + util::Format()( GetSomeDataTypes().e_ALA.GetIndex())
        + util::Format()( " GetIndex=") + util::Format()( GetSomeDataTypes().e_ALA.GetName())
      );

      // possible to make a copy
      SomeDataTypes::EnumType test_copy( GetSomeDataTypes().e_XXX);
      // a copy does not change the number of iterators
      BCL_MessageDbg( "GetEnumCount=" + util::Format()( GetSomeDataTypes().GetEnumCount()));

      // test the copied value
      BCL_MessageDbg( util::Format()( test_copy)
        + util::Format()( " Defined=") + util::Format()( test_copy.IsDefined())
        + util::Format()( " GetIndex=") + util::Format()( test_copy.GetIndex()));
      BCL_Example_Check
      (
        GetSomeDataTypes().e_XXX.IsDefined() == test_copy.IsDefined(),
        "XXX.IsDefined() == test_copy.IsDefined() does not return the expected value TRUE!"
      );
      BCL_Example_Check
      (
        GetSomeDataTypes().e_XXX.GetIndex() == test_copy.GetIndex(),
        "XXX.GetIndex() == test_copy.GetIndex() does not return the expected value TRUE!"
      );
      BCL_Example_Check
      (
        GetSomeDataTypes().e_XXX->GetValue() == test_copy->GetValue(),
        "XXX.GetData().GetValue() == test_copy.GetData().GetValue() does not return the expected value TRUE!"
      );

    ///////////////
    // operators //
    ///////////////

      // check that Is Defined works
      BCL_Example_Check
      (
        GetSomeDataTypes().e_XXX.IsDefined(),
        "XXX.IsDefined() does not return the expected value TRUE!"
      );

      // access data behind enum with operator '->' and check that the GetValue function returns the correct result
      BCL_Example_Check
      (
        GetSomeDataTypes().e_XXX->GetValue() == 1.1,
        "XXX.GetData().GetValue() does not return the expected value for the first element of 1.1 but" +
        util::Format()( GetSomeDataTypes().e_XXX->GetValue()) + "\"!"
      );
      BCL_MessageDbg
      (
        util::Format()( GetSomeDataTypes().e_ALA)
        + util::Format()( " Defined=") + util::Format()( GetSomeDataTypes().e_ALA.IsDefined())
        + util::Format()( " GetIndex=") + util::Format()( GetSomeDataTypes().e_ALA.GetIndex())
      );
      BCL_MessageDbg
      (
        util::Format()( GetSomeDataTypes().e_ARG)
        + util::Format()( " Defined=") + util::Format()( GetSomeDataTypes().e_ARG.IsDefined())
        + util::Format()( " GetIndex=") + util::Format()( GetSomeDataTypes().e_ARG.GetIndex())
      );

      // undefined enumerator value
      BCL_MessageDbg
      (
        util::Format()( util::GetUndefined< SomeDataTypes::EnumType>())
        + util::Format()( " Defined=") + util::Format()( util::GetUndefined< SomeDataTypes::EnumType>().IsDefined())
        + util::Format()( " GetIndex=") + util::Format()( util::GetUndefined< SomeDataTypes::EnumType>().GetIndex())
      );
      BCL_Example_Check
      (
        !util::GetUndefined< SomeDataTypes::EnumType>().IsDefined(),
        "util::GetUndefined< SomeDataTypes::EnumType>().IsDefined() does not return the expected value FALSE!"
      );

    //////////////////////////////
    // operations and operators //
    //////////////////////////////

      // GetEnumIteratorFromIndex
      SomeDataTypes::const_iterator found = GetSomeDataTypes().GetEnumIteratorFromIndex( 0);

      // test iterator, equal, unequal
      for( SomeDataTypes::const_iterator itr( GetSomeDataTypes().Begin()), itr_end( GetSomeDataTypes().End()); itr != itr_end; ++itr)
      {
        BCL_MessageDbg( "Iterator = " + util::Format()( *itr));
        BCL_MessageDbg
        (
          "Equal " + util::Format()( *found)
          + " == " + util::Format()( *itr)
          + " -> " + util::Format()( *found == *itr)
        );
        BCL_MessageDbg
        (
          "Equal " + util::Format()( *found)
          + " != " + util::Format()( *itr)
          + " -> " + util::Format()( *found != *itr)
        );
      }
      if( found == GetSomeDataTypes().End())
      {
        BCL_MessageDbg( "Index not found");
      }
      else
      {
        BCL_MessageDbg( "Enum=" + util::Format()( *found));
      }
      BCL_Example_Check
      (
        found != GetSomeDataTypes().End() && found->GetIndex() == 0, "No enumerator with value of 0 found!"
      );

      // GetEnumIteratorFromName
      found = GetSomeDataTypes().GetEnumIteratorFromName( std::string( "ASN"));
      if( found == GetSomeDataTypes().End())
      {
        BCL_MessageDbg( "Object not found");
      }
      else
      {
        BCL_MessageDbg( "Enum=" + util::Format()( *found));
      }

    //////////////////////
    // input and output //
    //////////////////////

      // write all enums to file
      io::OFStream write;
      std::string out_filename( AddExampleOutputPathToFilename( util::GetNamespaceIdentifier(), "somedata.enum"));
      BCL_ExampleMustOpenOutputFile( write, out_filename);
      write << GetSomeDataTypes();
      io::File::CloseClearFStream( write);

      // read additional enums to file
      io::IFStream read;
      std::string in_filename( AddExampleInputPathToFilename( e_Biology, "somedata_more.enum"));
      BCL_ExampleMustOpenInputFile( read, in_filename);

      read >> GetSomeDataTypes();

      io::File::CloseClearFStream( read);

      BCL_MessageStd
      (
        "changed enum XXX: " + util::Format()( SomeDataTypes::EnumDataType( SomeDataTypes::EnumType( "XXX")))
      );
      BCL_Example_Check
      (
        SomeDataTypes::EnumType( "XXX")->GetValue() == 1.5,
        "changed value for XXX should be 1.5 but is: " + util::Format()( SomeDataTypes::EnumType( "XXX")->GetValue())
      );
      BCL_MessageStd
      (
        "added enum ADD: " + util::Format()( SomeDataTypes::EnumDataType( SomeDataTypes::EnumType( "ADD")))
      );
      BCL_Example_Check
      (
        SomeDataTypes::EnumType( "ADD")->GetValue() == 9.9,
        "value for new enum ADD should be 9.9 but is: " + util::Format()( SomeDataTypes::EnumType( "ADD")->GetValue())
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAATypeData

  const ExampleClass::EnumType ExampleUtilEnum::s_Instance
  (
    GetExamples().AddEnum( ExampleUtilEnum())
  );

} // namespace bcl
