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
#include "sdf/bcl_sdf_mdl_entry_type_data.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_sdf_mdl_entry_type_data.cpp
  //!
  //! @author butkiem1, alexanns
  //! @date May 14, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleSdfMdlEntryTypeData :
    public ExampleInterface
  {
  public:

    ExampleSdfMdlEntryTypeData *Clone() const
    {
      return new ExampleSdfMdlEntryTypeData( *this);
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

      //! default constructor
      sdf::MdlEntryTypeData entry_type;

      //! constructor with parameters
      sdf::MdlEntryTypeData entry_type_param
      (
        sdf::e_HeaderLine,
        size_t( 0),
        size_t( 5),
        std::string( "99"),
        util::Format().R(),
        util::CPPDataTypes::e_SizeT
      );

      //! test clone
      util::ShPtr< sdf::MdlEntryTypeData> sp_entry_type( entry_type_param.Clone());

    /////////////////
    // data access //
    /////////////////

      // check GetMdlLineType
      BCL_ExampleCheck( entry_type_param.GetMdlLineType(), sdf::e_HeaderLine);

      // check GetStart
      BCL_ExampleCheck( entry_type_param.GetStart(), size_t( 0));

      // check GetLength
      BCL_ExampleCheck( entry_type_param.GetLength(), size_t( 5));

      // check GetDefault
      BCL_ExampleCheck( entry_type_param.GetDefault(), std::string( "99"));

      // check GetFormat
      BCL_ExampleCheck( entry_type_param.GetFormat()( size_t( 11)), util::Format().W( 5).R()( size_t( 11)));

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // write score_accessibility to file
      WriteBCLObject( entry_type_param);
      // read
      sdf::MdlEntryTypeData entry_type_read;
      ReadBCLObject( entry_type_read);

      // check GetMdlLineType
      BCL_ExampleCheck( entry_type_read.GetMdlLineType(), sdf::e_HeaderLine);

      // check GetStart
      BCL_ExampleCheck( entry_type_read.GetStart(), size_t( 0));

      // check GetLength
      BCL_ExampleCheck( entry_type_read.GetLength(), size_t( 5));

      // check GetDefault
      BCL_ExampleCheck( entry_type_read.GetDefault(), std::string( "99"));

      // check GetFormat
      BCL_ExampleCheck( entry_type_param.GetFormat()( size_t( 11)), util::Format().W( 5).R()( size_t( 11)));

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSdfMdlEntryTypeData

  const ExampleClass::EnumType ExampleSdfMdlEntryTypeData::s_Instance
  (
    GetExamples().AddEnum( ExampleSdfMdlEntryTypeData())
  );

} // namespace bcl
