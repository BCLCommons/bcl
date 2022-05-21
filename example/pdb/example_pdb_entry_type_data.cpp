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
#include "pdb/bcl_pdb_entry_type_data.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_pdb_entry_type_data.cpp
  //!
  //! @author riddeljs
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExamplePdbEntryTypeData :
    public ExampleInterface
  {
  public:

    ExamplePdbEntryTypeData *Clone() const
    {
      return new ExamplePdbEntryTypeData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
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
      pdb::EntryTypeData entry_type_data_a;
      BCL_Example_Check
      (
        entry_type_data_a.GetLineType() == pdb::GetLineTypes().e_Undefined
          && entry_type_data_a.GetStart() == util::GetUndefined< size_t>()
          && entry_type_data_a.GetLength() == util::GetUndefined< size_t>()
          && entry_type_data_a.GetDataType() == util::CPPDataTypes::e_Unknown,
        "default constructor of entry_type_data_a failed."
      );

      // test constructor given Line Type, Start, Length, and DataType
      pdb::EntryTypeData entry_type_data_b
      (
        pdb::GetLineTypes().LINK,
        0,
        5,
        util::GetUndefined< size_t>(),
        false,
        util::CPPDataTypes::e_Float
      );
      BCL_Example_Check
      (
        entry_type_data_b.GetLineType() == pdb::GetLineTypes().LINK
          && entry_type_data_b.GetStart() == 0
          && entry_type_data_b.GetLength() == 5
          && entry_type_data_b.GetDataType() == util::CPPDataTypes::e_Float,
        "constructing entry_type_data_b from line typek, start, length, datatype failed."
      );

      // test copy constructor
      pdb::EntryTypeData entry_type_data_c( entry_type_data_b);
      BCL_Example_Check
      (
        entry_type_data_c.GetLineType() == entry_type_data_b.GetLineType()
          && entry_type_data_c.GetStart() == entry_type_data_b.GetStart()
          && entry_type_data_c.GetLength() == entry_type_data_b.GetLength()
          && entry_type_data_c.GetDataType() == entry_type_data_b.GetDataType(),
        "copy constructor of entry_type_data_c failed."
      );

      // test clone
      util::ShPtr< pdb::EntryTypeData> sp_entry_type_data_d( entry_type_data_c.Clone());
      BCL_Example_Check
      (
        entry_type_data_c.GetLineType() == sp_entry_type_data_d->GetLineType()
          && entry_type_data_c.GetStart() == sp_entry_type_data_d->GetStart()
          && entry_type_data_c.GetLength() == sp_entry_type_data_d->GetLength()
          && entry_type_data_c.GetDataType() == sp_entry_type_data_d->GetDataType(),
        "clone constructor of sp_entry_type_data_d failed."
      );

    /////////////////
    // data access //
    /////////////////

      // Get methods don't need to be tested, as if the above checks succeeded so did the gets.
      BCL_MessageStd
      (
        " entry_type_data_b Line Type: " + util::Format()( entry_type_data_b.GetLineType())
      );
      BCL_MessageStd
      (
        " entry_type_data_b Start: " + util::Format()( entry_type_data_b.GetStart())
      );
      BCL_MessageStd
      (
        " entry_type_data_b Length: " + util::Format()( entry_type_data_b.GetLength())
      );
      BCL_MessageStd
      (
        " entry_type_data_b Data Type: " + util::Format()( entry_type_data_b.GetDataType())
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

      // write to file
      WriteBCLObject( entry_type_data_b);
      // test reading entry_type_data from file
      pdb::EntryTypeData entry_type_data_e;
      ReadBCLObject( entry_type_data_e);
      BCL_MessageStd( util::Format()( entry_type_data_b));
      BCL_MessageStd( util::Format()( entry_type_data_e));
      BCL_Example_Check
      (
        entry_type_data_e.GetLineType() == entry_type_data_b.GetLineType()
          && entry_type_data_e.GetStart() == entry_type_data_b.GetStart()
          && entry_type_data_e.GetLength() == entry_type_data_b.GetLength()
          && entry_type_data_e.GetDataType() == entry_type_data_b.GetDataType(),
        "reading and writing entry_type_datas to and from files did not work"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExamplePdbEntryTypeData

  const ExampleClass::EnumType ExamplePdbEntryTypeData::s_Instance
  (
    GetExamples().AddEnum( ExamplePdbEntryTypeData())
  );

} // namespace bcl

