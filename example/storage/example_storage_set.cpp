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
#include "storage/bcl_storage_set.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_storage_set.cpp
  //!
  //! @author heinzes1
  //! @date Nov 6, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleStorageSet :
    public ExampleInterface
  {
  public:

    ExampleStorageSet *Clone() const
    { return new ExampleStorageSet( *this);}

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
    ////////////////////////
    // test set functions //
    ////////////////////////

      // test default constructor
      storage::Set< int> def_const;

      // test InsertElement function with standard pair
      BCL_Example_Check
      (
        ( def_const.Insert( 2)).second, "insert failed"
      );

      // test InsertElement function with standard pair
      BCL_Example_Check
      (
        ( def_const.Insert( 1)).second, "insertion failed"
      );

      // test InsertElement function with standard pair
      BCL_Example_Check
      (
        ( def_const.Insert( 5)).second, "insertion failed"
      );

      // test InsertElement function with standard pair
      BCL_Example_Check
      (
        ( def_const.Insert( 3)).second, "insertion failed"
      );

      BCL_MessageStd( "def_const" + util::Format()( def_const));

      // test construct from single initial object
      storage::Set< int> single_const( 4);
      BCL_Example_Check
      (
        *single_const.Begin() == 4,
        "constructing from single initial object did not work properly, first element is "
        + util::Format()( *single_const.Begin()) + " but should be 4"
       );

      // test construct from two initial objects
      storage::Set< int> double_const( 5, 6);
      BCL_Example_Check
      (
        *double_const.Begin() == 5
        && *( ++( double_const.Begin())) == 6,
        "constructing from two initial objects did not work properly. First is "
        + util::Format()( *double_const.Begin()) + " and second is "
        + util::Format()( *( ++( double_const.Begin())))
      );

      // test copy constructor
      storage::Set< int> copy_const( def_const);
      BCL_MessageStd( "copy_const" + util::Format()( copy_const));
//      BCL_Example_Check
//      (
//        ExampleClass::ExampleResult::e_Trivial,
//        copy_const == def_const, "copy constructor did not properly construct copy"
//      );

      // test clone copy constructor
      util::ShPtr< storage::Set< int> > virtual_copy( def_const.Clone());
//      BCL_Example_Check
//      (
//        ExampleClass::ExampleResult::e_Trivial,
//        *virtual_copy == def_const, "clone copy constructor failed"
//      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleStorageSet

  const ExampleClass::EnumType ExampleStorageSet::s_Instance
  (
    GetExamples().AddEnum( ExampleStorageSet())
  );

} // namespace bcl
