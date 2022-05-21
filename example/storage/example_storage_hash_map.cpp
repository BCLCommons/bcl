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
#include "storage/bcl_storage_hash_map.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_storage_hash_map.cpp
  //!
  //! @author karakam
  //! @date Sep 08, 2010
  //! @remarks status complete
  //! @remarks reviewed by fischea on Dec 2, 2015
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleStorageHashMap :
    public ExampleInterface
  {
  public:

    ExampleStorageHashMap *Clone() const
    { return new ExampleStorageHashMap( *this);}

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
      // create "insert_success" for testing insertions
      bool insert_success;

      // create a hashmap that stores a size_t for a double
      storage::HashMap< size_t, double> map;

      BCL_MessageStd( "Inserting (1, 0.1)");
      std::pair< storage::HashMap< size_t, double>::iterator, bool> insert_result
      (
        map.Insert( storage::Pair< size_t, double>( 1, 0.1))
      );

      // insert key and value
      BCL_MessageStd
      (
        "Insert( 1, 0.1): " +
        util::Format()( insert_result.first->first) + " and " + util::Format()( insert_result.first->second)
      );
      BCL_Example_Check
      (
        map.Find( 1)->first == 1 && map.Find( 1)->second == 0.1, "should be ( 1, 0.1)"
      );

      // test should not overwrite previously stored entry with key 1 and therefore should fail to insert
      insert_success = map.Insert( storage::Pair< size_t, double>( 1, 0.2)).second;
      BCL_MessageStd( " Try to Insert( 1, 0.2): fails " + util::Format()( storage::Pair< size_t, double>( *map.Find( 1))));
      BCL_Example_Check
      (
        !insert_success, "should fail to insert second instance of key 1"
      );

      // insert key and value
      insert_success = map.Insert( storage::Pair< size_t, double>( 2, 0.4)).second;
      BCL_Example_Check
      (
        insert_success, "should insert key 2 and value 0.4 " + util::Format()( map)
      );

      // test should not overwrite previously store entry with key 2 and therefore should fail to insert
      insert_success = map.Insert( storage::Pair< size_t, double>( 2, 0.5)).second;
      BCL_Example_Check
      (
        !insert_success, "should fail to insert second instance of key 2 and value 0.5" + util::Format()( map)
      );

      // insert key and value
      insert_success = map.Insert( storage::Pair< size_t, double>( 3, 0.6)).second;
      BCL_Example_Check
      (
        insert_success, "should insert key 3 and value 0.6"
      );

      // insert a std:pair that contains the key and value
      insert_success = map.Insert( std::make_pair( 4, 0.7)).second;
      BCL_Example_Check
      (
        insert_success, "should insert key 4 and value 0.7"
      );

      // test constructor from two iterators
      storage::HashMap< size_t, double> constr_itr( map.Begin(), map.End());
//      BCL_Example_Check
//      (
//        ExampleClass::ExampleResult::e_Trivial,
//        constr_itr == map, "construction from two iterators did not work"
//      );

      // test copy constructor
      storage::HashMap< size_t, double> constr_copy( constr_itr);
//      BCL_Example_Check
//      (
//        ExampleClass::ExampleResult::e_Trivial,
//        constr_itr == constr_copy, "copy constructor did not work"
//      );

      // test clone constructor
      util::ShPtr< storage::HashMap< size_t, double> > virtual_copy( map.Clone());
//      BCL_Example_Check
//      (
//        ExampleClass::ExampleResult::e_Trivial,
//        *virtual_copy == map, "clone copy constructor failed"
//      );

      // call the GetSize() function
      BCL_MessageStd( " Size of hashmap: " + util::Format()( map.GetSize()));
      BCL_Example_Check
      (
        map.GetSize() == 4, "size of map does not match - should be 4 but it is " + util::Format()( map.GetSize())
      );

      // iterate over elements and print out size_t and stored double values using Begin(), End() functions
      for( storage::HashMap< size_t, double>::const_iterator itr( map.Begin()), itr_end( map.End()); itr != itr_end; ++itr)
      {
        BCL_MessageStd( util::Format()( itr->first) + " ==> " + util::Format()( itr->second));
      }

      // find the element with key 4 using  Find function
      BCL_MessageStd
      (
        " looking for value 0.7 stored with key 4 and found for 4 the value : " + util::Format()( map.Find( 4)->second)
      );
      BCL_Example_Check
      (
        map.Find( 4)->second == 0.7, "finding key four should give value of 0.7"
      );

      // erase the element with key 4, and try to find it again
      BCL_MessageStd( " erasing element 4 using key");
      map.Erase( 4);
      BCL_MessageStd
      (
        " search for key 4 should give iterator past the end ( 1 for yes, 0 for no):  " +
        util::Format()( map.Find( 4) == map.End())
      );
      BCL_Example_Check
      (
        map.Find( 4) == map.End(), "entry 4 was not deleted correctly"
      );

      // erase the element with key 3, using the other Erase function which takes in iterator
      BCL_MessageStd( " erasing element 3 using iterator");
      map.RemoveElement( map.Find( 3));
      BCL_MessageStd
      (
        " search for key 3 should give iterator past the end ( 1 for yes, 0 for no): " +
        util::Format()( map.Find( 3) == map.End())
      );
      BCL_Example_Check
      (
        map.Find( 3) == map.End(), "entry 3 was not deleted correctly"
      );

      // check [] operator: accessing map with key that does not exist creates key and assigned value
      BCL_MessageStd( " trying map[ 5] = 0.9 ");
      map[ 5] = 0.9;
      BCL_MessageStd( " now looking for key 5 " + util::Format()( storage::Pair< size_t, double>( *map.Find( 5))));
      BCL_Example_Check
      (
        map[ 5] == 0.9, "inserted 5 with value 0.9 did not work"
      );

      // test that [] operator can be used to overwrite a previously existant value associated with a given key
      BCL_MessageStd( " trying map[ 5] = 1.0, it should overwrite ");
      map[ 5] = 1.0;
      BCL_MessageStd( " now looking for key 5 " + util::Format()( storage::Pair< size_t, double>( *map.Find( 5))));
      BCL_Example_Check
      (
        map[ 5] == 1.0, "assigned 5 with value 1.0 did not work"
      );

      // test Append function
      BCL_MessageStd( " testing Append: map size before Append " + util::Format()( map.GetSize()));
      storage::HashMap< size_t, double> map_B;
      map_B.Insert( storage::Pair< size_t, double>( 11, 1.1));
      map_B.Insert( storage::Pair< size_t, double>( 12, 1.2));
      map_B.Insert( storage::Pair< size_t, double>( 13, 1.3));
      map.InsertElements( map_B.Begin(), map_B.End());
      // iterate over elements and print out size_t and stored double values of Appended map
      for( storage::HashMap< size_t, double>::const_iterator itr( map.Begin()), itr_end( map.End()); itr != itr_end; ++itr)
      {
        BCL_MessageStd
        (
          " map after Append: " + util::Format()( itr->first) + " ==> " + util::Format()( itr->second)
        );
      }

      // test reset function and IsEmpty Function
      BCL_MessageStd( "resetting the map so that it is empty");
      map.Reset();
      BCL_MessageStd( "the resulting size is: " + util::Format()( map.GetSize()));
      BCL_Example_Check
      (
        map.IsEmpty(), "resetting the map failed!"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // end ExampleStorageHashMap

  const ExampleClass::EnumType ExampleStorageHashMap::s_Instance
  (
    GetExamples().AddEnum( ExampleStorageHashMap())
  );

} // namespace bcl
