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
#include "storage/bcl_storage_map.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_storage_map.cpp
  //!
  //! @author alexanns
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleStorageMap :
    public ExampleInterface
  {
  public:

    ExampleStorageMap *Clone() const
    {
      return new ExampleStorageMap( *this);
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

    ////////////////////////
    // test map functions //
    ////////////////////////

      // test default constructor
      storage::Map< std::string, int> def_const;

      // test InsertElement function with standard pair
      BCL_Example_Check
      (
        ( def_const.Insert( storage::Pair< std::string, int>( "b", 2))).second, "insert failed"
      );

      // test InsertElement function with standard pair
      BCL_Example_Check
      (
        ( def_const.Insert( storage::Pair< std::string, int>( "d", 2))).second, "insertion failed"
      );

      // test InsertElement function with standard pair
      BCL_Example_Check
      (
        ( def_const.Insert( storage::Pair< std::string, int>( "a", 1))).second, "insertion failed"
      );

      // test InsertElement function with standard pair
      BCL_Example_Check
      (
        ( def_const.Insert( storage::Pair< std::string, int>( "c", 3))).second, "insertion failed"
      );

      BCL_MessageStd( "def_const" + util::Format()( def_const));

      // test copy constructor
      storage::Map< std::string, int> copy_const( def_const);
      BCL_MessageStd( "copy_const" + util::Format()( copy_const));
      BCL_Example_Check
      (
        copy_const == def_const, "copy constructor did not properly construct copy"
      );

      // test clone copy constructor
      util::ShPtr< storage::Map< std::string, int> > virtual_copy( def_const.Clone());
      BCL_Example_Check
      (
        *virtual_copy == def_const, "clone copy constructor failed"
      );

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

      // test [] operator
      BCL_Example_Check
      (
        def_const[ "c"] == 3, "[] operator did not return correct data"
      );

    ////////////////
    // operations //
    ////////////////

    /////////////////////////////////////////////////
    // test sorted associative container functions //
    /////////////////////////////////////////////////

      // test get iterator to reverse begin
      storage::Map< std::string, int>::reverse_iterator reverse_itr_begin( def_const.ReverseBegin());
      BCL_Example_Check
      (
        reverse_itr_begin->first == "d", "unable to get correct reversible begin iterator"
      );

      // test get const iterator to reverse begin
      storage::Map< std::string, int>::const_reverse_iterator reverse_const_itr_begin
      (
        def_const.ReverseBegin()
      );
      BCL_Example_Check
      (
        reverse_const_itr_begin->first == "d", "unable to get correct const reversible begin iterator"
      );

      // test get iterator to reverse end
      storage::Map< std::string, int>::reverse_iterator reverse_itr_end( def_const.ReverseEnd());
      BCL_Example_Check
      (
        ( --reverse_itr_end)->first == "a", "unable to get correct reversible end iterator"
      );

      // test get const iterator to reverse end
      storage::Map< std::string, int>::const_reverse_iterator reverse_const_itr_end( def_const.ReverseEnd());
      BCL_Example_Check
      (
        ( --reverse_const_itr_end)->first == "a", "unable to get correct reversible end iterator"
      );

      // test lowerbound function
      def_const.InsertElement( storage::Pair< std::string, int>( "f", 17));
      def_const.InsertElement( storage::Pair< std::string, int>( "i", 3));
      BCL_Example_Check
      (
        def_const.LowerBound( "e")->first == "f", "did not give correct lowerbound"
      );

      // test lowerbound function giving const iterator
      BCL_Example_Check
      (
        def_const.LowerBound( "e")->first == "f", "did not give correct const iterator lowerbound"
      );

      // test upperbound function
      BCL_Example_Check
      (
        def_const.UpperBound( "e")->first == "f", "did not give correct upper bound"
      );

      // test upperbound function giving const iterator
      BCL_Example_Check
      (
        def_const.UpperBound( "e")->first == "f", "did not give correct const iterator upper bound"
      );

    /////////////////////////////////////////////////
    // test unique associative container functions //
    /////////////////////////////////////////////////

      // test construction from iterators
      storage::Map<  std::string, int>::iterator itr_begin( def_const.Begin());
      storage::Map<  std::string, int>::iterator itr_end( def_const.End());
      storage::Map<  std::string, int> itr_range_map( itr_begin, itr_end);
      BCL_Example_Check
      (
        itr_range_map.GetSize() == def_const.GetSize(), "construction from range of iterators did not work"
      );
      BCL_Example_Check
      (
        *( --( itr_range_map.ReverseEnd())) == *( def_const.Begin()), "range of iterators constr did not work"
      );
      BCL_MessageStd( util::Format()( itr_range_map));

      // test insert a bcl pair
      BCL_Example_Check
      (
        ( def_const.Insert( storage::Pair< std::string, int>( "l", 7))).second, "insert failed"
      );

      // test insert elements with two iterators
      copy_const.InsertElements( ( ++( ++( ++( ++itr_begin)))), itr_end);
      BCL_MessageStd( util::Format()( copy_const));

      // test insert elements of another unique associative container
      storage::Map< std::string, int> another_map;
      another_map.Insert( storage::Pair< std::string, int>( "m", 13));
      another_map.Insert( storage::Pair< std::string, int>( "p", 6));
      BCL_MessageStd( "before insert elements of another map" + util::Format()( copy_const));
      copy_const.InsertElements( another_map);
      BCL_MessageStd( "after insert elements of another map" + util::Format()( copy_const));
      BCL_Example_Check
      (
        copy_const.GetSize() == 9 && copy_const[ "m"] == 13 && copy_const[ "p"] == 6, "Insert map is broken"
      );

    //////////////////////////////////////////
    // test associative container functions //
    //////////////////////////////////////////

      // test count function
      BCL_Example_Check
      (
        copy_const.Count( "m") == 1, "count function did not work properly"
      );

      // test Erase function
      def_const.Erase( "a");
      BCL_Example_Check
      (
        def_const.Count( "a") == 0 && def_const.GetSize() == 6, "erase function did not work"
      );

      // test remove element with iterator
      def_const.RemoveElement( def_const.Begin());
      BCL_Example_Check
      (
        def_const.Count( "b") == 0 && def_const.GetSize() == 5, "remove element did not work"
      );

      // test erase all elements between two iterators
      storage::Map<  std::string, int>::iterator itr_A( def_const.Begin());
      storage::Map<  std::string, int>::iterator itr_B( def_const.End());
      ( ++( ++itr_A));
      BCL_MessageStd( "before erase elements between two itrs" + util::Format()( def_const));
      def_const.Erase( itr_A, itr_B);
      BCL_MessageStd( "after erase elements between two itrs" + util::Format()( def_const));
      BCL_Example_Check
      (
        def_const.GetSize() == 2, "Erase elements between two iterators did not work"
      );

      // test clear function
      def_const.Reset();
      BCL_Example_Check
      (
        def_const.IsEmpty(), "Clear function did not work");

      // test EqualRange function
      BCL_Example_Check
      (
        copy_const.EqualRange( "m").first->first == "m", "equal range function did not work"
      );

      // test max size function
      BCL_MessageStd( "max size of map = " + util::Format()( copy_const.MaxSize()));

      // test IsEmpty function
      BCL_Example_Check
      (
        def_const.IsEmpty(), "IsEmpty function did not work properly"
      );

      // test Swap function
      storage::Map<  std::string, int> swap_A;
      storage::Map<  std::string, int> swap_B;
      swap_A.Insert( storage::Pair< std::string, int>( "a", 1));
      swap_A.Insert( storage::Pair< std::string, int>( "b", 2));
      swap_B.Insert( storage::Pair< std::string, int>( "c", 3));
      swap_B.Insert( storage::Pair< std::string, int>( "d", 4));
      swap_A.Swap( swap_B);
      BCL_Example_Check
      (
        swap_A[ "c"] == 3 && swap_A[ "d"] && swap_B[ "a"] == 1 && swap_B[ "b"] == 2, "swap function is broken"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // Write
      // test write function
      WriteBCLObject( swap_A);
      // test read function
      storage::Map< std::string, int> map_read;
      // difficulty testing Read because keytype is implicitly const and written out as so, but serialization can't
      // read it back in since nonconst object is expected

      //ReadBCLObject( map_read);
      //BCL_ExampleIndirectCheck( swap_A, map_read, "I/O");

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleStorageMap

  const ExampleClass::EnumType ExampleStorageMap::s_Instance
  (
    GetExamples().AddEnum( ExampleStorageMap())
  );

} // namespace bcl
