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
#include "storage/bcl_storage_list.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //! binary predicate for inverse sorting, used by Sort
  template< class T> struct more_than :
    public std::binary_function< T, T, bool>
  {
    bool operator()( const T &x, const T &y) const
    { return x > y;}
  };

  //! binary predicate for making all elements equal by Unique
  template< class T> struct all_equal :
    public std::binary_function< T, T, bool>
  {
    bool operator()( const T &x, const T &y) const
    { return true;}
  };

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_storage_list.cpp
  //! @details Tests ~90% of the functions in the list class
  //!
  //! @author alexanns
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by Sten on ?
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleStorageList :
    public ExampleInterface
  {
  public:

    ExampleStorageList *Clone() const
    { return new ExampleStorageList( *this);}

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    { return GetStaticClassName( *this);}

    int Run() const
    {
      // size of c1 array with doubles
      static const size_t s_c1_size( 15);

      // initialize double array
      static const double s_c1[ s_c1_size] =
      {
         0.0,  -1.0,   2.0,
        -3.0,   4.0,  -5.0,
         6.0,  -7.0,   8.0,
        -9.0,  10.0, -11.0,
        12.0, -13.0,  14.0
      };

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      storage::List< double> list_A;
      BCL_MessageStd( "list_A: should be empty ");
      for( storage::List< double>::const_iterator ptr( list_A.Begin()), ptr_end( list_A.End()); ptr != ptr_end; ++ptr)
      {
        BCL_MessageStd( "list_A: " + util::Format()( *ptr));
      }
      BCL_Example_Check
      (
        list_A.IsEmpty(), "list_A size should be 0 but is not"
      );

      // construct with given size
      storage::List< double> list_B( 4);
      BCL_MessageStd( "list_B: should be all 0 and have size 4");
      BCL_Example_Check
      (
        list_B.GetSize() == 4, "list_B size should be 4 but is not"
      );
      for( storage::List< double>::const_iterator ptr( list_B.Begin()), ptr_end( list_B.End()); ptr != ptr_end; ++ptr)
      {
        BCL_MessageStd( "list_B: " + util::Format()( *ptr));
        BCL_Example_Check
        (
          *ptr == 0, "upon construction all four elements of list_B should be 0"
        );
      }

      // construct with given size and default value
      storage::List< double> list_C( 4, 3.1415);
      BCL_MessageStd( "list_C: should be all 3.1415");
      BCL_Example_Check
      (
        list_B.GetSize() == 4, "list_B size should be 4 but is not"
      );
      for( storage::List< double>::const_iterator ptr( list_C.Begin()), ptr_end( list_C.End()); ptr != ptr_end; ++ptr)
      {
        BCL_MessageStd( "list_C: " + util::Format()( *ptr));
        BCL_Example_Check
        (
          *ptr == 3.1415, "upon construction all four elements of list_C should be 3.1415"
        );
      }

      // copy constructor
      storage::List< double> list_D( list_B);
      BCL_MessageStd( "list_D: should be copy of list_B");
      for
      (
        storage::List< double>::const_iterator ptr_D( list_D.Begin()), ptr_D_end( list_D.End()), ptr_B( list_B.Begin());
        ptr_D != ptr_D_end;
        ++ptr_D, ++ptr_B
      )
      {
        BCL_MessageStd( "list_D: " + util::Format()( *ptr_D));
        BCL_Example_Check
        (
          *ptr_D == *ptr_B, "elements of list_D should be the same as each element of list_B but is not"
        );
      }

      // construct from size and pointer to data
      storage::List< double> list_E( s_c1_size, s_c1);
      BCL_MessageStd( "list_E: should be copy of c1");
      size_t index( 0);
      for( storage::List< double>::const_iterator ptr( list_E.Begin()), ptr_end( list_E.End()); ptr != ptr_end; ++ptr)
      {
        BCL_MessageStd( "list_E: " + util::Format()( *ptr));
        BCL_Example_Check
        (
          *ptr == s_c1[ index], "elements of list_E should be the same as each element of c1 but is not"
        );
        ++index;
      }

      // construct from range denoted by iterators
      storage::List< double> list_F( ++list_E.Begin(), --( --list_E.End()));
      BCL_MessageStd( "list_F: should be range of list_E");
      BCL_Example_Check
      (
        list_F.GetSize() == list_E.GetSize() - 3,
        "size of list_F should be 3 less than size of list_E but is not list_E: " +
        util::Format()( list_E.GetSize()) + " list_F: " + util::Format()( list_F.GetSize())
      );
      for
      (
        storage::List< double>::const_iterator
         ptr_F( list_F.Begin()), ptr_F_end( list_F.End()), ptr_E( ++list_E.Begin());
        ptr_F != ptr_F_end;
        ++ptr_F, ++ptr_E
      )
      {
        BCL_MessageStd( "list_F: " + util::Format()( *ptr_F));
        BCL_Example_Check
        (
          *ptr_F == *ptr_E, "elements of list_F should be the same as range of list_E but are not"
        );
      }

      // test clone function
      list_A = *util::ShPtr< storage::List< double> >( list_C.Clone());
      BCL_MessageStd( "list_A: should be clone of list_C now");
      for
      (
        storage::List< double>::const_iterator ptr_A( list_A.Begin()), ptr_A_end( list_A.End()), ptr_C( list_C.Begin());
        ptr_A != ptr_A_end;
        ++ptr_A, ++ptr_C
      )
      {
        BCL_MessageStd( "list_A: " + util::Format()( *ptr_A));
        BCL_Example_Check
        (
          *ptr_A == *ptr_C, "list_a should be clone of list_C but is not"
        );
      }

    /////////////////
    // data access //
    /////////////////

      // test GetClassIdentifier function
      BCL_MessageStd( "list_A: GetClassIdentifier: " + list_A.GetClassIdentifier());

      // test GetSize function on list_A and list_B
      BCL_MessageStd( "list_A: Size: " + util::Format()( list_A.GetSize()));
      BCL_Example_Check
      (
        list_A.GetSize() == 4, "list_A size should be 4 but is not"
      );
      BCL_MessageStd( "list_B: Size: " + util::Format()( list_B.GetSize()));
      BCL_Example_Check
      (
        list_B.GetSize() == 4, "list_B size should be 4 but is not"
      );

      // check MaxSize
      BCL_MessageStd( "list_A: MaxSize: " + util::Format()( list_A.MaxSize()));
      BCL_MessageStd( "list_B: MaxSize: " + util::Format()( list_B.MaxSize()));

      // get InternalData (the std::list) of list_B and initialize std::list with it
      std::list< double> stdlist( list_B.InternalData());
      BCL_MessageStd( "list_B: stdlist = InternalData: ");
      for
      (
        std::list< double>::const_iterator stdptr( stdlist.begin()), stdptr_end( stdlist.end()), ptr_B( list_B.Begin());
        stdptr != stdptr_end;
        ++stdptr, ++ptr_B
      )
      {
        BCL_MessageStd( "stdlist: " + util::Format()( *stdptr));
        BCL_Example_Check
        (
          *stdptr == *ptr_B, "stdlist should be same as list_B but is not"
        );
      }

      // set InternalData of list_B using altered "stdlist"
      *stdlist.begin() = 1.41421;
      BCL_MessageStd( "stdlist: set first element");
      BCL_Example_Check
      (
        *stdlist.begin() == 1.41421, "first element of stdlist should be 0 but is not"
      );
      for
      (
        std::list< double>::const_iterator ptr( stdlist.begin()), ptr_end( stdlist.end());
        ptr != ptr_end; ++ptr
      )
      {
        BCL_MessageStd( "stdlist: " + util::Format()( *ptr));
      }
      list_B.InternalData() = stdlist;
      BCL_MessageStd( "list_B: set InternalData of list_B equal to stdlist : ");
      for
      (
        std::list< double>::const_iterator ptr_B( list_B.Begin()), ptr_B_end( list_B.End()), stdptr( stdlist.begin());
        ptr_B != ptr_B_end;
        ++ptr_B, ++stdptr
      )
      {
        BCL_MessageStd( "list_B: " + util::Format()( *ptr_B));
        BCL_Example_Check
        (
          *ptr_B == *stdptr, "list_B should be the same as stdlist but is not"
        );
      }

      // check iterators of ReverseBegin, ReverseEnd
      BCL_MessageStd( "list_B: reverse iteration: ");
      BCL_Example_Check
      (
        *list_B.ReverseBegin() == *list_B.Last(),
        "reverse begin iterators do not match the corresponding forward end iterators"
      );
      BCL_Example_Check
      (
        *( --( list_B.ReverseEnd())) == *list_B.Begin(),
        "reverse end iterators do not match the corresponding forward begin iterators"
      );
      for
      (
        storage::List< double>::reverse_iterator ptr( list_B.ReverseBegin()), ptr_end( list_B.ReverseEnd());
        ptr != ptr_end; ++ptr
      )
      {
        BCL_MessageStd( "list_B: " + util::Format()( *ptr));
      }

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // check FirstElement gives the first element of list_B
      BCL_MessageStd( "list_B: get FirstElement: " + util::Format()( list_B.FirstElement()));
      BCL_Example_Check
      (
        list_B.FirstElement() == 1.41421, "first element of list_B should be 1.41421 but is not"
      );

      // check that FirstElement function can be used to assign first element to a value
      BCL_MessageStd( "list_B: set FirstElement: ");
      list_B.FirstElement() = 0.707106;
      BCL_Example_Check
      (
        *list_B.Begin() == 0.707106, "first element of list_B should be 0.707106 but is not"
      );
      for( storage::List< double>::const_iterator ptr( list_B.Begin()), ptr_end( list_B.End()); ptr != ptr_end; ++ptr)
      {
        BCL_MessageStd( "list_B: " + util::Format()( *ptr));
      }

      // check LastElement function can be used to assign last element to a value
      BCL_MessageStd( "list_B: get LastElement: " + util::Format()( list_B.LastElement()));
      BCL_MessageStd( "list_B: set LastElement: ");
      list_B.LastElement() = 0.707106;
      BCL_Example_Check
      (
        *list_B.Last() == 0.707106, "last element should be 0.707106 but is not"
      );
      for( storage::List< double>::const_iterator ptr( list_B.Begin()), ptr_end( list_B.End()); ptr != ptr_end; ++ptr)
      {
        BCL_MessageStd( "list_B: " + util::Format()( *ptr));
      }

      // check IsEmpty and Reset functions
      BCL_MessageStd( "list_B: IsEmpty: " + util::Format()( list_B.IsEmpty()));
      BCL_Example_Check
      (
        !list_B.IsEmpty(), "IsEmpty should give false but does not"
      );
      BCL_MessageStd( "list_B: Reset: ");
      list_B.Reset();
      BCL_Example_Check
      (
        list_B.GetSize() == 0 && list_B.IsEmpty(), "Size should be 0 and IsEmpty true but are not"
      );
      BCL_MessageStd( "list_B: IsEmpty: " + util::Format()( list_B.IsEmpty()));

      // test InsertElement with t_DataType
      BCL_MessageStd( "list_D: InsertElement 2.71828");
      list_D.PushBack( 2.71828);
      BCL_Example_Check
      (
        list_D.LastElement() == 2.71828, "Last element of list_D should be 2.71828 but is not"
      );
      for( storage::List< double>::const_iterator ptr( list_D.Begin()), ptr_end( list_D.End()); ptr != ptr_end; ++ptr)
      {
        BCL_MessageStd( "list_D: " + util::Format()( *ptr));
      }

      // test InsertElement with iterator and t_DataType
      BCL_MessageStd( "list_D: InsertElement 0.707106 as first element");
      list_D.InsertElement( list_D.Begin(), 0.707106);
      BCL_Example_Check
      (
        *list_D.Begin() == 0.707106, "first element of list_D should be 0.707106 but is not"
      );
      for( storage::List< double>::const_iterator ptr( list_D.Begin()), ptr_end( list_D.End()); ptr != ptr_end; ++ptr)
      {
        BCL_MessageStd( "list_D: " + util::Format()( *ptr));
      }

      // test InsertElements with another list
      BCL_MessageStd( "list_D: InsertElements list_C");
      list_D.InsertElements( list_D.End(), list_C);
      for( storage::List< double>::const_iterator ptr( list_C.Begin()), ptr_end( list_C.End()); ptr != ptr_end; ++ptr)
      {
        BCL_MessageStd( "list_C: " + util::Format()( *ptr));
      }
      for( storage::List< double>::const_iterator ptr( list_D.Begin()), ptr_end( list_D.End()); ptr != ptr_end; ++ptr)
      {
        BCL_MessageStd( "list_D: " + util::Format()( *ptr));
      }

      // test InsertElements with a position and two iterators
      BCL_MessageStd( "test InsertElements with a position and two iterators");
      storage::List< double> pos_two_itr;
      pos_two_itr.InsertElements( pos_two_itr.End(), list_D.Begin(), list_D.End());
      BCL_Example_Check
      (
        pos_two_itr == list_D, "InsertElements with a position and two iterators did not work"
      );

      // test RemoveElement function with iterator
      BCL_MessageStd( "list_D: RemoveElement first element");
      list_D.RemoveElement( list_D.Begin());
      BCL_Example_Check
      (
        *list_D.Begin() == 0, "first element of list_D should be 0 but is not"
      );
      for( storage::List< double>::const_iterator ptr( list_D.Begin()), ptr_end( list_D.End()); ptr != ptr_end; ++ptr)
      {
        BCL_MessageStd( "list_D: " + util::Format()( *ptr));
      }

      // test Append function, which adds elements to the end of the List
      BCL_MessageStd( "list_D: Append element 0.707106");
      list_D.Append( 0.707106);
      BCL_Example_Check
      (
        list_D.LastElement() == 0.707106, "last element of list_D should be 0.707106 but is not"
      );
      BCL_MessageStd( "list_D: " + util::Format()( list_D));

      // test append function with adding another List to an existing list
      BCL_MessageStd( "list_D: Append list_C");
      list_D.Append( list_C);
      BCL_MessageStd( "list_C: " + util::Format()( list_C));
      BCL_MessageStd( "list_D: " + util::Format()( list_D));

      // test PushBack function which adds element to end of List
      BCL_MessageStd( "list_D: PushBack 0.707106");
      list_D.PushBack( 0.707106);
      BCL_Example_Check
      (
        list_D.LastElement() == 0.707106, "last element of list_D should be 0.707106 but is not"
      );
      BCL_MessageStd( "list_D: " + util::Format()( list_D));

      // test PushBack function with adding default element (which is 0)
      BCL_MessageStd( "list_D: PushBack default element");
      list_D.PushBack();
      BCL_Example_Check
      (
        list_D.LastElement() == 0, "last element of list_D should be 0 but is not"
      );
      BCL_MessageStd( "list_D: " + util::Format()( list_D));

      // test PushFront which adds element to beginning of List
      BCL_MessageStd( "list_D: PushFront 0.707106");
      list_D.PushFront( 0.707106);
      BCL_Example_Check
      (
        list_D.FirstElement() == 0.707106, "first element of list_D should be 0.707106 but is not"
      );
      BCL_MessageStd( "list_D: " + util::Format()( list_D));

      // test PushFront with adding default element (which is 0)
      BCL_MessageStd( "list_D: PushFront default element");
      list_D.PushFront();
      BCL_Example_Check
      (
        list_D.FirstElement() == 0, "first element of list_D should be 0 but is not"
      );
      BCL_MessageStd( "list_D: " + util::Format()( list_D));

      // test PopBack which removes the last element of a List
      BCL_MessageStd( "list_D: PopBack");
      list_D.PopBack();
      BCL_Example_Check
      (
        list_D.LastElement() == 0.707106, "last element of list_D should be 0.707106 but is not"
      );
      BCL_MessageStd( "list_D: " + util::Format()( list_D));

      // test PopFront which removes the first element of a List
      BCL_MessageStd( "list_D: PopFront");
      list_D.PopFront();
      BCL_Example_Check
      (
        list_D.FirstElement() == 0.707106, "first element of list_D should be 0.707106 but is not"
      );
      BCL_MessageStd( "list_D: " + util::Format()( list_D));

      // test Resize which gives the List a new size: in this case removes all but first two elements
      BCL_MessageStd( "list_D: Resize to length 2 with default element");
      list_D.Resize( 2);
      BCL_Example_Check
      (
        list_D.GetSize() == 2 && list_D.FirstElement() == 0.707106 && list_D.LastElement() == 0,
        "list_D size should be 2 but is not"
      );
      BCL_MessageStd( "list_D: " + util::Format()( list_D));

      // test Resize which in this case resizes the List to four and uses 0.707106 to fill new spaces
      BCL_MessageStd( "list_D: Resize to length 4 with element 0.707106");
      list_D.Resize( 4, 0.707106);
      BCL_Example_Check
      (
        list_D.GetSize() == 4 && list_D.LastElement() == 0.707106, "list_D size should be 4 but is not"
      );
      BCL_MessageStd( "list_D: " + util::Format()( list_D));

      // test Resize which in this case resizes the List to six and uses the default (0) to fill the new spaces
      BCL_MessageStd( "list_D: Resize to length 6 with default element");
      list_D.Resize( 6);
      BCL_Example_Check
      (
        list_D.GetSize() == 6 && list_D.LastElement() == 0, "list_D size should be 6 but is not"
      );
      BCL_MessageStd( "list_D: " + util::Format()( list_D));

      // test Swap which exchanges two elements
      BCL_MessageStd( "list_D: Swap first and second element");
      std::swap( *list_D.Begin(), *( ++list_D.Begin()));
            BCL_Example_Check
      (
        *list_D.Begin() == 0 && *++list_D.Begin() == 0.707106,
        "1st and 2nd elements of list_D were not swapped correctly"
      );
      BCL_MessageStd( "list_D: " + util::Format()( list_D));

      // test Splice which puts a list into another list at a position denoted by an iterator
      BCL_MessageStd( "list_D: Splice list_C list");
      for( storage::List< double>::const_iterator ptr( list_C.Begin()), ptr_end( list_C.End()); ptr != ptr_end; ++ptr)
      {
        BCL_MessageStd( "list_C: " + util::Format()( *ptr));
      }
      list_D.Splice( ++list_D.Begin(), list_C);
      for( storage::List< double>::const_iterator ptr( list_D.Begin()), ptr_end( list_D.End()); ptr != ptr_end; ++ptr)
      {
        BCL_MessageStd( "list_D: " + util::Format()( *ptr));
      }
      for
      (
        storage::List< double>::const_iterator ptr_C( list_C.Begin()), ptr_C_end( list_C.End()), ptr_D = ++list_D.Begin();
        ptr_C != ptr_C_end;
        ++ptr_C, ++ptr_D
      )
      {
        BCL_Example_Check
        (
          *ptr_C == *ptr_D, "list_C should be spliced into list_D but is not"
        );
      }

      // test Sort using default sorting algorithm (smallest to largest)
      BCL_MessageStd( "list_D: Sort");
      list_D.Sort( std::less< double>());
      double previous( 0);
      for( storage::List< double>::const_iterator ptr( list_D.Begin()), ptr_end( list_D.End()); ptr != ptr_end; ++ptr)
      {
        BCL_MessageStd( "list_D: " + util::Format()( *ptr));
        BCL_Example_Check
        (
          *ptr >= previous, "list_D should be sorted smallest to largest but is not"
        );
        previous = *ptr;
      }

      // test Sort using defined predicate (largest to smallest)
      BCL_MessageStd( "list_D: Sort with predicate");
      list_D.Sort( more_than< double>());
      previous = *list_D.Begin();
      for( storage::List< double>::const_iterator ptr( list_D.Begin()), ptr_end( list_D.End()); ptr != ptr_end; ++ptr)
      {
        BCL_MessageStd( "list_D: " + util::Format()( *ptr));
        BCL_Example_Check
        (
          *ptr <= previous, "list_D should be sorted largest to smallest but is not"
        );
        previous = *ptr;
      }

      // test Reverse which reverses elements in container
      BCL_MessageStd( "list_D: Reverse");
      list_D.Reverse();
      previous = *list_D.Begin();
      for( storage::List< double>::const_iterator ptr( list_D.Begin()), ptr_end( list_D.End()); ptr != ptr_end; ++ptr)
      {
        BCL_MessageStd( "list_D: " + util::Format()( *ptr));
        BCL_Example_Check
        (
          *ptr >= previous, "after reverse function list_D should be sorted smallest to largest but is not"
        );
        previous = *ptr;
      }

      // test Unique which removes all repetitive elements
      BCL_MessageStd( "list_D: Unique");
      list_D.Unique( std::equal_to< double>());
      for
      (
        storage::List< double>::const_iterator ptr( list_D.Begin()), ptr_end( list_D.End()), ptr_next( ++list_D.Begin());
        ptr_next != ptr_end; ++ptr, ++ptr_next
      )
      {
        BCL_MessageStd( "list_D: " + util::Format()( *ptr));
        BCL_Example_Check
        (
          *ptr != *( ptr_next), "all elements in list_D should be unique but are not"
        );
      }

      // test Unique with defined predicate which removes all elements not equal to default value of double (0)
      BCL_MessageStd( "list_D: Unique with predicate");
      list_D.Unique( all_equal< double>());
      for( storage::List< double>::const_iterator ptr( list_D.Begin()), ptr_end( list_D.End()); ptr != ptr_end; ++ptr)
      {
        BCL_MessageStd( "list_D: " + util::Format()( *ptr));
        BCL_Example_Check
        (
          *ptr == 0, "list_D should contain one unique element of value 0"
        );
      }

      // test Merge: prepare list_D
      BCL_MessageStd( "test Merge: prepare list_D: pushback + sort");
      list_D.PushBack( 0.707106);
      list_D.PushBack( 3.1415);
      list_D.Sort( std::less< double>());
      BCL_MessageStd( "list_D: " + util::Format()( list_D));
      // test Merge: prepare list_G
      BCL_MessageStd( "test Merge: prepare list_G: construct + pushback + sort");
      storage::List< double> list_G;
      list_G.PushBack( 2.71828);
      list_G.PushBack( 1.41421);
      list_G.Sort( std::less< double>());
      for( storage::List< double>::const_iterator ptr( list_G.Begin()), ptr_end( list_G.End()); ptr != ptr_end; ++ptr)
      {
        BCL_MessageStd( "list_G: " + util::Format()( *ptr));
      }
      // test Merge: do the merge
      BCL_MessageStd( "test Merge: list_D: Merge list_G");
      list_D.Merge( list_G, std::less< double>());
      BCL_Example_Check
      (
        list_D.GetSize() == 5 && list_G.IsEmpty(),
        "list_G should be merged into list_D but has not properly been done"
      );
      for( storage::List< double>::const_iterator ptr( list_D.Begin()), ptr_end( list_D.End()); ptr != ptr_end; ++ptr)
      {
        BCL_MessageStd( "list_D: " + util::Format()( *ptr));
      }

      // test Merge with predicate: prepare list_G
      BCL_MessageStd
      (
        "test Merge with predicate: prepare list_G: pushback + sort with predicate"
      );
      list_G.PushBack( 1.41421);
      list_G.PushBack( 2.71828);
      list_G.Sort( more_than< double>());
      for( storage::List< double>::const_iterator ptr( list_G.Begin()), ptr_end( list_G.End()); ptr != ptr_end; ++ptr)
      {
        BCL_MessageStd( "list_G: " + util::Format()( *ptr));
      }
      // test Merge with predicate: prepare list_D
      BCL_MessageStd
      (
        "test Merge with predicate: prepare list_D: empty + pushback + sort with predicate"
      );
      list_D.Reset();
      list_D.PushBack();
      list_D.PushBack( 0.707106);
      list_D.PushBack( 3.1415);
      list_D.Sort( more_than< double>());
      BCL_MessageStd( "list_D: " + util::Format()( list_D));

      // test Merge with predicate: do the merge
      BCL_MessageStd( "test Merge with predicate: list_D: Merge list_G");
      list_D.Merge( list_G, more_than< double>());
      BCL_Example_Check
      (
        list_D.GetSize() == 5 && list_G.IsEmpty(),
        "list_G should be merged with predicate into list_D but has not properly been done"
      );
      for( storage::List< double>::const_iterator ptr( list_D.Begin()), ptr_end( list_D.End()); ptr != ptr_end; ++ptr)
      {
        BCL_MessageStd( "list_D: " + util::Format()( *ptr));
      }
      for( storage::List< double>::const_iterator ptr( list_G.Begin()), ptr_end( list_G.End()); ptr != ptr_end; ++ptr)
      {
        BCL_MessageStd( "list_G: " + util::Format()( *ptr));
      }

      // test SetAllElements with given value
      BCL_MessageStd( "list_D: set all elements to 2.71828");
      list_D.SetAllElements( 2.71828);
      for( storage::List< double>::const_iterator ptr( list_D.Begin()), ptr_end( list_D.End()); ptr != ptr_end; ++ptr)
      {
        BCL_MessageStd( "list_D: " + util::Format()( *ptr));
        BCL_Example_Check
        (
          *ptr == 2.71828, "all elements of list_D should be 2.71828 but are not"
        );
      }

      // test SetAllElements with default value (0)
      BCL_MessageStd( "list_D: SetAllElements to 0");
      list_D.SetAllElements();
      for( storage::List< double>::const_iterator ptr( list_D.Begin()), ptr_end( list_D.End()); ptr != ptr_end; ++ptr)
      {
        BCL_MessageStd( "list_D: " + util::Format()( *ptr));
        BCL_Example_Check
        (
          *ptr == 0, "all elements of list_D should be 0 but are not"
        );
      }

      // test operator == when false
      BCL_MessageStd
      (
        "test operator==: list_B == list_D (should be no): " + util::Format()( list_B == list_D)
      );
      BCL_Example_Check
      (
        !( list_B == list_D), "operator == returns true when it should not"
      );

      // test operator == when true
       BCL_MessageStd
      (
        "test operator==: list_B == list_D (should be yes): first make the lists the same"
      );
      list_B.Resize( 4);
      list_B.SetAllElements( 0.707106);
      for( storage::List< double>::const_iterator ptr( list_B.Begin()), ptr_end( list_B.End()); ptr != ptr_end; ++ptr)
      {
        BCL_MessageStd( "list_B: " + util::Format()( *ptr));
      }
      list_D.Resize( 4);
      list_D.SetAllElements( 0.707106);
      for( storage::List< double>::const_iterator ptr( list_D.Begin()), ptr_end( list_D.End()); ptr != ptr_end; ++ptr)
      {
        BCL_MessageStd( "list_D: " + util::Format()( *ptr));
      }
      BCL_MessageStd
      (
        "test operator==: list_B == list_D (should be yes): " + util::Format()( list_B == list_D)
      );
      BCL_Example_Check
      (
        list_B == list_D, "operator == returns false when it should not"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // Write
      // test write function
      WriteBCLObject( list_D);
      // test read function
      storage::List< double> list_read;
      ReadBCLObject( list_read);
      BCL_ExampleIndirectCheck( list_D, list_read, "I/O");

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleStorageList

  const ExampleClass::EnumType ExampleStorageList::s_Instance
  (
    GetExamples().AddEnum( ExampleStorageList())
  );

} // namespace bcl
