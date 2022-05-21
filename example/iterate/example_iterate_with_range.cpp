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
#include "iterate/bcl_iterate_with_range.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //! @example example_iterator_generic_iterate_with_range.cpp
  //! @details this examples demonstrates how to use the generic iterator implementation
  //!
  //! @author mendenjl
  //! @date Jul 14, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleIterateWithRange :
    public ExampleInterface
  {
  public:

    ExampleIterateWithRange *Clone() const
    {
      return new ExampleIterateWithRange( *this);
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

      // construct iterator data for vectors
      storage::Vector< double> dozen_doubles( 12);
      iterate::WithRange< storage::Vector< double>::iterator, double>
        mutable_vect_itr( dozen_doubles.Begin(), dozen_doubles.End());
      iterate::WithRange< storage::Vector< double>::const_iterator, const double>
        const_vect_itr( dozen_doubles.Begin(), dozen_doubles.End());

      // construct iterator data for lists
      storage::List< double> ten_doubles( 10, 5.0);
      iterate::WithRange< storage::List< double>::iterator, double>
        mutable_list_itr( ten_doubles.Begin(), ten_doubles.End());
      iterate::WithRange< storage::List< double>::const_iterator, const double>
        const_list_itr( ten_doubles.Begin(), ten_doubles.End());

      // construct iterator for sets
      storage::Set< size_t> set_of_nums( storage::Set< size_t>::Create( 100, 1, 75, 92));
      iterate::WithRange< storage::Set< size_t>::iterator, const size_t>
        mutable_set_itr( set_of_nums.Begin(), set_of_nums.End());
      iterate::WithRange< storage::Set< size_t>::const_iterator, const size_t>
        const_set_itr( set_of_nums.Begin(), set_of_nums.End());

      // construct iterator data for a set of only one element
      storage::Set< size_t> set_of_one( 100);
      iterate::WithRange< storage::Set< size_t>::iterator, const size_t>
        mutable_set_of_one_itr( set_of_one.Begin(), set_of_one.End());
      iterate::WithRange< storage::Set< size_t>::const_iterator, const size_t>
        const_set_of_one_itr( set_of_one.Begin(), set_of_one.End());

    /////////////////
    // data access //
    /////////////////

      // test constructor
      BCL_ExampleIndirectCheck
      (
        mutable_list_itr == ten_doubles.Begin(), true,
        "constructor from range should start at the beginning"
      );

      // try coping the mutable list itr and seeing if it too is at the beginning
      iterate::WithRange< storage::List< double>::iterator, double> mutable_list_itr_copy( mutable_list_itr);
      BCL_ExampleIndirectCheck
      (
        mutable_list_itr_copy.NotAtEnd()
        && mutable_list_itr_copy == ten_doubles.Begin()
        && mutable_list_itr == mutable_list_itr_copy,
        true,
        "copy constructor"
      );

      // move the copy forward by 1, make sure it doesn't change mutable_list_itr
      BCL_ExampleCheck( ( ++mutable_list_itr_copy).GetPosition(), 1);
      BCL_ExampleIndirectCheck
      (
        mutable_list_itr_copy != mutable_list_itr && mutable_list_itr == ten_doubles.Begin(),
        true,
        "Copy constructed iterator operator ++ should not change iterator that was copied"
      );

      // try moving the iterator forward 10 times to see if it reaches the end
      for( size_t i( 0), size( ten_doubles.GetSize()); i < size; i++)
      {
        ++mutable_list_itr;
      }

      // check that we are indeed at the end
      BCL_ExampleIndirectCheck
      (
        mutable_list_itr.NotAtEnd(),
        false,
        "mutable_list_itr should be at the end after 10 increments"
      );

      // check get displacement
      BCL_ExampleCheck( mutable_list_itr.GetPosition(), 10);

      // try Restarting the iterator to the beginning
      mutable_list_itr.GotoBegin();
      BCL_ExampleIndirectCheck
      (
        mutable_list_itr == ten_doubles.Begin(),
        true,
        "mutable_list_itr.GotoBegin()"
      );

      mutable_list_itr.GotoEnd(); // Restart the iterator to the last element
      BCL_ExampleIndirectCheck
      (
        !mutable_list_itr.NotAtEnd() && mutable_list_itr == ten_doubles.End(),
        true,
        "mutable_list_itr.GotoEnd()"
      );

      mutable_list_itr.Restart();
      ++mutable_list_itr;

      storage::List< double>::iterator itr_list_one( ++ten_doubles.Begin());
      // check that the iterators are the same and return the same value when dereferenced
      BCL_ExampleIndirectCheck
      (
        mutable_list_itr == itr_list_one && *mutable_list_itr == *itr_list_one,
        true,
        "increment on list"
      );

      // check setting a value
      *mutable_list_itr = 5.0;
      BCL_ExampleIndirectCheck
      (
        *itr_list_one, 5.0, "setting a value through a generic iterator implemention"
      );

      BCL_ExampleIndirectCheck
      (
        const_set_of_one_itr == set_of_one.Begin(),
        true,
        "construction of generic iterator on a set with one element"
      );
      const_set_of_one_itr.Restart();
      BCL_ExampleIndirectCheck( const_set_of_one_itr == set_of_one.Begin(), true, "Restart");

      // check that moving a set iterator forward by 1 works
      BCL_ExampleCheck( ( ++const_set_of_one_itr).GetPosition(), 1);
      BCL_ExampleIndirectCheck
      (
        const_set_of_one_itr.NotAtEnd(),
        false,
        "Moving to position 1 on a set of size one reaches the end"
      );

      BCL_ExampleIndirectCheck
      (
        const_set_itr == set_of_nums.Begin(),
        true,
        "generic iterator implementation constructor on set"
      );

      // check that restarting the itr puts it back to the beginning
      ++const_set_itr;
      const_set_itr.Restart();
      BCL_ExampleIndirectCheck( const_set_itr == set_of_nums.Begin(), true, "( ++const_set_itr).Restart()");

      // check that moving a set iterator forward by 4 puts it at the end
      BCL_ExampleCheck( const_set_itr.GotoPosition( set_of_nums.GetSize()), false);
      BCL_ExampleIndirectCheck
      (
        const_set_itr == set_of_nums.End(),
        true,
        "const_set_itr.GotoPosition( set_of_nums.GetSize())"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleIterateWithRange

  const ExampleClass::EnumType ExampleIterateWithRange::s_Instance
  (
    GetExamples().AddEnum( ExampleIterateWithRange())
  );

} // namespace bcl
