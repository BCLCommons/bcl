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
#include "iterate/bcl_iterate_generic.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_plus_equals.h"
#include "util/bcl_util_si_ptr_vector.h"
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //! @example example_iterate_Generic.cpp
  //! @details this examples demonstrates how to use the Generic iterator and times them compared to
  //!          using siptrvectors to the interfaces
  //!
  //! @author mendenjl
  //! @date Jul 14, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleIterateGeneric :
    public ExampleInterface
  {
  public:

    ExampleIterateGeneric *Clone() const
    {
      return new ExampleIterateGeneric( *this);
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

    // enum for the functions that will be tested on various Generic iterators
    enum Test
    {
      e_None,
      e_Validity,
      e_Displacement,
      e_Increment,
      e_GetPositionStart,
      e_Decrement,
      e_GetSize,
      e_IncrementToEnd,
      e_GetPositionEnd,
      e_DecrementToStart,
      e_GetPositionRestart,
      e_DereferenceStart,
      e_DereferenceSecond,
      e_Pointer,
      e_GotoPositionSecond,
      e_GotoPositionStart,
      e_GotoRandomPosition,
      e_CountElementsThatSatisfy,
      s_NumberTests
    };

    const std::string &GetTestName( const Test &TEST) const
    {
      static const std::string s_tests[ int( s_NumberTests) + 1] =
      {
        "None",
        "Validity",
        "Displacement",
        "Increment",
        "GetPositionStart",
        "Decrement",
        "GetSize",
        "IncrementToEnd",
        "GetPositionEnd",
        "DecrementToStart",
        "GetPositionRestart",
        "DereferenceStart",
        "DereferenceSecond",
        "Pointer",
        "GotoPositionSecond",
        "GotoPositionStart",
        "GotoRandomPosition",
        "CountElementsThatSatisfy",
        "All"
      };
      return s_tests[ int( TEST)];
    }

    int Run() const
    {
      // this example runs over 100 checks (20 for each container)
      // Because example checks cannot have the same line number, however, we just check to see that the function
      // that tests out a Generic iterator on a given container returns that it passed all the tests

      // test Generic iterators on an empty vector
      BCL_ExampleCheck( GetTestName( GetHighestTestPassed( storage::Vector< char>())), GetTestName( s_NumberTests));

      // test Generic iterators on a vector containing only 'a'
      BCL_ExampleCheck
      (
        GetTestName( GetHighestTestPassed( storage::Vector< char>( 1, 'a'))),
        GetTestName( s_NumberTests)
      );

      // test Generic iterators on the set of vowels
      BCL_ExampleCheck
      (
        GetTestName( GetHighestTestPassed( storage::Set< char>::Create( 'a', 'e', 'i', 'o', 'u'))),
        GetTestName( s_NumberTests)
      );

      // test Generic iterators on a vector of chars
      BCL_ExampleCheck
      (
        GetTestName( GetHighestTestPassed( storage::Vector< char>( 9, "aeiouaaa"))),
        GetTestName( s_NumberTests)
      );

      // test Generic iterators on list of chars
      BCL_ExampleCheck
      (
        GetTestName
        (
          GetHighestTestPassed( storage::List< char>( std::string( "abracadabra").size(), "abracadabra"))
        ),
        GetTestName( s_NumberTests)
      );

      // compare how long it takes to access First in 1000 elements 1000 times using the
      // normal iterator vs the Generic iterator
      storage::List< storage::Pair< double, size_t> > pair_list( 1000);

      const size_t num_loops( 100000);

      storage::Pair< double, size_t> default_pair( 0.0, 0);

      util::Stopwatch native_itr_time_nocall
      (
        "100 million iterations with native iterator on list",
        util::Message::e_Standard,
        true
      );

      storage::List< storage::Pair< double, size_t> >::iterator native_itr( pair_list.Begin()), itr_end( pair_list.End());
      for( size_t i( 0); i < num_loops; ++i)
      {
        for( native_itr = pair_list.Begin(); native_itr != itr_end; ++native_itr)
        {
          *native_itr = default_pair;
        }
      }
      native_itr_time_nocall.Stop();

      util::Stopwatch generic_itr_time_nocall
      (
        "100 million iterations with Generic iterator on list",
        util::Message::e_Standard,
        true
      );

      iterate::Generic< storage::Pair< double, size_t> > generic_itr_pair( pair_list.Begin(), pair_list.End());
      for( size_t i( 0); i < num_loops; ++i)
      {
        for( generic_itr_pair.Restart(); generic_itr_pair.NotAtEnd(); ++generic_itr_pair)
        {
          *generic_itr_pair = default_pair;
        }
      }
      generic_itr_time_nocall.Stop();

      storage::Vector< math::PlusEquals< float> > plus_equals( 1000);
      iterate::Generic< const math::PlusEquals< float> > generic_itr_pe( plus_equals.Begin(), plus_equals.End());
      iterate::Generic< const math::AssignmentOperationInterface< float> >
        generic_itr_assignment( plus_equals.Begin(), plus_equals.End());

      util::Stopwatch generic_itr_time_base
      (
        "100 million iterations with Generic iterator to base on vector",
        util::Message::e_Standard,
        true
      );

      for( size_t i( 0); i < num_loops; ++i)
      {
        for( generic_itr_assignment.Restart(); generic_itr_assignment.NotAtEnd(); ++generic_itr_assignment)
        {
          i += 1;
        }
        i -= 1000;
      }
      generic_itr_time_base.Stop();

      util::Stopwatch generic_itr_time_der
      (
        "100 million iterations with Generic iterator to derived class on vector",
        util::Message::e_Standard,
        true
      );

      for( size_t i( 0); i < num_loops; ++i)
      {
        for( generic_itr_pe.Restart(); generic_itr_pe.NotAtEnd(); ++generic_itr_pe)
        {
          i += 1;
        }
        i -= 1000;
      }
      generic_itr_time_der.Stop();

      util::Stopwatch native_itr_time
      (
        "100 million iterations with native iterator on vector",
        util::Message::e_Standard,
        true
      );
      for( size_t i( 0); i < num_loops; ++i)
      {
        for
        (
          storage::Vector< math::PlusEquals< float> >::const_iterator
            itr( plus_equals.Begin()), itr_end( plus_equals.End());
          itr < itr_end;
          ++itr
        )
        {
          i += 1;
        }
        i -= 1000;
      }
      native_itr_time.Stop();

      util::Stopwatch native_itr_time_si_ptr
      (
        "100 million iterations with constructed simple pointer vector to interface",
        util::Message::e_Standard,
        true
      );
      for( size_t i( 0); i < num_loops; ++i)
      {
        util::SiPtrVector< const math::AssignmentOperationInterface< float> >
          si_ptr_vector( plus_equals.Begin(), plus_equals.End());
        for
        (
          util::SiPtrVector< const math::AssignmentOperationInterface< float> >::const_iterator
            itr( si_ptr_vector.Begin()), itr_end( si_ptr_vector.End());
          itr < itr_end;
          ++itr
        )
        {
          i += 1;
        }
        i -= 1000;
      }
      native_itr_time_si_ptr.Stop();

      // demonstrate that a generic iterator can be made using an existing container of pointers, with no need
      // to manually dereference the values
      util::SiPtrVector< const math::PlusEquals< float> >
        si_ptr_vector( plus_equals.Begin(), plus_equals.End());

      iterate::Generic< const math::AssignmentOperationInterface< float> >
        itr_si_ptr_vector( si_ptr_vector.Begin(), si_ptr_vector.End());

      // ensure that the returned types are equal
      const float to_add( 5.0);
      const float original_value( 0.1);
      float added_to( 0.1);
      itr_si_ptr_vector->operator ()( added_to, to_add);
      BCL_ExampleIndirectCheck
      (
        to_add + original_value,
        added_to,
        "Generic iterator can handle pointer vectors without additional dereferencing"
      );
      return 0;
    }

    //! @brief a function that tests a Generic iterator to const char &
    //! @param CONTAINER the container to make a Generic iterator of chars over
    //! @return the last test passed
    //! Note: this function will be run multiple times with Generic iterators
    //! that operate on different containers.  Since example checks must have unique line numbers, this function
    //! unfortunately cannot have example checks.
    //! If all goes well, this function returns s_NumberTests.
    //! Otherwise, it returns the highest test that was passed
    template< typename t_Container>
    Test GetHighestTestPassed( const t_Container &CONTAINER) const
    {
      // construct a new Generic iterator, which is a Generic version of the iterator returned by Begin() and End()
      // on the container
      iterate::Generic< const char> itr( CONTAINER.Begin(), CONTAINER.End());

      Test highest_test_passed( e_None);

      // test whether the iterator appears valid
      if( CONTAINER.GetSize() > 0)
      {
        if( !itr.NotAtEnd()) // iterator thought the container was invalid, but it has elements
        {
          return highest_test_passed;
        }
      }
      else if( itr.NotAtEnd()) // iterator thought the container was valid, but it was empty
      {
        return highest_test_passed;
      }
      highest_test_passed = e_Validity;

      // test whether the iterator knows its at the start
      if( itr.GetPosition() != 0)
      {
        return highest_test_passed;
      }
      highest_test_passed = e_GetPositionStart;

      ++itr;
      if( itr.GetPosition() != 1)
      {
        return highest_test_passed;
      }
      highest_test_passed = e_Increment;

      --itr;
      if( itr.GetPosition() != 0)
      {
        return highest_test_passed;
      }
      highest_test_passed = e_Decrement;

      if( itr.GetSize() != CONTAINER.GetSize())
      {
        return highest_test_passed;
      }
      highest_test_passed = e_GetSize;

      for( size_t i( 0); i < CONTAINER.GetSize(); ++i, ++itr)
      {
        if( !itr.NotAtEnd())
        {
          return highest_test_passed;
        }
      }
      if( itr.NotAtEnd())
      {
        return highest_test_passed;
      }
      highest_test_passed = e_IncrementToEnd;

      if
      (
        size_t( itr.GetPosition()) != CONTAINER.GetSize()
        || size_t( itr.GetPosition()) != itr.GetSize()
      )
      {
        return highest_test_passed;
      }
      highest_test_passed = e_GetPositionEnd;

      // check whether we can decrement to the start of the container without any problems
      size_t position( itr.GetPosition());
      while( position > 0)
      {
        --position;
        --itr;
      }
      if( !itr.NotAtEnd() && CONTAINER.GetSize() > 0)
      {
        return highest_test_passed;
      }
      highest_test_passed = e_DecrementToStart;

      // check whether current position looks right
      if( itr.GetPosition() != 0)
      {
        return highest_test_passed;
      }
      highest_test_passed = e_GetPositionRestart;

      // the remaining operations require at least one element
      if( CONTAINER.GetSize() > 1)
      {
        // now compare the results of dereferencing the iterators with the normal iterator
        typename t_Container::const_iterator normal_itr( CONTAINER.Begin());
        if( *itr != *normal_itr)
        {
          return highest_test_passed;
        }
        highest_test_passed = e_DereferenceStart;

        // what about on the second element?
        ++itr;
        ++normal_itr;
        if( *itr != *normal_itr)
        {
          return highest_test_passed;
        }
        highest_test_passed = e_DereferenceSecond;

        itr.Restart();
        itr.GotoPosition( 1);

        if( itr != normal_itr)
        {
          return highest_test_passed;
        }
        highest_test_passed = e_GotoPositionSecond;

        itr.GotoPosition( 0);
        if( itr != CONTAINER.Begin())
        {
          BCL_MessageStd
          (
            "position should have been 0 but was: " + util::Format()( itr.GetPosition())
          );
          return highest_test_passed;
        }
        highest_test_passed = e_GotoPositionStart;

        // check whether goto random position returns the correct position
        size_t new_position = itr.GotoRandomPosition();
        normal_itr = CONTAINER.Begin();
        std::advance( normal_itr, new_position);
        if( !itr.NotAtEnd() || itr != normal_itr)
        {
          return highest_test_passed;
        }
        highest_test_passed = e_GotoRandomPosition;

        // check whether we can count the elements that satisfy a particular condition
        if
        (
          itr.CountElementsThatSatisfy( std::bind2nd( std::equal_to< char>(), 'a'))
          != size_t( std::count_if( CONTAINER.Begin(), CONTAINER.End(), std::bind2nd( std::equal_to< char>(), 'a')))
        )
        {
          return highest_test_passed;
        }
        highest_test_passed = e_CountElementsThatSatisfy;
      }
      highest_test_passed = s_NumberTests;
      return highest_test_passed;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleIterateGeneric

  const ExampleClass::EnumType ExampleIterateGeneric::s_Instance
  (
    GetExamples().AddEnum( ExampleIterateGeneric())
  );

} // namespace bcl
