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
#include "storage/bcl_storage_vector.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_storage_vector.cpp
  //!
  //! @author heinzes1
  //! @date Nov 6, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleStorageVector :
    public ExampleInterface
  {
  public:

    ExampleStorageVector *Clone() const
    {
      return new ExampleStorageVector( *this);
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
      // initialize random storage::Vector
      const double c1[] =
      {
        0.0,
        -1.0,
        2.0,
        -3.0,
        4.0,
        -5.0,
        6.0,
        -7.0,
        8.0,
        -9.0,
        10.0,
        -11.0,
        12.0,
        -13.0,
        14.0
      };

      // construction of an empty storage vector
      storage::Vector< double> sv1;
      BCL_Example_Check
      (
        sv1.IsEmpty(), "default constructor of empty vector does not give empty vector"
      );

      // construction of an empty storage vector with 15 elements
      storage::Vector< double> sv2( 15);
      BCL_Example_Check
      (
        sv2.GetSize() == 15, "construction of empty storage vector with given size did not work"
      );

      // construction of a vector with 15 elements which are given with pointer to data of c1
      storage::Vector< double> sv3( 15, c1);
      BCL_Example_Check
      (
        sv3.GetSize() == 15 && *sv3.Begin() == 0.0 && *sv3.Last() == 14.0,
        "construction of vector with given size and pointer to data did not work correctly"
      );

      // test GetSize function
      BCL_MessageStd( "size of sv3: " + util::Format()( sv3.GetSize()));

      // ask vector whether it is empty ( is faster than asking size == 0) in false situation
      BCL_MessageStd( "is sv3 empty?: " + util::Format()( sv3.IsEmpty()));
      BCL_Example_Check
      (
        !( sv3.IsEmpty()), "IsEmpty function gives true when it should give false"
      );

      // ask vector whether it is empty ( is faster than asking size == 0) in true situation
      BCL_MessageStd( "is sv1 empty?: " + util::Format()( sv1.IsEmpty()));
      BCL_Example_Check
      (
        sv1.IsEmpty(), "IsEmpty function gives false when it should give true"
      );

      // test output of entire vector
      BCL_MessageStd( "test output of all elements of sv3");
      BCL_MessageStd( util::Format()( sv3));

      // test output and selection of specified element (numbering begins with 0)
      BCL_MessageStd( "element 2 of sv3: " + util::Format()( sv3( 1)));
      BCL_Example_Check
      (
        sv3( 1) == -1.0, "element accession did not work properly"
      );

      // test assignment of selected element
      BCL_MessageStd( "element 2 of sv3 gets new value: ");
      sv3( 1) = 4;
      BCL_MessageStd( util::Format()( sv3( 1)));
      BCL_Example_Check
      (
        sv3( 1) == 4.0, "element value change did not work properly"
      );

      // output of first and last element using indexes
      BCL_MessageStd( "1st element of sv3: " + util::Format()( sv3( 0)));
      BCL_MessageStd( "last element of sv3: " + util::Format()( sv3( sv3.GetSize() - 1)));
      BCL_Example_Check
      (
        sv3.GetSize() - 1 == 14.0 && sv3( 0) == 0.0, "accessing elements by index did not work properly"
      );

      // test SetAllElements to 1.1
      sv2.SetAllElements( 1.1);
      BCL_MessageStd( "all elements of sv2 get value 1.1 ");
      for( size_t index( 0); index < sv2.GetSize(); ++index)
      {
        BCL_Example_Check
        (
          sv2( index) == 1.1, "SetAllElements function did not work properly"
        );
      }

      // test Reset function to delete all elements of sv2
      BCL_MessageStd( "Deletion of all elements of sv2");
      sv2.Reset();
      BCL_Example_Check
      (
        sv2.IsEmpty(), "Reset function did not delete all elements of the vector"
      );

      // test RemoveElements function
      BCL_MessageStd( "test RemoveElements function: ");
      BCL_MessageStd( "all elements of sv3: " + util::Format()( sv3));
      BCL_MessageStd
      (
        "erase five elements of sv3 starting with 5th element (count from 0)"
      );
      sv3.RemoveElements( 5, 5);
      BCL_MessageStd( "output of all elements of sv3 after erasing");
      BCL_MessageStd( util::Format()( sv3));
      BCL_Example_Check
      (
        sv3.GetSize() == 10 && sv3( 4) == 4 && sv3( 5) == 10, "Remove elements function malfunctioned"
      );

      // test copy constructor: sv5 of sv3
      BCL_MessageStd
      (
        "test copy constructor: new storagevector v5 constructed from sv3"
      );
      storage::Vector< double> sv5( sv3);
      BCL_MessageStd( util::Format()( sv5));
      BCL_Example_Check
      (
        sv5 == sv3, "copy constructor did not work properly"
      );

      // test PushBack which adds element to the last position
      BCL_MessageStd( "test PushBack: add element to the last position of sv5");
      sv5.PushBack( 2.0);
      BCL_MessageStd( util::Format()( sv5));
      BCL_Example_Check
      (
        *sv5.Last() == 2.0, "PushBack did not properly add element"
      );

      // test PopBack which removes element from the last position
      BCL_MessageStd( "test PopBack: remove element from the last position of sv5");
      sv5.PopBack();
      BCL_MessageStd( util::Format()( sv5));
      BCL_Example_Check
      (
        *sv5.Last() == 14.0, "PopBack did not properly remove element"
      );

      // test Resize by adding 2 elements to the end of value 2.0
      BCL_MessageStd( "test Resize: size of sv5: " + util::Format()( sv5.GetSize()));
      BCL_MessageStd( "resize sv5 -> add 2 elements");
      sv5.Resize( 12, 2.0);
      BCL_MessageStd( "size of sv5: " + util::Format()( sv5.GetSize()));
      BCL_MessageStd( util::Format()( sv5));
      BCL_Example_Check
      (
        sv5.GetSize() == 12 && sv5( 10) == 2.0 && sv5( 11) == 2.0, "Resize did not resize the vector correctly"
      );

      // test construct storage vector sv6 from five elements of sv5 starting at Pos 4
      BCL_MessageStd
      (
        "test: construct sv6 from five elements of sv5 starting at POS 4: "
      );
      storage::Vector< double> sv6( sv5, 3, 5);
      BCL_MessageStd( "size of sv6: " + util::Format()( sv6.GetSize()));
      BCL_MessageStd( util::Format()( sv6));

      // test InsertElements which inserts storage vector elements into another storage vector at POS
      BCL_MessageStd
      (
        "test InsertElements: insert elements of sv5 to sv4 before POS 3 "
      );
      storage::Vector< double> sv4( 15, c1);
      BCL_MessageStd( "new vector: sv4: " + util::Format()( sv4));
      BCL_MessageStd( "insert sv5 into sv4 before position 3");
      sv4.InsertElements( 2, sv5);
      BCL_MessageStd( "sv4: " + util::Format()( sv4));
      BCL_Example_Check
      (
        sv4.GetSize() == 27 && sv4( 2) == *sv5.Begin() && sv4( 13) == *sv5.Last(),
        "InsertElements did not work"
      );

      // insert 5 elements in other storage vector at POS
      BCL_MessageStd
      (
        "test InsertElements: insert 5.0 three times into sv5 before POS 1 "
      );
      BCL_MessageStd( "size of sv4: " + util::Format()( sv4.GetSize()));
      sv4.InsertElements( 1, 5.0, 3);
      BCL_MessageStd( "size of sv4: " + util::Format()( sv4.GetSize()));
      BCL_MessageStd( util::Format()( sv4));
      BCL_Example_Check
      (
        sv4.GetSize() == 30 && sv4( 1) == 5 && sv4( 2) == 5 && sv4( 3) == 5, "InsertElements did not work"
      );

      // test writing vector to file
      BCL_MessageStd( "write storage::Vector sv6 to file");
      WriteBCLObject( sv6);
      // test reading vector from file
      BCL_MessageStd( "read storage::Vector sv7 from file");
      storage::Vector< double> sv7;
      ReadBCLObject( sv7);
      BCL_MessageStd( util::Format()( sv7));
      BCL_Example_Check
      (
        sv7 == sv6, "reading and writing vectors to and from files did not work"
      );

      // test Swap two elements POS1 and POS2
      BCL_MessageStd( "test Swap function");
      std::swap( sv7( 1), sv7( 2));
      BCL_MessageStd( util::Format()( sv7));
      BCL_Example_Check
      (
        sv7( 1) == 10 && sv7( 2) == 4, "Swap did not correctly switch elements"
      );

      // test SetAllElements
      BCL_MessageStd( "test SetAllElements function: set sv7 to hold all zeros");
      sv7.SetAllElements( double( 0.0));
      BCL_MessageStd( util::Format()( sv7));
      for( size_t index( 0); index < sv7.GetSize(); ++index)
      {
        BCL_Example_Check
        (
          sv7( index) == 0, "test of SetAllElements did not work"
        );
      }

      // test construction from size and pointer to data
      BCL_MessageStd( "test construction from size and pointer to data");
      std::string s[ 2] =
      {
        "C",
        "CB"
      };

      storage::Vector< std::string> ssv1( 2, s);
      BCL_MessageStd( util::Format()( ssv1));
      BCL_Example_Check
      (
        ssv1( 0) == "C" && ssv1( 1) == "CB", "construction from size and pointer to data did not work"
      );

      // test Reset function
      BCL_MessageStd( "test Reset function");
      ssv1.Reset();
      BCL_MessageStd
      (
        "clear storagevector - is empty? : " + util::Format()( ssv1.IsEmpty())
      );
      BCL_MessageStd( "clear storagevector: Size : " + util::Format()( ssv1.GetSize()));
      BCL_Example_Check
      (
        ssv1.IsEmpty(), "Reset function did not create empty vector"
      );
      BCL_Example_Check
      (
        ssv1.IsEmpty(), "Reset function did not create vector of size 0"
      );

      storage::Vector< size_t> sv9( 64);
      for( size_t i = 0; i < sv9.GetSize(); i++)
      {
        sv9( i) = i;
      }

      storage::Vector< size_t> sv10( sv9);
      BCL_MessageStd
      (
        "test shuffle function: set sv9 to hold all integers between 0 and "
        + util::Format()( sv9.GetSize())
      );
      sv9.Shuffle();
      BCL_Example_Check
      (
        !( sv9 == sv10),
        "Shuffle function did not change anything (this may occur by chance with likelihood 2^-"
        + util::Format()( sv9.GetSize()) + ")"
      );
      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;
  }; //end Examplestorage::Vector

  const ExampleClass::EnumType ExampleStorageVector::s_Instance
  (
    GetExamples().AddEnum( ExampleStorageVector())
  );
} // namespace bcl

