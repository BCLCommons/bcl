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
#include "align/bcl_align_assignment.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_back_bone.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_align_assignment.cpp
  //!
  //! @author heinzes1
  //! @date 2010/11/01
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAlignAssignment :
    public ExampleInterface
  {
  public:

    ExampleAlignAssignment *Clone() const
    { 
      return new ExampleAlignAssignment( *this);
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
      storage::List< size_t> assignment_element;
      util::SiPtr< const storage::List< size_t> > sp_assignment_element( assignment_element);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      align::Assignment< storage::List< size_t> > default_assignment;

      // test constructor taking one member
      align::Assignment< storage::List< size_t> > one_member_assignment( sp_assignment_element);

      // test constructor taking two members
      align::Assignment< storage::List< size_t> > two_member_assignment( sp_assignment_element, sp_assignment_element);

      // test clone and copy constructor together
      util::ShPtr< align::Assignment< storage::List< size_t> > > assignment_copy( two_member_assignment.Clone());

      // create some more for IsIdentical/IsEqual tests
      biol::AABackBone alanine( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().ALA)));
      util::SiPtr< const biol::AABase> sp_alanine( alanine);
      align::Assignment< biol::AABase> test_assignment( sp_alanine);

      biol::AABackBone alanine_copy( alanine);
      util::SiPtr< const biol::AABase> sp_alanine_copy( alanine_copy);
      align::Assignment< biol::AABase> test_assignment_equal( sp_alanine_copy);

      biol::AABackBone tyrosine( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().TYR)));
      util::SiPtr< const biol::AABase> sp_tyrosine( tyrosine);
      align::Assignment< biol::AABase> test_assignment_different( sp_tyrosine);

    /////////////////
    // data access //
    /////////////////

      // test GetMembers and their size on default_constructed_assignment
      BCL_MessageStd( "1: Assignment()");
      BCL_ExampleCheck( default_assignment.GetMembers().GetSize(), 0);

      // test one_member_constructed_assignment
      BCL_MessageStd( "2: Assignment( sp_assignment)");
      BCL_ExampleCheck( one_member_assignment.GetMembers().GetSize(), 1);

      // test one_member_constructed_assignment
      BCL_MessageStd( "3: Assignment( sp_assignment, sp_assignment)");
      BCL_ExampleCheck( two_member_assignment.GetMembers().GetSize(), 2);

      // test one_member_constructed_assignment
      BCL_MessageStd( "4: Clone()");
      BCL_ExampleCheck( assignment_copy->GetMembers().GetSize(), 2);

    ////////////////
    // operations //
    ////////////////

      // test Append
      default_assignment.Append( sp_assignment_element);
      BCL_MessageStd( "5: default_assignment.Append()");
      BCL_ExampleCheck( default_assignment.GetMembers().GetSize(), 1);

      // test IsIdentical
      BCL_MessageStd( "6: test_assignment.IsIdentical( test_assignment)");
      BCL_ExampleCheck( test_assignment.IsIdentical( test_assignment), true);
      BCL_MessageStd( "7: test_assignment.IsIdentical( test_assignment_equal)");
      BCL_ExampleCheck( test_assignment.IsIdentical( test_assignment_equal), false);
      BCL_MessageStd( "8: test_assignment.IsIdentical( test_assignment_different)");
      BCL_ExampleCheck( test_assignment.IsIdentical( test_assignment_different), false);

      // test IsEqual
      BCL_MessageStd( "9: test_assignment.IsEqual( test_assignment)");
      BCL_ExampleCheck( test_assignment.IsEqual( test_assignment), true);
      BCL_MessageStd( "10: test_assignment.IsEqual( test_assignment_equal)");
      BCL_ExampleCheck( test_assignment.IsEqual( test_assignment_equal), true);
      BCL_MessageStd( "11: test_assignment.IsEqual( test_assignment_different)");
      BCL_ExampleCheck( test_assignment.IsEqual( test_assignment_different), false);

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

      // test ToString
      BCL_MessageStd
      (
        "12: test_assignment.ToString() = " + util::Format()( test_assignment.ToString())
      );
      std::string expected_result( "A");
      BCL_ExampleCheck( test_assignment.ToString(), expected_result);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRestraintAssignment

  const ExampleClass::EnumType ExampleAlignAssignment::s_Instance( GetExamples().AddEnum( ExampleAlignAssignment()));

} // namespace bcl
