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
#include "restraint/bcl_restraint_group_collection.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_group_collection.cpp
  //!
  //! @author heinzes1
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRestraintGroupCollection :
    public ExampleInterface
  {
  public:

    ExampleRestraintGroupCollection *Clone() const
    { return new ExampleRestraintGroupCollection( *this);}

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
      // test default constructor
      BCL_MessageStd( "test default const");
      restraint::GroupCollection< std::string, biol::Atom> def_const;
      BCL_Example_Check
      (
        def_const.IsEmpty(), "default constructor did not properly construct copy"
      );

      // test construct from an initial group identifier and Group
      BCL_MessageStd( "test construct from an initial group identifier and Group");
      biol::Atom cb( biol::GetAtomTypes().CB);
      restraint::GroupCollection< std::string, biol::Atom> initial_group_constr
      (
        "non-backbone",
        restraint::Group< biol::Atom>( 3, util::SiPtr< const biol::Atom>( cb))
      );
      BCL_Example_Check
      (
        initial_group_constr[ "non-backbone"].GetSize() == 3, "test construct from an initial group identifier and Group gave group wrong size"
      );
      for
      (
        restraint::Group< biol::Atom>::const_iterator
          itr( initial_group_constr[ "non-backbone"].Begin()), itr_end( initial_group_constr[ "non-backbone"].End());
        itr != itr_end;
        ++itr
      )
      {
        BCL_Example_Check
        (
          ( *itr)->GetType() == biol::GetAtomTypes().CB,
          "constructing from a Group and a number of empty members did not add Group properly"
        );
      }

      // test Insert taking a t_GroupIdentifier and a Group
      BCL_MessageStd( "test Insert taking a t_GroupIdentifier and a Group");
      biol::Atom ca( biol::GetAtomTypes().CA);
      BCL_Example_Check
      (
        initial_group_constr.Insert
        (
          "backbone",
          restraint::Group< biol::Atom>( 3, util::SiPtr< const biol::Atom>( ca))
        ).second,
        "Insert taking a t_GroupIdentifier and a Group failed"
      );

      // test Insert function with standard pair
      BCL_MessageStd( "test Insert with std pair");
      BCL_Example_Check
      (
        (
          def_const.Insert
          (
            std::pair< std::string, restraint::Group< biol::Atom> >
            (
              "b", restraint::Group< biol::Atom>( 3, util::SiPtr< const biol::Atom>( cb))
            )
          )
        ).second, "insert failed"
      );

      // test Insert function with bcl Pair
      BCL_MessageStd( "test Insert with bcl Pair");
      BCL_Example_Check
      (
        (
          def_const.Insert
          (
            storage::Pair< std::string, restraint::Group< biol::Atom> >
            (
              "a", restraint::Group< biol::Atom>( 3, util::SiPtr< const biol::Atom>( ca))
            )
          )
        ).second, "insert failed"
      );

      // print the GroupCollection
      BCL_MessageStd( "def_const" + util::Format()( def_const));

      // test copy constructor
      BCL_MessageStd( "test copy const");
      restraint::GroupCollection< std::string, biol::Atom> copy_const( def_const);
      BCL_MessageStd( "copy_const" + util::Format()( copy_const));
      BCL_Example_Check
      (
        copy_const.GetSize() == 2, "copy constructor did not properly construct copy"
      );

      // test clone copy constructor
      BCL_MessageStd( "test clone copy const");
      util::ShPtr< restraint::GroupCollection< std::string, biol::Atom> > virtual_copy( def_const.Clone());
      BCL_Example_Check
      (
        virtual_copy->GetSize() == 2, "clone copy constructor failed"
      );

      // test [] operator
      BCL_MessageStd( "test [] operator" + util::Format()( def_const[ "a"]));
      BCL_Example_Check
      (
        ( *def_const[ "a"].Begin())->GetType() == biol::GetAtomTypes().CA, "[] operator did not return correct data"
      );

      // test insert elements with two iterators
      BCL_MessageStd( "test insert elements using two iterators");
      biol::Atom o( biol::GetAtomTypes().O);
      BCL_Example_Check
      (
        (
          copy_const.Insert
          (
            storage::Pair< std::string, restraint::Group< biol::Atom> >
            (
              "d", restraint::Group< biol::Atom>( 2, util::SiPtr< const biol::Atom>( o))
            )
          )
        ).second, "insert failed test insert elements with two iterators"
      );
      def_const.InsertElements( copy_const.Begin(), copy_const.End());
      BCL_MessageStd( "def_const after insert between iterators" + util::Format()( def_const));
      BCL_Example_Check
      (
        ( *( --def_const.End())->second.Begin())->GetType() == biol::GetAtomTypes().O, "insert iterators is broken"
      );

      // test RemoveGroup function taking a Group identifier
      BCL_MessageStd( "test RemoveGroup function taking a group identifier");
      def_const.Erase( "a");
      BCL_Example_Check
      (
        def_const.Find( "a") == def_const.End(), "RemoveGroup function taking an identifier did not work"
      );

      // test RemoveGroup function taking an iterator
      BCL_MessageStd( "test RemoveGroup function taking an iterator");
      def_const.RemoveElement( def_const.Begin());
      BCL_Example_Check
      (
        def_const.Find( "b") == def_const.End(), "RemoveGroup function taking an iterator did not work"
      );

      // test Swap function
      BCL_MessageStd( "test Swap function");
      def_const.Swap( copy_const);
      BCL_Example_Check
      (
        def_const.GetSize() == 3 && copy_const.GetSize() == 1, "Swap function did not work"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRestraintGroupCollection

  const ExampleClass::EnumType ExampleRestraintGroupCollection::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintGroupCollection())
  );

} // namespace bcl
