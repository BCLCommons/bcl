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
#include "restraint/bcl_restraint_group.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_wrapper_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  //! binary predicate for inverse sorting, used by Sort
  template< typename t_DataType> struct more_than :
    public std::binary_function< util::SiPtr< const t_DataType>, util::SiPtr< const t_DataType>, bool>
  {
    bool operator()( const util::SiPtr< const t_DataType> &DATA_A, const util::SiPtr< const t_DataType> &DATA_B) const
    { return *DATA_A > *DATA_B;}
  };

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_group.cpp
  //!
  //! @author heinzes1
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRestraintGroup :
    public ExampleInterface
  {
  public:

    ExampleRestraintGroup *Clone() const
    { return new ExampleRestraintGroup( *this);}

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
      // create array to t_DataType for inserting into SiPtrLists
      double data[ 4] =
      {
        double( 0.5),
        double( 0.2),
        double( 0.7),
        double( 0.1)
      };

      // test default constructor
      BCL_MessageStd( "test default constructor");
      restraint::Group< double> def_constr;
      BCL_Example_Check
      (
        def_constr.IsEmpty(), "default constructor did not work"
      );

      // test constructing from size and giving default element
      BCL_MessageStd( "test constructing from size and giving default element");
      double one_point_eight( 1.8);
      restraint::Group< double> size_elem_constr
      (
        6, //< size
        util::ToSiPtr( one_point_eight) //< default element
      );
      BCL_Example_Check
      (
        size_elem_constr.GetSize() == 6, "constructing from size and giving default element wrong size"
      );
      for
      (
        restraint::Group< double>::const_iterator
          itr( size_elem_constr.Begin()), itr_end( size_elem_constr.End());
        itr != itr_end;
        ++itr
      )
      {
        BCL_Example_Check
        (
          **itr == one_point_eight, "constructing from size and giving default element wrong element"
        );
      }

      // test construct from size and pointer to data
      BCL_MessageStd( "test construct from size and pointer to data");
      restraint::Group< double> size_pointer_constr( 4, data);
      double *ptr( data), *ptr_end( data + 4);
      for
      (
        restraint::Group< double>::const_iterator
          itr( size_pointer_constr.Begin()), itr_end( size_pointer_constr.End());
        itr != itr_end && ptr != ptr_end;
        ++itr, ++ptr
      )
      {
        BCL_Example_Check
        (
          **itr == *ptr, "test construct from size and pointer to data wrong "
          + util::Format()( **itr) + " != " + util::Format()( *ptr)
        );
      }

      // test construct from SiPtrList
      BCL_MessageStd( "test construct SiPtrList");
      double two_point_three( 2.3);
      util::SiPtrList< const double> si_ptr_list
      (
        3, util::SiPtr< const double>( two_point_three)
      );
      restraint::Group< double> si_ptr_list_constr( si_ptr_list);
      BCL_Example_Check
      (
        si_ptr_list == si_ptr_list_constr, "construction from SiPtrList did not work"
      );

      // test construct from iterator range
      BCL_MessageStd( "test construct from iterator range");
      restraint::Group< double> itr_range_constr
      (
        si_ptr_list_constr.Begin(), si_ptr_list_constr.End()
      );
      BCL_Example_Check
      (
        si_ptr_list_constr == itr_range_constr, "construct from range of iterators did not work"
      );

      // test construct from two existing Groups
      BCL_MessageStd( "test construct from two existing Groups");
      restraint::Group< double> two_groups_constr( size_elem_constr, size_elem_constr);
      BCL_Example_Check
      (
        two_groups_constr.GetSize() == 12, "construct from two existing Groups gave wrong size"
      );
      for
      (
        restraint::Group< double>::const_iterator
          itr( two_groups_constr.Begin()), itr_end( two_groups_constr.End());
        itr != itr_end;
        ++itr
      )
      {
        BCL_Example_Check
        (
          **itr == one_point_eight, "constructing from two existing Groups did not work"
        );
      }

      // test construct from a Group and a number of empty members
      BCL_MessageStd( "test construct from a Group and a number of empty members");
      restraint::Group< double> group_empty_constr( si_ptr_list_constr, 1);
      BCL_Example_Check
      (
        group_empty_constr.GetSize() == 4, "construct from a Group and a number of empty members gave wrong size"
      );
      BCL_Example_Check
      (
        !( ( --group_empty_constr.End())->IsDefined()), "constructing from a Group and a number of empty members did not add empty members properly"
      );
      for
      (
        restraint::Group< double>::const_iterator
          itr( group_empty_constr.Begin()), itr_end( --group_empty_constr.End());
        itr != itr_end;
        ++itr
      )
      {
        BCL_Example_Check
        (
          **itr == two_point_three, "constructing from a Group and a number of empty members did not add Group properly"
        );
      }

      // test construct from a number of empty members and a Group
      BCL_MessageStd( "test construct from a number of empty members and a Group");
      restraint::Group< double> empty_group_constr( 1, si_ptr_list_constr);
      BCL_Example_Check
      (
        group_empty_constr.GetSize() == 4, "construct from a number of empty members and a Group gave wrong size"
      );
      BCL_Example_Check
      (
        !( empty_group_constr.Begin()->IsDefined()), "constructing from a number of empty members and a Group did not add empty members properly"
      );
      for
      (
        restraint::Group< double>::const_iterator
          itr( ++empty_group_constr.Begin()), itr_end( empty_group_constr.End());
        itr != itr_end;
        ++itr
      )
      {
        BCL_Example_Check
        (
          **itr == two_point_three, "constructing from a number of empty members and a Group did not add Group properly"
        );
      }

      // test copy constructor
      BCL_MessageStd( "test copy constructor");
      restraint::Group< double> copy_constr( size_pointer_constr);
      BCL_Example_Check
      (
        copy_constr == size_pointer_constr, "construct from range of iterators did not work"
      );

      // test virtual copy constructor
      BCL_MessageStd( "test virtual copy constructor");
      BCL_Example_Check
      (
        *util::ShPtr< restraint::Group< double> >( copy_constr.Clone()) == copy_constr,
        "clone did not work"
      );

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRestraintGroup

  const ExampleClass::EnumType ExampleRestraintGroup::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintGroup())
  );

} // namespace bcl
