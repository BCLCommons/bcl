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
#include "util/bcl_util_sh_ptr_list.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_sh_ptr_vector.h"
#include "util/bcl_util_wrapper_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //! @example example_util_sh_ptr_list.cpp
  //!
  //! @author alexanns
  //! @date Nov 06, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleUtilShPtrList :
    public ExampleInterface
  {
  public:

    ExampleUtilShPtrList *Clone() const
    {
      return new ExampleUtilShPtrList( *this);
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
      // create array to t_DataType for inserting into ShPtrLists
      double data[ 4] =
      {
        double( 0.5),
        double( 0.2),
        double( 0.7),
        double( 0.1)
      };

      // create shared pointer vector to be used for adding elements to ShPtrLists
      util::ShPtrVector< double> sh_ptr_vector;
      sh_ptr_vector.PushBack( util::ShPtr< double >( new double( data[ 0])));
      sh_ptr_vector.PushBack( util::ShPtr< double >( new double( data[ 1])));
      sh_ptr_vector.PushBack( util::ShPtr< double >( new double( data[ 2])));
      sh_ptr_vector.PushBack( util::ShPtr< double >( new double( data[ 3])));
      BCL_MessageStd( util::Format()( sh_ptr_vector));

      // test default constructor with default arguments
      BCL_MessageStd( "test default constructor with default arguments");
      util::ShPtrList< double> def_const;
      //BCL_Example_Check
      //(
      //  ExampleClass::ExampleResult::e_Trivial,
      //  def_const.GetSize() == 0, "default constructor did not work"
      //);

      // test default constructor with passing size and default element argument
      BCL_MessageStd( "test default constructor passing size and default element argument");
      double nine_point_three( 9.3);
      util::ShPtr< double> default_element( new double( nine_point_three));
      util::ShPtrList< double> arg_const( 5, default_element);
      for
      (
        storage::List< util::ShPtr< double> >::iterator ptr_begin( arg_const.Begin()),
          ptr_end( arg_const.End());
        ptr_begin != ptr_end;
        ++ptr_begin
      )
      {
        BCL_Example_Check
        (
          *( *ptr_begin) == 9.3,
          "default constructor passing size and default element failed:wrong element"
        );
      }
      BCL_Example_Check
      (
        arg_const.GetSize() == 5,
        "default constructor passing size and default element failed: size wrong"
      );
      BCL_MessageStd( util::Format()( arg_const));

      // test construct from number and pointer to data
      BCL_MessageStd( "test construct from number and pointer to data");
      util::ShPtrList< double> num_and_data( 4, data);
      BCL_Example_Check
      (
        num_and_data.GetSize() == 4,
        "size is not correct in construction from number and pointer to data"
      );
      BCL_MessageStd( util::Format()( num_and_data));
      size_t index( 0);
      for
      (
        storage::List< util::ShPtr< double> >::iterator ptr_begin( num_and_data.Begin()),
          ptr_end( num_and_data.End());
        ptr_begin != ptr_end;
        ++ptr_begin
      )
      {
        BCL_Example_Check
        (
          *( *ptr_begin) == *sh_ptr_vector( index++),
          "default constructor passing size and default element failed:wrong element"
        );
      }

      // test construct from range of iterators
      BCL_MessageStd( "test construct from range of iterators");
      util::ShPtrList< double> itr_range_constr( sh_ptr_vector.Begin(), sh_ptr_vector.End());
      //BCL_Example_Check
      //(
      //  ExampleClass::ExampleResult::e_Trivial,
      //  itr_range_constr == sh_ptr_vector, "construct from range of iterators is broken"
      //);

      // test copy constructor
      BCL_MessageStd( "test copy constructor");
      util::ShPtrList< double> copy_constr( arg_const);
      BCL_Example_Check
      (
        copy_constr == arg_const, "copy constructor did not work properly"
      );

      // test clone constructor
      BCL_MessageStd( "test clone constructor");
      BCL_Example_Check
      (
        *( arg_const.Clone()) == arg_const, "clone constructor did not work"
      );

      // test InsertElement function
      BCL_MessageStd( "test InsertElement function taking a ShPtr<t_DataType>");
      def_const.PushBack( default_element);
      BCL_Example_Check
      (
        *( *( --def_const.ReverseEnd())) == 9.3,
        "InsertElement taking ShPtr<t_DataType> is broken"
      );

      // test HardCopy function
      BCL_MessageStd( "test HardCopy function");
      util::ShPtrList< double> hard_copy( def_const.HardCopy());
      for
      (
        storage::List< util::ShPtr< double> >::iterator
          itr_def_const( def_const.Begin()),
          itr_end_def_const( def_const.End()),
          itr_hard_copy( hard_copy.Begin()),
          itr_end_hard_copy( hard_copy.End());
        itr_def_const != itr_end_def_const && itr_hard_copy != itr_end_hard_copy;
        ++itr_def_const, ++itr_hard_copy
      )
      {
        BCL_Example_Check
        (
          ( *itr_def_const != *itr_hard_copy) && ( **itr_def_const == **itr_hard_copy),
          "HardCopy function did not work"
        );
      }

      // test GetSiPtrList function
      BCL_MessageStd( "test GetSiPtrList function");
      util::SiPtrList< double> si_ptr_list( hard_copy);
//      BCL_Example_Check
//      (
//        ExampleClass::ExampleResult::e_Trivial,
//        hard_copy == si_ptr_list,
//        "conversion from ShPtrList to SiPtrList function did not work"
//      );

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleUtilShPtrList

  const ExampleClass::EnumType ExampleUtilShPtrList::s_Instance
  (
    GetExamples().AddEnum( ExampleUtilShPtrList())
  );

} // namespace bcl
