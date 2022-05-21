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
#include "util/bcl_util_si_ptr_list.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //! @example example_util_si_ptr_list.cpp
  //!
  //! @author alexanns
  //! @date Nov 06, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleUtilSiPtrList :
    public ExampleInterface
  {
  public:

    ExampleUtilSiPtrList *Clone() const
    {
      return new ExampleUtilSiPtrList( *this);
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
      // create array to t_DataType for inserting into SiPtrLists
      double data[ 4] =
      {
        double( 0.5),
        double( 0.2),
        double( 0.7),
        double( 0.1)
      };

      // create simple pointer vector to be used for adding elements to SiPtrLists
      util::SiPtrVector< double> si_ptr_vector;
      si_ptr_vector.PushBack( util::ToSiPtrNonConst( data[ 0]));
      si_ptr_vector.PushBack( util::ToSiPtrNonConst( data[ 1]));
      si_ptr_vector.PushBack( util::ToSiPtrNonConst( data[ 2]));
      si_ptr_vector.PushBack( util::ToSiPtrNonConst( data[ 3]));
      BCL_MessageStd( util::Format()( si_ptr_vector));

      // test default constructor with default arguments
      BCL_MessageStd( "test default constructor with default arguments");
      util::SiPtrList< double> def_const;
      BCL_Example_Check
      (
        def_const.GetSize() == 0, "default constructor did not work"
      );

      // test default constructor with passing size and default element argument
      BCL_MessageStd( "test default constructor passing size and default element argument");
      double nine_point_three( 9.3);
      util::SiPtr< double> default_element( &nine_point_three);
      util::SiPtrList< double> arg_const( 5, default_element);
      for
      (
        storage::List< util::SiPtr< double> >::iterator ptr_begin( arg_const.Begin()),
          ptr_end( arg_const.End());
        ptr_begin != ptr_end;
        ++ptr_begin
      )
      {
        BCL_Example_Check
        (
          *( *ptr_begin) == 9.3, "default constructor passing size and default element failed:wrong element"
        );
      }
      BCL_Example_Check
      (
        arg_const.GetSize() == 5, "default constructor passing size and default element failed: size wrong"
      );
      BCL_MessageStd( util::Format()( arg_const));

      // test construct from number and pointer to data
      BCL_MessageStd( "test construct from number and pointer to data");
      util::SiPtrList< double> num_and_data( 4, data);
      BCL_Example_Check
      (
        num_and_data.GetSize() == 4, "size is not correct in construction from number and pointer to data"
      );
      BCL_MessageStd( util::Format()( num_and_data));
      //BCL_Example_Check
      //(
      //  ExampleClass::ExampleResult::e_Trivial,
      //  num_and_data == si_ptr_vector, "construct from number and pointer to data does not work"
      //);

      // test construct from range of iterators
      BCL_MessageStd( "test construct from range of iterators");
      util::SiPtrList< double> itr_range_constr( si_ptr_vector.Begin(), si_ptr_vector.End());
      //BCL_Example_Check
      //(
      //  ExampleClass::ExampleResult::e_Trivial,
      //  itr_range_constr == si_ptr_vector, "construct from range of iterators is broken"
      //);

      // test copy constructor
      BCL_MessageStd( "test copy constructor");
      util::SiPtrList< double> copy_constr( arg_const);
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

      // test PushFront with argument being pointer to datatype
      BCL_MessageStd( "test PushFront with argument being pointer to datatype");
      double five_point_two( 5.2);
      arg_const.PushFront( &five_point_two);
      BCL_Example_Check
      (
        *( *( arg_const.Begin())) == 5.2, "PushFront did not work"
      );

      // test PushBack with argument being pointer to datatype
      BCL_MessageStd( "test PushBack with argument being pointer to datatype");
      double four_point_size( 4.6);
      arg_const.PushBack( &four_point_size);
      BCL_Example_Check
      (
        *( *( arg_const.Last())) == 4.6, "PushFront did not work"
      );

      // test ConvertToSiPtrList function
      BCL_MessageStd( "test ConvertToSiPtrList function");
      storage::Vector< double> vector;
      vector.PushBack( double( double( 7.5)));
      vector.PushBack( double( double( 8.2)));
      vector.PushBack( double( double( 2.7)));
      vector.PushBack( double( double( 10.1)));
      util::SiPtrList< double> convert_to_si_ptr_list
      (
        util::ConvertToSiPtrList< double>( data, data + 4)
      );
      util::SiPtrList< const double> convert_to_const_si_ptr_list
      (
        util::ConvertToSiPtrList< const double>( data, data + 4)
      );
      util::SiPtrList< const double> convert_to_const_si_ptr_list_from_vector
      (
        util::ConvertToSiPtrList< const double>( vector.Begin(), vector.End())
      );

      // test equal comparison operator with SiPtrList to a SiPtrVector
      //BCL_Example_Check
      //(
      //  ExampleClass::ExampleResult::e_Trivial,
      //  convert_to_si_ptr_list == si_ptr_vector,
      //  "SiPtr list and SiPtrVector should be identical - when comparing addresses they are pointing to"
      //);

      // test equal comparison operator with a SiPtr and a SiPtr of const data
      //BCL_Example_Check
      //(
      //  ExampleClass::ExampleResult::e_Trivial,
      //  convert_to_si_ptr_list == convert_to_const_si_ptr_list,
      //  "SiPtr list and SiPtr const list should be identical - when comparing addresses they are pointing to"
      //);

      // test equal comparison operator with SiPtrList to SiPtrVector where the addresses differ
      //BCL_Example_Check
      //  ExampleClass::ExampleResult::e_Trivial,
      //(
      //  !( convert_to_const_si_ptr_list == convert_to_const_si_ptr_list_from_vector),
      //  "SiPtr list and SiPtrVector should be different - when comparing adresses they are pointing to"
      //);

      // test InsertElement function
      BCL_MessageStd( "test InsertElement function taking a SiPtr<t_DataType>");
      def_const.PushBack( default_element);
      BCL_Example_Check
      (
        *( *( --def_const.ReverseEnd())) == 9.3, "InsertElement taking SiPtr<t_DataType> is broken"
      );

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleUtilSiPtrList

  const ExampleClass::EnumType ExampleUtilSiPtrList::s_Instance
  (
    GetExamples().AddEnum( ExampleUtilSiPtrList())
  );

} // namespace bcl
