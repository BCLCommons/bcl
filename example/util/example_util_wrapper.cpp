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
#include "util/bcl_util_wrapper.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_util_wrapper.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleUtilWrapper :
    public ExampleInterface
  {
  public:

    ExampleUtilWrapper *Clone() const
    {
      return new ExampleUtilWrapper( *this);
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

      // instantiate the s_Instance so that we can read/write from ShPtr's
      util::Wrapper< std::string>::s_Instance.IsDefined();

      // default constructor
      util::Wrapper< std::string> wrapped_string_empty;
      BCL_Example_Check
      (
        wrapped_string_empty.empty(),
        "the wrapped std::string default constructor does not create an empty string"
      );

      //wrap a string as a bcl class
      util::Wrapper< std::string> wrapped_string( std::string( "hello"));
      BCL_MessageStd( "wrapped std::string: " + util::Format()( wrapped_string));

      BCL_Example_Check
      (
        wrapped_string == std::string( "hello"),
        "the wrapped std::string is not initialized to \"hello\""
      );

      // copy constructor
      util::Wrapper< std::string> wrapped_string_copy( wrapped_string);
      BCL_Example_Check
      (
        wrapped_string_copy == wrapped_string && wrapped_string_copy == std::string( "hello"),
        "the copy of the wrapped std::string is not initialized to \"hello\""
      );

      // clone
      util::Wrapper< std::string> *ptr( wrapped_string.Clone());

    /////////////////
    // data access //
    /////////////////

      BCL_Example_Check
      (
        GetStaticClassName< util::Wrapper< std::string> >() == ptr->GetClassIdentifier(),
        "incorrect class identifier"
      );

      BCL_MessageStd
      (
        "name of Wrapped std::string: " + GetStaticClassName< util::Wrapper< std::string> >()
      );

      // direct access is granted by access function gedata (const and non-const)
      BCL_MessageStd
      (
        "content of wrapped std::string: " + wrapped_string.GetData()
      );
      BCL_Example_Check
      (
        wrapped_string.GetData() == "hello",
        "GetData() does not return hello but: " + wrapped_string.GetData()
      );

      // change data
      wrapped_string_copy.GetData() += " bcl user";
      BCL_Example_Check
      (
        wrapped_string_copy.GetData() == "hello bcl user",
        "GetData() does not return hello bcl user but: " + wrapped_string_copy.GetData()
      );

    ////////////////
    // operations //
    ////////////////

      // operations of the wrapped class are still accessible
      BCL_MessageStd( "1st three chars of wrapped string: " + wrapped_string.substr( 0, 3));
      BCL_Example_Check
      (
        wrapped_string.substr( 0, 3) == std::string( "hel"),
        "the wrapped std::string substr( 0, 3) does not return \"hel\""
      );

      // create a util::ShPtr on a wrapped string
      // the following code would not be possible:
      // util::ShPtr< std::string> sp_string( new std::string( "hello"));
      // sp_string.HardCopy();
      // but this works with a wrapper
      util::ShPtr< util::Wrapper< std::string> > sp_wrapped_string( wrapped_string.Clone());
      BCL_MessageStd
      (
        "this is a clone of the wrapped string: " + util::Format()( sp_wrapped_string)
      );

      // convert the wrapper back into its original type
      std::string original( wrapped_string);
      // origianl datatype behaves like wrapped datatype
      BCL_Example_Check
      (
        wrapped_string == original,
        "the wrapped std::string and original data type are not the same"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // write object
      WriteBCLObject( wrapped_string);
      // read object
      ReadBCLObject( wrapped_string_empty);

      BCL_Example_Check
      (
        wrapped_string_empty == "hello" && wrapped_string_empty == wrapped_string,
        "writing and reading of wrapped std::string did not work"
      );

      // cleanup
      delete ptr;

      //successful
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // end ExampleUtilWrapper

  const ExampleClass::EnumType ExampleUtilWrapper::s_Instance
  (
    GetExamples().AddEnum( ExampleUtilWrapper())
    );

} // namespace bcl
