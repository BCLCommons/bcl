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
#include "math/bcl_math_const_function.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_const_function.cpp
  //!
  //! @author butkiem1
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathConstFunction :
    public ExampleInterface
  {
  public:

    ExampleMathConstFunction *Clone() const
    {
      return new ExampleMathConstFunction( *this);
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

      // default constructor
      math::ConstFunction< size_t, size_t> const_function_default;

      // constructor with parameter
      math::ConstFunction< size_t, size_t> const_function_init( size_t( 26));

      // test clone method
      util::ShPtr< math::ConstFunction< size_t, size_t> > const_function_clone( const_function_init.Clone());

    /////////////////
    // data access //
    /////////////////

      // check class identifier
      BCL_Example_Check
      (
        ( GetStaticClassName< math::ConstFunction< size_t, size_t> >() == "bcl::math::ConstFunction<size_t,size_t>"),
        "incorrect static class name"
      );

      BCL_Example_Check
      (
        ( GetStaticClassName< math::ConstFunction< size_t, size_t> >() == const_function_clone->GetClassIdentifier()),
        "incorrect class identifier"
      );

      BCL_Example_Check
      (
        const_function_default.GetScheme() == "ConstFunction",
        "incorrect scheme identifier"
      );

    ///////////////
    // operators //
    ///////////////

      // check operator
      BCL_Example_Check
      (
        const_function_init( size_t( 5)) == size_t( 26),
        "ConstFunction operator does not return argument value! is" +
        util::Format()( const_function_init( size_t( 555))) + " should be " + util::Format()( size_t( 26))
      );

      // check operator
      BCL_Example_Check
      (
        const_function_init( size_t( 0)) == size_t( 26),
        "ConstFunction operator does not return argument value! is" +
        util::Format()( const_function_init( size_t( 0))) + " should be " + util::Format()( size_t( 26))
      );

      // check operator
      BCL_Example_Check
      (
        const_function_init( size_t( -99)) == size_t( 26),
        "ConstFunction operator does not return argument value! is" +
        util::Format()( const_function_init( size_t( -99))) + " should be " + util::Format()( size_t( 26))
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // write ConstFunction to file
      WriteBCLObject( const_function_init);

      // ConstFunction object for read from bcl file
      math::ConstFunction< size_t, size_t> const_function_read_in;

      // read ConstFunction
      ReadBCLObject( const_function_read_in);

      // check class identifier
      BCL_Example_Check
      (
        const_function_init( size_t( 7)) == size_t( 26),
        "ConstFunction operator does not return argument value! is" +
        util::Format()( const_function_init( size_t( 7))) + " should be " + util::Format()( size_t( 26))
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMathConstFunction

  const ExampleClass::EnumType ExampleMathConstFunction::s_Instance
  (
    GetExamples().AddEnum( ExampleMathConstFunction())
  );

} // namespace bcl
