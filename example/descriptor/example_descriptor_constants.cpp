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
#include "descriptor/bcl_descriptor_constants.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_example_string_sequence.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_constants.cpp
  //!
  //! @author mendenjl
  //! @date Jan 30, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorConstants :
    public ExampleInterface
  {
  public:

    ExampleDescriptorConstants *Clone() const
    {
      return new ExampleDescriptorConstants( *this);
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
      // form the storage::Vector< float> to be returned
      linal::Vector< float> list( 2);
      list( 0) = 1.0;
      list( 1) = 2.0;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      descriptor::Constants< char, float> constant_default;

      // constructor from a description
      descriptor::Constants< char, float> constant( list);

      // clone
      util::ShPtr< descriptor::Constants< char, float> > sp_constant( constant.Clone());

    /////////////////
    // data access //
    /////////////////

      // test the GetStaticClassName
      BCL_ExampleCheck( constant.GetClassIdentifier(), GetStaticClassName( constant));

      // test the GetLength function
      BCL_ExampleCheck( constant_default.GetNormalSizeOfFeatures(), 0);

      // test the GetLength function
      BCL_ExampleCheck( constant.GetNormalSizeOfFeatures(), 2);

    ///////////////
    // operators //
    ///////////////

      // initialize strings to be passed
      const descriptor::StringSequence string_a( "asdf"), string_b( "fff");

      // test the operators

      // check the result of default constructor
      BCL_ExampleCheck( constant_default( string_a), linal::Vector< float>());

      // check the result constructor from a list
      BCL_ExampleCheck( constant( string_a), list);

      // check that same result is given also with string b
      BCL_ExampleCheck( constant( string_a), constant( string_b));

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities");
      WriteBCLObject( constant);
      descriptor::Constants< char, float> constant_read;
      ReadBCLObject( constant_read);

      // check that same result is given also with object constructed with clone
      BCL_ExampleIndirectCheck( constant( string_a), constant_read( string_a), "I/O");

      // Check that we can also create the constants via implementation
      util::Implementation< descriptor::Base< char, float> > impl_constants;
      BCL_ExampleIndirectAssert
      (
        impl_constants.TryRead( util::ObjectDataLabel( "Constant(1.0,2.0)"), util::GetLogger()),
        true,
        "Reading into an implentation"
      );
      impl_constants->SetObject( string_a);
      descriptor::Iterator< char> itr( impl_constants->GetType(), string_a);
      BCL_ExampleIndirectCheck( ( *impl_constants)( itr), constant( string_a), "Implementation");
      BCL_ExampleIndirectAssert
      (
        impl_constants.TryRead( util::ObjectDataLabel( "Constant(nan)"), util::GetLogger()),
        true,
        "Reading into an implentation with undefined values"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorConstants

  const ExampleClass::EnumType ExampleDescriptorConstants::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorConstants())
  );

} // namespace bcl
