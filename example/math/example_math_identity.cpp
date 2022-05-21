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
#include "math/bcl_math_identity.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_identity.cpp
  //!
  //! @author butkiem1
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathIdentity :
    public ExampleInterface
  {
  public:

    ExampleMathIdentity *Clone() const
    {
      return new ExampleMathIdentity( *this);
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
      math::Identity< size_t> identity_default;

      // test clone method
      util::ShPtr< math::Identity< size_t> > identity_clone( identity_default.Clone());

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck( GetStaticClassName< math::Identity< size_t> >(), identity_clone->GetClassIdentifier());
      BCL_ExampleCheck( identity_default.GetScheme(), "identity");

    ///////////////
    // operators //
    ///////////////

      BCL_ExampleCheck( identity_default( size_t( 5)), 5);
      BCL_ExampleCheck( identity_default( size_t( -99)), -99);

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // write identity to file
      WriteBCLObject( identity_default);

      // identity object for read from bcl file
      math::Identity< size_t> identity_read_in;

      // read identity
      ReadBCLObject( identity_read_in);

      // check class identifier
      BCL_ExampleIndirectCheck( identity_read_in( size_t( 7)), 7, "Read/Write");

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMathIdentity

  const ExampleClass::EnumType ExampleMathIdentity::s_Instance
  (
    GetExamples().AddEnum( ExampleMathIdentity())
  );

} // namespace bcl
