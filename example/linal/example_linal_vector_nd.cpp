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
#include "linal/bcl_linal_vector_nd.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_linal_vector_nd.cpp
  //!
  //! @author woetzen
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleLinalVectorND :
    public ExampleInterface
  {
  public:

    ExampleLinalVectorND *Clone() const
    {
      return new ExampleLinalVectorND( *this);
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
      linal::VectorND< double, 3> vector_doub_3_default;

      // construct from an element
      linal::VectorND< double, 4> vector_doub_4_default( 1.0);

    /////////////////
    // data access //
    /////////////////

      // size
      BCL_ExampleCheck( ( linal::VectorND< double, 3>().GetSize()), 3);

      // element access
      BCL_ExampleCheck( ( linal::VectorND< double, 4>( 1.0))( 0), 1.0);

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////
    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleLinalVectorND

  const ExampleClass::EnumType ExampleLinalVectorND::s_Instance
  (
    GetExamples().AddEnum( ExampleLinalVectorND())
  );

} // namespace bcl
