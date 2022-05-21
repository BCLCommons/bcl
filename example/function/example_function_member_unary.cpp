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
#include "function/bcl_function_member_unary.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_function_member_unary.cpp
  //!
  //! @author woetzen
  //! @date Nov 22, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFunctionMemberUnary :
    public ExampleInterface
  {
  public:

    ExampleFunctionMemberUnary *Clone() const
    {
      return new ExampleFunctionMemberUnary( *this);
    }

    struct Person
    {

      std::string m_Name;

      Person( const std::string &NAME) :
        m_Name( NAME)
      {
      }

      const std::string &GetName() const
      {
        return m_Name;
      }

      void SetName( const std::string &NAME)
      {
        m_Name = NAME;
      }

      bool SetNameSuccess( const std::string &NAME)
      {
        m_Name = NAME;
        return true;
      }

    }; // class Person

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
      const std::string joes_name( "Joe");

      // make a person
      Person joe( joes_name);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct from ptr to member function
      function::MemberUnary< Person, const std::string, void> person_set_name( &Person::SetName);
      function::MemberUnary< Person, const std::string, bool> person_set_name_success( &Person::SetNameSuccess);

      // MakeFunction
      util::ShPtr< function::BinaryInterface< Person, const std::string, void> > binary_function_ptr( function::MemberFunction( &Person::SetName).Clone());
      util::ShPtr< function::BinaryInterface< Person, const std::string, bool> > binary_function_ptr_success( function::MemberFunction( &Person::SetNameSuccess).Clone());

    /////////////////
    // data access //
    /////////////////

      // check class identifier
      BCL_ExampleCheck( person_set_name.GetClassIdentifier(), GetStaticClassName( person_set_name));

    ///////////////
    // operators //
    ///////////////

      BCL_ExampleCheck( joe.GetName(), joes_name);
      person_set_name( joe, "Joey");
      BCL_ExampleIndirectCheck( joe.GetName() == joes_name, false, "test of SetName function through functor");

      // test through pointer
      binary_function_ptr->operator()( joe, joes_name);
      BCL_ExampleIndirectCheck( joe.GetName() == joes_name, true, "test of SetName function through functor ptr");

    //////////////////////
    // input and output //
    //////////////////////

      // write object
      WriteBCLObject( person_set_name);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFunctionMemberUnary

  const ExampleClass::EnumType ExampleFunctionMemberUnary::s_Instance
  (
    GetExamples().AddEnum( ExampleFunctionMemberUnary())
  );

} // namespace bcl
