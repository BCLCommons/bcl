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
#include "function/bcl_function_member_const.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_function_member_const.cpp
  //!
  //! @author woetzen
  //! @date Nov 26, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFunctionMemberConst :
    public ExampleInterface
  {
  public:

    ExampleFunctionMemberConst *Clone() const
    {
      return new ExampleFunctionMemberConst( *this);
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
      const Person joe( joes_name);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct from ptr to member function
      function::MemberConst< Person, const std::string &> person_get_name( &Person::GetName);

      // MakeFunction
      util::ShPtr< function::UnaryInterface< const Person, const std::string &> > unary_function_ptr( function::MemberFunction( &Person::GetName).Clone());

    /////////////////
    // data access //
    /////////////////

      // check class identifier
      BCL_ExampleCheck( person_get_name.GetClassIdentifier(), GetStaticClassName( person_get_name));

    ///////////////
    // operators //
    ///////////////

      BCL_ExampleCheck( joe.GetName(), joes_name);
      BCL_ExampleCheck( person_get_name( joe), joes_name);

      // test through pointer
      BCL_ExampleCheck( unary_function_ptr->operator ()( joe), joes_name);

    //////////////////////
    // input and output //
    //////////////////////

      // write object
      WriteBCLObject( person_get_name);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFunctionMemberConst

  const ExampleClass::EnumType ExampleFunctionMemberConst::s_Instance
  (
    GetExamples().AddEnum( ExampleFunctionMemberConst())
  );

} // namespace bcl
