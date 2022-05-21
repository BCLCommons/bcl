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
#include "util/bcl_util_object_interface.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //! @example example_util_object_interface.cpp
  //!
  //! @author woetzen
  //! @date
  //! @remarks status empty
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleUtilObjectInterface :
    public ExampleInterface
  {
  public:

    ExampleUtilObjectInterface *Clone() const
    { return new ExampleUtilObjectInterface( *this);}

    class TestClass1 :
      public util::ObjectInterface
    {
    public:
      TestClass1 *Clone() const
      {
        return new TestClass1( *this);
      }

      std::string const &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      virtual std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    };

    template< typename T1, typename T2>
    class TestClass2 :
      public util::ObjectInterface
    {
    public:
      TestClass2< T1, T2> *Clone() const
      {
        return new TestClass2< T1, T2>( *this);
      }

      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      virtual std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    };

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
      util::ShPtr< util::ObjectInterface> test1( new TestClass1());
      util::ShPtr< util::ObjectInterface> test2( new TestClass2< TestClass1, int>());

      util::ShPtr< util::ObjectInterface> test_copy1( test1->Clone());
      util::ShPtr< util::ObjectInterface> test_copy2( test2->Clone());

      BCL_MessageStd( "ClassName of test1: |" + test1->GetClassIdentifier() + "|");
      BCL_MessageStd( "ClassName of test2: |" + test2->GetClassIdentifier() + "|");

      BCL_MessageStd( "ClassName of test_copy1: |" + test_copy1->GetClassIdentifier() + "|");
      BCL_MessageStd( "ClassName of test_copy2: |" + test_copy2->GetClassIdentifier() + "|");

//      BCL_MessageStd( "Template list of TestClass2: |" + ( TestClass2< double, int>().GetTemplteList()) + "|");

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // end ExampleBCL

  const ExampleClass::EnumType ExampleUtilObjectInterface::s_Instance
  (
    GetExamples().AddEnum( ExampleUtilObjectInterface())
  );

} // namespace bcl
