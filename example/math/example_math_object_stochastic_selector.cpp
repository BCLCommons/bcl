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
#include "math/bcl_math_object_stochastic_selector.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_object_stochastic_selector.cpp
  //! @brief this example tests the implementation of ObjectStochasticSelector
  //!
  //! @author geanesar
  //! @date Nov 1, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathObjectStochasticSelector :
    public ExampleInterface
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleMathObjectStochasticSelector
    ExampleMathObjectStochasticSelector *Clone() const
    {
      return new ExampleMathObjectStochasticSelector( *this);
    }

  //////////
  // data //
  //////////

    //! @brief a test function
    class TestFunction :
      public util::SerializableInterface
    {
    private:

      //! static instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

      //! @brief default constructor
      TestFunction()
      {
      }

      //! @brief Clone function
      //! @return a copy of this class
      TestFunction *Clone() const
      {
        return new TestFunction( *this);
      }

      //! @brief get class name
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( this);
      }

      //! @brief get class name
      //! @return the class name
      const std::string &GetAlias() const
      {
        static const std::string name( "Test");
        return name;
      }

      //! @brief get information about member data
      //! @return serializer containing member data information
      io::Serializer GetSerializer() const
      {
        return io::Serializer();
      }

    };

      //! @brief a second test function
    class TestFunctionTwo :
      public util::SerializableInterface
    {
    public:

      //! @brief default constructor
      TestFunctionTwo()
      {
      }

      //! @brief Clone function
      //! @return a copy of this class
      TestFunctionTwo *Clone() const
      {
        return new TestFunctionTwo( *this);
      }

      //! @brief get class name
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( this);
      }

      //! @brief get information about member data
      //! @return serializer containing member data information
      io::Serializer GetSerializer() const
      {
        return io::Serializer();
      }

    };

  /////////////////
  // data access //
  /////////////////

      //! @brief returns the class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

  ////////////////
  // operations //
  ////////////////

      //! @brief run routine
      //! this is performing the execution of the example
      int Run() const
      {

        // Construction
        math::ObjectStochasticSelector< util::SerializableInterface> selector;

        // Serialized construction
        std::string constructor_string( "(options=(Test),probabilities=(1.0))");
        math::ObjectStochasticSelector< util::SerializableInterface> serial_selector;
        BCL_ExampleCheck( serial_selector.TryRead( constructor_string, util::GetLogger()), true);

        // Class name
        BCL_ExampleCheck( selector.GetAlias(), "StochasticSelector");

        // Implementation of test function
        TestFunction test_fxn;
        util::Implementation< util::SerializableInterface> impl( test_fxn);

        // Add
        selector.AddImplementation( impl, 0.7);

        // Test the implementation was added
        BCL_ExampleCheck( selector.GetSize(), 1);

        // Create a map for construction of another selector
        storage::Map< util::Implementation< util::SerializableInterface>, double> obj_map;

        obj_map[ impl] = 1.0;

        math::ObjectStochasticSelector< util::SerializableInterface> map_selector( obj_map);

        BCL_ExampleCheck( map_selector.GetSize(), 1);

        const util::SerializableInterface &test( selector.SelectRandomCase());

        // Check selection works
        BCL_ExampleIndirectCheck
        ( 
          test.GetClassIdentifier(), 
          TestFunction().GetClassIdentifier(),
          "checking that selecting a case works"
        );

        TestFunctionTwo test_fxn_two;

        //! Add a second implementation using a reference
        selector.AddImplementation( test_fxn_two, 0.3);

        BCL_ExampleIndirectCheck
        (
          selector.GetSize(),
          2,
          "checking if addition of another implementation works"
        );

        return 0;
      }

      static const ExampleClass::EnumType s_Instance;

  }; // class ExampleMathObjectStochasticSelector

  const ExampleClass::EnumType ExampleMathObjectStochasticSelector::s_Instance
  (
    GetExamples().AddEnum( ExampleMathObjectStochasticSelector())
  );

  const util::SiPtr< const util::ObjectInterface> ExampleMathObjectStochasticSelector::TestFunction::s_Instance
  (
    util::Enumerated< util::SerializableInterface>::AddInstance
    (
      new ExampleMathObjectStochasticSelector::TestFunction()
    )
  );

} // namespace bcl
