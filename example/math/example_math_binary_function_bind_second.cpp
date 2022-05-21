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
#include "math/bcl_math_binary_function_bind_second.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_binary_function_bind_second.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathBinaryFunctionBindSecond :
    public ExampleInterface
  {
  public:

    ExampleMathBinaryFunctionBindSecond *Clone() const
    { return new ExampleMathBinaryFunctionBindSecond( *this);}

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    // class Less
    class Less :
      public util::BinaryFunctionInterface< double, double, bool>
    {
      static const util::SiPtr< const util::ObjectInterface> s_Instance;
    public:
      Less *Clone() const
      {
        return new Less( *this);
      }

      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      bool operator()( const double &LHS, const double &RHS) const
      {
        return LHS < RHS;
      }
      std::istream &Read( std::istream &ISTREAM)
      {
        // return
        return ISTREAM;
      }

      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // end
        return OSTREAM;
      }
    };

    // class Larger
    class Larger :
      public util::BinaryFunctionInterface< double, double, bool>
    {
      static const util::SiPtr< const util::ObjectInterface> s_Instance;
    public:
      Larger *Clone() const
      {
        return new Larger( *this);
      }

      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      bool operator()( const double &LHS, const double &RHS) const
      {
        return LHS > RHS;
      }
      std::istream &Read( std::istream &ISTREAM)
      {
        // return
        return ISTREAM;
      }

      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // end
        return OSTREAM;
      }
    };

    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      math::BinaryFunctionBindSecond< double, double, bool> binary_bind_second_default;

      // construct from reference to function interface and an argument to bind to as second
      math::BinaryFunctionBindSecond< double, double, bool> binary_bind_second_less_than_1( Less(), double( 1.0));

      // construct from sharedpointer to function interface and an argument to bind to second
      math::BinaryFunctionBindSecond< double, double, bool> binary_bind_second_greater_than_m1( util::ShPtr< util::BinaryFunctionInterface< double, double, bool> >( new Larger()), double( -1.0));

      // copy constructor
      math::BinaryFunctionBindSecond< double, double, bool> binary_bind_second_greater_than_m1_copy( binary_bind_second_greater_than_m1);

      // clone
      util::ShPtr< math::FunctionInterfaceSerializable< double, bool> > binder_cloned( binary_bind_second_greater_than_m1.Clone());

    /////////////////
    // data access //
    /////////////////

      // class identifiers
      BCL_MessageStd( "class name: " + binder_cloned->GetClassIdentifier());
      BCL_Example_Check
      (
        ( GetStaticClassName< math::BinaryFunctionBindSecond< double, double, bool> >() == "bcl::math::BinaryFunctionBindSecond<double,double,bool>")
        && ( binder_cloned->GetClassIdentifier() == GetStaticClassName< math::BinaryFunctionBindSecond< double, double, bool> >()),
          "incorrect class name: static class name: " + ( GetStaticClassName< math::BinaryFunctionBindSecond< double, double, bool> >())
        + " class identifier: " + binder_cloned->GetClassIdentifier()
      );

    ///////////////
    // operators //
    ///////////////

      // operator() on binary_bind_second_less_than_1
      BCL_Example_Check
      (
        binary_bind_second_less_than_1( double( 0.0)),
        "0.0 should be less than 1.0"
      );
      BCL_Example_Check
      (
        !binary_bind_second_less_than_1( double( 1.0)),
        "1.0 should not be less than 1.0"
      );
      BCL_Example_Check
      (
        !binary_bind_second_less_than_1( double( 2.0)),
        "2.0 should not be less than 1.0"
      );

      // operator() on binary_bind_second_greater_than_m1
      BCL_Example_Check
      (
        binary_bind_second_greater_than_m1( double( 0.0)),
        "0.0 should be greater than -1.0"
      );
      BCL_Example_Check
      (
        !binary_bind_second_greater_than_m1( double( -1.0)),
        "-1.0 should not be greater than -1.0"
      );
      BCL_Example_Check
      (
        !binary_bind_second_greater_than_m1( double( -2.0)),
        "-2.0 should not be greater than -1.0"
      );

      BCL_Example_Check
      (
        ( binary_bind_second_greater_than_m1( double( 0.0)) == binary_bind_second_greater_than_m1_copy( double( 0.0)))
        && ( binary_bind_second_greater_than_m1( double(  0.0)) == binder_cloned->operator()( double(  0.0)))
        && ( binary_bind_second_greater_than_m1( double( -1.0)) == binary_bind_second_greater_than_m1_copy( double( -1.0)))
        && ( binary_bind_second_greater_than_m1( double( -1.0)) == binder_cloned->operator()( double( -1.0)))
        && ( binary_bind_second_greater_than_m1( double( -2.0)) == binary_bind_second_greater_than_m1_copy( double( -2.0)))
        && ( binary_bind_second_greater_than_m1( double( -2.0)) == binder_cloned->operator()( double( -2.0))),
        "cloned and copied greater -1 should gave same results as original"
      );
      
    //////////////////////
    // input and output //
    //////////////////////

      // write to file
      WriteBCLObject( binary_bind_second_greater_than_m1);
      // read
      ReadBCLObject( binary_bind_second_default);

      BCL_Example_Check
      (
        binary_bind_second_greater_than_m1( double( 0.0)) == binary_bind_second_default( double( 0.0))
        && binary_bind_second_greater_than_m1( double( -1.0)) == binary_bind_second_default( double( -1.0))
        && binary_bind_second_greater_than_m1( double( -2.0)) == binary_bind_second_default( double( -2.0)),
        "read greater -1 should gave same results as original"
      );
    
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // end ExampleMathBinaryFunctionBindSecond

  const util::SiPtr< const util::ObjectInterface> ExampleMathBinaryFunctionBindSecond::Less::s_Instance
  (
    GetObjectInstances().AddInstance( new ExampleMathBinaryFunctionBindSecond::Less())
  );

  const util::SiPtr< const util::ObjectInterface> ExampleMathBinaryFunctionBindSecond::Larger::s_Instance
  (
    GetObjectInstances().AddInstance( new ExampleMathBinaryFunctionBindSecond::Larger())
  );

  const ExampleClass::EnumType ExampleMathBinaryFunctionBindSecond::s_Instance
  (
    GetExamples().AddEnum( ExampleMathBinaryFunctionBindSecond())
  );

} // namespace bcl
