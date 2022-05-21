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
#include "util/bcl_util_implementation.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_util_implementation.cpp
  //!
  //! @author mendenjl
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class DoubleFuncInterface :
    public util::SerializableInterface
  {

  public:

    //! @brief Clone function
    //! @return pointer to new class
    virtual DoubleFuncInterface *Clone() const = 0;

    //! @brief virtual operator () to perform some operation with a double
    //! @return the value of the operation
    virtual double operator()( const double &VAL) const = 0;

  //////////////////////
  // input and output //
  //////////////////////

  protected:

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  };

  class SquareDouble :
    public DoubleFuncInterface
  {
  public:
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    static const util::SiPtr< const util::ObjectInterface> s_Instance; //!< True if/when this class has been registered with Enumerated

    //! @brief Clone function
    //! @return pointer to new class_name
    SquareDouble *Clone() const
    {
      return new SquareDouble( *this);
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

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &GetAlias() const
    {
      static const std::string s_Alias( "Square");
      return s_Alias;
    }

  ///////////////
  // operators //
  ///////////////

    double operator()( const double &VAL) const
    {
      return VAL * VAL;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Squares a double");
      return parameters;
    }

  };

  // add the interface to the set of known implementations
  const util::SiPtr< const util::ObjectInterface> SquareDouble::s_Instance( util::Enumerated< DoubleFuncInterface>::AddInstance( new SquareDouble()));

  class MultiplyDouble :
    public DoubleFuncInterface
  {
  private:

    double m_Coefficient; //!< coefficient to multiply by
    storage::Vector< double> m_Coefficients;
    storage::Vector< storage::Vector< double> > m_Coefficients2D;

  public:

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    static const util::SiPtr< const util::ObjectInterface> s_Instance; //!< True if/when this class has been registered with Enumerated

    //! @brief default constructor
    MultiplyDouble() :
      m_Coefficient( 1.0)
    {
    }

    //! @brief Clone function
    //! @return pointer to new class_name
    MultiplyDouble *Clone() const
    {
      return new MultiplyDouble( *this);
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

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &GetAlias() const
    {
      static const std::string s_Alias( "Multiply");
      return s_Alias;
    }

  ///////////////
  // operators //
  ///////////////

    double operator()( const double &VAL) const
    {
      return m_Coefficient * VAL;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Multiplies a double by a coefficient");

      parameters.AddInitializer
      (
        "coefficient",
        "value to multiply other values by",
        io::Serialization::GetAgent( &m_Coefficient),
        "1.0"
      );

      return parameters;
    }
  };

  const util::SiPtr< const util::ObjectInterface> MultiplyDouble::s_Instance( util::Enumerated< DoubleFuncInterface>::AddInstance( new MultiplyDouble()));

  // this class holds a vector of coefficients, and returns their sum times the passed in number
  class VectorScalarProduct :
    public DoubleFuncInterface
  {
  private:

    linal::Vector< double> m_Coefficients;

  public:

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    static const util::SiPtr< const util::ObjectInterface> s_Instance; //!< True if/when this class has been registered with Enumerated

    //! @brief Clone function
    //! @return pointer to new class_name
    VectorScalarProduct *Clone() const
    {
      return new VectorScalarProduct( *this);
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

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &GetAlias() const
    {
      static const std::string s_Alias( "VectorScalarProduct");
      return s_Alias;
    }

  ///////////////
  // operators //
  ///////////////

    double operator()( const double &VAL) const
    {
      return m_Coefficients.Sum() * VAL;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Multiplies a double by a vector of #'s, adding them up");

      parameters.AddInitializer
      (
        "",
        "values to multiply other values by",
        io::Serialization::GetAgent( &m_Coefficients)
      );

      return parameters;
    }

  };

  const util::SiPtr< const util::ObjectInterface> VectorScalarProduct::s_Instance( util::Enumerated< DoubleFuncInterface>::AddInstance( new VectorScalarProduct()));

  class ExampleUtilImplementation :
    public ExampleInterface
  {
  public:

    ExampleUtilImplementation *Clone() const
    {
      return new ExampleUtilImplementation( *this);
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
      util::Implementation< DoubleFuncInterface> default_obj_from_data_label;
      util::Implementation< DoubleFuncInterface> square_double( "Square");
      util::Implementation< DoubleFuncInterface> multiply_double_identity( "Multiply");
      util::Implementation< DoubleFuncInterface> multiply_double_by_3( "Multiply(coefficient=3)");
      util::Implementation< DoubleFuncInterface> scalar_product_5_9( "VectorScalarProduct(9.0,5)");

      // check that the default object behaves as expected
      BCL_ExampleCheck( default_obj_from_data_label.IsDefined(), false);
      BCL_ExampleCheck( default_obj_from_data_label.GetAlias(), GetStaticClassName< DoubleFuncInterface>());
      BCL_ExampleCheck( default_obj_from_data_label.GetString(), GetStaticClassName< DoubleFuncInterface>());

      // ensure that the other object from data labels were defined; if not, there is nothing more to test
      BCL_ExampleAssert( square_double.IsDefined(), true);
      BCL_ExampleAssert( multiply_double_identity.IsDefined(), true);
      BCL_ExampleAssert( multiply_double_by_3.IsDefined(), true);
      BCL_ExampleAssert( scalar_product_5_9.IsDefined(), true);

      // check the aliases
      BCL_ExampleCheck( square_double.GetAlias(), "Square");
      BCL_ExampleCheck( multiply_double_identity.GetAlias(), "Multiply");
      BCL_ExampleCheck( multiply_double_by_3.GetAlias(), "Multiply");
      BCL_ExampleCheck( scalar_product_5_9.GetAlias(), "VectorScalarProduct");

      // check the data strings
      BCL_ExampleCheck( square_double.GetString(), "Square"); // no arguments -> no parenthesis
      BCL_ExampleCheck( multiply_double_identity.GetString(), "Multiply(coefficient=1)");
      BCL_ExampleCheck( multiply_double_by_3.GetString(), "Multiply(coefficient=3)");
      BCL_ExampleCheck( multiply_double_by_3->GetString(), "3");
      BCL_ExampleCheck( scalar_product_5_9.GetString(), "VectorScalarProduct(9,5)");
      BCL_ExampleCheck( scalar_product_5_9->GetString(), "(9,5)");

      // check that the operators work as expected
      BCL_ExampleCheck( ( *square_double)( 5.0), 25.0);
      BCL_ExampleCheck( ( *multiply_double_identity)( 5.0), 5.0 * 1.0);
      BCL_ExampleCheck( multiply_double_by_3->operator()( 5.0), 5.0 * 3.0);
      BCL_ExampleCheck( scalar_product_5_9->operator()( 3.0), 3.0 * ( 5.0 + 9.0));

      // write out help for these implementations
      if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Standard))
      {
        io::FixedLineWidthWriter writer( 0, 120);
        writer << "This is the help for the double function interface:\n";
        util::Implementation< DoubleFuncInterface>().WriteInstancesHelp( writer);
        util::GetLogger() << writer.String();
      }

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleUtilImplementation

  const ExampleClass::EnumType ExampleUtilImplementation::s_Instance
  (
    GetExamples().AddEnum( ExampleUtilImplementation())
  );

} // namespace bcl
