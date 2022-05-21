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

#ifndef BCL_MATH_SUM_FUNCTION_INTERFACE_H_
#define BCL_MATH_SUM_FUNCTION_INTERFACE_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_math_function_interface_serializable.h"
#include "io/bcl_io_serialization.h"
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SumFunctionInterface
    //! @brief This class is a template for t_ResultType = Function( t_ArgumentType&) allowing addition
    //! of functions with coefficients.
    //! @details Assuming y being element of t_ResultType, x being element of t_ArgumentType and s being scalar
    //! this class can represent the following functions:
    //! y1 = s1 * g1(x)
    //! y2 = s2 * g2(x)
    //!
    //! f(x) = y0 + y1 + y2 = y0 + s1 * g1(x) + s2 * g2(x)
    //!
    //! operators provide the functionality to divide the sumfunction, add an absolute to y0, add coefficients
    //! s3 * g3(x), add another sumfunction and to create sumfunctions from function interfaces and coefficients.
    //!
    //! @remarks example unnecessary
    //! @author meilerj, woetzen
    //! @date 27.08.2004
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Interface>
    class SumFunctionInterface :
      public t_Interface
      // public FunctionInterfaceSerializable< typename t_BaseClass::argument_type, typename t_BaseClass::return_type>
    {
    public:

      typedef typename t_Interface::argument_type t_ArgumentType;
      typedef typename t_Interface::return_type   t_ResultType;

    ///////////
    // types //
    ///////////

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class Term
      //! @brief storage for each weight and interface
      //! @remarks example unnecessary
      //! @author mendenjl
      //! @date Sep 17, 2016
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      class Term :
        public storage::Pair< double, util::Implementation< t_Interface> >,
        public util::SerializableInterface
      {
      private:

        typedef storage::Pair< double, util::Implementation< t_Interface> > Base;

      public:

        //! @brief default constructor
        Term() :
          Base()
        {
        }

        //! @brief construct Pair from two values
        //! @param FIRST the object which will be m_First
        //! @param SECOND the object which will be m_Second
        Term
        (
          const double &FIRST,
          const util::Implementation< t_Interface> &SECOND
        ) :
          Base( FIRST, SECOND)
        {
        }

        //! @brief construct Pair from std::pair
        //! @param PAIR which will be the Pair
        Term( const Base &PAIR) :
          Base( PAIR)
        {
        }

        //! copy constructor
        Term( const Term &PAIR) :
          Base( PAIR)
        {
        }

        //! virtual copy constructor (Clone)
        Term *Clone() const
        {
          return new Term( *this);
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

        //! @brief alias
        //! @return an alias for the class
        const std::string &GetAlias() const
        {
          static const std::string s_empty;
          return s_empty;
        }

        //! @brief alias
        //! @return an alias for the class
        const std::string &GetScheme() const
        {
          static const std::string s_empty;
          return s_empty;
        }

        //! @brief return parameters for member data that are set up from the labels
        //! @return parameters for member data that are set up from the labels
        io::Serializer GetSerializer() const
        {
          io::Serializer parameters;
          parameters.SetClassDescription( "weight and term");
          parameters.AddInitializer
          (
            "weight",
            "weight for the term",
            io::Serialization::GetAgent( &this->First())
          );
          parameters.AddInitializer
          (
            "",
            "the term of interest",
            io::Serialization::GetAgent( &this->Second())
          );
          return parameters;
        }

      };

    //////////
    // data //
    //////////

      //! @brief return the value table vertical column names
      static const storage::Vector< std::string> &GetValueTableVerticalColumnNames()
      {
        // initialize column names
        static const storage::Vector< std::string> s_column_names
        (
          storage::Vector< std::string>::Create( "weight", "value", "weighted_value")
        );

        // return
        return s_column_names;
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief return the absolute value
      //! @return the absolute value
      virtual const t_ResultType &GetAbsolute() const = 0;

      //! @brief access to term with SCHEME
      //! @param SCHEME the scheme of the function in the term to access
      //! @return reference to term - if invalid scheme, term has undefined scalar and function pointer
      virtual const Term &GetTerm( const std::string &SCHEME) const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief set the coefficient for a function with the scheme
      //! @param SCHEME the scheme of the function
      //! @param COEFFICIENT s like in s * g(x)
      //! @return true, if successful, false id there is no such function with that scheme
      virtual bool SetCoefficient( const std::string &SCHEME, const double &COEFFICIENT) = 0;

      //! @brief add an new operand with an optional coefficient to the object s1 * g1(x)
      //! @param SP_FUNCTION function g(x) as Implementation to Function
      //! @param COEFFICIENT s like in s1 * g1(x) - default is 1
      virtual void NewOperand
      (
        const util::Implementation< t_Interface> &SP_FUNCTION,
        const double &COEFFICIENT
      ) = 0;

      //! @brief returns a vector of string that has schemes of the individual functions
      //! @return the vector of strings that has schemes of the individual functions
      virtual storage::Vector< std::string> GetFunctionSchemes() const = 0;

      //! @brief creates a table with individual function values, weights and weighted values
      //! the table has scores as the columns
      //! @param ARGUMENT Argument to be used for calculating the function value
      //! @return a table with individual functions and their weighted sums
      virtual storage::Table< t_ResultType> CreateValueTableHorizontal( const t_ArgumentType &ARGUMENT) const = 0;

      //! @brief creates a table with individual function values and their weighted sums for the given argument
      //! the table has scores as rows
      //! @param ARGUMENT Argument to be used for calculating the function value
      //! @return a table with individual functions and their weighted sums
      virtual storage::Table< t_ResultType> CreateValueTableVertical( const t_ArgumentType &ARGUMENT) const = 0;

    }; // template class SumFunctionInterface

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_SUM_FUNCTION_INTERFACE_H_
