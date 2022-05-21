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

#ifndef BCL_MATH_SUM_FUNCTION_MIXIN_H_
#define BCL_MATH_SUM_FUNCTION_MIXIN_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_sum_function_interface.h"
#include "io/bcl_io_serialization.h"
#include "storage/bcl_storage_table.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically

// Uncomment the next line to output timings for each score
// This can get excessive for inner-loop optimizations, so it is not on by default
//#define PROFILE_SUMMED_FUNCTIONS

namespace bcl
{
  namespace math
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SumFunctionMixin
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
    //! @see @link example_math_sum_function_mixin.cpp @endlink
    //! @author mendenjl
    //! @date Sep 19, 2016
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Interface>
    class SumFunctionMixin :
      public SumFunctionInterface< t_Interface>
    {
    public:

      //! @typedefs
      typedef typename SumFunctionInterface< t_Interface>::t_ArgumentType t_ArgumentType;
      typedef typename SumFunctionInterface< t_Interface>::t_ResultType   t_ResultType;
      typedef typename SumFunctionInterface< t_Interface>::Term Term;

    protected:

    //////////
    // data //
    //////////

      //! absolute value y0
      t_ResultType m_Absolute;

      //! @brief vector of pairs of coefficients and pointer on functions - representing s1 * g1(x) + s2 * g2(x)
      //! where s * g(x) is the Pair, and s the scalar and g(x) the ShPtr< FunctionInterface>
      storage::Vector< Term> m_Functions;

#ifdef PROFILE_SUMMED_FUNCTIONS
      //! @brief vector of stopwatches
      mutable storage::Vector< util::Stopwatch> m_Stopwatches;
#endif

      typedef typename storage::Vector< Term>::const_iterator const_iterator;
      typedef typename storage::Vector< Term>::iterator iterator;

      //! scheme
      std::string m_Scheme;

    public:

    //////////
    // data //
    //////////

      //! @brief return the value table vertical column names
      static const storage::Vector< std::string> &GetValueTableVerticalColumnNames()
      {
        // initialize column names
        static const storage::Vector< std::string> s_column_names
        (
          storage::Vector< std::string>::Create( "weight", "value", "weighted_val")
        );

        // return
        return s_column_names;
      }

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SumFunctionMixin( const std::string &SCHEME = "") :
        m_Absolute(),
        m_Functions(),
        m_Scheme( SCHEME)
      {
      }

      //! @brief construct from ABS_CONST
      //! @param ABS_CONST y0 for the sumfunction
      explicit SumFunctionMixin( const t_ResultType &ABS_CONST) :
        m_Absolute( ABS_CONST),
        m_Functions(),
        m_Scheme()
      {
      }

      //! @brief construct from ShPtr on Function, optional coefficient, and absolute value
      //! y = y0 + s1*g1(x)
      //! @param SP_FUNCTION g1(x) as part of sumfunction as ShPtr to FunctionInterfaceSerializable
      //! @param COEFFICIENT s1(x) as coefficient for g1(x) -> s1 * g1(x) - default 1
      //! @param ABS_CONST    y0 as absolute in sumfunction - default 0
      SumFunctionMixin
      (
        const util::Implementation< t_Interface> &SP_FUNCTION,
        const double &COEFFICIENT,
        const t_ResultType &ABS_CONST
      ) :
        m_Absolute( ABS_CONST),
        m_Functions(),
        m_Scheme()
      {
        NewOperand( SP_FUNCTION, COEFFICIENT);
      }

      //! @brief construct from Function, optional coefficient, and optional absolute value
      //! y = y0 + s1*g1(x)
      //! @param FUNCTION    g1(x) as part of sumfunction as reference to FunctionInterfaceSerializable
      //! @param COEFFICIENT s1(x) as coefficient for g1(x) -> s1 * g1(x) - default 1
      //! @param ABS_CONST    y0 as absolute in sumfunction - default 0
      SumFunctionMixin
      (
        const t_Interface &FUNCTION,
        const double &COEFFICIENT,
        const t_ResultType &ABS_CONST
      ) :
        m_Absolute( ABS_CONST),
        m_Functions(),
        m_Scheme()
      {
        NewOperand( util::Implementation< t_Interface>( FUNCTION.Clone()), COEFFICIENT);
      }

      //! @brief construct from two Functions, optional coefficients, and optional absolute value
      //! y = y0 + s1*g1(x) + s2*g2(x)
      //! @param SP_FUNCTION_A g1(x) as part of sumfunction as ShPtr to FunctionInterfaceSerializable
      //! @param SP_FUNCTION_B g2(x) as part of sumfunction as ShPtr to FunctionInterfaceSerializable
      //! @param COEFFICIENT_A s1(x) as coefficient for g1(x) -> s1 * g1(x) - default 1
      //! @param COEFFICIENT_B s2(x) as coefficient for g2(x) -> s2 * g2(x) - default 1
      //! @param ABS_CONST    y0 as absolute in sumfunction - default 0
      SumFunctionMixin
      (
        const util::Implementation< t_Interface> &SP_FUNCTION_A,
        const util::Implementation< t_Interface> &SP_FUNCTION_B,
        const double &COEFFICIENT_A,
        const double &COEFFICIENT_B,
        const t_ResultType &ABS_CONST
      ) :
        m_Absolute( ABS_CONST),
        m_Functions(),
        m_Scheme()
      {
        NewOperand( SP_FUNCTION_A, COEFFICIENT_A);
        NewOperand( SP_FUNCTION_B, COEFFICIENT_B);
      }

      //! @brief construct from two Functions, optional coefficients, and optional absolute value
      //! y = y0 + s1*g1(x) + s2*g2(x)
      //! @param FUNCTION_A    g1(x) as part of sumfunction as reference to FunctionInterfaceSerializable
      //! @param FUNCTION_B    g2(x) as part of sumfunction as reference to FunctionInterfaceSerializable
      //! @param COEFFICIENT_A s1(x) as coefficient for g1(x) -> s1 * g1(x) - default 1
      //! @param COEFFICIENT_B s2(x) as coefficient for g2(x) -> s2 * g2(x) - default 1
      //! @param ABS_CONST    y0 as absolute in sumfunction - default 0
      SumFunctionMixin
      (
        const t_Interface &FUNCTION_A,
        const t_Interface &FUNCTION_B,
        const double &COEFFICIENT_A,
        const double &COEFFICIENT_B,
        const t_ResultType &ABS_CONST
      ) :
        m_Absolute( ABS_CONST),
        m_Functions(),
        m_Scheme()
      {
        NewOperand( util::Implementation< t_Interface>( FUNCTION_A.Clone()), COEFFICIENT_A);
        NewOperand( util::Implementation< t_Interface>( FUNCTION_B.Clone()), COEFFICIENT_B);
      }

      //! @brief copy constructor
      //! @return pointer to new SumFunctionMixin which is a copy of this
      SumFunctionMixin *Clone() const
      {
        return new SumFunctionMixin( *this);
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

      //! @brief return the absolute value
      //! @return the absolute value
      const t_ResultType &GetAbsolute() const
      {
        return m_Absolute;
      }

      //! @brief returns the vector of functions and the associated weights
      //! @return the vector of functions and the associated weights
      const storage::Vector< Term> &
      GetFunction() const
      {
        return m_Functions;
      }

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetAlias() const
      {
        static const std::string s_name( "WeightedSum");
        return s_name;
      }

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

      //! @brief access to term with SCHEME
      //! @param SCHEME the scheme of the function in the term to access
      //! @return reference to term - if invalid scheme, term has undefined scalar and function pointer
      const Term &GetTerm( const std::string &SCHEME) const
      {
        // undefined term
        static const Term s_undefined
        (
          util::GetUndefined< double>(),
          util::Implementation< t_Interface>()
        );

        // iterate over terms
        for( const_iterator itr( m_Functions.Begin()), itr_end( m_Functions.End()); itr != itr_end; ++itr)
        {
          if( itr->Second()->GetScheme() == SCHEME)
          {
            return *itr;
          }
        }

        // not function with SCHEME found
        return s_undefined;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief assign Function
      //! @param FUNCTION g1( x) with s1 = 1
      //! @return reference to this SumFunctionMixin
      SumFunctionMixin &operator =( const t_Interface &FUNCTION)
      {
        // reset all funtions in sum
        m_Functions.Reset();

        // add this function as new operand
        NewOperand( FUNCTION, t_ResultType( 1));

        // reset m_Absoulte to 0
        m_Absolute = t_ResultType( 0);

        // return this
        return *this;
      }

      //! @brief assignment operator
      //! @param SUM_FUNCTION sum function this should be assigned to
      //! @return reference to this SumFunctionMixin
      SumFunctionMixin &operator =( const SumFunctionMixin &SUM_FUNCTION)
      {
        // assign data member
        m_Absolute = SUM_FUNCTION.m_Absolute;
        m_Functions = SUM_FUNCTION.m_Functions;
        #ifdef PROFILE_SUMMED_FUNCTIONS
        m_Stopwatches = SUM_FUNCTION.m_Stopwatches;
        #endif

        // return this
        return *this;
      }

      //! @brief actual f(x) operator taking ARGUMENT and returning the sum of all function values
      //! @brief ARGUMENT x to f(x)
      //! @return the value of f(x) = sum
      t_ResultType operator()( const t_ArgumentType &ARGUMENT) const;

    ////////////////
    // operations //
    ////////////////

      //! @brief set the coefficient for a function with the scheme
      //! @param SCHEME the scheme of the function
      //! @param COEFFICIENT s like in s * g(x)
      //! @return true, if successful, false id there is no such function with that scheme
      bool SetCoefficient( const std::string &SCHEME, const double &COEFFICIENT)
      {
        // iterate over terms
        for( iterator itr( m_Functions.Begin()), itr_end( m_Functions.End()); itr != itr_end; ++itr)
        {
          if( itr->Second()->GetScheme() == SCHEME)
          {
            itr->First() = COEFFICIENT;
            return true;
          }
        }

        // not function with SCHEME found
        return false;
      }

      //! @brief add an new operand with an optional coefficient to the object s1 * g1(x)
      //! @param SP_FUNCTION function g(x) as ShPtr to Function
      //! @param COEFFICIENT s like in s1 * g1(x) - default is 1
      void NewOperand
      (
        const util::Implementation< t_Interface> &SP_FUNCTION,
        const double &COEFFICIENT = double( 1)
      )
      {
        this->NewOperand( *SP_FUNCTION, COEFFICIENT);
      }

      //! @brief add an new operand with an optional coefficient to the object s1 * g1(x)
      //! @param SP_FUNCTION function g(x) as ShPtr to Function
      //! @param COEFFICIENT s like in s1 * g1(x) - default is 1
      void NewOperand
      (
        const t_Interface &FUNCTION,
        const double &COEFFICIENT = double( 1)
      )
      {
        // insert pair of function and coefficient
        m_Functions.PushBack( Term( COEFFICIENT, FUNCTION));

        // update the identifier
        m_Scheme += FUNCTION.GetScheme();

        #ifdef PROFILE_SUMMED_FUNCTIONS
        m_Stopwatches.PushBack
        (
          util::Stopwatch
          (
            "Time to compute " + FUNCTION.GetScheme(),
            util::Time( 10000, 0.0),
            util::Message::e_Standard,
            true,
            false
          )
        );
        #endif
      }

      //! @brief returns a vector of string that has schemes of the individual functions
      //! @return the vector of strings that has schemes of the individual functions
      storage::Vector< std::string> GetFunctionSchemes() const;

      //! @brief creates a table with individual function values, weights and weighted values
      //! the table has scores as the columns
      //! @param ARGUMENT Argument to be used for calculating the function value
      //! @return a table with individual functions and their weighted sums
      storage::Table< t_ResultType> CreateValueTableHorizontal( const t_ArgumentType &ARGUMENT) const;

      //! @brief creates a table with individual function values and their weighted sums for the given argument
      //! the table has scores as rows
      //! @param ARGUMENT Argument to be used for calculating the function value
      //! @return a table with individual functions and their weighted sums
      storage::Table< t_ResultType> CreateValueTableVertical( const t_ArgumentType &ARGUMENT) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator += Function
      //! @param FUNCTION function g(x) to be added to the sumfunction as reference - with coefficient 1
      //! @return reference to this altered sum function
      SumFunctionMixin operator -() const
      {
        return SumFunctionMixin( *this).operator *=( -1.0);
      }

      //! @brief operator += Function
      //! @param FUNCTION function g(x) to be added to the sumfunction as reference - with coefficient 1
      //! @return reference to this altered sum function
      SumFunctionMixin &operator +=( const t_Interface &FUNCTION)
      {
        NewOperand( util::Implementation< t_Interface>( FUNCTION.Clone()));
        return *this;
      }

      //! @brief operator -= Function
      //! @param FUNCTION function g(x) to be subtracted from the sumfunction as reference - with coefficient 1
      //! @return reference to this altered sum function
      SumFunctionMixin &operator -=( const t_Interface &FUNCTION)
      {
        NewOperand( util::Implementation< t_Interface>( FUNCTION.Clone()), t_ResultType( -1));

        return *this;
      }

      //! @brief operator += SumFunctionMixin
      //! @param SUM_FUNCTION sum function to be added to the sum function as reference
      //! @return reference to this altered sum function
      SumFunctionMixin &
      operator +=
      (
        const SumFunctionMixin &SUM_FUNCTION
      )
      {
        // append functions in SUM_FUNCTION to this sum function
        m_Functions.Append( SUM_FUNCTION.m_Functions);

        // add the absolute
        m_Absolute += SUM_FUNCTION.m_Absolute;

        #ifdef PROFILE_SUMMED_FUNCTIONS
        // add stopwatches, if profiling
        m_Stopwatches.Append( SUM_FUNCTION.m_Stopwatches);
        #endif

        // return
        return *this;
      }

      //! @brief operator -= SumFunctionMixin
      //! @param SUM_FUNCTION sum function to be subtracted from the sum function as reference
      //! @return reference to this altered sum function
      SumFunctionMixin &
      operator -=
      (
        const SumFunctionMixin &SUM_FUNCTION
      )
      {
        return operator +=( -SUM_FUNCTION);
      }

      //! @brief operator += VALUE adds y to y0
      //! @param VALUE y to be added to y0
      SumFunctionMixin &
      operator +=
      (
        const t_ResultType &VALUE
      )
      {
        // add to absolute
        m_Absolute += VALUE;

        // return
        return *this;
      }

      //! @brief operator -= VALUE subtract y from y0
      //! @param VALUE y to be subtracted from y0
      SumFunctionMixin &
      operator -=
      (
        const t_ResultType &VALUE
      )
      {
        // subtract from absolute
        m_Absolute -= VALUE;

        // return
        return *this;
      }

      //! @brief operator *= SCALAR
      //! mutlipy absolute and each coefficient with the given scalar: s*f(x) = s*y0 + s*s1*g1(x) + s*s2*g2(x)...
      //! @param SCALAR scalar to be used to multiply every term
      //! @return reference to altered sum function
      SumFunctionMixin &
      operator *=
      (
        const double &SCALAR
      )
      {
        //multiply the absolute value
        m_Absolute *= SCALAR;

        // multiply each coefficient by the given SCALAR
        //iterate over all functions, calculate their function value and add them multiplied with their coefficient
        for
        (
          iterator pair_itr( m_Functions.Begin()), pair_itr_end( m_Functions.End());
          pair_itr != pair_itr_end;
          ++pair_itr
        )
        {
          pair_itr->First() *= SCALAR;
        }

        // return this
        return *this;
      }

      //! @brief operator /= SCALAR
      //! devide absolute and each coefficient with the given scalar: f(x)/s = y0/s + s1/s*g1(x) + s2/s*g2(x)...
      //! @param SCALAR scalar to be used to devide every term
      //! @return reference to altered sum function
      SumFunctionMixin &
      operator /=
      (
        const double &SCALAR
      )
      {
        // check that scalr is not 0
        BCL_Assert( SCALAR != double( 0), "attempt division by zero!");

        // call *= with 1/s
        return operator *=( double( 1) / SCALAR);
      }

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write detailed scheme and values to OSTREAM
      //! @param ARGUMENT Argument to be used to evaluate the function
      //! @param OSTREAM std::ostream to be written to
      //! @return std::ostream which was written to
      std::ostream &WriteDetailedSchemeAndValues
      (
        const t_ArgumentType &ARGUMENT,
        std::ostream &OSTREAM
      ) const;

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer parameters;
        parameters.SetClassDescription( "Weighted sum of " + GetStaticClassName< t_ArgumentType>());
        parameters.AddInitializer
        (
          "offset",
          "constant to add to the weighted sum",
          io::Serialization::GetAgent( &m_Absolute),
          "0"
        );
        parameters.AddInitializer
        (
          "terms",
          "weights and functions to sum",
          io::Serialization::GetAgent( &m_Functions)
        );
        return parameters;
      }

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
      {
        #ifdef PROFILE_SUMMED_FUNCTIONS
        // add stopwatches, if profiling
        m_Stopwatches.Reset();
        for
        (
          typename storage::Vector< Term>::const_iterator itr( m_Functions.Begin()), itr_end( m_Functions.End());
          itr != itr_end;
          ++itr
        )
        {
          m_Stopwatches.PushBack
          (
            util::Stopwatch
            (
              "Time to compute " + itr->GetScheme(),
              util::Time( 10000, 0.0),
              util::Message::e_Standard,
              true,
              false
            )
          );
        }
        #endif
        return true;
      }

    }; // template class SumFunctionMixin

  //////////
  // data //
  //////////

    // instantiate s_Instance
    template< typename t_Interface>
    const util::SiPtr< const util::ObjectInterface> SumFunctionMixin< t_Interface>::s_Instance
    (
      util::Enumerated< t_Interface>::AddInstance( new SumFunctionMixin< t_Interface>())
    );

  /////////////////
  // data access //
  /////////////////

    //! @brief implement operator()
    template< typename t_Interface>
    inline
    typename SumFunctionMixin< t_Interface>::t_ResultType
    SumFunctionMixin< t_Interface>::operator()( const t_ArgumentType &ARGUMENT) const
    {
//      static util::Stopwatch timer;

      // initialize sum with the absolute
      t_ResultType sum( m_Absolute);

      #ifdef PROFILE_SUMMED_FUNCTIONS
      storage::Vector< util::Stopwatch>::iterator itr_stop( m_Stopwatches.Begin());
      #endif
      // iterate over all functions, calculate their function value and add them multiplied with their coefficient
      for
      (
        const_iterator pair_itr( m_Functions.Begin()), pair_itr_end( m_Functions.End());
        pair_itr != pair_itr_end;
        #ifdef PROFILE_SUMMED_FUNCTIONS
        ++pair_itr, ++itr_stop
        #else
        ++pair_itr
        #endif
      )
      {
//        timer.Reset();
//        timer.Start();
          #ifdef PROFILE_SUMMED_FUNCTIONS
          itr_stop->Start();
          #endif
          if( pair_itr->First())
          {
            sum += pair_itr->First() * pair_itr->Second()->operator()( ARGUMENT);
          }
          #ifdef PROFILE_SUMMED_FUNCTIONS
          itr_stop->Stop();
          #endif
//        const util::Time this_time( timer.GetProcessDuration());
//        std::cout << pair_itr->Second()->GetScheme() << '\t' << this_time.GetTotalMicroseconds() << '\t';
//        timer.Reset();
      }

      // return sum f(x) = y
      return sum;
    }

    //! @brief returns a vector of string that has schemes of the individual functions
    //! @return the vector of strings that has schemes of the individual functions
    template< typename t_Interface>
    inline
    storage::Vector< std::string>
    SumFunctionMixin< t_Interface>::GetFunctionSchemes() const
    {
      // initialize vector to store the function names
      storage::Vector< std::string> function_names;

      // iterate again this time to collect the scores
      for
      (
        const_iterator function_itr( m_Functions.Begin()), function_itr_end( m_Functions.End());
        function_itr != function_itr_end;
        ++function_itr
      )
      {
        // insert the name for this function
        function_names.PushBack( function_itr->Second()->GetScheme());
      }

      // push sum column
      function_names.PushBack( "sum");

      // end
      return function_names;
    }

    //! @brief creates a table with individual function values, weights and weighted values
    //! the table has scores as the columns
    //! @param ARGUMENT Argument to be used for calculating the function value
    //! @return a table with individual functions and their weighted sums
    template< typename t_Interface>
    inline
    storage::Table< typename SumFunctionMixin< t_Interface>::t_ResultType>
    SumFunctionMixin< t_Interface>::CreateValueTableHorizontal( const t_ArgumentType &ARGUMENT) const
    {
      // initialize vectors to store values and allocate memory
      storage::Vector< t_ResultType> weights_vector, values_vector, weighted_values_vector;
      weights_vector.AllocateMemory( m_Functions.GetSize() + 1);
      values_vector.AllocateMemory( m_Functions.GetSize() + 1);
      weighted_values_vector.AllocateMemory( m_Functions.GetSize() + 1);

      // initialize the sums
      double sum_weights( 0.0), sum_values( 0.0), sum_weighted_values( 0.0);

      // initialize vector to store the function names
      storage::Vector< std::string> function_names;

      // iterate again this time to collect the scores
      for
      (
        const_iterator function_itr( m_Functions.Begin()), function_itr_end( m_Functions.End());
        function_itr != function_itr_end;
        ++function_itr
      )
      {
        // insert the name for this function
        function_names.PushBack( function_itr->Second()->GetScheme());

        // calculate the value for the argument using this function
        const t_ResultType this_weight( function_itr->First());
        const t_ResultType this_value( function_itr->Second()->operator()( ARGUMENT));
        const t_ResultType this_weighted_value( this_value * function_itr->First());

        // record the value, weight and the weighted value
        weights_vector.PushBack( this_weight);
        values_vector.PushBack( this_value);
        weighted_values_vector.PushBack( this_weighted_value);

        // update the sums
        sum_weights += this_weight;
        sum_values += this_value;
        sum_weighted_values += this_weighted_value;
      }
      // add the sum column
      function_names.PushBack( "sum");

      // add the summed values
      weights_vector.PushBack( sum_weights);
      values_vector.PushBack( sum_values);
      weighted_values_vector.PushBack( sum_weighted_values);

      // Insert the rows
      util::ShPtr< storage::TableHeader> sp_table_header( new storage::TableHeader( function_names));
      storage::Table< t_ResultType> function_table( sp_table_header);
      function_table.InsertRow( "weights", weights_vector);
      function_table.InsertRow( "value", values_vector);
      function_table.InsertRow( "weighted_val", weighted_values_vector);

      // end
      return function_table;
    }

    //! @brief creates a table with individual function values and their weighted sums for the given argument
    //! the table has scores as rows
    //! @param ARGUMENT Argument to be used for calculating the function value
    //! @return a table with individual functions and their weighted sums
    template< typename t_Interface>
    inline
    storage::Table< typename SumFunctionMixin< t_Interface>::t_ResultType>
    SumFunctionMixin< t_Interface>::CreateValueTableVertical( const t_ArgumentType &ARGUMENT) const
    {
      // initialize table header and a table
      util::ShPtr< storage::TableHeader> sp_table_header( new storage::TableHeader( GetValueTableVerticalColumnNames()));
      storage::Table< t_ResultType> function_table( sp_table_header);

      // initialize the sums
      storage::Vector< t_ResultType> sums_vector( 3, t_ResultType( 0));

      // iterate again this time to collect the scores
      for
      (
        const_iterator function_itr( m_Functions.Begin()), function_itr_end( m_Functions.End());
        function_itr != function_itr_end;
        ++function_itr
      )
      {
        // calculate the value for the argument using this function
        const t_ResultType this_weight( function_itr->First());
        const t_ResultType this_value( function_itr->Second()->operator()( ARGUMENT));
        const t_ResultType this_weighted_value( this_value * function_itr->First());

        // create and insert the row
        function_table.InsertRow
        (
          function_itr->Second()->GetScheme(),
          storage::Vector< t_ResultType>::Create( this_weight, this_value, this_weighted_value)
        );

        // update the sums
        sums_vector( 0) += this_weight;
        sums_vector( 1) += this_value;
        sums_vector( 2) += this_weighted_value;
      }
      // add the sum row
      function_table.InsertRow( "sum", sums_vector);

      // end
      return function_table;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write detailed scheme and values to OSTREAM
    //! @param ARGUMENT Argument to be used to evaluate the function
    //! @param OSTREAM std::ostream to be written to
    //! @return std::ostream which was written to
    template< typename t_Interface>
    inline
    std::ostream &
    SumFunctionMixin< t_Interface>::WriteDetailedSchemeAndValues
    (
      const t_ArgumentType &ARGUMENT,
      std::ostream &OSTREAM
    ) const
    {
      for
      (
        const_iterator pair_itr( m_Functions.Begin()), pair_itr_end( m_Functions.End());
        pair_itr != pair_itr_end;
        ++pair_itr
      )
      {
        OSTREAM << pair_itr->Second()->GetScheme();
        OSTREAM << '\t' << ( pair_itr->First() * pair_itr->Second()->operator()( ARGUMENT))
                << "                " << '\n';

        pair_itr->Second()->WriteDetailedSchemeAndValues( ARGUMENT, OSTREAM) << '\n';
      }

      //end
      return OSTREAM;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator SumFunctionMixin + SumFunctionMixin
    //! @param FUNCTION_A ya = y0a + s1a*g1a(x) as reference to a sum function
    //! @param FUNCTION_B yb = y0b + s1b*g1b(x) as reference to a sum function
    //! @return new SumFunctionMixin y = ya + yb = y0a + y0b + s1a*g1a(x) + s1b*g1b(x)
    template< typename t_Interface>
    inline SumFunctionMixin< t_Interface> operator +
    (
      const SumFunctionMixin< t_Interface> &FUNCTION_A,
      const SumFunctionMixin< t_Interface> &FUNCTION_B
    )
    {
      return SumFunctionMixin< t_Interface>( FUNCTION_A).operator +=( FUNCTION_B);
    }

    //! @brief operator SumFunctionMixin - SumFunctionMixin
    //! @param FUNCTION_A ya = y0a + s1a*g1a(x) as reference to a sum function
    //! @param FUNCTION_B yb = y0b + s1b*g1b(x) as reference to a sum function
    //! @return new SumFunctionMixin y = ya - yb = y0a - y0b + s1a*g1a(x) - s1b*g1b(x)
    template< typename t_Interface>
    inline SumFunctionMixin< t_Interface> operator -
    (
      const SumFunctionMixin< t_Interface> &FUNCTION_A,
      const SumFunctionMixin< t_Interface> &FUNCTION_B
    )
    {
      return SumFunctionMixin< t_Interface>( FUNCTION_A).operator -=( FUNCTION_B);
    }

    //! @brief operator FunctionInterfaceSerializable + SumFunctionMixin
    //! @param FUNCTION_A     ya = g1a(x) as reference to a function
    //! @param SUM_FUNCTION_B yb = y0b + s1b*g1b(x) as reference to a sum function
    //! @return new SumFunctionMixin y = ya + yb = y0b + g1a(x) + s1b*g1b(x)
    template< typename t_Interface>
    inline SumFunctionMixin< t_Interface> operator +
    (
      const t_Interface     &FUNCTION_A,
      const SumFunctionMixin< t_Interface> &SUM_FUNCTION_B
    )
    {
      return SumFunctionMixin< t_Interface>( SUM_FUNCTION_B).operator +=( FUNCTION_A);
    }

    //! @brief operator FunctionInterfaceSerializable - SumFunctionMixin
    //! @param FUNCTION_A     ya = g1a(x) as reference to a function
    //! @param SUM_FUNCTION_B yb = y0b - s1b*g1b(x) as reference to a sum function
    //! @return new SumFunctionMixin y = ya - yb = - y0b + g1a(x) - s1b*g1b(x)
    template< typename t_Interface>
    inline SumFunctionMixin< t_Interface> operator -
    (
      const t_Interface &FUNCTION_A,
      const SumFunctionMixin< t_Interface> &SUM_FUNCTION_B
    )
    {
      return SumFunctionMixin< t_Interface>( -SUM_FUNCTION_B).operator +=( FUNCTION_A);
    }

    //! @brief operator SumFunctionMixin + FunctionInterfaceSerializable
    //! @param SUM_FUNCTION_A ya = y0a + s1a*g1a(x) as reference to a sum function
    //! @param FUNCTION_B     yb = g1b(x) as reference to a function
    //! @return new SumFunctionMixin y = ya + yb = y0a + s1a*g1a(x) + g1b(x)
    template< typename t_Interface>
    inline SumFunctionMixin< t_Interface> operator +
    (
      const SumFunctionMixin< t_Interface> &SUM_FUNCTION_A,
      const t_Interface &FUNCTION_B
    )
    {
      return SumFunctionMixin< t_Interface>( SUM_FUNCTION_A).operator +=( FUNCTION_B);
    }

    //! @brief operator SumFunctionMixin - FunctionInterfaceSerializable
    //! @param SUM_FUNCTION_A ya = y0a + s1a*g1a(x) as reference to a sum function
    //! @param FUNCTION_B     yb = g1b(x) as reference to a function
    //! @return new SumFunctionMixin y = ya - yb = y0a + s1a*g1a(x) - g1b(x)
    template< typename t_Interface>
    inline SumFunctionMixin< t_Interface> operator -
    (
      const SumFunctionMixin< t_Interface> &SUM_FUNCTION_A,
      const t_Interface &FUNCTION_B
    )
    {
      SumFunctionMixin< t_Interface>( SUM_FUNCTION_A).operator -=( FUNCTION_B);
    }

    //! @brief operator SumFunctionMixin + VALUE
    //! create a function y = y0 + y0' + 1*g1(x)
    //! @param VALUE y0'
    //! @param SUM_FUNCTION function y = y0 + s1*g1(x) as reference to a SumFunctionMixin
    //! @return SumFunctionMixin f(x) = y0 + y0' + s1*g1(x)
    template< typename t_Interface>
    inline SumFunctionMixin< t_Interface> operator +
    (
      const SumFunctionMixin< t_Interface> &SUM_FUNCTION,
      const typename SumFunctionMixin< t_Interface>::t_ResultType &VALUE
    )
    {
      return SumFunctionMixin< t_Interface>( SUM_FUNCTION).operator +=( VALUE);
    }

    //! @brief operator VALUE + SumFunctionMixin
    //! create a function y = y0 + y0' + 1*g1(x)
    //! @param VALUE y0'
    //! @param SUM_FUNCTION function y = y0 + s1*g1(x) as reference to a SumFunctionMixin
    //! @return SumFunctionMixin f(x) = y0 + y0' + s1*g1(x)
    template< typename t_Interface>
    inline SumFunctionMixin< t_Interface> operator +
    (
      const typename SumFunctionMixin< t_Interface>::t_ResultType &VALUE,
      const SumFunctionMixin< t_Interface> &SUM_FUNCTION
    )
    {
      return SumFunctionMixin< t_Interface>( SUM_FUNCTION).operator +=( VALUE);
    }

    //! @brief operator SumFunctionMixin - VALUE
    //! create a function y = y0 - y0' + 1*g1(x)
    //! @param VALUE y0'
    //! @param SUM_FUNCTION function y = y0 + s1*g1(x) as reference to a SumFunctionMixin
    //! @return SumFunctionMixin f(x) = y0 - y0' + s1*g1(x)
    template< typename t_Interface>
    inline SumFunctionMixin< t_Interface> operator -
    (
      const SumFunctionMixin< t_Interface> &SUM_FUNCTION,
      const typename SumFunctionMixin< t_Interface>::t_ResultType &VALUE
    )
    {
      return SumFunctionMixin< t_Interface>( SUM_FUNCTION).operator -=( VALUE);
    }

    //! @brief operator VALUE - SumFunctionMixin
    //! create a function y = y0 - y0' + 1*g1(x)
    //! @param VALUE y0'
    //! @param SUM_FUNCTION function y = y0 + s1*g1(x) as reference to a SumFunctionMixin
    //! @return SumFunctionMixin f(x) = y0 - y0' + s1*g1(x)
    template< typename t_Interface>
    inline SumFunctionMixin< t_Interface> operator -
    (
      const typename SumFunctionMixin< t_Interface>::t_ResultType &VALUE,
      const SumFunctionMixin< t_Interface> &SUM_FUNCTION
    )
    {
      return SumFunctionMixin< t_Interface>( -SUM_FUNCTION).operator +=( VALUE);
    }

    //! @brief operator SumFunctionMixin * SCALAR
    //! create function f(x) = s*y = s*y0 + s*s1*g1(x)
    //! @param SUM_FUNCTION y = y0 + s1*g1(x) as reference to SumFunctionMixin
    //! @param SCALAR scalar s
    //! @return reference to new sum function f(x) = s*y0 + s*s1*g1(x)
    template< typename t_Interface>
    inline SumFunctionMixin< t_Interface> operator *
    (
      const SumFunctionMixin< t_Interface> &SUM_FUNCTION,
      const double &SCALAR
    )
    {
      return SumFunctionMixin< t_Interface>( SUM_FUNCTION).operator *=( SCALAR);
    }

    //! @brief operator SCALAR * SumFunctionMixin
    //! create function f(x) = s*y = s*y0 + s*s1*g1(x)
    //! @param SUM_FUNCTION y = y0 + s1*g1(x) as reference to SumFunctionMixin
    //! @param SCALAR scalar s
    //! @return reference to new sum function f(x) = s*y0 + s*s1*g1(x)
    template< typename t_Interface>
    inline SumFunctionMixin< t_Interface> operator *
    (
      const double &SCALAR,
      const SumFunctionMixin< t_Interface> &SUM_FUNCTION
    )
    {
      return SumFunctionMixin< t_Interface>( SUM_FUNCTION).operator *=( SCALAR);
    }

    //! @brief operator SumFunctionMixin / SCALAR
    //! create function f(x) = 1/s*y0 + 1/s*g(x)
    //! @param SUM_FUNCTION f(x) = y0 + s1*g1(x) as reference to SumFunctionMixin
    //! @param SCALAR scalar s
    //! @return reference to new sum function f(x) = 1/s*y0 + 1/s*g1(x)
    template< typename t_Interface>
    inline SumFunctionMixin< t_Interface> operator /
    (
      const SumFunctionMixin< t_Interface> &SUM_FUNCTION,
      const double &SCALAR
    )
    {
      return SumFunctionMixin< t_Interface>( SUM_FUNCTION).operator /=( SCALAR);
    }

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_SUM_FUNCTION_MIXIN_H_
