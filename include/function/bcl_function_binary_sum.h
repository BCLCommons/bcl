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

#ifndef BCL_FUNCTION_BINARY_SUM_H_
#define BCL_FUNCTION_BINARY_SUM_H_

// include the namespace header
#include "bcl_function.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_function_binary_interface.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace function
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BinarySum
    //! @brief is a binary function version of BinarySum
    //!
    //! @tparam t_ArgumentType1 Type of the first argument to the function
    //! @tparam t_ArgumentType2 Type of the second argument to the function
    //! @tparam t_ResulType     Type of the result of the function
    //!
    //! @see @link example_function_binary_sum.cpp @endlink
    //! @author karakam
    //! @date Jun 6, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    class BinarySum :
      public BinaryInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType>
    {
    public:

    ///////////
    // types //
    ///////////

      //! @typedef type of terms
      typedef BinaryInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> FunctionType;

    private:

    //////////
    // data //
    //////////

      //! absolute value y0
      t_ResultType m_Absolute;

      //! @brief vector of pairs of coefficients and pointer on functions - representing s1 * g1(x1,x2) + s2 * g2(x1,x2)
      //! where s * g(x1,x2) is the Pair, and s the scalar and g(x1,x2) the ShPtr< FunctionInterface>
      storage::Vector< storage::Pair< double, util::ShPtr< FunctionType> > > m_Function;

      //! scheme
      std::string m_Scheme;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      BinarySum() :
        m_Absolute(),
        m_Function(),
        m_Scheme()
      {
      }

      //! @brief construct from ABSOLUTE
      //! @param ABSOLUTE y0 for the sumfunction
      BinarySum( const t_ResultType &ABSOLUTE) :
        m_Absolute( ABSOLUTE),
        m_Function(),
        m_Scheme()
      {
      }

      //! @brief construct from ShPtr on Function, optional coefficient, and optional absolute value
      //! y = y0 + s1*g1(x1,x2)
      //! @param SP_FUNCTION g1(x1,x2) as part of sumfunction as ShPtr to FunctionInterface
      //! @param COEFFICIENT s1(x1,x2) as coefficient for g1(x1,x2) -> s1 * g1(x1,x2) - default 1
      //! @param ABSOLUTE    y0 as absolute in sumfunction - default 0
      BinarySum
      (
        const util::ShPtr< FunctionType> &SP_FUNCTION,
        const double &COEFFICIENT = double( 1),
        const t_ResultType &ABSOLUTE = t_ResultType()
      ) :
        m_Absolute( ABSOLUTE),
        m_Function(),
        m_Scheme()
      {
        NewOperand( SP_FUNCTION, COEFFICIENT);
      }

      //! @brief construct from two Functions, optional coefficients, and optional absolute value
      //! y = y0 + s1*g1(x1,x2) + s2*g2(x1,x2)
      //! @param SP_FUNCTION_A g1(x1,x2) as part of sumfunction as ShPtr to FunctionInterface
      //! @param SP_FUNCTION_B g2(x1,x2) as part of sumfunction as ShPtr to FunctionInterface
      //! @param COEFFICIENT_A s1(x1,x2) as coefficient for g1(x1,x2) -> s1 * g1(x1,x2) - default 1
      //! @param COEFFICIENT_B s2(x1,x2) as coefficient for g2(x1,x2) -> s2 * g2(x1,x2) - default 1
      //! @param ABSOLUTE    y0 as absolute in sumfunction - default 0
      BinarySum
      (
        const util::ShPtr< FunctionType> &SP_FUNCTION_A,
        const util::ShPtr< FunctionType> &SP_FUNCTION_B,
        const double &COEFFICIENT_A = double( 1),
        const double &COEFFICIENT_B = double( 1),
        const t_ResultType &ABSOLUTE = t_ResultType()
      ) :
        m_Absolute( ABSOLUTE),
        m_Function(),
        m_Scheme()
      {
        NewOperand( SP_FUNCTION_A, COEFFICIENT_A);
        NewOperand( SP_FUNCTION_B, COEFFICIENT_B);
      }

      //! @brief copy constructor
      //! @return pointer to new BinarySum which is a copy of this
      BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType> *Clone() const
      {
        return new BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>( *this);
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

      //! @brief returns the vector of functions and the associated weights
      //! @return the vector of functions and the associated weights
      const storage::Vector< storage::Pair< double, util::ShPtr< FunctionType> > > &
      GetFunction() const
      {
        return m_Function;
      }

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief assign Function
      //! @param FUNCTION g1( x) with s1 = 1
      //! @return reference to this BinarySum
      BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType> &
      operator =( const FunctionType &FUNCTION)
      {
        // reset all funtions in sum
        m_Function.Reset();

        // add this function as new operand
        NewOperand( util::CloneToShPtr( FUNCTION), double( 1));

        // reset m_Absoulte to 0
        m_Absolute = t_ResultType( 0);

        // return this
        return *this;
      }

      //! @brief actual f(x1,x2) operator taking two arguments and returning the sum of all function values
      //! @param ARGUMENT_1 x1 to f(x1,x2)
      //! @param  ARGUMENT_2 x2 to f(x1,x2)
      //! @return the value of f(x1,x2) = sum
      t_ResultType operator()
      (
        t_ArgumentType1 &ARGUMENT_1,
        t_ArgumentType2 &ARGUMENT_2
      ) const
      {
        // initialize sum with the absolute
        t_ResultType sum( m_Absolute);

        // iterate over all functions, calculate their function value and add them multiplied with their coefficient
        for
        (
          typename storage::Vector< storage::Pair< t_ResultType, util::ShPtr< FunctionType> > >::const_iterator
            pair_itr( m_Function.Begin()), pair_itr_end( m_Function.End());
          pair_itr != pair_itr_end;
          ++pair_itr
        )
        {
          sum += pair_itr->First() * pair_itr->Second()->operator()( ARGUMENT_1, ARGUMENT_2);
        }

        // return sum f(x1,x2) = y
        return sum;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief add an new operand with an optional coefficient to the object s1 * g1(x1,x2)
      //! @param SP_FUNCTION function g(x1,x2) as ShPtr to Function
      //! @param COEFFICIENT s like in s1 * g1(x1,x2) - default is 1
      void NewOperand
      (
        const util::ShPtr< FunctionType> &SP_FUNCTION,
        const double &COEFFICIENT = double( 1)
      )
      {
        // insert pair of function and coefficient
        m_Function.PushBack( storage::Pair< double, util::ShPtr< FunctionType> >( COEFFICIENT, SP_FUNCTION));

        // update the scheme
        m_Scheme += SP_FUNCTION->GetScheme();
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator -
      //! negate BinarySum y = y0 + s1*g1(x1,x2) becomes y = -y0 - s1*g1(x1,x2)
      //! @param SUM_FUNCTION y = y0 + s1*g1(x1,x2) as reference to sum function
      //! @return new negated BinarySum
      BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>
      operator -() const
      {
        return BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>( *this).operator *=( t_ResultType( -1));
      }

      //! @brief operator += Function
      //! @param FUNCTION function g(x1,x2) to be added to the sumfunction as reference - with coefficient 1
      //! @return reference to this altered sum function
      BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType> &
      operator +=( const FunctionType &FUNCTION)
      {
        NewOperand( util::CloneToShPtr( FUNCTION));
        return *this;
      }

      //! @brief operator -= Function
      //! @param FUNCTION function g(x1,x2) to be subtracted from the sumfunction as reference - with coefficient 1
      //! @return reference to this altered sum function
      BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType> &
      operator -=( const FunctionType &FUNCTION)
      {
        NewOperand( util::CloneToShPtr( FUNCTION), double( -1));

        return *this;
      }

      //! @brief operator += BinarySum
      //! @param SUM_FUNCTION sum function to be added to the sum function as reference
      //! @return reference to this altered sum function
      BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType> &
      operator +=( const BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType> &SUM_FUNCTION)
      {
        // append functions in SUM_FUNCTION to this sum function
        m_Function.Append( SUM_FUNCTION.m_Function);

        // add the absolute
        m_Absolute += SUM_FUNCTION.m_Absolute;

        // return
        return *this;
      }

      //! @brief operator -= BinarySum
      //! @param SUM_FUNCTION sum function to be subtracted from the sum function as reference
      //! @return reference to this altered sum function
      BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType> &
      operator -=( const BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType> &SUM_FUNCTION)
      {
        return operator +=( -SUM_FUNCTION);
      }

      //! @brief operator += VALUE adds y to y0
      //! @param VALUE y to be added to y0
      BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType> &
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
      BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType> &
      operator -=( const t_ResultType &VALUE)
      {
        // subtract from absolute
        m_Absolute -= VALUE;

        // return
        return *this;
      }

      //! @brief operator *= SCALAR
      //! mutlipy absolute and each coefficient with the given scalar: s*f(x1,x2) = s*y0 + s*s1*g1(x1,x2) + s*s2*g2(x1,x2)...
      //! @param SCALAR scalar to be used to multiply every term
      //! @return reference to altered sum function
      BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType> &
      operator *=( const double &SCALAR)
      {
        //multiply the absolute value
        m_Absolute *= SCALAR;

        // multiply each coefficient by the given SCALAR
        //iterate over all functions, calculate their function value and add them multiplied with their coefficient
        for
        (
          typename std::vector< storage::Pair< double, util::ShPtr< FunctionType> > >::iterator
            pair_itr( m_Function.Begin()), pair_itr_end( m_Function.End());
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
      //! devide absolute and each coefficient with the given scalar: f(x1,x2)/s = y0/s + s1/s*g1(x1,x2) + s2/s*g2(x1,x2)...
      //! @param SCALAR scalar to be used to devide every term
      //! @return reference to altered sum function
      BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType> &
      operator /=( const double &SCALAR)
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
      //! @param ARGUMENT_1 First argument to be used to evaluate the function
      //! @param ARGUMENT_2 Second argument to be used to evaluate the function
      //! @param OSTREAM std::ostream to be written to
      //! @return std::ostream which was written to
      std::ostream &WriteDetailedSchemeAndValues
      (
        t_ArgumentType1 &ARGUMENT_1,
        t_ArgumentType2 &ARGUMENT_2,
        std::ostream &OSTREAM
      ) const
      {
        for
        (
          typename std::vector< storage::Pair< t_ResultType, util::ShPtr< FunctionType> > >::const_iterator
            pair_itr( m_Function.Begin()), pair_itr_end( m_Function.End());
          pair_itr != pair_itr_end;
          ++pair_itr
        )
        {
          OSTREAM << pair_itr->Second()->GetScheme();
          OSTREAM << '\t' << ( pair_itr->First() * pair_itr->Second()->operator()( ARGUMENT_1, ARGUMENT_2))
                  << "                " << '\n';

          pair_itr->Second()->WriteDetailedSchemeAndValues( ARGUMENT_1, ARGUMENT_2, OSTREAM) << '\n';
        }

        //end
        return OSTREAM;
      }

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // read members
        io::Serialize::Read( m_Absolute, ISTREAM);
        io::Serialize::Read( m_Function, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_Absolute, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Function, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // template class BinarySum

    // instantiate s_Instance
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    const util::SiPtr< const util::ObjectInterface> BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>::s_Instance
    (
      GetObjectInstances().AddInstance( new BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>())
    );

  ///////////////
  // operators //
  ///////////////

    //! @brief operator BinarySum + BinarySum
    //! @param FUNCTION_A ya = y0a + s1a*g1a(x1,x2) as reference to a sum function
    //! @param FUNCTION_B yb = y0b + s1b*g1b(x1,x2) as reference to a sum function
    //! @return new BinarySum y = ya + yb = y0a + y0b + s1a*g1a(x1,x2) + s1b*g1b(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    inline
    BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>
    operator +
    (
      const BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION_A,
      const BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION_B
    )
    {
      return BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>( FUNCTION_A).operator +=( FUNCTION_B);
    }

    //! @brief operator FunctionInterface + BinarySum
    //! @param FUNCTION_A     ya = g1a(x1,x2) as reference to a function
    //! @param SUM_FUNCTION_B yb = y0b + s1b*g1b(x1,x2) as reference to a sum function
    //! @return new BinarySum y = ya + yb = y0b + g1a(x1,x2) + s1b*g1b(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    inline
    BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>
    operator +
    (
      const BinaryInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType>     &FUNCTION_A,
      const BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType> &SUM_FUNCTION_B
    )
    {
      return BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>( SUM_FUNCTION_B).operator +=( FUNCTION_A);
    }

    //! @brief operator BinarySum + FunctionInterface
    //! @param SUM_FUNCTION_A ya = y0a + s1a*g1a(x1,x2) as reference to a sum function
    //! @param FUNCTION_B     yb = g1b(x1,x2) as reference to a function
    //! @return new BinarySum y = ya + yb = y0a + s1a*g1a(x1,x2) + g1b(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    inline
    BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>
    operator +
    (
      const BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType> &SUM_FUNCTION_A,
      const BinaryInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION_B
    )
    {
      return BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>( SUM_FUNCTION_A).operator +=( FUNCTION_B);
    }

    //! @brief operator FunctionInterface + FunctionInterface
    //! @param FUNCTION_A ya = g1a(x1,x2) as reference to a function
    //! @param FUNCTION_B yb = g1b(x1,x2) as reference to a function
    //! @return new BinarySum y = ya + yb = 0 + g1a(x1,x2) + g1b(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    inline
    BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>
    operator +
    (
      const BinaryInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION_A,
      const BinaryInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION_B
    )
    {
      return BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>( FUNCTION_A, FUNCTION_B);
    }

    //! @brief operator BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> + VALUE
    //! create a function y = y0 + 1*g1(x1,x2)
    //! @param FUNCTION function g1(x1,x2) as reference to a Function
    //! @param VALUE y0
    //! @return BinarySum f(x1,x2) = y0 + 1*g1(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    inline
    BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>
    operator +
    (
      const BinaryInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION,
      const t_ResultType &VALUE
    )
    {
      return BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>( FUNCTION, double( 1), VALUE);
    }

    //! @brief operator BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> + VALUE
    //! create a function y = y0 + 1*g1(x1,x2)
    //! @param VALUE y0
    //! @param FUNCTION function g1(x1,x2) as reference to a Function
    //! @return BinarySum f(x1,x2) = y0 + 1*g1(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    inline
    BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>
    operator +
    (
      const t_ResultType &VALUE,
      const BinaryInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION
    )
    {
      return BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>( FUNCTION, double( 1), VALUE);
    }

    //! @brief operator BinarySum + VALUE
    //! create a function y = y0 + y0' + 1*g1(x1,x2)
    //! @param VALUE y0'
    //! @param SUM_FUNCTION function y = y0 + s1*g1(x1,x2) as reference to a BinarySum
    //! @return BinarySum f(x1,x2) = y0 + y0' + s1*g1(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    inline
    BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>
    operator +
    (
      const BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType> &SUM_FUNCTION,
      const t_ResultType &VALUE
    )
    {
      return BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>( SUM_FUNCTION).operator +=( VALUE);
    }

    //! @brief operator VALUE + BinarySum
    //! create a function y = y0 + y0' + 1*g1(x1,x2)
    //! @param VALUE y0'
    //! @param SUM_FUNCTION function y = y0 + s1*g1(x1,x2) as reference to a BinarySum
    //! @return BinarySum f(x1,x2) = y0 + y0' + s1*g1(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    inline
    BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>
    operator +
    (
      const t_ResultType &VALUE,
      const BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType> &SUM_FUNCTION
    )
    {
      return BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>( SUM_FUNCTION).operator +=( VALUE);
    }

    //! @brief operator BinarySum - BinarySum
    //! @param FUNCTION_A ya = y0a + s1a*g1a(x1,x2) as reference to a sum function
    //! @param FUNCTION_B yb = y0b + s1b*g1b(x1,x2) as reference to a sum function
    //! @return new BinarySum y = ya - yb = y0a - y0b + s1a*g1a(x1,x2) - s1b*g1b(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    inline
    BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>
    operator -
    (
      const BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION_A,
      const BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION_B
    )
    {
      return BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>( FUNCTION_A).operator -=( FUNCTION_B);
    }

    //! @brief operator FunctionInterface - BinarySum
    //! @param FUNCTION_A     ya = g1a(x1,x2) as reference to a function
    //! @param SUM_FUNCTION_B yb = y0b - s1b*g1b(x1,x2) as reference to a sum function
    //! @return new BinarySum y = ya - yb = - y0b + g1a(x1,x2) - s1b*g1b(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    inline
    BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>
    operator -
    (
      const BinaryInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION_A,
      const BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType> &SUM_FUNCTION_B
    )
    {
      return BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>( -SUM_FUNCTION_B).operator +=( FUNCTION_A);
    }

    //! @brief operator BinarySum - FunctionInterface
    //! @param SUM_FUNCTION_A ya = y0a + s1a*g1a(x1,x2) as reference to a sum function
    //! @param FUNCTION_B     yb = g1b(x1,x2) as reference to a function
    //! @return new BinarySum y = ya - yb = y0a + s1a*g1a(x1,x2) - g1b(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    inline
    BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>
    operator -
    (
      const BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType> &SUM_FUNCTION_A,
      const BinaryInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION_B
    )
    {
      BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>( SUM_FUNCTION_A).operator -=( FUNCTION_B);
    }

    //! @brief operator FunctionInterface - FunctionInterface
    //! @param FUNCTION_A ya = g1a(x1,x2) as reference to a function
    //! @param FUNCTION_B yb = g1b(x1,x2) as reference to a function
    //! @return new BinarySum y = ya - yb = 0 + g1a(x1,x2) - g1b(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    inline
    BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>
    operator -
    (
      const BinaryInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION_A,
      const BinaryInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION_B
    )
    {
      return BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>( FUNCTION_A, FUNCTION_B, t_ResultType( 1), t_ArgumentType1( -1));
    }

    //! @brief operator BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> - VALUE
    //! create a function y = -y0 + 1*g1(x1,x2)
    //! @param VALUE y0
    //! @param FUNCTION function g1(x1,x2) as reference to a Function
    //! @return BinarySum f(x1,x2) = -y0 + 1*g1(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    inline
    BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>
    operator -
    (
      const BinaryInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION,
      const t_ResultType &VALUE
    )
    {
      return BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>( util::CloneToShPtr( FUNCTION), double( 1), -VALUE);
    }

    //! @brief operator  VALUE - BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType>
    //! create a function y = y0 - 1*g1(x1,x2)
    //! @param VALUE y0
    //! @param FUNCTION function g1(x1,x2) as reference to a Function
    //! @return BinarySum f(x1,x2) = y0 - 1*g1(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    inline
    BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>
    operator -
    (
      const t_ResultType &VALUE,
      const BinaryInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION
    )
    {
      return BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>( FUNCTION, double( -1), VALUE);
    }

    //! @brief operator BinarySum - VALUE
    //! create a function y = y0 - y0' + 1*g1(x1,x2)
    //! @param VALUE y0'
    //! @param SUM_FUNCTION function y = y0 + s1*g1(x1,x2) as reference to a BinarySum
    //! @return BinarySum f(x1,x2) = y0 - y0' + s1*g1(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    inline
    BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>
    operator -
    (
      const BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType> &SUM_FUNCTION,
      const t_ResultType &VALUE
    )
    {
      return BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>( SUM_FUNCTION).operator -=( VALUE);
    }

    //! @brief operator VALUE - BinarySum
    //! create a function y = y0 - y0' + 1*g1(x1,x2)
    //! @param VALUE y0'
    //! @param SUM_FUNCTION function y = y0 + s1*g1(x1,x2) as reference to a BinarySum
    //! @return BinarySum f(x1,x2) = y0 - y0' + s1*g1(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    inline
    BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>
    operator -
    (
      const t_ResultType &VALUE,
      const BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType> &SUM_FUNCTION
    )
    {
      return BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>( -SUM_FUNCTION).operator +=( VALUE);
    }

    //! @brief operator Function * SCALAR
    //! create function f(x1,x2) = 0 + s*g(x1,x2)
    //! @param FUNCTION g(x1,x2) as reference to Function
    //! @param SCALAR scalar s
    //! @return reference to new sum function f(x1,x2) = 0 + s*g(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    inline
    BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>
    operator *
    (
      const BinaryInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION,
      const double &SCALAR
    )
    {
      return BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>( FUNCTION, SCALAR);
    }

    //! @brief operator SCALAR * Function
    //! create function f(x1,x2) = 0 + s*g(x1,x2)
    //! @param FUNCTION g(x1,x2) as reference to Function
    //! @param SCALAR scalar s
    //! @return reference to new sum function f(x1,x2) = 0 + s*g(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    inline
    BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>
    operator *
    (
      const double &SCALAR,
      const BinaryInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION
    )
    {
      return BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>( util::CloneToShPtr( FUNCTION), SCALAR);
    }

    //! @brief operator BinarySum * SCALAR
    //! create function f(x1,x2) = s*y = s*y0 + s*s1*g1(x1,x2)
    //! @param SUM_FUNCTION y = y0 + s1*g1(x1,x2) as reference to BinarySum
    //! @param SCALAR scalar s
    //! @return reference to new sum function f(x1,x2) = s*y0 + s*s1*g1(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    inline
    BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>
    operator *
    (
      const BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType> &SUM_FUNCTION,
      const double &SCALAR
    )
    {
      return BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>( SUM_FUNCTION).operator *=( SCALAR);
    }

    //! @brief operator SCALAR * BinarySum
    //! create function f(x1,x2) = s*y = s*y0 + s*s1*g1(x1,x2)
    //! @param SUM_FUNCTION y = y0 + s1*g1(x1,x2) as reference to BinarySum
    //! @param SCALAR scalar s
    //! @return reference to new sum function f(x1,x2) = s*y0 + s*s1*g1(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    inline
    BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>
    operator *
    (
      const double &SCALAR,
      const BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType> &SUM_FUNCTION
    )
    {
      return BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>( SUM_FUNCTION).operator *=( SCALAR);
    }

    //! @brief operator Function / SCALAR
    //! create function f(x1,x2) = 0 + 1/s*g(x1,x2)
    //! @param FUNCTION g(x1,x2) as reference to Function
    //! @param SCALAR scalar s
    //! @return reference to new sum function f(x1,x2) = 0 + 1/s*g(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    inline
    BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>
    operator /
    (
      const BinaryInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION,
      const double &SCALAR
    )
    {
      return BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>( util::CloneToShPtr( FUNCTION), double( 1) / SCALAR);
    }

    //! @brief operator BinarySum / SCALAR
    //! create function f(x1,x2) = 1/s*y0 + 1/s*g(x1,x2)
    //! @param SUM_FUNCTION f(x1,x2) = y0 + s1*g1(x1,x2) as reference to BinarySum
    //! @param SCALAR scalar s
    //! @return reference to new sum function f(x1,x2) = 1/s*y0 + 1/s*g1(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    inline
    BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>
    operator /
    (
      const BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType> &SUM_FUNCTION,
      const double &SCALAR
    )
    {
      return BinarySum< t_ArgumentType1, t_ArgumentType2, t_ResultType>( SUM_FUNCTION).operator /=( SCALAR);
    }

  } // namespace function
} // namespace bcl

#endif // BCL_FUNCTION_BINARY_SUM_H_
