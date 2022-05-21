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

#ifndef BCL_MATH_BINARY_SUM_FUNCTION_H_
#define BCL_MATH_BINARY_SUM_FUNCTION_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_binary_function_interface.h"

// external includes - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serializer.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_vector.h"

namespace bcl
{
  namespace math
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BinarySumFunction
    //! @brief is a binary function version of BinarySumFunction
    //!
    //! @tparam t_ArgumentType1 Type of the first argument to the function
    //! @tparam t_ArgumentType2 Type of the second argument to the function
    //! @tparam t_ResulType     Type of the result of the function
    //!
    //! @see @link example_math_binary_sum_function.cpp @endlink
    //! @author karakam
    //! @date Jun 6, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType, typename t_Scalar>
    class BinarySumFunction :
      public BinaryFunctionInterfaceSerializable< t_ArgumentType1, t_ArgumentType2, t_ResultType>
    {
    private:

    //////////
    // data //
    //////////

      //! absolute value y0
      t_ResultType m_Absolute;

      //! @brief vector of pairs of coefficients and pointer on functions - representing s1 * g1(x1,x2) + s2 * g2(x1,x2)
      //! where s * g(x1,x2) is the Pair, and s the scalar and g(x1,x2) the ShPtr< FunctionInterface>
      storage::Vector
      <
        storage::Pair
        <
          t_Scalar,
          util::ShPtr< BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> >
        >
      >
      m_Function;

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
      BinarySumFunction
      (
      ) :
        m_Absolute(),
        m_Function(),
        m_Scheme()
      {
      }

      //! @brief construct from ABS_CONST
      //! @param ABS_CONST y0 for the sumfunction
      BinarySumFunction
      (
        const t_ResultType &ABS_CONST
      ) :
        m_Absolute( ABS_CONST),
        m_Function(),
        m_Scheme()
      {
      }

      //! @brief copy constructor
      //! @param SUM_FUNCTION the sum_function which should be copied
      BinarySumFunction
      (
        const BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &SUM_FUNCTION
      ) :
        m_Absolute( SUM_FUNCTION.m_Absolute),
        m_Function( SUM_FUNCTION.m_Function),
        m_Scheme( SUM_FUNCTION.m_Scheme)
      {
      }

      //! @brief construct from ShPtr on Function, optional coefficient, and optional absolute value
      //! y = y0 + s1*g1(x1,x2)
      //! @param SP_FUNCTION g1(x1,x2) as part of sumfunction as ShPtr to FunctionInterface
      //! @param COEFFICIENT s1(x1,x2) as coefficient for g1(x1,x2) -> s1 * g1(x1,x2) - default 1
      //! @param ABS_CONST    y0 as absolute in sumfunction - default 0
      BinarySumFunction
      (
        const util::ShPtr< BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> > &SP_FUNCTION,
        const t_Scalar &COEFFICIENT = t_Scalar( 1),
        const t_ResultType &ABS_CONST = t_ResultType()
      ) :
        m_Absolute( ABS_CONST),
        m_Function(),
        m_Scheme()
      {
        NewOperand( SP_FUNCTION, COEFFICIENT);
      }

      //! @brief construct from Function, optional coefficient, and optional absolute value
      //! y = y0 + s1*g1(x1,x2)
      //! @param FUNCTION    g1(x1,x2) as part of sumfunction as reference to FunctionInterface
      //! @param COEFFICIENT s1(x1,x2) as coefficient for g1(x1,x2) -> s1 * g1(x1,x2) - default 1
      //! @param ABS_CONST    y0 as absolute in sumfunction - default 0
      BinarySumFunction
      (
        const BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION,
        const t_Scalar &COEFFICIENT = t_Scalar( 1),
        const t_ResultType &ABS_CONST = t_ResultType()
      ) :
        m_Absolute( ABS_CONST),
        m_Function(),
        m_Scheme()
      {
        NewOperand
        (
          util::ShPtr< BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> >( FUNCTION.Clone()),
          COEFFICIENT
        );
      }

      //! @brief construct from two Functions, optional coefficients, and optional absolute value
      //! y = y0 + s1*g1(x1,x2) + s2*g2(x1,x2)
      //! @param SP_FUNCTION_A g1(x1,x2) as part of sumfunction as ShPtr to FunctionInterface
      //! @param SP_FUNCTION_B g2(x1,x2) as part of sumfunction as ShPtr to FunctionInterface
      //! @param COEFFICIENT_A s1(x1,x2) as coefficient for g1(x1,x2) -> s1 * g1(x1,x2) - default 1
      //! @param COEFFICIENT_B s2(x1,x2) as coefficient for g2(x1,x2) -> s2 * g2(x1,x2) - default 1
      //! @param ABS_CONST    y0 as absolute in sumfunction - default 0
      BinarySumFunction
      (
        const util::ShPtr< BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> > &SP_FUNCTION_A,
        const util::ShPtr< BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> > &SP_FUNCTION_B,
        const t_Scalar &COEFFICIENT_A = t_Scalar( 1),
        const t_Scalar &COEFFICIENT_B = t_Scalar( 1),
        const t_ResultType &ABS_CONST = t_ResultType()
      ) :
        m_Absolute( ABS_CONST),
        m_Function(),
        m_Scheme()
      {
        NewOperand( SP_FUNCTION_A, COEFFICIENT_A);
        NewOperand( SP_FUNCTION_B, COEFFICIENT_B);
      }

      //! @brief construct from two Functions, optional coefficients, and optional absolute value
      //! y = y0 + s1*g1(x1,x2) + s2*g2(x1,x2)
      //! @param FUNCTION_A    g1(x1,x2) as part of sumfunction as reference to FunctionInterface
      //! @param FUNCTION_B    g2(x1,x2) as part of sumfunction as reference to FunctionInterface
      //! @param COEFFICIENT_A s1(x1,x2) as coefficient for g1(x1,x2) -> s1 * g1(x1,x2) - default 1
      //! @param COEFFICIENT_B s2(x1,x2) as coefficient for g2(x1,x2) -> s2 * g2(x1,x2) - default 1
      //! @param ABS_CONST    y0 as absolute in sumfunction - default 0
      BinarySumFunction
      (
        const BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION_A,
        const BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION_B,
        const t_Scalar &COEFFICIENT_A = t_Scalar( 1),
        const t_Scalar &COEFFICIENT_B = t_Scalar( 1),
        const t_ResultType &ABS_CONST = t_ResultType()
      ) :
        m_Absolute( ABS_CONST),
        m_Function(),
        m_Scheme()
      {
        NewOperand
        (
          util::ShPtr< BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> >( FUNCTION_A.Clone()),
          COEFFICIENT_A
        );
        NewOperand
        (
          util::ShPtr< BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> >( FUNCTION_B.Clone()),
          COEFFICIENT_B
        );
      }

      //! @brief copy constructor
      //! @return pointer to new BinarySumFunction which is a copy of this
      BinarySumFunction *Clone() const
      {
        return new BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      virtual const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief returns the vector of functions and the associated weights
      //! @return the vector of functions and the associated weights
      const storage::Vector
      <
        storage::Pair< t_Scalar, util::ShPtr< BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> > >
      > &
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
      //! @return reference to this BinarySumFunction
      BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &
      operator =
      (
        const BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION
      )
      {
        // reset all funtions in sum
        m_Function.Reset();

        // add this function as new operand
        NewOperand
        (
          util::ShPtr< BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> >( FUNCTION.Clone()),
          t_ResultType( 1)
        );

        // reset m_Absoulte to 0
        m_Absolute = t_ResultType( 0);

        // return this
        return *this;
      }

      //! @brief assignment operator
      //! @param SUM_FUNCTION sum function this should be assigned to
      //! @return reference to this BinarySumFunction
      BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &
      operator  =
      (
        const BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &SUM_FUNCTION
      )
      {
        // assign data member
        m_Absolute = SUM_FUNCTION.m_Absolute;
        m_Function = SUM_FUNCTION.m_Function;

        // return this
        return *this;
      }

      //! @brief actual f(x1,x2) operator taking two arguments and returning the sum of all function values
      //! @param ARGUMENT_1 x1 to f(x1,x2)
      //! @param  ARGUMENT_2 x2 to f(x1,x2)
      //! @return the value of f(x1,x2) = sum
      t_ResultType
      operator()
      (
        const t_ArgumentType1 &ARGUMENT_1,
        const t_ArgumentType2 &ARGUMENT_2
      ) const;

    ////////////////
    // operations //
    ////////////////

      //! @brief add an new operand with an optional coefficient to the object s1 * g1(x1,x2)
      //! @param SP_FUNCTION function g(x1,x2) as ShPtr to Function
      //! @param COEFFICIENT s like in s1 * g1(x1,x2) - default is 1
      void NewOperand
      (
        const util::ShPtr< BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> > &SP_FUNCTION,
        const t_Scalar &COEFFICIENT = t_Scalar( 1)
      )
      {
        // insert pair of function and coefficient
        m_Function.PushBack
        (
          storage::Pair
          <
            t_Scalar, util::ShPtr< BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> >
          >
          (
            COEFFICIENT, SP_FUNCTION
          )
        );

        // update the identifier
        m_Scheme += SP_FUNCTION->GetScheme();
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator += Function
      //! @param FUNCTION function g(x1,x2) to be added to the sumfunction as reference - with coefficient 1
      //! @return reference to this altered sum function
      BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &
      operator +=
      (
        const BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION
      )
      {
        NewOperand
        (
          util::ShPtr< BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> >( FUNCTION.Clone())
        );
        return *this;
      }

      //! @brief operator -= Function
      //! @param FUNCTION function g(x1,x2) to be subtracted from the sumfunction as reference - with coefficient 1
      //! @return reference to this altered sum function
      BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &
      operator -=
      (
        const BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION
      )
      {
        NewOperand
        (
          util::ShPtr< BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> >( FUNCTION.Clone()),
          t_ResultType( -1)
        );

        return *this;
      }

      //! @brief operator += BinarySumFunction
      //! @param SUM_FUNCTION sum function to be added to the sum function as reference
      //! @return reference to this altered sum function
      BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &
      operator +=
      (
        const BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &SUM_FUNCTION
      )
      {
        // append functions in SUM_FUNCTION to this sum function
        m_Function.Append( SUM_FUNCTION.m_Function);

        // add the absolute
        m_Absolute += SUM_FUNCTION.m_Absolute;

        // return
        return *this;
      }

      //! @brief operator -= BinarySumFunction
      //! @param SUM_FUNCTION sum function to be subtracted from the sum function as reference
      //! @return reference to this altered sum function
      BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &
      operator -=
      (
        const BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &SUM_FUNCTION
      )
      {
        return operator +=( -SUM_FUNCTION);
      }

      //! @brief operator += VALUE adds y to y0
      //! @param VALUE y to be added to y0
      BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &
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
      BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &
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
      //! mutlipy absolute and each coefficient with the given scalar: s*f(x1,x2) = s*y0 + s*s1*g1(x1,x2) + s*s2*g2(x1,x2)...
      //! @param SCALAR scalar to be used to multiply every term
      //! @return reference to altered sum function
      BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &
      operator *=
      (
        const t_Scalar &SCALAR
      )
      {
        //multiply the absolute value
        m_Absolute *= SCALAR;

        // multiply each coefficient by the given SCALAR
        //iterate over all functions, calculate their function value and add them multiplied with their coefficient
        for( typename std::vector< storage::Pair< t_Scalar, util::ShPtr< BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> > > >::iterator
               pair_itr( m_Function.Begin()), pair_itr_end( m_Function.End());
             pair_itr != pair_itr_end;
             ++pair_itr)
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
      BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &
      operator /=
      (
        const t_Scalar &SCALAR
      )
      {
        // check that scalr is not 0
        BCL_Assert( SCALAR != t_Scalar( 0), "attempt division by zero!");

        // call *= with 1/s
        return operator *=( t_Scalar( 1) / SCALAR);
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
        const t_ArgumentType1 &ARGUMENT_1,
        const t_ArgumentType2 &ARGUMENT_2,
        std::ostream &OSTREAM
      ) const;

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

    }; // class BinarySumFunction

    // instantiate s_Instance
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType, typename t_Scalar>
    const util::SiPtr< const util::ObjectInterface> BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>::s_Instance
    (
      GetObjectInstances().AddInstance
      (
        new BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>()
      )
    );

  /////////////////
  // data access //
  /////////////////

    //! @brief implement operator()
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType, typename t_Scalar>
    t_ResultType
    BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>::operator()
    (
      const t_ArgumentType1 &ARGUMENT_1,
      const t_ArgumentType2 &ARGUMENT_2
    ) const
    {
      // initialize sum with the absolute
      t_ResultType sum( m_Absolute);

      // iterate over all functions, calculate their function value and add them multiplied with their coefficient
      for
      (
        typename storage::Vector
                 <
                   storage::Pair< t_ResultType, util::ShPtr< BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> > >
                 >::const_iterator
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

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write detailed scheme and values to OSTREAM
    //! @param ARGUMENT_1 First argument to be used to evaluate the function
    //! @param ARGUMENT_2 First argument to be used to evaluate the function
    //! @param OSTREAM std::ostream to be written to
    //! @return std::ostream which was written to
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType, typename t_Scalar>
    inline
    std::ostream &
    BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>::WriteDetailedSchemeAndValues
    (
      const t_ArgumentType1 &ARGUMENT_1,
      const t_ArgumentType2 &ARGUMENT_2,
      std::ostream &OSTREAM
    ) const
    {
      for
      (
        typename std::vector< storage::Pair< t_ResultType, util::ShPtr< BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> > > >::const_iterator
          pair_itr( m_Function.Begin()), pair_itr_end( m_Function.End());
          pair_itr != pair_itr_end; ++pair_itr
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

  ///////////////
  // operators //
  ///////////////

    //! @brief operator BinarySumFunction + BinarySumFunction
    //! @param FUNCTION_A ya = y0a + s1a*g1a(x1,x2) as reference to a sum function
    //! @param FUNCTION_B yb = y0b + s1b*g1b(x1,x2) as reference to a sum function
    //! @return new BinarySumFunction y = ya + yb = y0a + y0b + s1a*g1a(x1,x2) + s1b*g1b(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType, typename t_Scalar>
    inline
    BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>
    operator +
    (
      const BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &FUNCTION_A,
      const BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &FUNCTION_B
    )
    {
      return BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>( FUNCTION_A).operator +=( FUNCTION_B);
    }

    //! @brief operator BinarySumFunction - BinarySumFunction
    //! @param FUNCTION_A ya = y0a + s1a*g1a(x1,x2) as reference to a sum function
    //! @param FUNCTION_B yb = y0b + s1b*g1b(x1,x2) as reference to a sum function
    //! @return new BinarySumFunction y = ya - yb = y0a - y0b + s1a*g1a(x1,x2) - s1b*g1b(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType, typename t_Scalar>
    inline
    BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>
    operator -
    (
      const BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &FUNCTION_A,
      const BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &FUNCTION_B
    )
    {
      return BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>( FUNCTION_A).operator -=( FUNCTION_B);
    }

    //! @brief operator FunctionInterface + BinarySumFunction
    //! @param FUNCTION_A     ya = g1a(x1,x2) as reference to a function
    //! @param SUM_FUNCTION_B yb = y0b + s1b*g1b(x1,x2) as reference to a sum function
    //! @return new BinarySumFunction y = ya + yb = y0b + g1a(x1,x2) + s1b*g1b(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType, typename t_Scalar>
    inline
    BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>
    operator +
    (
      const BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType>     &FUNCTION_A,
      const BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &SUM_FUNCTION_B
    )
    {
      return BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>( SUM_FUNCTION_B).operator +=( FUNCTION_A);
    }

    //! @brief operator FunctionInterface - BinarySumFunction
    //! @param FUNCTION_A     ya = g1a(x1,x2) as reference to a function
    //! @param SUM_FUNCTION_B yb = y0b - s1b*g1b(x1,x2) as reference to a sum function
    //! @return new BinarySumFunction y = ya - yb = - y0b + g1a(x1,x2) - s1b*g1b(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType, typename t_Scalar>
    inline
    BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>
    operator -
    (
      const BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION_A,
      const BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &SUM_FUNCTION_B
    )
    {
      return BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>( -SUM_FUNCTION_B).operator +=( FUNCTION_A);
    }

    //! @brief operator BinarySumFunction + FunctionInterface
    //! @param SUM_FUNCTION_A ya = y0a + s1a*g1a(x1,x2) as reference to a sum function
    //! @param FUNCTION_B     yb = g1b(x1,x2) as reference to a function
    //! @return new BinarySumFunction y = ya + yb = y0a + s1a*g1a(x1,x2) + g1b(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType, typename t_Scalar>
    inline
    BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>
    operator +
    (
      const BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &SUM_FUNCTION_A,
      const BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION_B
    )
    {
      return BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>( SUM_FUNCTION_A).operator +=( FUNCTION_B);
    }

    //! @brief operator BinarySumFunction - FunctionInterface
    //! @param SUM_FUNCTION_A ya = y0a + s1a*g1a(x1,x2) as reference to a sum function
    //! @param FUNCTION_B     yb = g1b(x1,x2) as reference to a function
    //! @return new BinarySumFunction y = ya - yb = y0a + s1a*g1a(x1,x2) - g1b(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType, typename t_Scalar>
    inline
    BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>
    operator -
    (
      const BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &SUM_FUNCTION_A,
      const BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION_B
    )
    {
      BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>( SUM_FUNCTION_A).operator -=( FUNCTION_B);
    }

    //! @brief operator FunctionInterface + FunctionInterface
    //! @param FUNCTION_A ya = g1a(x1,x2) as reference to a function
    //! @param FUNCTION_B yb = g1b(x1,x2) as reference to a function
    //! @return new BinarySumFunction y = ya + yb = 0 + g1a(x1,x2) + g1b(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    inline
    BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, double>
    operator +
    (
      const BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION_A,
      const BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION_B
    )
    {
      return BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, double>( FUNCTION_A, FUNCTION_B);
    }

    //! @brief operator FunctionInterface - FunctionInterface
    //! @param FUNCTION_A ya = g1a(x1,x2) as reference to a function
    //! @param FUNCTION_B yb = g1b(x1,x2) as reference to a function
    //! @return new BinarySumFunction y = ya - yb = 0 + g1a(x1,x2) - g1b(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    inline
    BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, double>
    operator -
    (
      const BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION_A,
      const BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION_B
    )
    {
      return BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, double>( FUNCTION_A, FUNCTION_B, t_ResultType( 1), t_ArgumentType1( -1));
    }

    //! @brief operator - BinarySumFunction
    //! negate BinarySumFunction y = y0 + s1*g1(x1,x2) becomes y = -y0 - s1*g1(x1,x2)
    //! @param SUM_FUNCTION y = y0 + s1*g1(x1,x2) as reference to sum function
    //! @return new negated BinarySumFunction
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType, typename t_Scalar>
    inline
    BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>
    operator -
    (
      const BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &SUM_FUNCTION
    )
    {
      return BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>( SUM_FUNCTION).operator *=( t_ResultType( -1));
    }

    //! @brief operator BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> + VALUE
    //! create a function y = y0 + 1*g1(x1,x2)
    //! @param FUNCTION function g1(x1,x2) as reference to a Function
    //! @param VALUE y0
    //! @return BinarySumFunction f(x1,x2) = y0 + 1*g1(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    inline
    BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, double>
    operator +
    (
      const BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION,
      const t_ResultType &VALUE
    )
    {
      return BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, double>( FUNCTION, double( 1), VALUE);
    }

    //! @brief operator BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> + VALUE
    //! create a function y = y0 + 1*g1(x1,x2)
    //! @param VALUE y0
    //! @param FUNCTION function g1(x1,x2) as reference to a Function
    //! @return BinarySumFunction f(x1,x2) = y0 + 1*g1(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    inline
    BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, double>
    operator +
    (
      const t_ResultType &VALUE,
      const BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION
    )
    {
      return BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, double>( FUNCTION, double( 1), VALUE);
    }

    //! @brief operator BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> - VALUE
    //! create a function y = -y0 + 1*g1(x1,x2)
    //! @param VALUE y0
    //! @param FUNCTION function g1(x1,x2) as reference to a Function
    //! @return BinarySumFunction f(x1,x2) = -y0 + 1*g1(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    inline
    BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, double>
    operator -
    (
      const BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION,
      const t_ResultType &VALUE
    )
    {
      return BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, double>( FUNCTION, double( 1), -VALUE);
    }

    //! @brief operator  VALUE - BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType>
    //! create a function y = y0 - 1*g1(x1,x2)
    //! @param VALUE y0
    //! @param FUNCTION function g1(x1,x2) as reference to a Function
    //! @return BinarySumFunction f(x1,x2) = y0 - 1*g1(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    inline
    BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, double>
    operator -
    (
      const t_ResultType &VALUE,
      const BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION
    )
    {
      return BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, double>( FUNCTION, double( -1), VALUE);
    }

    //! @brief operator Function * SCALAR
    //! create function f(x1,x2) = 0 + s*g(x1,x2)
    //! @param FUNCTION g(x1,x2) as reference to Function
    //! @param SCALAR scalar s
    //! @return reference to new sum function f(x1,x2) = 0 + s*g(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType, typename t_Scalar>
    inline
    BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>
    operator *
    (
      const BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION,
      const t_Scalar &SCALAR
    )
    {
      return BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>( FUNCTION, SCALAR);
    }

    //! @brief operator SCALAR * Function
    //! create function f(x1,x2) = 0 + s*g(x1,x2)
    //! @param FUNCTION g(x1,x2) as reference to Function
    //! @param SCALAR scalar s
    //! @return reference to new sum function f(x1,x2) = 0 + s*g(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType, typename t_Scalar>
    inline
    BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>
    operator *
    (
      const t_Scalar &SCALAR,
      const BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION
    )
    {
      return BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>( FUNCTION, SCALAR);
    }

    //! @brief operator Function / SCALAR
    //! create function f(x1,x2) = 0 + 1/s*g(x1,x2)
    //! @param FUNCTION g(x1,x2) as reference to Function
    //! @param SCALAR scalar s
    //! @return reference to new sum function f(x1,x2) = 0 + 1/s*g(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType, typename t_Scalar>
    inline
    BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>
    operator /
    (
      const BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION,
      const t_Scalar &SCALAR
    )
    {
      return BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>( FUNCTION, SCALAR( 1) / SCALAR);
    }

    //! @brief operator BinarySumFunction + VALUE
    //! create a function y = y0 + y0' + 1*g1(x1,x2)
    //! @param VALUE y0'
    //! @param SUM_FUNCTION function y = y0 + s1*g1(x1,x2) as reference to a BinarySumFunction
    //! @return BinarySumFunction f(x1,x2) = y0 + y0' + s1*g1(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType, typename t_Scalar>
    inline
    BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>
    operator +
    (
      const BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &SUM_FUNCTION,
      const t_ResultType &VALUE
    )
    {
      return BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>( SUM_FUNCTION).operator +=( VALUE);
    }

    //! @brief operator VALUE + BinarySumFunction
    //! create a function y = y0 + y0' + 1*g1(x1,x2)
    //! @param VALUE y0'
    //! @param SUM_FUNCTION function y = y0 + s1*g1(x1,x2) as reference to a BinarySumFunction
    //! @return BinarySumFunction f(x1,x2) = y0 + y0' + s1*g1(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType, typename t_Scalar>
    inline
    BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>
    operator +
    (
      const t_ResultType &VALUE,
      const BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &SUM_FUNCTION
    )
    {
      return BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>( SUM_FUNCTION).operator +=( VALUE);
    }

    //! @brief operator BinarySumFunction - VALUE
    //! create a function y = y0 - y0' + 1*g1(x1,x2)
    //! @param VALUE y0'
    //! @param SUM_FUNCTION function y = y0 + s1*g1(x1,x2) as reference to a BinarySumFunction
    //! @return BinarySumFunction f(x1,x2) = y0 - y0' + s1*g1(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType, typename t_Scalar>
    inline
    BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>
    operator -
    (
      const BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &SUM_FUNCTION,
      const t_ResultType &VALUE
    )
    {
      return BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>( SUM_FUNCTION).operator -=( VALUE);
    }

    //! @brief operator VALUE - BinarySumFunction
    //! create a function y = y0 - y0' + 1*g1(x1,x2)
    //! @param VALUE y0'
    //! @param SUM_FUNCTION function y = y0 + s1*g1(x1,x2) as reference to a BinarySumFunction
    //! @return BinarySumFunction f(x1,x2) = y0 - y0' + s1*g1(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType, typename t_Scalar>
    inline
    BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>
    operator -
    (
      const t_ResultType &VALUE,
      const BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &SUM_FUNCTION
    )
    {
      return BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>( -SUM_FUNCTION).operator +=( VALUE);
    }

    //! @brief operator BinarySumFunction * SCALAR
    //! create function f(x1,x2) = s*y = s*y0 + s*s1*g1(x1,x2)
    //! @param SUM_FUNCTION y = y0 + s1*g1(x1,x2) as reference to BinarySumFunction
    //! @param SCALAR scalar s
    //! @return reference to new sum function f(x1,x2) = s*y0 + s*s1*g1(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType, typename t_Scalar>
    inline
    BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>
    operator *
    (
      const BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &SUM_FUNCTION,
      const t_Scalar &SCALAR
    )
    {
      return BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>( SUM_FUNCTION).operator *=( SCALAR);
    }

    //! @brief operator SCALAR * BinarySumFunction
    //! create function f(x1,x2) = s*y = s*y0 + s*s1*g1(x1,x2)
    //! @param SUM_FUNCTION y = y0 + s1*g1(x1,x2) as reference to BinarySumFunction
    //! @param SCALAR scalar s
    //! @return reference to new sum function f(x1,x2) = s*y0 + s*s1*g1(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType, typename t_Scalar>
    inline
    BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>
    operator *
    (
      const t_Scalar &SCALAR,
      const BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &SUM_FUNCTION
    )
    {
      return BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>( SUM_FUNCTION).operator *=( SCALAR);
    }

    //! @brief operator BinarySumFunction / SCALAR
    //! create function f(x1,x2) = 1/s*y0 + 1/s*g(x1,x2)
    //! @param SUM_FUNCTION f(x1,x2) = y0 + s1*g1(x1,x2) as reference to BinarySumFunction
    //! @param SCALAR scalar s
    //! @return reference to new sum function f(x1,x2) = 1/s*y0 + 1/s*g1(x1,x2)
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType, typename t_Scalar>
    inline
    BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>
    operator /
    (
      const BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar> &SUM_FUNCTION,
      const t_Scalar &SCALAR
    )
    {
      return BinarySumFunction< t_ArgumentType1, t_ArgumentType2, t_ResultType, t_Scalar>( SUM_FUNCTION).operator /=( SCALAR);
    }

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_BINARY_SUM_FUNCTION_H_
