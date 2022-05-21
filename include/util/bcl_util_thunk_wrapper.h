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

#ifndef BCL_UTIL_THUNK_WRAPPER_H_
#define BCL_UTIL_THUNK_WRAPPER_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_function_interface.h"
#include "bcl_util_function_interface_nonconst.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ThunkWrapper
    //! @brief This class wraps functions with no input and one output (thunks)
    //! @details Use this class to prepare to pass in a function to create a job for parallelization.
    //! It's assumed the non-static function belongs to typename t_FunctionClass
    //! Reason for separate wrapper:  thunks require FunctionInterface not BinaryFunctionInterface
    //!
    //! @see @link example_util_thunk_wrapper.cpp @endlink
    //! @author riddeljs
    //! @date 11.17.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_FunctionClass, typename t_ResultType>
    class ThunkWrapper :
      public FunctionInterfaceNonConst< t_FunctionClass, t_ResultType>
    {
    public:

      typedef t_ResultType ( t_FunctionClass::*PtrToMemberFunction)();
      typedef t_ResultType ( t_FunctionClass::*PtrToConstMemberFunction)() const;

    private:

      //! pointer to the function, belonging to class FunctionClass, we call for a job
      // Function, query, or static function?
      PtrToMemberFunction      m_Function;
      PtrToConstMemberFunction m_Query;
      t_ResultType ( *m_StaticFun)();
      size_t m_FunctionType;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! given a thunk static function to submit to the job, initialize instance of the wrapper.
      ThunkWrapper( t_ResultType ( *FUNCTION_TO_SUBMIT)()) :
        m_StaticFun( FUNCTION_TO_SUBMIT),
        m_FunctionType( 0)
      {
      }

      //! given a thunk query to submit to the job, initialize instance of the wrapper.
      ThunkWrapper( PtrToConstMemberFunction QUERY_TO_SUBMIT) :
        m_Query( QUERY_TO_SUBMIT),
        m_FunctionType( 1)
      {
      }

      //! given a thunk to submit to the job, initialize instance of the wrapper.
      ThunkWrapper( PtrToMemberFunction FUNCTION_TO_SUBMIT) :
        m_Function( FUNCTION_TO_SUBMIT),
        m_FunctionType( 2)
      {
      }

      //! virtual copy constructor
      ThunkWrapper *Clone() const
      {
        return new ThunkWrapper( *this);
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

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that implements the thunk returning a result
      //! @return the Result of the operation, even if result is void
      t_ResultType operator()( t_FunctionClass &OBJECT) const
      {
        switch( m_FunctionType)
        {
          case 0:
            return ( *m_StaticFun)();
          case 1:
            return ( OBJECT.*m_Query)();
          default:
            return ( OBJECT.*m_Function)();
        }
      }

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
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    }; // template class ThunkWrapper

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! Specialization of ThunkWrapper for const classes
    //!
    //! @author mendenjl
    //! @date   06/25/10
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_FunctionClass, typename t_ResultType>
    class ThunkWrapper< const t_FunctionClass, t_ResultType> :
      public FunctionInterface< t_FunctionClass, t_ResultType>
    {
    private:

      typedef t_ResultType ( t_FunctionClass::*PtrToConstMemberFunction)() const;

      //! pointer to the function, belonging to class FunctionClass, we call for a job
      PtrToConstMemberFunction m_Query;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! given a thunk query to submit to the job, initialize instance of the wrapper.
      ThunkWrapper( PtrToConstMemberFunction QUERY_TO_SUBMIT) :
        m_Query( QUERY_TO_SUBMIT)
      {
      }

      //! virtual copy constructor
      ThunkWrapper *Clone() const
      {
        return new ThunkWrapper( *this);
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

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that implements the thunk returning a result
      //! @return the Result of the operation, even if result is void
      t_ResultType operator()( const t_FunctionClass &OBJECT) const
      {
        return ( OBJECT.*m_Query)();
      }

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
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    }; // template class ThunkWrapper

  } // namespace util
} // namespace bcl

#endif //BCL_UTIL_THUNK_WRAPPER_H_

