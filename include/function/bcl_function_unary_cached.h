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

#ifndef BCL_FUNCTION_UNARY_CACHED_H_
#define BCL_FUNCTION_UNARY_CACHED_H_

// include the namespace header
#include "bcl_function.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_function_unary_interface.h"
#include "signal/bcl_signal_signal.h"
#include "signal/bcl_signal_slots.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically
#include <list>
#include <map>

namespace bcl
{
  namespace function
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class UnaryCached
    //! @brief caches results for Arguments, if they have not been changed, but the function
    //! is evaluated multiple times.
    //! @details The UnaryCached class wraps a around a normal UnaryInterface derived class, and stores results
    //! if they are evaluated for the first time, and returns them every time a reevaluation is invoked.
    //!
    //! @tparam t_ArgumentType1 the argument type for the operator the function is invoked for. It is expected
    //!         that it has a GetDestructorSignal() that gives a signal handler, that emits on object destruction, so
    //!         that results can be removed from the cache if the argument objects do not exist anymore
    //! @tparam t_ResultType the type of the result, the function returns
    //!
    //! @see @link example_function_unary_cached.cpp @endlink
    //! @author woetzen
    //! @date Jun 4, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class UnaryCached :
      public UnaryInterface< t_ArgumentType, t_ResultType>,
      public signal::Slots
    {

    private:

    //////////
    // data //
    //////////

      //! @brief the actual function
      util::ShPtr< UnaryInterface< t_ArgumentType, t_ResultType> > m_Function;

      //! @brief the cache that stores for each Argument the result
      //! the key is the pointer (address) of the argument, the data is a pointer to a result, NULL if the result does
      //! not exist (because the argument has changed)
      typedef std::map< t_ArgumentType *, const t_ResultType *> CacheType;
      mutable CacheType m_Cache;

      //! @brief typedef for signal connector functions
      typedef signal::Signal1< const t_ArgumentType &> &( t_ArgumentType::*SignalHandlerFunctionPtr)() const;

      //! destructor signal handler which UnaryCached connects to on the Argument
      SignalHandlerFunctionPtr m_DestructorHandlerToConnectTo;

      //! set of signal handlers wo which UnaryCached connects to on the Argument
      std::list< SignalHandlerFunctionPtr> m_HandlersToConnectTo;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from function and a destructor signal handler
      //! @param SP_FUNCTION function that determines the result for an argument
      //! @param DESTRUCTOR_SIGNAL_HANDLER Destructor signal handler
      UnaryCached
      (
        const util::ShPtr< UnaryInterface< t_ArgumentType, t_ResultType> > &SP_FUNCTION,
        SignalHandlerFunctionPtr DESTRUCTOR_SIGNAL_HANDLER
      ) :
        m_Function( SP_FUNCTION),
        m_DestructorHandlerToConnectTo( DESTRUCTOR_SIGNAL_HANDLER)
      {
      }

      //! @brief copy constructor
      //! will not copy the actual cache and signal connection to already seen arguments
      //! @param FUNCTION_CACHED_RHS cache function to copy
      UnaryCached
      (
        const UnaryCached< t_ArgumentType, t_ResultType> &FUNCTION_CACHED_RHS
      ) :
        signal::Slots(),
        m_Function( FUNCTION_CACHED_RHS.m_Function),
        m_Cache(),
        m_DestructorHandlerToConnectTo( FUNCTION_CACHED_RHS.m_DestructorHandlerToConnectTo),
        m_HandlersToConnectTo( FUNCTION_CACHED_RHS.m_HandlersToConnectTo)
      {
      }

      //! @brief Clone function
      //! @return pointer to new UnaryCached< t_ArgumentType, t_ResultType>
      UnaryCached< t_ArgumentType, t_ResultType> *Clone() const
      {
        return new UnaryCached< t_ArgumentType, t_ResultType>( *this);
      }

      //! @brief destructor
      ~UnaryCached()
      {
        // copy the cache, so that deletion of results does not trigger any of the member functions that removes results
        // on the cache
        CacheType copy;
        m_Cache.swap( copy);

        // iterate over cache and delete all result pointer
        for( typename CacheType::iterator itr( copy.begin()), itr_end( copy.end()); itr != itr_end; ++itr)
        {
          // delete only if there is an actual object
          if( itr->second != NULL)
          {
            delete itr->second;
          }
        }
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

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const
      {
        return m_Function->GetScheme();
      }

      //! @brief access to function
      const util::ShPtr< UnaryInterface< t_ArgumentType, t_ResultType> > &GetFunction() const
      {
        return m_Function;
      }

      //! @brief cache size
      //! @return number of entries in cache
      size_t GetCacheSize() const
      {
        return m_Cache.size();
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief add a SignalHandler of t_ArgumentType to connect to, that triggers the result invalidation
      //! @param SIGNAL_HANDLER pointer to t_ArgumentType member function with signature "signal::Signal1< const t_ArgumentType&> &t_ArgumentType::SignalHandlerFunction() const"
      //! @return true, if successful, false if error occurs (e.i. the handler was added already)
      bool AddSignalHandlerForArgument( SignalHandlerFunctionPtr SIGNAL_HANDLER)
      {
        // check if this signal handler member function is already added
        if
        (
          std::find( m_HandlersToConnectTo.begin(), m_HandlersToConnectTo.end(), SIGNAL_HANDLER) ==
            m_HandlersToConnectTo.end()
        )
        {
          m_HandlersToConnectTo.push_back( SIGNAL_HANDLER);
          return true;
        }
        else
        {
          BCL_MessageCrt( "given Signal Handler was already added to cache function")
          return false;
        }
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator taking an ARGUMENT and returning a t_ResultType object
      //! If there is a result stored in the cache, return that
      //! @param ARGUMENT Argument to be used to evaluate the function
      //! @return function value of the given argument
      t_ResultType operator()( const t_ArgumentType &ARGUMENT) const
      {
        // search in cache
        typename CacheType::iterator itr( m_Cache.find( typename CacheType::key_type( &ARGUMENT)));

        // no result found
        if( itr == m_Cache.end())
        {
          // try to insert
          const std::pair< typename CacheType::iterator, bool> insert_pair
          (
            m_Cache.insert
            (
              typename CacheType::value_type
              (
                typename CacheType::key_type( &ARGUMENT),
                typename CacheType::mapped_type( new t_ResultType( m_Function->operator()( ARGUMENT)))
              )
            )
          );

          // check that insert was successful
          BCL_Assert( insert_pair.second, "could not insert the result!");

          // connect this to all signal handlers on ARGUMENT
          ConnectToAllSignalHandlers( ARGUMENT);

          // return the result
          return *insert_pair.first->second;
        }
        // result was found
        else
        {
          // check if result exists
          if( itr->second == NULL)
          {
            itr->second = new t_ResultType( m_Function->operator()( ARGUMENT));
          }

          // return the result
          return *itr->second;
        }
      }

      //! @brief assignment operator
      //! will not copy the actual cache and signal connection to already seen arguments
      //! @param FUNCTION_CACHED_RHS cache function to copy
      //! @return reference to this UnaryCached< t_ArgumentType, t_ResultType> object
      UnaryCached< t_ArgumentType, t_ResultType> &operator =
      (
        const UnaryCached< t_ArgumentType, t_ResultType> &FUNCTION_CACHED_RHS
      )
      {
        // only copy if argument is different than this
        if( this != &FUNCTION_CACHED_RHS)
        {
          // only assign the function and the HandlerToConnectListTo
          m_Function = FUNCTION_CACHED_RHS.m_Function;
          m_DestructorHandlerToConnectTo = FUNCTION_CACHED_RHS.m_DestructorHandlerToConnectTo;
          m_HandlersToConnectTo = FUNCTION_CACHED_RHS.m_HandlersToConnectTo;
        }

        return *this;
      }

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write detailed scheme and values to OSTREAM
      //! @param ARGUMENT Argument to be used to evaluate the function
      //! @param OSTREAM std::ostream to be written to
      //! @param FORMAT Format object to be used in output
      //! @return std::ostream which was written to
      std::ostream &WriteDetailedSchemeAndValues
      (
        const t_ArgumentType &ARGUMENT,
        std::ostream &OSTREAM,
        const util::Format &FORMAT = util::Format()
      ) const
      {
        return m_Function->WriteDetailedSchemeAndValues( ARGUMENT, OSTREAM, FORMAT);
      }

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // read member
        io::Serialize::Read( m_Function, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write member
        io::Serialize::Write( m_Function, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief remove results for this object from the cache
      //! @param ARGUMENT argument that should be removed
      void RemoveResultFromCache( const t_ArgumentType &ARGUMENT)
      {
        // find result for that Argument
        typename CacheType::iterator itr( m_Cache.find( typename CacheType::key_type( &ARGUMENT)));

        // skip if there is no result cached
        if( itr == m_Cache.end())
        {
          return;
        }

        // copy of result - if result is deleted before the itr is removed from cache, if could be that the deletion of
        // the result triggers another call to this function
        const t_ResultType *result( itr->second);

        // erase from cache
        m_Cache.erase( itr);

        // delete result
        delete result;
      }

      //! @brief delete a result if argument changes and function has to be reevaluated
      //! @param ARGUMENT the argument for which the result has to be forgotten
      void InvalidateResult( const t_ArgumentType &ARGUMENT)
      {
        // find result for that Argument
        typename CacheType::iterator itr( m_Cache.find( typename CacheType::key_type( &ARGUMENT)));

        // delete if there is a result cached
        if( itr != m_Cache.end() && itr->second != NULL)
        {
          // delete result
          const t_ResultType *result( itr->second);
          itr->second = NULL;
          delete result;
        }
      }

      //! @brief connect to all signals on argument
      //! @param ARGUMENT argument to which to connect to
      void ConnectToAllSignalHandlers( const t_ArgumentType &ARGUMENT) const
      {
        UnaryCached< t_ArgumentType, t_ResultType> *non_const_this_ptr
        (
          const_cast< UnaryCached< t_ArgumentType, t_ResultType> *>( this)
        );

        // register this to the destructor signal with RemoveResultFromCache function
        ( ARGUMENT.*m_DestructorHandlerToConnectTo)().Connect
        (
          non_const_this_ptr, &UnaryCached< t_ArgumentType, t_ResultType>::RemoveResultFromCache
        );

        // iterate over handlers cache has to connect to with InvalidateResult function
        for
        (
          typename std::list< SignalHandlerFunctionPtr>::const_iterator
            func_itr( m_HandlersToConnectTo.begin()), func_itr_end( m_HandlersToConnectTo.end());
          func_itr != func_itr_end;
          ++func_itr
        )
        {
          ( ARGUMENT.**func_itr)().Connect
          (
            non_const_this_ptr, &UnaryCached< t_ArgumentType, t_ResultType>::InvalidateResult
          );
        }
      }

    }; // template class UnaryCached

    // instantiate s_Instance
    template< typename t_ArgumentType, typename t_ResultType>
    const util::SiPtr< const util::ObjectInterface> UnaryCached< t_ArgumentType, t_ResultType>::s_Instance
    (
      GetObjectInstances().AddInstance
      (
        UnaryCached< t_ArgumentType, t_ResultType>( util::ShPtr< UnaryInterface< t_ArgumentType, t_ResultType> >())
      )
    );

  } // namespace function
} // namespace bcl

#endif // BCL_FUNCTION_UNARY_CACHED_H_ 
