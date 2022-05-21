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

#ifndef BCL_FUNCTION_BINARY_CACHED_H_
#define BCL_FUNCTION_BINARY_CACHED_H_

// include the namespace header
#include "bcl_function.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_function_binary_interface.h"
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
    //! @class BinaryCached
    //! @brief caches results for Arguments, if they have not been changed, but the function
    //! is evaluated multiple times.
    //! @details The BinaryCached class wraps a around a normal BinaryInterface derived class, and
    //! stores results if they are evaluated for the first time, and returns them every time a reevaluation is invoked.
    //!
    //! @tparam t_ArgumentType1 the first argument type for the operator the function is invoked for. It is expected
    //!         that it has a GetDestructorSignal() that gives a signal handler, that emits on object destruction, so
    //!         that results can be removed from the cache if the argument objects do not exist anymore
    //! @tparam t_ArgumentType2 the second argument type for the operator the function is invoked for. It is expected
    //!         that it has a GetDestructorSignal() that gives a signal handler, that emits on object destruction, so
    //!         that results can be removed from the cache if the argument objects do not exist anymore
    //! @tparam t_ResultType the type of the result, the function returns
    //!
    //! @see @link example_function_binary_cached.cpp @endlink
    //! @author woetzen
    //! @date Jun 4, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    class BinaryCached :
      public BinaryInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType>,
      public signal::Slots
    {

    private:

    //////////
    // data //
    //////////

      //! @brief the actual function
      util::ShPtr< BinaryInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> > m_Function;

      //! @brief the cache that stores for each Argument the result
      //! the key is the pointer (address) of the argument, the data is a pointer to a result, NULL if the result does
      //! not exist (because the argument has changed)
      typedef std::map< std::pair< const t_ArgumentType1 *, const t_ArgumentType2 *>, const t_ResultType *> CacheType;
      mutable CacheType m_Cache;

      //! @brief typedef for signal connector functions for t_ArgumentType1
      typedef signal::Signal1< const t_ArgumentType1 &> &( t_ArgumentType1::*SignalHandlerFunctionPtr1)() const;

      //! @brief typedef for signal connector functions for t_ArgumentType2
      typedef signal::Signal1< const t_ArgumentType2 &> &( t_ArgumentType2::*SignalHandlerFunctionPtr2)() const;

      //! destructor signal handler which BinaryCached connects to on the t_ArgumentType1
      SignalHandlerFunctionPtr1 m_DestructorHandlerToConnectTo1;

      //! destructor signal handler which BinaryCached connects to on the t_ArgumentType2
      SignalHandlerFunctionPtr2 m_DestructorHandlerToConnectTo2;

      //! set of signal handlers  which BinaryCached connects to on the t_ArgumentType1
      std::list< SignalHandlerFunctionPtr1> m_HandlersToConnectTo1;

      //! set of signal handlers which BinaryCached connects to on the t_ArgumentType2
      std::list< SignalHandlerFunctionPtr2> m_HandlersToConnectTo2;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from function
      //! @param SP_FUNCTION function that determines the result for an argument
      //! @param DESTRUCTOR_SIGNAL_1 Destructor signal handler for first argument type
      //! @param DESTRUCTOR_SIGNAL_2 Destructor signal handler for second argument type
      BinaryCached
      (
        const util::ShPtr< BinaryInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> > &SP_FUNCTION,
        SignalHandlerFunctionPtr1 DESTRUCTOR_SIGNAL_1,
        SignalHandlerFunctionPtr2 DESTRUCTOR_SIGNAL_2
      ) :
        m_Function( SP_FUNCTION),
        m_DestructorHandlerToConnectTo1( DESTRUCTOR_SIGNAL_1),
        m_DestructorHandlerToConnectTo2( DESTRUCTOR_SIGNAL_2)
      {
      }

      //! @brief copy constructor
      //! will not copy the actual cache and signal connection to already seen arguments
      //! @param FUNCTION_CACHED_RHS cache function to copy
      BinaryCached
      (
        const BinaryCached< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION_CACHED_RHS
      ) :
        signal::Slots(),
        m_Function( FUNCTION_CACHED_RHS.m_Function),
        m_Cache(),
        m_DestructorHandlerToConnectTo1( FUNCTION_CACHED_RHS.m_DestructorHandlerToConnectTo1),
        m_DestructorHandlerToConnectTo2( FUNCTION_CACHED_RHS.m_DestructorHandlerToConnectTo2),
        m_HandlersToConnectTo1( FUNCTION_CACHED_RHS.m_HandlersToConnectTo1),
        m_HandlersToConnectTo2( FUNCTION_CACHED_RHS.m_HandlersToConnectTo2)
      {
      }

      //! @brief Clone function
      //! @return pointer to new BinaryCached< t_ArgumentType1, t_ArgumentType2, t_ResultType>
      BinaryCached< t_ArgumentType1, t_ArgumentType2, t_ResultType> *Clone() const
      {
        return new BinaryCached< t_ArgumentType1, t_ArgumentType2, t_ResultType>( *this);
      }

      //! @brief destructor
      ~BinaryCached()
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
      const util::ShPtr< BinaryInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> > &GetFunction() const
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

      //! @brief add a SignalHandler of t_ArgumentType1 to connect to, that triggers the result invalidation
      //! @param SIGNAL_HANDLER1 pointer to t_ArgumentType1 member function with signature "signal::Signal1< const t_ArgumentType1&> &t_ArgumentType1::SignalHandlerFunction() const"
      //! @return true, if successful, false if error occurs (e.i. the handler was added already)
      bool AddSignalHandlerForArgument( SignalHandlerFunctionPtr1 SIGNAL_HANDLER1)
      {
        // check if this signal handler member function is already added
        if
        (
          std::find( m_HandlersToConnectTo1.begin(), m_HandlersToConnectTo1.end(), SIGNAL_HANDLER1) ==
            m_HandlersToConnectTo1.end()
        )
        {
          m_HandlersToConnectTo1.push_back( SIGNAL_HANDLER1);
          return true;
        }
        else
        {
          return false;
        }
      }

      //! @brief add a SignalHandler of t_ArgumentType2 to connect to, that triggers the result invalidation
      //! @param SIGNAL_HANDLER2 pointer to t_ArgumentType2 member function with signature "signal::Signal1< const t_ArgumentType2&> &t_ArgumentType2::SignalHandlerFunction() const"
      //! @return true, if successful, false if error occurs (e.i. the handler was added already)
      bool AddSignalHandlerForArgument( SignalHandlerFunctionPtr2 SIGNAL_HANDLER2)
      {
        // check if this signal handler member function is already added
        if
        (
          std::find( m_HandlersToConnectTo2.begin(), m_HandlersToConnectTo2.end(), SIGNAL_HANDLER2) ==
            m_HandlersToConnectTo2.end()
        )
        {
          m_HandlersToConnectTo2.push_back( SIGNAL_HANDLER2);
          return true;
        }
        else
        {
          return false;
        }
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator taking an ARGUMENT and returning a t_ResultType object
      //! If there is a result stored in the cache, return that
      //! @param ARGUMENT1 Argument1 to be used to evaluate the function
      //! @param ARGUMENT2 Argument2 to be used to evaluate the function
      //! @return function value of the given argument
      t_ResultType operator()( t_ArgumentType1 &ARGUMENT1, t_ArgumentType2 &ARGUMENT2) const
      {
        // search in cache
        typename CacheType::iterator itr( m_Cache.find( typename CacheType::key_type( &ARGUMENT1, &ARGUMENT2)));

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
                std::make_pair( &ARGUMENT1, &ARGUMENT2),
                new t_ResultType( m_Function->operator()( ARGUMENT1, ARGUMENT2))
              )
            )
          );

          // check that insert was successful
          BCL_Assert( insert_pair.second, "could not insert the result!");

          // connect this to all signal handlers on ARGUMENT1 and ARGUMENT 2
          ConnectToAllSignalHandlers1( ARGUMENT1);
          ConnectToAllSignalHandlers2( ARGUMENT2);

          // return the result
          return *insert_pair.first->second;
        }
        // result was found
        else
        {
          // check if result exists
          if( itr->second == NULL)
          {
            itr->second = new t_ResultType( m_Function->operator()( ARGUMENT1, ARGUMENT2));
          }

          // return the result
          return *itr->second;
        }
      }

      //! @brief assignment operator
      //! will not copy the actual cache and signal connection to already seen arguments
      //! @param FUNCTION_CACHED_RHS cache function to copy
      //! @return reference to this BinaryCached< t_ArgumentType, t_ResultType> object
      BinaryCached< t_ArgumentType1, t_ArgumentType2, t_ResultType> &operator =
      (
        const BinaryCached< t_ArgumentType1, t_ArgumentType2, t_ResultType> &FUNCTION_CACHED_RHS
      )
      {
        // only copy if argument is different than this
        if( this != &FUNCTION_CACHED_RHS)
        {
          // only assign the function and the HandlerToConnectListTo
          m_Function = FUNCTION_CACHED_RHS.m_Function;
          m_DestructorHandlerToConnectTo1 = FUNCTION_CACHED_RHS.m_DestructorHandlerToConnectTo1;
          m_HandlersToConnectTo1 = FUNCTION_CACHED_RHS.m_HandlersToConnectTo1;
          m_DestructorHandlerToConnectTo2 = FUNCTION_CACHED_RHS.m_DestructorHandlerToConnectTo2;
          m_HandlersToConnectTo2 = FUNCTION_CACHED_RHS.m_HandlersToConnectTo2;
        }

        return *this;
      }

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write detailed scheme and values to OSTREAM
      //! @param ARGUMENT1 First argument to be used to evaluate the function
      //! @param ARGUMENT2 Second argument to be used to evaluate the function
      //! @param OSTREAM std::ostream to be written to
      //! @param FORMAT Format object to be used in output
      //! @return std::ostream which was written to
      std::ostream &WriteDetailedSchemeAndValues
      (
        const t_ArgumentType1 &ARGUMENT1,
        const t_ArgumentType2 &ARGUMENT2,
        std::ostream &OSTREAM,
        const util::Format &FORMAT = util::Format()
      ) const
      {
        return m_Function->WriteDetailedSchemeAndValues( ARGUMENT1, ARGUMENT2, OSTREAM, FORMAT);
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
      //! @param ARGUMENT of type t_ArgumentType1 which was connected to through the GetDestructorSignal()
      void RemoveResultFromCache1( const t_ArgumentType1 &ARGUMENT)
      {
        // cast to t_ArgumentType1 assuming that DESTRUCTOr object is
        const t_ArgumentType1 *ptr1( &ARGUMENT);

        // collect all results to delete - otherwise deleting a result could trigger another call to this function
        // causing the iterator on the map to be invalidated
        std::list< typename CacheType::mapped_type> results_to_delete;

        // find result for that Argument
        for( typename CacheType::iterator itr( m_Cache.begin()), itr_end( m_Cache.end()); itr != itr_end;)
        {
          // skip results that are not associated with that ptr1
          if( itr->first.first != ptr1)
          {
            ++itr;
            continue;
          }

          // mark result for erase and remove iterator from cache
          typename CacheType::iterator itr_to_erase( itr);
          ++itr;

          // copy of result - if result is deleted before the itr is removed from cache, if could be that the deletion of
          // the result triggers another call to this function
          results_to_delete.push_back( itr_to_erase->second);

          // erase from cache
          m_Cache.erase( itr_to_erase);
        }

        // delete all results that need deletion
        for
        (
          typename std::list< typename CacheType::mapped_type>::iterator
            itr( results_to_delete.begin()), itr_end( results_to_delete.end());
          itr != itr_end;
          ++itr
        )
        {
          delete *itr;
        }
      }

      //! @brief remove results for this object from the cache
      //! @param ARGUMENT of type t_ArgumentType2 which was connected to through the GetDestructorSignal()
      void RemoveResultFromCache2( const t_ArgumentType2 &ARGUMENT)
      {
        // cast to t_ArgumentType2 assuming that DESTRUCTOr object is
        const t_ArgumentType2 *ptr2( &ARGUMENT);

        // collect all results to delete - otherwise deleting a result could trigger another call to this function
        // causing the iterator on the map to be invalidated
        std::list< typename CacheType::mapped_type> results_to_delete;

        // find result for that Argument
        for( typename CacheType::iterator itr( m_Cache.begin()), itr_end( m_Cache.end()); itr != itr_end;)
        {
          // skip results that are not associated with that ptr2
          if( itr->first.second != ptr2)
          {
            ++itr;
            continue;
          }

          // erase from cache
          typename CacheType::iterator itr_to_erase( itr);
          ++itr;

          // copy of result - if result is deleted before the itr is removed from cache, if could be that the deletion of
          // the result triggers another call to this function
          results_to_delete.push_back( itr_to_erase->second);

          // delete entry from map
          m_Cache.erase( itr_to_erase);
        }

        // delete all results that need deletion
        for
        (
          typename std::list< typename CacheType::mapped_type>::iterator
            itr( results_to_delete.begin()), itr_end( results_to_delete.end());
          itr != itr_end;
          ++itr
        )
        {
          delete *itr;
        }
      }

      //! @brief delete a result if argument changes and function has to be reevaluated
      //! @param ARGUMENT the argument for which the result has to be forgotten
      void InvalidateResult1( const t_ArgumentType1 &ARGUMENT)
      {
        // cast to t_ArgumentType1 assuming that DESTRUCTOr object is
        const t_ArgumentType1 *ptr1( &ARGUMENT);

        // find result for that Argument
        for( typename CacheType::iterator itr( m_Cache.begin()), itr_end( m_Cache.end()); itr != itr_end; ++itr)
        {
          // skip results that are not associated with that argument
          if( itr->first.first == ptr1 && itr->second != NULL)
          {
            // delete result
            const t_ResultType *result( itr->second);
            itr->second = NULL;
            delete result;
          }
        }
      }

      //! @brief delete a result if argument changes and function has to be reevaluated
      //! @param ARGUMENT the argument for which the result has to be forgotten
      void InvalidateResult2( const t_ArgumentType2 &ARGUMENT)
      {
        // cast to t_ArgumentType2 assuming that destructor object is
        const t_ArgumentType2 *ptr2( &ARGUMENT);

        // find result for that Argument
        for( typename CacheType::iterator itr( m_Cache.begin()), itr_end( m_Cache.end()); itr != itr_end; ++itr)
        {
          // skip results that are not associated with that argument
          if( itr->first.second == ptr2 && itr->second != NULL)
          {
            // delete result
            const t_ResultType *result( itr->second);
            itr->second = NULL;
            delete result;
          }
        }
      }

      //! @brief connect to all signals on argument
      //! @param ARGUMENT argument to which to connect to
      void ConnectToAllSignalHandlers1( const t_ArgumentType1 &ARGUMENT) const
      {
        BinaryCached< t_ArgumentType1, t_ArgumentType2, t_ResultType> *non_const_this_ptr
        (
          const_cast< BinaryCached< t_ArgumentType1, t_ArgumentType2, t_ResultType> *>( this)
        );

        // register this to the desructor signal with RemoveResultFromCache function
        ( ARGUMENT.*m_DestructorHandlerToConnectTo1)().Connect
        (
          non_const_this_ptr, &BinaryCached< t_ArgumentType1, t_ArgumentType2, t_ResultType>::RemoveResultFromCache1
        );

        // iterate over handlers cache has to conncet to with InvalidateResult function
        for
        (
          typename std::list< SignalHandlerFunctionPtr1>::const_iterator
            func_itr( m_HandlersToConnectTo1.begin()), func_itr_end( m_HandlersToConnectTo1.end());
          func_itr != func_itr_end;
          ++func_itr
        )
        {
          (ARGUMENT.**func_itr)().Connect
          (
            non_const_this_ptr, &BinaryCached< t_ArgumentType1, t_ArgumentType2, t_ResultType>::InvalidateResult1
          );
        }
      }

      //! @brief connect to all signals on argument
      //! @param ARGUMENT argument to which to connect to
      void ConnectToAllSignalHandlers2( const t_ArgumentType2 &ARGUMENT) const
      {
        BinaryCached< t_ArgumentType1, t_ArgumentType2, t_ResultType> *non_const_this_ptr
        (
          const_cast< BinaryCached< t_ArgumentType1, t_ArgumentType2, t_ResultType> *>( this)
        );

        // register this to the desructor signal with RemoveResultFromCache function
        ( ARGUMENT.*m_DestructorHandlerToConnectTo2)().Connect
        (
          non_const_this_ptr, &BinaryCached< t_ArgumentType1, t_ArgumentType2, t_ResultType>::RemoveResultFromCache2
        );

        // iterate over handlers cache has to conncet to with InvalidateResult function
        for
        (
          typename std::list< SignalHandlerFunctionPtr2>::const_iterator
            func_itr( m_HandlersToConnectTo2.begin()), func_itr_end( m_HandlersToConnectTo2.end());
          func_itr != func_itr_end;
          ++func_itr
        )
        {
          (ARGUMENT.**func_itr)().Connect
          (
            non_const_this_ptr, &BinaryCached< t_ArgumentType1, t_ArgumentType2, t_ResultType>::InvalidateResult2
          );
        }
      }

    }; // template class BinaryCached

    // instantiate s_Instance
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    const util::SiPtr< const util::ObjectInterface> BinaryCached< t_ArgumentType1, t_ArgumentType2, t_ResultType>::s_Instance
    (
      GetObjectInstances().AddInstance( new BinaryCached< t_ArgumentType1, t_ArgumentType2, t_ResultType>( util::ShPtr< BinaryInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> >()))
    );

    //! @brief partial specialization of BinaryCached for the special case t_ArgumentType1 = t_ArgumentType2
    template< typename t_ArgumentType, typename t_ResultType>
    class BinaryCached< t_ArgumentType, t_ArgumentType, t_ResultType> :
      public BinaryInterface< t_ArgumentType, t_ArgumentType, t_ResultType>,
      public signal::Slots
    {

    private:

    //////////
    // data //
    //////////

      //! @brief the actual function
      util::ShPtr< BinaryInterface< t_ArgumentType, t_ArgumentType, t_ResultType> > m_Function;

      //! @brief is the function symmetric - true if Arg1 and Arg2 can be swapped without changing the result
      bool m_Symmetric;

      //! @brief the cache that stores for each Argument the result
      //! the key is the pointer (address) of the argument, the data is a pointer to a result, NULL if the result does
      //! not exist (because the argument has changed)
      typedef std::map< std::pair< const t_ArgumentType *, const t_ArgumentType *>, const t_ResultType *> CacheType;
      mutable CacheType m_Cache;

      //! @brief typedef for signal connector functions for t_ArgumentType
      typedef signal::Signal1< const t_ArgumentType &> &( t_ArgumentType::*SignalHandlerFunctionPtr)() const;

      //! destructor signal handler which BinaryCached connects to on the t_ArgumentType
      SignalHandlerFunctionPtr m_DestructorHandlerToConnectTo;

      //! set of signal handlers which BinaryCached connects to on the t_ArgumentType
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

      //! @brief construct from function
      //! @param SP_FUNCTION function that determines the result for an argument
      //! @param DESTRUCTOR_SIGNAL_HANDLER Destructor signal handler function pointer
      //! @param SYMMETRIC true, if argument 1 and 2 to the operator can be exchanged without changing the result
      BinaryCached
      (
        const util::ShPtr< BinaryInterface< t_ArgumentType, t_ArgumentType, t_ResultType> > &SP_FUNCTION,
        SignalHandlerFunctionPtr DESTRUCTOR_SIGNAL_HANDLER,
        const bool SYMMETRIC
      ) :
        m_Function( SP_FUNCTION),
        m_Symmetric( SYMMETRIC),
        m_DestructorHandlerToConnectTo( DESTRUCTOR_SIGNAL_HANDLER)
      {
      }

      //! @brief copy constructor
      //! will not copy the actual cache and signal connection to already seen arguments
      //! @param FUNCTION_CACHED_RHS cache function to copy
      BinaryCached
      (
        const BinaryCached< t_ArgumentType, t_ArgumentType, t_ResultType> &FUNCTION_CACHED_RHS
      ) :
        signal::Slots(),
        m_Function( FUNCTION_CACHED_RHS.m_Function),
        m_Symmetric( FUNCTION_CACHED_RHS.m_Symmetric),
        m_Cache(),
        m_DestructorHandlerToConnectTo( FUNCTION_CACHED_RHS.m_DestructorHandlerToConnectTo),
        m_HandlersToConnectTo( FUNCTION_CACHED_RHS.m_HandlersToConnectTo)
      {
      }

      //! @brief Clone function
      //! @return pointer to new BinaryCached< t_ArgumentType, t_ArgumentType2, t_ResultType>
      BinaryCached< t_ArgumentType, t_ArgumentType, t_ResultType> *Clone() const
      {
        return new BinaryCached< t_ArgumentType, t_ArgumentType, t_ResultType>( *this);
      }

      //! @brief destructor
      ~BinaryCached()
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
      const util::ShPtr< BinaryInterface< t_ArgumentType, t_ArgumentType, t_ResultType> > &GetFunction() const
      {
        return m_Function;
      }

      //! @brief cache size
      //! @return number of entries in cache
      size_t GetCacheSize() const
      {
        return m_Cache.size();
      }

      //! @brief is symmetric
      //! @return if given function is symmetric
      bool IsSymmetric() const
      {
        return m_Symmetric;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief add a SignalHandler of t_ArgumentType to connect to, that triggers the result invalidation
      //! @param SIGNAL_HANDLER1 pointer to t_ArgumentType member function with signature "signal::Signal1< const t_ArgumentType&> &t_ArgumentType::SignalHandlerFunction() const"
      //! @return true, if successful, false if error occurs (e.i. the handler was added already)
      bool AddSignalHandlerForArgument( SignalHandlerFunctionPtr SIGNAL_HANDLER1)
      {
        // check if this signal handler member function is already added
        if
        (
          std::find( m_HandlersToConnectTo.begin(), m_HandlersToConnectTo.end(), SIGNAL_HANDLER1) ==
            m_HandlersToConnectTo.end()
        )
        {
          m_HandlersToConnectTo.push_back( SIGNAL_HANDLER1);
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
      //! @param ARGUMENT1 Argument1 to be used to evaluate the function
      //! @param ARGUMENT2 Argument2 to be used to evaluate the function
      //! @return function value of the given argument
      t_ResultType operator()( t_ArgumentType &ARGUMENT1, t_ArgumentType &ARGUMENT2) const
      {
        // pair of argument pointers
        typename CacheType::key_type arg_pair( &ARGUMENT1, &ARGUMENT2);

        // order pair by address, if symmetric
        if( m_Symmetric)
        {
          if( arg_pair.first > arg_pair.second)
          {
            std::swap( arg_pair.first, arg_pair.second);
          }
        }

        // search in cache
        typename CacheType::iterator itr( m_Cache.find( arg_pair));

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
                arg_pair,
                new t_ResultType( m_Function->operator()( ARGUMENT1, ARGUMENT2))
              )
            )
          );

          // check that insert was successful
          BCL_Assert( insert_pair.second, "could not insert the result!");

          // connect this to all signal handlers on ARGUMENT1 and ARGUMENT2
          ConnectToAllSignalHandlers( ARGUMENT1);
          if( &ARGUMENT1 != &ARGUMENT2)
          {
            ConnectToAllSignalHandlers( ARGUMENT2);
          }

          // return the result
          return *insert_pair.first->second;
        }
        // result was found
        else
        {
          // check if result exists
          if( itr->second == NULL)
          {
            itr->second = new t_ResultType( m_Function->operator()( ARGUMENT1, ARGUMENT2));
          }

          // return the result
          return *itr->second;
        }
      }

      //! @brief assignment operator
      //! will not copy the actual cache and signal connection to already seen arguments
      //! @param FUNCTION_CACHED_RHS cache function to copy
      //! @return reference to this BinaryCached< t_ArgumentType, t_ResultType> object
      BinaryCached< t_ArgumentType, t_ArgumentType, t_ResultType> &operator =
      (
        const BinaryCached< t_ArgumentType, t_ArgumentType, t_ResultType> &FUNCTION_CACHED_RHS
      )
      {
        // only copy if argument is different than this
        if( this != &FUNCTION_CACHED_RHS)
        {
          // only assign the function and the HandlerToConnectListTo
          m_Function = FUNCTION_CACHED_RHS.m_Function;
          m_HandlersToConnectTo = FUNCTION_CACHED_RHS.m_HandlersToConnectTo;
        }

        return *this;
      }

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write detailed scheme and values to OSTREAM
      //! @param ARGUMENT1 First argument to be used to evaluate the function
      //! @param ARGUMENT2 Second argument to be used to evaluate the function
      //! @param OSTREAM std::ostream to be written to
      //! @param FORMAT Format object to be used in output
      //! @return std::ostream which was written to
      std::ostream &WriteDetailedSchemeAndValues
      (
        const t_ArgumentType &ARGUMENT1,
        const t_ArgumentType &ARGUMENT2,
        std::ostream &OSTREAM,
        const util::Format &FORMAT = util::Format()
      ) const
      {
        return m_Function->WriteDetailedSchemeAndValues( ARGUMENT1, ARGUMENT2, OSTREAM, FORMAT);
      }

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // read member
        io::Serialize::Read( m_Function, ISTREAM);
        io::Serialize::Read( m_Symmetric, ISTREAM);

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
        io::Serialize::Write( m_Function, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Symmetric, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief remove results for this object from the cache
      //! @param ARGUMENT of type t_ArgumentType which was connected to through the GetDestructorSignal()
      void RemoveResultFromCache( const t_ArgumentType &ARGUMENT)
      {
        // cast to t_ArgumentType
        const t_ArgumentType *ptr( &ARGUMENT);

        // collect all results to delete - otherwise deleting a result could trigger another call to this function
        // causing the iterator on the map to be invalidated
        std::list< typename CacheType::mapped_type> results_to_delete;

        // find result for that Argument
        for( typename CacheType::iterator itr( m_Cache.begin()), itr_end( m_Cache.end()); itr != itr_end;)
        {
          // address pair
          const typename CacheType::key_type &address_pair( itr->first);

          // skip results that are not associated with that ptr
          if( address_pair.first != ptr && address_pair.second != ptr)
          {
            ++itr;
            continue;
          }

          // erase from cache, first make copy, increment the current iterator, and delete the copy
          typename CacheType::iterator itr_to_erase( itr);
          ++itr;

          // copy of result - if result is deleted before the itr is removed from cache, if could be that the deletion of
          // the result triggers another call to this function
          results_to_delete.push_back( itr_to_erase->second);

          // erase from cache
          m_Cache.erase( itr_to_erase);
        }

        // delete all results that need deletion
        for
        (
          typename std::list< typename CacheType::mapped_type>::iterator
            itr( results_to_delete.begin()), itr_end( results_to_delete.end());
          itr != itr_end;
          ++itr
        )
        {
          delete *itr;
        }
      }

      //! @brief delete a result if argument changes and function has to be reevaluated
      //! @param ARGUMENT the argument for which the result has to be forgotten
      void InvalidateResult( const t_ArgumentType &ARGUMENT)
      {
        // ptr to argument
        const t_ArgumentType *ptr( &ARGUMENT);

        // find result for that Argument
        for( typename CacheType::iterator itr( m_Cache.begin()), itr_end( m_Cache.end()); itr != itr_end; ++itr)
        {
          // skip results that are not associated with that argument
          if( ( itr->first.first == ptr || itr->first.second == ptr) && itr->second != NULL)
          {
            // delete result
            const t_ResultType *result( itr->second);
            itr->second = NULL;
            delete result;
          }
        }
      }

      //! @brief connect to all signals on argument
      //! @param ARGUMENT argument to which to connect to
      void ConnectToAllSignalHandlers( const t_ArgumentType &ARGUMENT) const
      {
        BinaryCached< t_ArgumentType, t_ArgumentType, t_ResultType> *non_const_this_ptr
        (
          const_cast< BinaryCached< t_ArgumentType, t_ArgumentType, t_ResultType> *>( this)
        );

        // register this to the desructor signal with RemoveResultFromCache function
        ( ARGUMENT.*m_DestructorHandlerToConnectTo)().Connect
        (
          non_const_this_ptr,
          &BinaryCached< t_ArgumentType, t_ArgumentType, t_ResultType>::RemoveResultFromCache
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
          (ARGUMENT.**func_itr)().Connect
          (
            non_const_this_ptr, &BinaryCached< t_ArgumentType, t_ArgumentType, t_ResultType>::InvalidateResult
          );
        }
      }

    }; // template class BinaryCached

  } // namespace function
} // namespace bcl

#endif // BCL_FUNCTION_BINARY_CACHED_H_ 
