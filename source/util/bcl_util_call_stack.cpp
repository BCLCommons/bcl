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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "util/bcl_util_call_stack.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically
#include <sstream>

#if __GNUC__ && !defined(_WIN32) && !defined(__APPLE__) && !defined(__ANDROID__) && !defined(NOCALLSTACK)
#include <cxxabi.h>
#include <dlfcn.h>
#include <execinfo.h>
#include <stdio.h>
#include <stdlib.h>
#endif // __GNUC__

#define MAX_DEPTH 32

namespace bcl
{
  namespace util
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor
    CallStack::Entry::Entry( const std::string &FILE_NAME, const size_t LINE_NUMBER, const std::string &FUNCTION) :
      m_File( FILE_NAME),
      m_LineNumber( LINE_NUMBER),
      m_Function( FUNCTION)
    {
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief return entry as string
    std::string CallStack::Entry::String() const
    {
      std::ostringstream os;
      os << m_File;
      if( IsDefined( m_LineNumber))
      {
        os << " (" << m_LineNumber << ")";
      }
      os << ": " << m_Function;

      return os.str();
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

#if __GNUC__ && !defined(_WIN32) && !defined(__APPLE__) && !defined(__ANDROID__) && !defined(NOCALLSTACK)

    //! @brief constructor
    //! @param DISCARD number of stack entries to discard at the top
    CallStack::CallStack( const size_t DISCARD)
    {
      // retrieve call-stack
      void *trace[ MAX_DEPTH];
      int stack_depth = backtrace( trace, MAX_DEPTH);

      for( int i = DISCARD + 1; i < stack_depth; i++)
      {
        Dl_info dlinfo;
        if( !dladdr( trace[ i], &dlinfo))
        {
          break;
        }

        const char *symname = dlinfo.dli_sname;

        int    status;
        char *demangled = abi::__cxa_demangle( symname, NULL, 0, &status);
        if( status == 0 && demangled)
        {
          symname = demangled;
        }

        // store entry to stack
        if( dlinfo.dli_fname && symname)
        {
          m_Stack.push_back
          (
            Entry
            (
              dlinfo.dli_fname,
              GetUndefined< size_t>(), // unsupported
              symname
            )
          );
        }
        else
        {
          // skip last entries below main
          continue;
        }

        if( demangled)
        {
          free( demangled);
        }
      }
    }
#else

    //! @brief constructor
    //! @param DISCARD number of stack entries to discard at the top
    CallStack::CallStack( const size_t DISCARD)
    {
    }

#endif // __GNUC__

  ////////////////
  // operations //
  ////////////////

    //! @brief convert the stacktrace to a string
    std::string CallStack::String() const
    {
      std::stringstream os;
#if __GNUC__ && !defined(_WIN32) && !defined(__APPLE__) && !defined(__ANDROID__) && !defined(NOCALLSTACK)
      if( m_Stack.empty())
      {
        os << "Link with '-rdynamic' for call stack information!\n";
      }
      else
      {
        os << "Call Stack\n";
      }
      for( std::vector< Entry>::const_iterator itr( m_Stack.begin()), itr_end( m_Stack.end()); itr != itr_end; ++itr)
      {
        os << itr->String() << '\n';
      }
#endif

      return os.str();
    }

  } // namespace util
} // namespace bcl
#undef MAX_DEPTH
