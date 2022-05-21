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
#include "util/bcl_util_assert.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_call_stack.h"
#include "util/bcl_util_loggers.h"
#include "util/bcl_util_runtime_environment_interface.h"

// external includes - sorted alphabetically
#include <stdlib.h>

namespace bcl
{
  namespace util
  {

    //! an error message is output to the logger and program is terminated
    //! @param ERRORMSG message describing the error
    //! @param ERROR_CODE the error code that should be passed to exit
    //! @param FILE_NAME the filename the assert was raised
    //! @param LINE_NUMBER the line number in the file
    //! @param FUNCTION_NAME the name of the function, the assert was raised
    //! @param DISPLAY_CALLSTACK whether to display the call stack, if possible
    void Assert::Exit
    (
      const std::string &ERRORMSG,
      const int ERROR_CODE,
      const char *FILE_NAME,
      const int LINE_NUMBER,
      const char *FUNCTION_NAME,
      const bool &DISPLAY_CALLSTACK
    )
    {
      // remove path from filename
      std::string file_name( FILE_NAME);
      file_name = file_name.substr( file_name.find_last_of( PATH_SEPARATOR) + 1);

      // if exit code does not indicate an error
      if( !DISPLAY_CALLSTACK)
      {
        if( !ERRORMSG.empty())
        {
          GetLogger() << "===> " << ERRORMSG << std::endl;
        }
      }
      else
      {
        std::cerr << "BCL Error: "   << ERROR_CODE
                  << " | function: " << FUNCTION_NAME
                  << " | file: "     << file_name
                  << " | line: "     << LINE_NUMBER
                  << "\n===> "       << ERRORMSG << std::endl;

        // stack till here
        const CallStack stack( 1);

        // output error
        std::cerr                    << stack.String();

        // if the logger is not the default logger, write the error out to the logger as well
        if( &GetLogger() != &**GetLoggers().e_Default)
        {
          // output error also to logger, in case std::cerr is not being monitored
          GetLogger() << "BCL Error: "   << ERROR_CODE
                      << " | function: " << FUNCTION_NAME
                      << " | file: "     << file_name
                      << " | line: "     << LINE_NUMBER
                      << "\n===> "       << ERRORMSG << std::endl;
        }
      }

      // abort program
      GetRuntimeEnvironment().Finalize( ERROR_CODE);
      exit( ERROR_CODE);
    }

  } // namespace util
} // namespace bcl
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
#include "util/bcl_util_class_descriptor.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_name_standardizer.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    //! @brief Standardizes the output of __PRETTY_FUNCTION__ between compilers / architectures
    //! This function should only be called from a (non-member) template function (like GetStaticClassName)
    //! @param PRETTY_FUNCTION_NAME should be passed from the  '__PRETTY_FUNCTION__' macro
    //! @return the standardized class, struct, union, or enum name as std::string
    std::string ClassNameViaTemplatedFunction( const std::string &PRETTY_FUNCTION_NAME)
    {
      return ClassNameStandardizer::StandardizeTemplateFunctionParameter( PRETTY_FUNCTION_NAME);
    }

    //! @brief extracts the class name
    //! The '__PRETTY_FUNCTION__' macro provided by the compiler should be used.
    //! If the Line looks like "virtual const std::string& bcl::Test< T>::GetName() const [with T = double]" the
    //! string "double" will be extracted
    //! @param PRETTY_FUNCTION_NAME should be passed from the  '__PRETTY_FUNCTION__' macro
    //! @return the Class name as std::string
    std::string StandardizeClassName( const std::string &PRETTY_FUNCTION_NAME)
    {
      return ClassNameStandardizer::Standardize( PRETTY_FUNCTION_NAME);
    }

    //! @brief extracts the namespace name before the "::"
    //! The '__PRETTY_FUNCTION__' macro provided by the compiler should be used.
    //! Example: ExtractNamespaceIdentifier( "virtual const std::string& bcl::Test::GetName() const) " the string "bcl::Test" will be extracted
    //! Algorithm: Finds the first '(' after the first
    //! @param PRETTY_FUNCTION_NAME should be passed from the  '__PRETTY_FUNCTION__' macro
    //! @return the Class name as std::string
    const std::string ExtractNamespaceIdentifier( const std::string &PRETTY_FUNCTION_NAME)
    {
      // standardize name first
      const std::string standardized_name( StandardizeClassName( PRETTY_FUNCTION_NAME));

      // find the first :: in the string
      const size_t first_scope( standardized_name.find( "::"));
      if( first_scope == std::string::npos)
      {
        return standardized_name; // no scope; the user probably passed a string that was already a namespace
      }

      // find the first of (, and < after the first :: in the string.  These cannot be in a namespace name
      const size_t function_parens( standardized_name.find_first_of( "(<", first_scope));
      if( function_parens == std::string::npos)
      {
        return standardized_name; // no function start, the string was likely already a namespace
      }

      const size_t last_pos( standardized_name.rfind( "::", function_parens));
      if( function_parens == std::string::npos)
      {
        return std::string(); // no scope, just a function name, so no namespace
      }

      // find the first non-alphanumeric character other than ':' and '_', starting with the character before (
      size_t first_pos( function_parens - 1);
      while
      (
        first_pos > 0
        &&
        (
          isalnum( standardized_name[ first_pos - 1])
          || standardized_name[ first_pos - 1] == '_'
          || standardized_name[ first_pos - 1] == ':'
        )
      )
      {
        --first_pos;
      }

      //return the substring between begin and end of class name
      return standardized_name.substr( first_pos, last_pos - first_pos);
    }

  } // namespace util
} // namespace bcl

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
#include "util/bcl_util_class_name_standardizer.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_string_functions.h"
#include "util/bcl_util_string_replacement.h"

// external includes - sorted alphabetically
#include <climits>
#include <map>

namespace bcl
{
  namespace util
  {

    //! @brief Standardize a class name that came from a function template (e.g. GetStaticClassName<>)
    //! @param NAME the name from the GetStaticClassName<>
    //! @return the standardized class name
    std::string ClassNameStandardizer::StandardizeTemplateFunctionParameter( const std::string &NAME)
    {
      return Standardize( ExtractTemplateParameter( NAME));
    }

    //! @brief extracts the template list
    //! The '__PRETTY_FUNCTION__' macro provided by the compiler should be used.
    //! If the Line looks like "virtual const std::string& bcl::Test< T>::GetName() const [with T = double]" the
    //! string "double" will be extracted
    //! @param PRETTY_FUNCTION_NAME should be passed from the  '__PRETTY_FUNCTION__' macro
    //! @return the template parameter in PRETTY_FUNCTION_NAME
    std::string ClassNameStandardizer::ExtractTemplateParameter( const std::string &PRETTY_FUNCTION_NAME)
    {
      //find the last "[with" or "[ with"
      size_t first_pos( PRETTY_FUNCTION_NAME.rfind( "[with"));
      size_t last_pos = std::string::npos;
      if( first_pos == std::string::npos) // also allow a space between [ and with
      {
        first_pos = PRETTY_FUNCTION_NAME.rfind( "[ with");
        if( first_pos == std::string::npos)
        {
          first_pos = PRETTY_FUNCTION_NAME.rfind( "[t_Class");
        }
      }

      if( first_pos != std::string::npos) // if either one exists, find the matching bracket
      {
        last_pos = PRETTY_FUNCTION_NAME.find( ']', first_pos);
        if( last_pos != std::string::npos)
        {
          first_pos = PRETTY_FUNCTION_NAME.find( "= ", first_pos); // go past the equals to get to the first parameter
          if( first_pos != std::string::npos)
          {
            first_pos += 2;
          }
        }
      }
      else
      {
        last_pos = PRETTY_FUNCTION_NAME.rfind( '>');
        if( last_pos != std::string::npos)
        {
          first_pos = RFindMatchingBracket( PRETTY_FUNCTION_NAME, last_pos, '<', '>');
          if( first_pos != std::string::npos)
          {
            first_pos++;
          }
        }
      }

      // if neither exists, assume the whole string is the name
      if( first_pos == std::string::npos || last_pos == std::string::npos)
      {
        first_pos = 0;
        last_pos  = PRETTY_FUNCTION_NAME.size();
      }

      const size_t semicolon_pos( PRETTY_FUNCTION_NAME.rfind( ';'));
      if( semicolon_pos != std::string::npos && semicolon_pos > first_pos)
      {
        last_pos = semicolon_pos;
      }

      // extract the relevant part of the string
      std::string parameter
      (
        TrimString
        (
          PRETTY_FUNCTION_NAME.size() && first_pos != last_pos
          ? PRETTY_FUNCTION_NAME.substr( first_pos, last_pos - first_pos)
            : std::string()
        )
      );

      return parameter;
    }

    //! @brief a helper function for standardizing a class name
    //! @param NAME should be the name passed only from any of the standardize class name functions
    void ClassNameStandardizer::AbbreviateEnums( std::string &NAME)
    {
      static StringReplacement finder( StringReplacement::e_Prefix, "Enum<");

      // For each position where Enum< is found
      for
      (
        // start looking at position 0
        size_t position( finder.FindNextMatch( NAME, 0));
        position != std::string::npos;
        position = finder.FindNextMatch( NAME, position + 1)
      )
      {
        size_t namespace_length( 0);
        bool namespace_okay( true);
        if( position >= 6 && NAME[ position - 1] == ':') // do we have "util::Enum<"?
        {
          // make sure the name really is prefixed by util:: and not some other namespace
          if( NAME.substr( position - 6, 6) == "util::" && ( position == 6 || !isalnum( NAME[ position - 7])))
          {
            namespace_length = 6;
            if( position > 6 && NAME[ position - 7] == ':')
            {
              // since there appears to be a namespace before the "util::", better make sure it is bcl, with no
              // other namespaces before it
              if
              (
                position >= 11 && NAME.substr( position - 11, 5) == "bcl::"
                &&
                ( position == 11 || !isalnum( NAME[ position - 12]))
              )
              {
                namespace_length = 11;
              }
              else
              {
                namespace_okay = false;
              }
            }
          }
          else
          {
            namespace_okay = false;
          }
        }

        if( !namespace_okay)
        {
          position++;
          continue;
        }

        if( namespace_length)
        {
          position -= namespace_length;
          // erase the "bcl::util::" scope of the Enum to make the Enum have the same scope as the derived object
          NAME.erase( position, namespace_length);
        }

        // now that the namespace component has been erased, NAME[ position] is the 'E' in Enum<; position + 4 is '<'
        size_t open_angle_pos( position + 4);

        size_t pos( open_angle_pos + 1);             // index currently under examination
        size_t param_two_start_pos( open_angle_pos); // the position of the second template parameter
        size_t angle_depth( 1);                      // number of unclosed template < (angles) at index pos

        // find the comma and the end angle
        for( ; angle_depth > 0; pos++)
        {
          if( NAME[ pos] == ',')
          {
            if( angle_depth == 1)
            {
              param_two_start_pos = pos + 1;
            }
          }
          else if( NAME[ pos] == '<')
          {
            ++angle_depth;
          }
          else if( NAME[ pos] == '>')
          {
            --angle_depth;
          }
        }
        size_t close_angle_pos( pos - 1);

        std::string derived_class_name
        (
          StringReplacement::SafeSubstr( NAME, param_two_start_pos, close_angle_pos - param_two_start_pos)
        );

        // remove everything in between the templates
        NAME.erase( open_angle_pos, close_angle_pos - open_angle_pos + 1);

        // insert the derived classes name
        NAME.insert( position, derived_class_name + "::");
      }
    }

    //! @brief a helper function for standardizing a class name
    //! @param NAME should be the name passed only from Standardize class name
    //! @param WORD should be the word to migrate forward in the type declaration to the last sequence character (*>,)
    void ClassNameStandardizer::MigrateWordForward( std::string &NAME, const std::string &WORD)
    {
      StringReplacement finder( StringReplacement::e_Word, WORD);

      // For each position of WORD
      for
      (
        // start looking at position 1 ( a const at the start cannot move forward anymore)
        size_t position = finder.FindNextMatch( NAME, 1);
        position != std::string::npos;
        position = finder.FindNextMatch( NAME, position + 1)
      )
      {
        size_t new_position = position;

        // Put WORD immediately after the nearest of the A. start of the string and B. any of "*,<"
        // If any of ">)]" are found before one of the above characters, then find the end of the brace context first
        {
          size_t first_sequencing_char_pos( NAME.find_last_of( "*,<", position));
          size_t first_rfence_char_pos( NAME.find_last_of( ">)]", position));
          if( first_sequencing_char_pos == std::string::npos)
          {
            // Then the const must be moved to the front of the string
            new_position = 0;
          }
          else if( first_rfence_char_pos == std::string::npos || first_sequencing_char_pos > first_rfence_char_pos)
          {
            // Then the const must be moved right after the sequencing character
            new_position = first_sequencing_char_pos + 1;
          }
          else
          {
            static std::map< char, char> s_matching_fence; // matching fence, e.g. ( <-> )
            if( !s_matching_fence.size()) // need to initialize the map
            {
              s_matching_fence[ '<'] = '>';
              s_matching_fence[ '>'] = '<';
              s_matching_fence[ '('] = ')';
              s_matching_fence[ ')'] = '(';
              s_matching_fence[ '['] = ']';
              s_matching_fence[ ']'] = '[';
            }

            while( first_rfence_char_pos > first_sequencing_char_pos)
            {
              // have to skip over the fence
              size_t matching_lfence_char_pos
              (
                RFindMatchingBracket
                (
                  NAME,
                  first_rfence_char_pos,
                  s_matching_fence[ NAME[ first_rfence_char_pos]],
                  NAME[ first_rfence_char_pos]
                )
              );

              first_sequencing_char_pos = NAME.find_last_of( "*,<", matching_lfence_char_pos - 1);
              first_rfence_char_pos     = NAME.find_last_of( ">)]", matching_lfence_char_pos - 1);
              if( first_sequencing_char_pos == std::string::npos)
              {
                // "const" must be moved to the front of the string
                new_position = 0;
                break;
              }
              else if( first_rfence_char_pos == std::string::npos || first_sequencing_char_pos > first_rfence_char_pos)
              {
                // "const" must be moved right after the sequencing character
                new_position = first_sequencing_char_pos + 1;
                break;
              }
            }
          }
        }

        if( position == new_position) // no need to change anything
        {
          continue;
        }

        // new_position < position
        // now erase WORD or -WORD.
        if( NAME[ position - 1] == '-')
        {
          NAME.erase( position - 1, WORD.size() + 1);
        }
        else if( position + 1 < NAME.size() && NAME[ position + 1] == '-')
        {
          // If WORD is not preceded by dash but a dash follows WORD, remove the dash too
          NAME.erase( position, WORD.size() + 1);
        }
        else
        {
          // just remove the WORD
          NAME.erase( position, WORD.size());
        }

        // add WORD back where it should be
        NAME.insert( new_position, WORD + "-");
      }
    }

    //! @brief Remove struct/class/union/enum artifacts put into __PRETTY_NAME__ by the MSVS compiler
    //! MSVS adds a space between class/struct/union/enum and the type-name, but not for a pointer to a type.
    //! So some class names start with "class " and others with "class" (no space).
    //! so for example, with class class_type {}, GetStaticClassName< class_type> would
    //! return class-class_type, but GetStaticClassName< class_type*> returns classclass_type
    //! in MSVS.
    //! @param MSVS_NAME the output of the __PRETTY_NAME__
    //! @return MSVS_NAME without the MSVS artifacts.  Should not be used if the string came from GCC, as it could
    //!         potentially alter the class name
    std::string ClassNameStandardizer::HandleMSVSDataTypeString( const std::string MSVS_NAME)
    {
      std::string parameter( MSVS_NAME);

      //! To standardize the behavior of this function across GCC and MSVS then, any word that is prepended by
      //! Here we want to remove only the first class, enum, struct, or union prefixes, since the second may be part of
      //! the type name, e.g.
      //! class class_type {}; should come back as class_type, not _type.
      std::vector< std::string> type_names( 4);
      type_names[ 0] = "class";
      type_names[ 1] = "struct";
      type_names[ 2] = "union";
      type_names[ 3] = "enum";

      // start by initializing a set with the start and end positions of each type_name.  The end position is the first
      // letter of the next word following any type_name
      std::set< size_t> type_start_positions;
      std::set< size_t> type_end_positions;
      for( size_t type_num = 0; type_num < type_names.size(); type_num++)
      {
        std::string &sought = type_names[ type_num];
        for
        (
          size_t index = parameter.find( sought);
          index != std::string::npos;
          index = parameter.find( sought, index + 1)
        )
        {
          if
          (
            (
              !index || ( ispunct( parameter[ index - 1]) && parameter[ index - 1] != '_')
            )
            &&
            index + sought.size() < parameter.size()
          )
          {
            if( parameter[ index + sought.size()] == '-')
            {
              type_start_positions.insert( index);
              type_end_positions.insert( index + sought.size() + 1);
            }
          }
        }
      }

      for( size_t type_num = 0; type_num < type_names.size(); type_num++)
      {
        std::string &sought = type_names[ type_num];
        for
        (
          size_t index = parameter.find( sought);
          index != std::string::npos;
          index = parameter.find( sought, index + 1)
        )
        {
          if
          (
            (
              !index || ( ispunct( parameter[ index - 1]) && parameter[ index - 1] != '_')
            )
            &&
            index + sought.size() < parameter.size()
          )
          {
            if( parameter[ index + sought.size()] != '-' && type_end_positions.count( index) == 0)
            {
              type_start_positions.insert( index);
              type_end_positions.insert( index + sought.size());
            }
          }
        }
      }

      // remove all the type names that were found
      for
      (
        std::set< size_t>::reverse_iterator
          end_of_type( type_end_positions.rbegin()),
          start_of_type( type_start_positions.rbegin()),
          end_of_set( type_start_positions.rend());
        start_of_type != end_of_set;
        end_of_type++, start_of_type++
      )
      {
        parameter = StringReplacement::SafeReplaceAt( parameter, *start_of_type, *end_of_type - *start_of_type, "");
      }

      return parameter;
    }

    //! @brief remove any ul ui, l, i, or other suffix from numerical template parameters
    void ClassNameStandardizer::StandardizeNumericTemplateParameters( std::string &NAME)
    {
      size_t next_number_pos( NAME.find_first_of( "0123456789", 0));
      while( next_number_pos != std::string::npos)
      {
        if( StringReplacement::IsNonVariableCharacter( NAME[ next_number_pos - 1]))
        {
          size_t next_non_number_pos( NAME.find_first_not_of( "0123456789", next_number_pos + 1));
          size_t next_sequencing_pos( NAME.find_first_of( ",>", next_number_pos + 1));
          if( next_non_number_pos != next_sequencing_pos)
          {
            NAME.erase( next_non_number_pos, next_sequencing_pos - next_non_number_pos);
          }
          next_number_pos = next_non_number_pos - 1;
        }

        next_number_pos = NAME.find_first_of( "0123456789", next_number_pos + 1);
      }
    }

    //////////////////////////////////////////////////////////////////////////////////////
    //! @class TypesAreEqual
    //! @brief allows determination, at compile time, whether types A and B are the same
    //! @tparam A, B the two types involved
    //////////////////////////////////////////////////////////////////////////////////////
    template< typename A, typename B>
    struct TypesAreEqual
    {
      // base case, types are different
      enum { equal = 0, not_equal = 1};
    };

    template< typename A>
    struct TypesAreEqual< A, A>
    {
      // case when A and B are identical
      enum { equal = 1, not_equal = 0};
    };

    //! @brief Standardize a class name
    //! @param NAME should be passed from ClassNameViaTemplatedFunction or from a file
    //! @return the standardized class, struct, union, or enum name as std::string
    std::string ClassNameStandardizer::Standardize( const std::string &NAME)
    {
      std::string parameter( TrimString( NAME));

      // removes spaces after punctuation (other than '_', which may be used in class names) or other spaces
      for
      (
        size_t index = parameter.find( ' ');
        index != std::string::npos;
        index = parameter.find( ' ', index + 1)
      )
      {
        // no need to check whether index == 0 or index == parameter.size() - 1 because string was just trimmed
        if
        (
          ( isalnum( parameter[ index - 1]) || parameter[ index - 1] == '_') // was the last letter part of a type name?
          &&
          ( isalnum( parameter[ index + 1]) || parameter[ index + 1] == '_') // is the next letter part of a type name?
        )
        {
          // Then change spaces between words/class names into underscores
          parameter[ index] = '-';
        }
      }

      // remove all remaining spaces
      const size_t number_of_spaces( std::count( parameter.begin(), parameter.end(), ' '));
      if( number_of_spaces > 0)
      {
        std::remove( parameter.begin(), parameter.end(), ' ');
        parameter.resize( parameter.size() - number_of_spaces);
      }

#if defined (_MSC_VER)
      parameter = HandleMSVSDataTypeString( parameter);
#endif

      // move unsigned as far ahead in the type name as possible
      MigrateWordForward( parameter, "unsigned");

      // move volatile as far ahead in the type name as possible
      MigrateWordForward( parameter, "volatile");

      // move const as far ahead in the type name as possible
      MigrateWordForward( parameter, "const");

      // abbreviate Enum<ARG1,ARG2> to ARG2::Enum
      AbbreviateEnums( parameter);

      // change VectorND<2ul,double> to VectorND<2,double>
      StandardizeNumericTemplateParameters( parameter);

      // replace all long-ints with long
      StringReplacement( StringReplacement::e_Word, "long-int", "long").ReplaceEachIn( parameter);

      // replace all short-ints with short
      StringReplacement( StringReplacement::e_Word, "short-int", "short").ReplaceEachIn( parameter);

      // typedef std::string
      StringReplacement
      (
        StringReplacement::e_Prefix,
        "std::__cxx11::",
        "std::"
      ).ReplaceEachIn( parameter);

      // typedef std::string
      StringReplacement
      (
        StringReplacement::e_Prefix,
        "std::__1::",
        "std::"
      ).ReplaceEachIn( parameter);

      // typedef std::string
      StringReplacement
      (
        StringReplacement::e_Word,
        "std::basic_string<char,std::char_traits<char>,std::allocator<char>>",
        "std::string"
      ).ReplaceEachIn( parameter);

      // mingw (and perhaps other compilers) typedef std::basic_string<char> as std::string
      StringReplacement
      (
        StringReplacement::e_Word,
        "std::basic_string<char>",
        "std::string"
      ).ReplaceEachIn( parameter);

      // remove the typedef for string from the template parameter if it was included
      StringReplacement
      (
        StringReplacement::e_Any,
        "std::string=std::string,",
        ""
      ).ReplaceEachIn( parameter);
      StringReplacement
      (
        StringReplacement::e_Any,
        ",std::string=std::string",
        ""
      ).ReplaceEachIn( parameter);

      // MSVS 2008 converts long-longs and unsigned long longs to __intXXX where x is the number of bits
      static std::string long_long_str( Format()( sizeof( long long) * CHAR_BIT));
      StringReplacement
      (
        StringReplacement::e_Word,
        "__int" + long_long_str,
        "long-long"
      ).ReplaceEachIn( parameter);

      StringReplacement
      (
        StringReplacement::e_Word,
        "__uint" + long_long_str,
        "unsigned-long-long"
      ).ReplaceEachIn( parameter);

      // use the common names for size_t
      if( TypesAreEqual< unsigned long long, size_t>::equal == 1)
      {
        StringReplacement
        (
          StringReplacement::e_Word,
          "unsigned-long-long",
          "size_t"
        ).ReplaceEachIn( parameter);
      }
      if( TypesAreEqual< unsigned long, size_t>::equal == 1)
      {
        // make a list of 1 element, which is a string replacement that matches "unsigned-long-long"
        // which we want to leave unchanged
        static const storage::List< StringReplacement>
              unsigned_long_long
              (
                1,
                StringReplacement( StringReplacement::e_Word, "unsigned-long-long")
              );

        // replace unsigned-long with size_t
        StringReplacement
        (
          StringReplacement::e_Word,
          "unsigned-long",
          "size_t"
        ).ReplaceEachWithExclusions( parameter, unsigned_long_long);
      }
      if( TypesAreEqual< unsigned int, size_t>::equal == 1)
      {
        StringReplacement
        (
          StringReplacement::e_Word,
          "unsigned-int",
          "size_t"
        ).ReplaceEachIn( parameter);
      }

      // remove any '-' next to punctuation
      StringReplacement( StringReplacement::e_Prefix, "-").ReplaceAllIn( parameter);
      StringReplacement( StringReplacement::e_Suffix, "-").ReplaceAllIn( parameter);

      // and make the other replacements necessary to standardize the pretty-function macro
      return parameter;
    } // Standardize

  } // namespace util
} // namespace bcl
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
#include "util/bcl_util_cleanable_interface.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_call_stack.h"
#include "util/bcl_util_format.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically
#include <map>
#include <signal.h>
#include <stdlib.h>

#ifndef _SIGHANDLER_T_DEFINED
#define _SIGHANDLER_T_DEFINED
//! @brief signal_handler_t
typedef void( *sighandler_t)( int);
#endif

namespace bcl
{
  namespace util
  {

    //! @brief access to singelton map that stores previous signal handlers
    //! @return map that has the signal handler function pointer for the bcl overwritten signal handler
    static std::map< int, sighandler_t> &GetSignalHandlerMap()
    {
      static std::map< int, sighandler_t> *s_signal_handler_map( new std::map< int, sighandler_t>());

      return *s_signal_handler_map;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief private constructor
    //! registers the overall cleanup function
    CleanableInterface::Cleanables::Cleanables()
    {
      atexit( &AtExitClean);

      #ifndef BCL_NO_OS_SIGNAL_HANDLING
      GetSignalHandlerMap()[ SIGILL ] = ::signal( SIGILL , &SignalHandler);
      GetSignalHandlerMap()[ SIGINT ] = ::signal( SIGINT , &SignalHandler);
      GetSignalHandlerMap()[ SIGFPE ] = ::signal( SIGFPE , &SignalHandler);
      GetSignalHandlerMap()[ SIGABRT] = ::signal( SIGABRT, &SignalHandler);
      GetSignalHandlerMap()[ SIGSEGV] = ::signal( SIGSEGV, &SignalHandler);
      GetSignalHandlerMap()[ SIGTERM] = ::signal( SIGTERM, &SignalHandler);
      #endif
    }

  //////////
  // data //
  //////////

    //! @brief access to the only instance of Cleanables
    CleanableInterface::Cleanables &CleanableInterface::Cleanables::GetCleanables()
    {
      // only instance of cleanables - needs to exist to he absolute end, so it is a new pointer, that will never be
      // deleted and only cleaned up by the operating system
      static Cleanables *s_cleanables( new Cleanables());

      // return
      return *s_cleanables;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief register the CLEANABLE
    void CleanableInterface::Cleanables::Register( CleanableInterface *CLEANABLE)
    {
      m_Mutex.Lock();
      m_Cleanables.insert( CLEANABLE);
      m_Mutex.Unlock();
    }

    //! @brief unregister the CLEANABLE
    void CleanableInterface::Cleanables::UnRegister( CleanableInterface *CLEANABLE)
    {
      m_Mutex.Lock();
      m_Cleanables.erase( CLEANABLE);
      m_Mutex.Unlock();
    }

    //! @brief clean all cleanables that are registered
    void CleanableInterface::Cleanables::CleanAllCleanables()
    {
      // iterate over all cleanables and call clean
      while( !m_Cleanables.empty())
      {
        CleanableInterface *current( *m_Cleanables.begin());
        m_Cleanables.erase( m_Cleanables.begin());
        current->CleanUp();
      }
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Cleanup all cleanables
    void CleanableInterface::Cleanables::AtExitClean()
    {
      GetCleanables().CleanAllCleanables();
    }

    //! @brief signal handler, for non exit calls
    //! @param SIGNAL
    void CleanableInterface::Cleanables::SignalHandler( int SIGNAL)
    {
      static volatile sig_atomic_t flag( 0);
      if( flag != 0)
      {
        return;
      }

      // set flag, so that signal is not processed twice
      flag = 1;
      BCL_MessageTop
      (
        "caught signal: " + util::Format()( SIGNAL) + " cleaning! " + util::CallStack().String()
      );

      // clean all the things that needed to be cleaned
      GetCleanables().CleanAllCleanables();
      std::map< int, sighandler_t>::iterator itr( GetSignalHandlerMap().find( SIGNAL));

      // reset the signal handler
      ::signal( SIGNAL, itr->second);

      // remove the itr
      GetSignalHandlerMap().erase( itr);

      // raise the signal again
      raise( SIGNAL);
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor
    //! registers this with the cleanables
    CleanableInterface::CleanableInterface()
    {
      Cleanables::GetCleanables().Register( this);
    }

    //! @brief virtual destructor
    CleanableInterface::~CleanableInterface()
    {
      Cleanables::GetCleanables().UnRegister( this);
    }

  } // namespace util
} // namespace bcl
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
#include "util/bcl_util_color_gradient.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "util/bcl_util_colors.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const SiPtr< const ObjectInterface> ColorGradient::s_Instance
    (
      Enumerated< FunctionInterfaceSerializable< double, linal::Vector3D> >::AddInstance( new ColorGradient())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ColorGradient::ColorGradient() :
      m_Range( 0.0, 1.0),
      m_GradientColors( GetColors().GetRainbow())
    {
    }

    //! @brief constructor from a range and vector of gradient colors
    //! @param RANGE range of values
    //! @param GRADIENT_POINTS all of the colors which will be used in the color gradient
    ColorGradient::ColorGradient
    (
      const math::Range< double> &RANGE,
      const storage::Vector< Color> &GRADIENT_POINTS
    ) :
      m_Range( RANGE),
      m_GradientColors( GRADIENT_POINTS)
    {
      BCL_Assert( GRADIENT_POINTS.GetSize() > 1, "The gradient needs to have more than one color");
    }

    //! @brief virtual copy constructor
    //! @return new copy of this class
    ColorGradient *ColorGradient::Clone() const
    {
      return new ColorGradient( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ColorGradient::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &ColorGradient::GetAlias() const
    {
      static const std::string s_Name( "ColorGradient");
      return s_Name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking a value and returning the corresponding color point
    //! @param VALUE the value which will be converted into a color
    //! @return returns the corresponding color
    linal::Vector3D ColorGradient::operator()( const double &VALUE) const
    {
      // if value larger than the upper bound
      if( VALUE >= m_Range.GetMax())
      {
        BCL_MessageCrt
        (
          Format()( VALUE) + " is not within range " + m_Range.GetString() + " returning max color"
        );

        // return the color for the largest expected value
        return m_GradientColors.LastElement();
      }
      // if value is smaller than the lower bound
      else if( VALUE < m_Range.GetMin())
      {
        BCL_MessageCrt
        (
          Format()( VALUE) + " is not within range " + m_Range.GetString() + " returning min color"
        );

        // return the color for the smallest expected value
        return m_GradientColors.FirstElement();
      }

      // rescale the value according to a new range
      const double position( ( VALUE - m_Range.GetMin()) / m_Range.GetWidth() * ( m_GradientColors.GetSize() - 1));

      // calculate the index of the color in the gradient this VALUE corresponds to
      const size_t interval_index( ( size_t)( position));

      // calculate the difference between start and end points of this gradient
      const linal::Vector3D color_difference( *m_GradientColors( interval_index + 1) - *m_GradientColors( interval_index));

      // calculate the new color and return it
      return *m_GradientColors( interval_index) + color_difference * ( position - double( interval_index));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ColorGradient::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Range, ISTREAM);
      io::Serialize::Read( m_GradientColors, ISTREAM);
      BCL_Assert( m_GradientColors.GetSize() > 1, "The gradient needs at least 2 colors!");

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &ColorGradient::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Range, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_GradientColors, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ColorGradient::GetSerializer() const
    {
      io::Serializer member_data;
      member_data.SetClassDescription
      (
        "For creating a taking a value and providing a corresponding color."
        "It takes a double value (e.g. score) and returns the color point at which the value falls within the defined"
        "gradient based on min and max values. You provide an arbitrary number of color points which you want to"
        "transition between."
        "The gradient between colors is linear. The range between the min and max values is split up evenly so that each"
        "color transition takes place over 1/(number color transitions) of the min to max value range."
        "The color points are typically thought of as Red, Green, Blue format, although there"
        "is nothing in this class hard coded to force this. The color transitions take place in the order in which"
        "they are provided. The double is the value which will be converted into color points."
        "The color is an instance Colors enumerator."
      );

      member_data.AddInitializer
      (
        "range",
        "the range of values over which a color code will be applied (must be quoted if range string contains () or ,)",
        io::Serialization::GetAgent( &m_Range),
        "\"(0.0,0.0)\""
      );

      member_data.AddInitializer
      (
        "color_list",
        "The list of colors that the color gradient will use. First color is for range minima, last is for range max.",
        io::Serialization::GetAgentWithSizeLimits( &m_GradientColors, 2),
        ObjectDataLabel
        (
          storage::Vector< ObjectDataLabel>::Create
          (
            ObjectDataLabel( GetColors().e_Black->GetName()),
            ObjectDataLabel( GetColors().e_White->GetName())
          )
        ).ToString()
      );

      return member_data;
    }

  } // namespace util
} // namespace bcl

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
#include "util/bcl_util_colors.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_range.h"
#include "math/bcl_math_running_min_max.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct all AAClasses
    Colors::Colors() :
      e_White(   AddEnum( "White",   linal::Vector3D( 1.00, 1.00, 1.00))),
      e_Black(   AddEnum( "Black",   linal::Vector3D( 0.00, 0.00, 0.00))),
      e_Red(     AddEnum( "Red",     linal::Vector3D( 1.00, 0.00, 0.00))),
      e_Green(   AddEnum( "Green",   linal::Vector3D( 0.00, 1.00, 0.00))),
      e_Blue(    AddEnum( "Blue",    linal::Vector3D( 0.00, 0.00, 1.00))),
      e_Yellow(  AddEnum( "Yellow",  linal::Vector3D( 1.00, 1.00, 0.00))),
      e_Cyan(    AddEnum( "Cyan",    linal::Vector3D( 0.00, 1.00, 1.00))),
      e_Magenta( AddEnum( "Magenta", linal::Vector3D( 1.00, 0.00, 1.00)))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Colors::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief return vector of rainbow colors
    //! @return vector of rainbow
    const storage::Vector< Color> &Colors::GetRainbow() const
    {
      // initialize static vector
      static const storage::Vector< Color> s_rainbow
      (
        storage::Vector< Color>::Create( e_Blue, e_Green, e_Yellow, e_Red)
      );

      // end
      return s_rainbow;
    }

    //! @brief construct on access function for all Colors
    //! @return reference to only instances of Colors
    const Colors &GetColors()
    {
      return Colors::GetEnums();
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief convert an rgd color to a hsv color
    //! @param RGB_COLOR the color as rgb values 0-1
    //! @return hsv color
    linal::Vector3D Colors::ConvertRGBToHSV( const linal::Vector3D &RGB_COLOR)
    {
      math::RunningMinMax< double> min_max;
      min_max += RGB_COLOR.X();
      min_max += RGB_COLOR.Y();
      min_max += RGB_COLOR.Z();

      const math::Range< double> range( min_max.GetMin(), min_max.GetMax());

      double s(  0);
      double h( -1);

      // s
      if( range.GetMax() != 0)
      {
        s = range.GetWidth() / range.GetMax();
      }

      if( RGB_COLOR.X() == range.GetMax())
      {
        h = ( RGB_COLOR.Y() - RGB_COLOR.Z() ) / range.GetWidth();   // between yellow & magenta
      }
      else if( RGB_COLOR.Y() == range.GetMax())
      {
        h = 2 + ( RGB_COLOR.Z() - RGB_COLOR.X()) / range.GetWidth(); // between cyan & yellow
      }
      else
      {
        h = 4 + ( RGB_COLOR.X() - RGB_COLOR.Y()) / range.GetWidth(); // between magenta & cyan
      }

      h /= 6.0; // normalize to [0,1]
      if( h < 0)
      {
        h += 1.0;
      }

      // end
      return linal::Vector3D( h, s, range.GetMax());
    }

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< linal::Vector3D, Colors>;

  } // namespace util
} // namespace bcl
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
#include "util/bcl_util.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

    //! @brief identifier for the name space
    //! @return the name of the namespace
    const std::string &GetNamespaceIdentifier()
    {
      static const std::string *s_namespace_name( new std::string( ExtractNamespaceIdentifier( __PRETTY_FUNCTION__)));
      return *s_namespace_name;
    }

  } // namespace util
} // namespace bcl
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
#include "util/bcl_util_cpp_data_types.h"

// includes from bcl - sorted alphabetically
#include "type/bcl_type_compare.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    //! @brief conversion to a string from a Types
    //! @param DATA_TYPE the data type to get a string for
    //! @return a string representing that data type
    const std::string &CPPDataTypes::GetCPPDatatypeName( const CPPDataTypes::Types &DATA_TYPE)
    {
      static const std::string s_descriptors[] =
      {
        "DT_UNKNOWN",
        "DT_CHAR",
        "DT_DOUBLE",
        "DT_FLOAT",
        "DT_INT",
        "DT_UNSIGNED_INT",
        "DT_SIZE_T",
        "DT_STRING",
        "DT_BOOL",
        GetStaticClassName< CPPDataTypes::Types>()
      };
      return s_descriptors[ size_t( DATA_TYPE)];
    }

    //! @brief c++ string for datatype
    //! @param DATA_TYPE the data type to get a string for
    //! @return the string representing that data type as it would be written in the source code
    const std::string &CPPDataTypes::GetCPPString( const CPPDataTypes::Types &DATA_TYPE)
    {
      static const std::string s_strings[] =
      {
        "",
        "char",
        "double",
        "float",
        "int",
        "unsigned int",
        "size_t",
        "std::string",
        "bool",
        ""
      };
      return s_strings[ size_t( DATA_TYPE)];
    }

    //! @brief Types from template parameter for char
    //! @return Types for given t_DataType
    template<>
    BCL_API CPPDataTypes::TypeEnum CPPDataTypes::DataTypeFromTemplate< char>()
    {
      return e_Char;
    }

    //! @brief Types from template parameter for double
    //! @return Types for given t_DataType
    template<>
    BCL_API CPPDataTypes::TypeEnum CPPDataTypes::DataTypeFromTemplate< double>()
    {
      return e_Double;
    }

    //! @brief Types from template parameter for float
    //! @return Types for given t_DataType
    template<>
    BCL_API CPPDataTypes::TypeEnum CPPDataTypes::DataTypeFromTemplate< float>()
    {
      return e_Float;
    }

    //! @brief Types from template parameter for int
    //! @return Types for given t_DataType
    template<>
    BCL_API CPPDataTypes::TypeEnum CPPDataTypes::DataTypeFromTemplate< int>()
    {
      return e_Int;
    }

    //! @brief Types from template parameter for unsigned ints
    //! @return Types for given t_DataType
    template<>
    BCL_API CPPDataTypes::TypeEnum CPPDataTypes::DataTypeFromTemplate< unsigned int>()
    {
      // Prefer the real type (unsigned int) to the typedef-ed size_t,
      //   unsigned int is often the same as size_t on 32-bit architectures
      return e_UnsignedInt;
    }

    //! @brief Types from template parameter for unsigned longs
    //! @return Types for given t_DataType
    template<>
    BCL_API CPPDataTypes::TypeEnum CPPDataTypes::DataTypeFromTemplate< unsigned long>()
    {
      return type::Compare< unsigned long, size_t>::e_Same ? e_SizeT : e_UnsignedInt;
    }

    //! @brief Types from template parameter for string
    //! @return Types for given t_DataType
    template<>
    BCL_API CPPDataTypes::TypeEnum CPPDataTypes::DataTypeFromTemplate< std::string>()
    {
      return e_String;
    }

    //! @brief Types from template parameter for bool
    //! @return Types for given t_DataType
    template<>
    BCL_API CPPDataTypes::TypeEnum CPPDataTypes::DataTypeFromTemplate< bool>()
    {
      return e_Bool;
    }

  } // namespace util

} // namespace bcl
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
#include "util/bcl_util_data_type.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"

// external includes - sorted alphabetically
#include <limits>

namespace bcl
{
  namespace util
  {

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief get the type name
    //! @param TYPE the type name
    //! @return string describing the type
    const std::string &DataType::GetTypeName( const DataType::Type &TYPE)
    {
      static const std::string s_names[ s_NumberTypes + 1] =
      {
        "Bool",
        "Char",
        "Int",
        "UnsignedInt",
        "Float",
        "String",
        "Scalar",
        "Sequence",
        "Map",
        "Set",
        "StaticObject",
        "DynamicObject",
        "Pointer",
        GetStaticClassName< DataType::Type>()
      };
      return s_names[ TYPE];
    }

    //! @brief detect whether a given type must be scalar, e.g. not have any sub-arguments
    //! @param TYPE the type
    //! @return bool - true if the type must be scalar
    bool DataType::TestMustBeScalar( const DataType::Type &TYPE)
    {
      return TYPE <= s_NumberScalarTypes;
    }

    //! @brief test whether the underlying cpp type for a given object is known
    //! @param TYPE the type
    //! @return bool - true if the type is known
    bool DataType::IsUnderlyingTypeKnown( const Type &TYPE)
    {
      return TYPE < e_DynamicObject;
    }

    //! @brief convenience function for containers to write out their size requirements
    //! @param OSTREAM stream to write out size requirements to
    //! @param MIN_SIZE the minimum size of the container
    //! @param MAX_SIZE the maximum size of the container
    void DataType::WriteSizeRequirements
    (
      std::ostream &OSTREAM,
      const size_t &MIN_SIZE,
      const size_t &MAX_SIZE
    )
    {
      if( MIN_SIZE > size_t( 0))
      {
        OSTREAM << " with ";
        if( MAX_SIZE == std::numeric_limits< size_t>::max())
        {
          OSTREAM << "at least " << MIN_SIZE;
        }
        else if( MIN_SIZE == MAX_SIZE)
        {
          OSTREAM << MIN_SIZE;
        }
        else
        {
          OSTREAM << "between " << MIN_SIZE << " and " << MAX_SIZE;
        }
      }
      else if( MAX_SIZE < std::numeric_limits< size_t>::max())
      {
        OSTREAM << " with at most " << MAX_SIZE;
      }
      OSTREAM << ' ';
    }

  } // namespace util
} // namespace bcl
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

// include forward header of this class
#include "util/bcl_util_enumerated.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

  } // namespace util
} // namespace bcl
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
#include "util/bcl_util_enums_instances.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_parameter.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

  //////////
  // data //
  //////////

    //! @brief GetFlagEnumsFiles gives command line flag to pass custom enums definitions
    //! @return ShPtr to a FlagInterface which is used to pass filenames for custom enum definitions
    ShPtr< command::FlagInterface> &EnumsInstances::GetFlagEnumsFiles()
    {
      static ShPtr< command::FlagInterface> s_flag_enums_files
      (
        new command::FlagDynamic
        (
          "enums_files",
          "files for enum data that adds enums or overrides data of existing enum data",
          command::Parameter
          (
            "enum_file",
            "file that is similar to a written Enums derived class"
          ),
          0,
          EnumsInstances::GetEnumsInstances().GetSize(),
          &EnumsInstances::AddEnumsFromCommandline
        )
      );

      return s_flag_enums_files;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    EnumsInstances::EnumsInstances() :
      m_EnumsInstances()
    {
    }

    //! @brief return the number of instances
    size_t EnumsInstances::GetSize() const
    {
      return m_EnumsInstances.size();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief add enums from files that are given in the commandline
    void EnumsInstances::AddEnumsFromCommandline()
    {
      // for each filename given in commandline
      for
      (
        ShPtrVector< command::ParameterInterface>::const_iterator
          itr_param( GetFlagEnumsFiles()->GetParameterList().Begin()),
          itr_param_end( GetFlagEnumsFiles()->GetParameterList().End());
        itr_param != itr_param_end;
        ++itr_param
      )
      {
        // create and open stream to given file
        io::IFStream read;
        io::File::MustOpenIFStream( read, ( *itr_param)->GetValue());

        const std::string enums_name( ObjectInterface::ExtractIdentifier( read));

        // find Enums derived class with given identifier
        iterator itr_enum( GetEnumsInstances().m_EnumsInstances.find( enums_name));

        if( itr_enum == GetEnumsInstances().m_EnumsInstances.end())
        {
          BCL_MessageCrt
          (
            "Enums instance with name \"" + enums_name + "\" given in file \"" + ( *itr_param)->GetValue() +
            " \" could not be found as instance. It might be that this Enums derived class is not instantiated since"
            " Singleton constructor was not called yet!"
          );
        }
        else
        {
          // read enums to that
          itr_enum->second->Read( read);
        }

        // close and clear stream
        io::File::CloseClearFStream( read);
      }
    }

    //! @brief the one and only instance of this class
    EnumsInstances &EnumsInstances::GetEnumsInstances()
    {
      static EnumsInstances s_enums_insances;

      return s_enums_insances;
    }

    //! @brief insert an instance of an Enums derived class
    void EnumsInstances::Insert( ObjectInterface &ENUMS_OBJECT, const std::string &ENUMS_NAME)
    {
      // check if Enums derived class only has one instance
      const_iterator itr( m_EnumsInstances.find( ENUMS_NAME));

      if( itr != m_EnumsInstances.end())
      {
        BCL_Exit( "more than one instance of Enums derived class: " + ENUMS_NAME + " is not supported!", -1);
      }

      // insert that instance
      m_EnumsInstances[ ENUMS_NAME] = &ENUMS_OBJECT;

      GetObjectInstances().AddInstanceWithName( &ENUMS_OBJECT, ENUMS_NAME);
    }

  } // namespace util
} // namespace bcl
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
#include "util/bcl_util_format.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialize.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Format::Format() :
      m_Precision( GetUndefined< size_t>()),
      m_Width( GetUndefined< size_t>()),
      m_SingleSpace( false),
      m_ForceWidth( false),
      m_Fill( ' '),
      m_Format(),
      m_Initialized( false)
    {
    }

    //! copy function
    Format *Format::Clone() const
    {
      return new Format( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
    const std::string &Format::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief turn SingleSpace on and off
    //! @return Format reference to this
    Format &Format::S()
    {
      m_SingleSpace = !m_SingleSpace;
      return *this;
    }

    //! @brief turn fixed width on and off
    //! @return Format reference to this
    Format &Format::ForceW()
    {
      m_ForceWidth = !m_ForceWidth;
      return *this;
    }

    //! @brief set left aligned output
    //! @return Format reference to this
    Format &Format::L()
    {
      m_Format |= std::ios::left;
      m_Initialized = true;
      return *this;
    }

    //! @brief set right aligned output
    //! @return Format reference to this
    Format &Format::R()
    {
      m_Format |= std::ios::right;
      m_Initialized = true;
      return *this;
    }

    //! @brief set width
    //! @param WIDTH width of formatted numerical value
    //! @return Format reference to this
    Format &Format::W( const size_t WIDTH)
    {
      m_Width = WIDTH;
      return *this;
    }

    //! @brief set scientific format, gets width and precision as input
    //! @param PREC the precision - number of digits after '.'
    //! @return Format reference to this
    Format &Format::SFP( const size_t PREC)
    {
      m_Precision = PREC;
      m_Format |= std::ios::scientific;
      m_Initialized = true;
      return *this;
    }

    //! @brief set floating point precision, gets precision as input
    //! @param PREC the precision - number of digits after '.'
    //! @return Format reference to this
    Format &Format::FFP( const size_t PREC)
    {
      m_Precision = PREC;
      m_Format |= std::ios::fixed;
      m_Initialized = true;
      return *this;
    }

    //! @brief set filling character - will be used to fill unused width
    //! @param FILL the character to fill with - default ' '
    //! @return Format reference to this
    Format &Format::Fill( const char FILL)
    {
      m_Fill = FILL;
      if( FILL == '0')
      {
        m_Format |= std::ios_base::internal;
      }
      m_Initialized = true;
      return *this;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief apply flags on stream
    //! @param IOSTREAM the stream to set the flags on
    //! @return std::ios reference to the IOSTREAM
    std::ios &Format::SetFlags( std::ios &IOSTREAM) const
    {
      if( IsDefined( m_Precision))
      {
        IOSTREAM.precision( m_Precision);
      }

      if( IsDefined( m_Width))
      {
        IOSTREAM.width( m_Width);
      }

      if( m_Fill != ' ')
      {
        IOSTREAM.fill( m_Fill);
      }

      if( m_Initialized)
      {
        IOSTREAM.setf( m_Format);
      }

      return IOSTREAM;
    }

    //! @brief remove flags on stream
    //! @param IOSTREAM the stream to disable the flags on
    //! @return std::ios reference to the IOSTREAM
    std::ios &Format::UnsetFlags( std::ios &IOSTREAM) const
    {
      // reset all formats
      IOSTREAM.unsetf( std::ios::floatfield);

      // end
      return IOSTREAM;
    }

    //! output operator() for float type (checks for nans and writes them accordingly)
    std::string Format::operator ()( const float &VALUE) const
    {
      std::stringstream stream; // make a stream

      SetupStream( stream); // set it up with any flags that are set

      if( !IsDefined( VALUE))
      {
        static std::string s_nan( "nan");
        // write nan, since the actual undefined value is machine-dependent
        stream << s_nan; // write out the value
      }
      else
      {
        stream << VALUE;
      }

      // post-process the string in the stream and return it
      return Postprocess( stream.str());
    }

    //! output operator() for double type (checks for nans and writes them accordingly)
    std::string Format::operator ()( const double &VALUE) const
    {
      std::stringstream stream; // make a stream

      SetupStream( stream); // set it up with any flags that are set

      if( !IsDefined( VALUE))
      {
        static std::string s_nan( "nan");
        // write nan, since the actual undefined value is machine-dependent
        stream << s_nan; // write out the value
      }
      else
      {
        stream << VALUE;
      }

      // post-process the string in the stream and return it
      return Postprocess( stream.str());
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Format::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Precision, ISTREAM);
      io::Serialize::Read( m_Width, ISTREAM);
      io::Serialize::Read( m_SingleSpace, ISTREAM);
      io::Serialize::Read( m_ForceWidth, ISTREAM);
      io::Serialize::Read( m_Fill, ISTREAM);
      io::Serialize::Read( m_Initialized, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &Format::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Precision,   OSTREAM, INDENT) << ' ';
      io::Serialize::Write( m_Width,       OSTREAM, INDENT) << ' ';
      io::Serialize::Write( m_SingleSpace, OSTREAM, INDENT) << ' ';
      io::Serialize::Write( m_ForceWidth,  OSTREAM, INDENT) << ' ';
      io::Serialize::Write( m_Fill,        OSTREAM, INDENT) << ' ';
      io::Serialize::Write( m_Initialized, OSTREAM, INDENT) << ' ';

      //end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief setup the stream with the given flags and add the single space, if necessary
    //! @param STREAM the stream to set properties on
    void Format::SetupStream( std::stringstream &STREAM) const
    {
      // output stream
      if( m_SingleSpace)
      {
        STREAM << ' ';
      }

      SetFlags( STREAM);
    }

    //! @brief truncates the string if force width is set
    //! @details will truncate the string if width was forced
    //! @param STRING the string to process
    //! @return the string in the stream after it was post-processed
    std::string Format::Postprocess( std::string STRING) const
    {
      // fix length
      if( m_ForceWidth)
      {
        // resize the stream to width
        STRING.resize( m_Width + size_t( m_SingleSpace), m_Fill);
      }

      return STRING;
    }

  } // namespace util
} // namespace bcl
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
#include "util/bcl_util_functional_type.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief get the type name
    //! @param TYPE the type name
    //! @return string describing the type
    const std::string &FunctionalType::GetTypeName( const FunctionalType::Type &TYPE)
    {
      static const std::string s_names[ s_NumberTypes + 1] =
      {
        "Basic Implementations",        // No user-definable parameters (e.g. Atom_Identity)
        "Customizable Implementations", // User-definable parameters (e.g. DihedralBins(bin size=30))
        "Operations",                   // Takes an object of the same interface type, no other parameters (e.g. Sum)
        "Customizable Operations",      // Takes an object of the same interface type and other parameters (e.g. 3DA)
        GetStaticClassName< FunctionalType::Type>()
      };
      return s_names[ TYPE];
    }

    //! @brief get the type from the number of wrappers and whether there are other parameters
    //! @param N_WRAPPERS the number of wrappers
    //! @param OTHER_PARAMS whether there were other parameters
    //! @return string describing the type
    FunctionalType::Type FunctionalType::GetType( const size_t &N_WRAPPERS, const bool &OTHER_PARAMS)
    {
      return Type( ( N_WRAPPERS ? 2 : 0) + ( OTHER_PARAMS ? 1 : 0));
    }

    //! @brief detect whether a given type has non-decorator parameters
    //! @param TYPE the type
    //! @return bool - true if the type is parameterized
    bool FunctionalType::TestIsParameterized( const Type &TYPE)
    {
      return size_t( TYPE) & 1;
    }

    //! @brief detect the number of decorator parameters (3 means any number other than 0-2), e.g. 2 for Add, 1 for Sum
    //! @param TYPE the type
    //! @return number of user-selected classes with the same interface
    size_t FunctionalType::GetWrapperSize( const Type &TYPE)
    {
      return size_t( TYPE) / 2;
    }

  } // namespace util
} // namespace bcl
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

// include forward header of this class

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

  } // namespace util
} // namespace bcl
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

// include forward header of this class

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

  } // namespace util
} // namespace bcl
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

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "util/bcl_util_cleanable_interface.h"
#include "util/bcl_util_loggers.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LoggerFile
    //! @brief LoggerInterface implementation printing to a file
    //! @details implementation of the LoggerInterface printing all messages to the commandline-set file
    //!
    //! @author heinzes1
    //! @date 01/12/2009
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LoggerFile :
      public LoggerInterface,
      public CleanableInterface
    {

    private:

    //////////
    // data //
    //////////

      //! file stream
      io::OFStream m_LoggerStream;

      //! message identifier
      std::string m_Filename;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LoggerFile() :
        LoggerInterface( std::cout.rdbuf()),
        m_Filename()
      {
      }

      //! @brief overwritten copy constructor
      //! @param LOGGER_FILE reference to a LoggerFile object
      LoggerFile( const LoggerFile &LOGGER_FILE) :
        LoggerInterface( LOGGER_FILE.rdbuf()),
        m_LoggerStream(),
        m_Filename()
      {
        SetLogIdentifier( LOGGER_FILE.m_Filename);
      }

      //! @brief Clone function
      //! @return pointer to new LoggerDefault
      virtual LoggerFile *Clone() const
      {
        return new LoggerFile( *this);
      }

      //! @brief virtual destructor
      virtual ~LoggerFile()
      {
        CleanUp();
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      virtual const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief virtual function that will be called if cleanable is registered
      virtual void CleanUp()
      {
        if( m_LoggerStream.is_open())
        {
          m_LoggerStream.close();
          m_LoggerStream.clear();
        }
      }

      //! @brief get the ideal maximum line width for the given logger
      //! @return the ideal maximum line width for the given logger
      size_t GetMaxLineWidth() const
      {
        return LoggerInterface::GetDefaultMaxLineWidth();
      }

      //! @brief set an identifier string for the messages
      //! @param IDENTIFIER string that appears in front of each message
      //! @return true on success
      bool SetLogIdentifier( const std::string &IDENTIFIER)
      {
        // was already initialized
        if( m_Filename == IDENTIFIER)
        {
          return true;
        }

        // stream with different filename is open - needs to be closed
        if( m_LoggerStream.is_open())
        {
          m_LoggerStream.close();
          m_LoggerStream.clear();
        }

        // set new filename and open file
        m_Filename = IDENTIFIER;
        io::File::MustOpenOFStream( m_LoggerStream, m_Filename);

        // set the buffer to the files buffer
        std::ostream::rdbuf( m_LoggerStream.rdbuf());

        // end
        return true;
      }

      //! @brief get the identifier
      //! @return the identifier
      const std::string &GetLogIdentifier() const
      {
        return m_Filename;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief pass message to the log
      //! @param MESSAGE reference to message string
      void LogMessage( const std::string &MESSAGE)
      {
        *this << "=" << MESSAGE << std::endl;
      }

      //! @brief pass error message to the log
      //! @param ERR_STRING reference to error message string
      void LogStatus( const std::string &STATUS)
      {
      }

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      virtual std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! instance of the file logger added to the enum
      static Logger s_LoggerFileInstance;

    }; // class LoggerFile

    Logger LoggerFile::s_LoggerFileInstance
    (
      GetLoggers().AddEnum( "File", ShPtr< LoggerInterface>( new LoggerFile()))
    );

  } // namespace util
} // namespace bcl
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
#include "util/bcl_util_logger_interface.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_loggers.h"

// external includes - sorted alphabetically

#if defined(_WIN32)
  #include <windows.h>
#else
  // posix window size from ioctl
  #include <errno.h>
  #include <sys/ioctl.h>
#endif

namespace bcl
{
  namespace util
  {

    //! @brief get the default maximum line width to be used for a particular logger
    //! this value should be used whenever the output is not to a terminal or the terminal window size is unavailable
    size_t LoggerInterface::GetDefaultMaxLineWidth()
    {
      return 120;
    }

    //! @brief get the line width of the terminal
    //! this value should be used whenever the output is to a terminal
    size_t LoggerInterface::GetTerminalLineWidth()
    {
      if( GetIgnoreTerminalLineWidth())
      {
        return GetDefaultMaxLineWidth();
      }
      size_t actual_cols( 0);
#if defined(_WIN32)
      // Windows
      CONSOLE_SCREEN_BUFFER_INFO buffer_rows_columns;

      if( !GetConsoleScreenBufferInfo( GetStdHandle( STD_OUTPUT_HANDLE), &buffer_rows_columns))
      {
        // probably not running in a terminal
        actual_cols = LoggerInterface::GetDefaultMaxLineWidth();
      }
      else
      {
        // running in a terminal; get the actual width of the window by subtracting the column positions
        // NOTE: when calling this function on an executable run with Wine, be aware that wine does not pass linux console
        //       information to the application.  Likewise, this value will always be 80 unless run on A. a windows box or
        //       B. running it with a proper windows console using wineconsole --backend=user.  If doing B, you may want
        //       to set wineconsole up to not exit at program completion, this setting is in ~/.wine/user.reg under
        //       ExitOnDie, which should be set to 0 to keep the console open
        actual_cols = buffer_rows_columns.srWindow.Right - buffer_rows_columns.srWindow.Left + 1;
      }
#else
      // Linux / Apple
      struct winsize window_rows_columns;

      // cache the old error number and reset it to detect errors in ioctl,
      // since its return value is non-standardized (implementations may return 0 on success, others -1, etc.)
      const int old_errno( errno);
      errno = 0;
      ioctl( 0, TIOCGWINSZ, &window_rows_columns);
      if( !errno)
      {
        // get the column width, if it could be retrieved successfully
        actual_cols = window_rows_columns.ws_col;
      }
      else
      {
        // probably not running in a terminal
        actual_cols = LoggerInterface::GetDefaultMaxLineWidth();
      }
      errno = old_errno;
#endif
      static const size_t s_min_window_size( 60);
      if( actual_cols < s_min_window_size)
      {
        // if the terminal size is really small (say, smaller than 60), than most well formatted output will look bad
        // anyway, so scale up the output
        actual_cols = std::max( actual_cols, size_t( 1));

        // try to still keep actual columns as an integer multiple of its real value
        actual_cols *= ( ( s_min_window_size - 1) / actual_cols + 1);
      }

      return actual_cols;
    }

    //! @brief set a flag to ignore the terminal's line width when reporting the default max line width
    //! @param IGNORE_TERMINAL_LINE_WIDTH true to ignore the terminal's line width when reporting GetTerminalLineWidth
    void LoggerInterface::SetIgnoreTerminalLineWidth( const bool &IGNORE_TERMINAL_LINE_WIDTH)
    {
      GetIgnoreTerminalLineWidth() = IGNORE_TERMINAL_LINE_WIDTH;
    }

    //! @brief get currently used Logger
    LoggerInterface &GetLogger()
    {
      return Loggers::GetEnums().GetCurrentLogger();
    }

  } // namespace util
} // namespace bcl
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
#include "util/bcl_util_loggers.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "util/bcl_util_enumerate.hpp"
#include "util/bcl_util_logger_default.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

  //////////
  // data //
  //////////

    //! command line flag to be used to set Logger over the command line
    ShPtr< command::FlagInterface> &Loggers::GetFlagLogger()
    {
      static ShPtr< command::FlagInterface> s_logger_flag
      (
        new command::FlagStatic
        (
          "logger",
          "change the logger this executable uses",
          ShPtrVector< command::ParameterInterface>::Create
          (
            ShPtr< command::ParameterInterface>
            (
              new command::Parameter
              (
                "",
                "",
                command::ParameterCheckEnumerate< Loggers>(),
                GetLoggers().e_Default.GetName()
              )
            ),
            ShPtr< command::ParameterInterface>
            (
              new command::Parameter
              (
                "identifier",
                "define a logger identifier - for file, it is the filename to be opened",
                ""
              )
            )
          ),
          &Loggers::UpdateCurrentLoggerFromCommandLineFlag
        )
      );

      return s_logger_flag;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Default constructor
    Loggers::Loggers() :
      Enumerate< ShPtr< LoggerInterface>, Loggers>( false),
      e_Default( AddEnum( "Default", ShPtr< LoggerInterface>( new LoggerDefault()))),
      m_CurrentLogger( new LoggerDefault()),
      m_CurrentLoggerEnum( e_Default)
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Loggers::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief Initialize the logger from the command line flag
    void Loggers::UpdateCurrentLoggerFromCommandLineFlag()
    {
      const Logger new_logger( GetFlagLogger()->GetFirstParameter()->GetValue());
      const std::string &new_identifier( GetFlagLogger()->GetParameterList()( 1)->GetValue());

      Loggers &loggers( GetEnums());
      if( new_logger == loggers.m_CurrentLoggerEnum && new_identifier == loggers.m_CurrentLogger->GetLogIdentifier())
      {
        // no need to update the current logger
        return;
      }

      // switching the current logger
      BCL_MessageStd( "change Logger to: " + new_logger.GetName() + " with identifier: " + new_identifier);
      loggers.m_CurrentLogger = ShPtr< LoggerInterface>( ( *new_logger)->Clone());
      loggers.m_CurrentLoggerEnum = new_logger;

      // set the identifier of the current logger
      BCL_Assert
      (
        loggers.m_CurrentLogger->SetLogIdentifier( new_identifier),
        "unable to initialize logger: " + new_logger.GetName() + " with identifier: " + new_identifier
      );

      // switch was successful
      BCL_MessageStd( "Logger was changed to: " + new_logger.GetName() + " with identifier: " + new_identifier);
    }

    //! @brief Initialize and return the current logger
    //! @return one and only reference to one of the loggers
    LoggerInterface &Loggers::GetCurrentLogger()
    {
      // return the logger
      return *m_CurrentLogger;
    }

    //! @brief get enumerated list of Loggers
    Loggers &GetLoggers()
    {
      return Loggers::GetEnums();
    }

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< LoggerInterface>, Loggers>;

  } // namespace util
} // namespace bcl
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
#include "util/bcl_util_memory_usage.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically
#include <unistd.h>

#if defined(_WIN32)
  #if defined( psapi_FOUND)
    // to gain memory tracking on windows, link with the psapi library
    #include <windows.h>
    #define PSAPI_VERSION 1
    #include <psapi.h>
  #endif
#elif defined(__APPLE__)
  #include <iostream>
  #include <mach/mach.h>
  #include <stdint.h>
  #include <sys/sysctl.h>
  #include <unistd.h>
#endif

namespace bcl
{
  namespace util
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor, determines values if they are available
    MemoryUsage::MemoryUsage() :
      m_PeakVirtual( GetUndefined< size_t>()),
      m_Virtual( GetUndefined< size_t>()),
      m_PeakRAM( GetUndefined< size_t>()),
      m_RAM( GetUndefined< size_t>())
    {
      const size_t bytes_per_megabyte( 1 << 20);
  #if defined(_WIN32)
    #if defined( psapi_FOUND)
      // windows uses PROCESS_MEMORY_COUNTERS; which require linking with psapi.lib
      PROCESS_MEMORY_COUNTERS mem_counter;
      if( GetProcessMemoryInfo( GetCurrentProcess(), &mem_counter, sizeof( PROCESS_MEMORY_COUNTERS)))
      {
        m_Virtual     = mem_counter.PagefileUsage / bytes_per_megabyte;
        m_PeakVirtual = mem_counter.PeakPagefileUsage / bytes_per_megabyte;
        m_PeakRAM     = mem_counter.PeakWorkingSetSize / bytes_per_megabyte;
        m_RAM         = mem_counter.WorkingSetSize / bytes_per_megabyte;
      }
    #endif
  #elif defined(__APPLE__)
      // on apple, getrusage returns the peak memory
      rusage usage;
      getrusage( RUSAGE_SELF, &usage);
      m_PeakRAM = usage.ru_maxrss / bytes_per_megabyte;

      // task_info must be used to get the RAM and virtual memory size
      struct task_basic_info_64 task_information;
      mach_msg_type_number_t count( TASK_BASIC_INFO_64_COUNT);

      task_info
      (
        mach_task_self(),
        TASK_BASIC_INFO_64,
        task_info_t( &task_information),
        &count
      );

      m_RAM = task_information.resident_size / bytes_per_megabyte;

      // the virtual memory size on apple is obtained from task_information.virtual_size
      // note that it includes a large shared block (usually around 630MB)
      m_Virtual = task_information.virtual_size / bytes_per_megabyte;

  #else
      // on unix systems, there are several ways to get the memory consumption of a process
      // mallinfo has the downside that it does not work above 4GB and is generally deprecated
      // rusage/getrusage has the downside that Linux does not insert the memory values into the rusage struct
      // malloc_stats has the downside that it's a compiler-specific extension and can only be used to
      //   write the statistics to std::cerr, because it's format and level of information are likely to change
      // the way chosen here is to parse the same file parsed by top, htop, ps, and similar commands; which
      // live at /proc/{pid}/statm, /proc/{pid}/stat, /proc/{pid}/status.
      // /proc/{pid}/status was selected as the best choice because it has labeled fields, so new fields can be added
      // without changing the code below, and because it has the peak virtual and RAM levels
      io::IFStream input;

      const std::string filename( "/proc/" + util::Format()( getpid()) + "/status");
      // try to open the file at /proc/process-id/status.  This contains all the information desired
      if( !io::File::TryOpenIFStream( input, filename))
      {
        // could not open file, just return
        return;
      }
      // load in all the lines
      storage::Vector< std::string> lines( StringLineListFromIStream( input));
      // close the file stream
      io::File::CloseClearFStream( input);

      // set the field to value map to undefined for all the fields needed for this object
      storage::Map< std::string, size_t> field_to_value;
      field_to_value[ "VmPeak"] = util::GetUndefined< size_t>();
      field_to_value[ "VmSize"] = util::GetUndefined< size_t>();
      field_to_value[ "VmHWM"]  = util::GetUndefined< size_t>();
      field_to_value[ "VmRSS"]  = util::GetUndefined< size_t>();

      for
      (
        storage::Vector< std::string>::const_iterator itr( lines.Begin()), itr_end( lines.End());
        itr != itr_end;
        ++itr
      )
      {
        // find the position of the first ':'
        const size_t field_delimiter_pos( itr->find( ':'));

        // if the field delimiter was not found, skip this line
        if( field_delimiter_pos >= itr->size())
        {
          continue;
        }

        // find the first non-space character on the line
        const size_t field_start( itr->find_first_not_of( " \t"));

        // get the field
        const std::string field( itr->substr( field_start, field_delimiter_pos - field_start));

        // check if the field is needed
        if( !field_to_value.Has( field))
        {
          // unused field, skip to the next line
          continue;
        }

        // field was desired, parse the value
        std::string memory_value_with_unit( TrimString( itr->substr( field_delimiter_pos + 1)));

        // get the length of the memory size
        const size_t memory_size_num_chars( LengthOfUnsignedIntegerType( memory_value_with_unit));

        // load the memory size
        size_t memory_size
        (
          ConvertStringToNumericalValue< size_t>( memory_value_with_unit.substr( 0, memory_size_num_chars))
        );

        // get the unit, if there was one; otherwise, assume it is bytes
        if( memory_size_num_chars != memory_value_with_unit.size())
        {
          // look at the first letter after the number
          const size_t first_unit_letter_index( memory_value_with_unit.find_first_not_of( " \t", memory_size_num_chars));

          // cast it to lower case
          const char first_unit_letter( tolower( int( memory_value_with_unit[ first_unit_letter_index])));

          // convert the units into mb
          switch( first_unit_letter)
          {
            case 'b':
              // bytes, convert into mb
              memory_size /= bytes_per_megabyte;
              break;
            case 'k':
              // kilobytes
              memory_size /= size_t( 1024);
              break;
            case 'm':
              // megabytes
              break;
            case 'g':
              memory_size *= size_t( 1024);
              // gigabytes
              break;
            case 't':
              memory_size *= bytes_per_megabyte;
              // terabytes
              break;
            default:
              // unknown units
              BCL_MessageCrt
              (
                "Cannot read memory value: " + memory_value_with_unit
              );
              continue;
              break;
          }
        }
        else
        {
          // assume the value was in bytes, divide by bytes_per_megabyte to convert to mb
          memory_size /= bytes_per_megabyte;
        }
        field_to_value[ field] = memory_size;
      }

      // load the values for the desired fields into the appropriate size_t
      m_PeakVirtual = field_to_value[ "VmPeak"];
      m_Virtual = field_to_value[ "VmSize"];
      m_PeakRAM = field_to_value[ "VmHWM"];
      m_RAM = field_to_value[ "VmRSS"];
  #endif
    }

    //! @brief Clone function
    //! @return pointer to new MemoryUsage
    MemoryUsage *MemoryUsage::Clone() const
    {
      return new MemoryUsage( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MemoryUsage::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the peak virtual memory used, in megabytes
    //! @return the peak virtual memory used, in megabytes
    size_t MemoryUsage::GetPeakVirtualMemoryUsed() const
    {
      return m_PeakVirtual;
    }

    //! @brief get the virtual memory in use, in megabytes
    //! @return the virtual memory in use, in megabytes
    size_t MemoryUsage::GetVirtualMemoryInUse() const
    {
      return m_Virtual;
    }

    //! @brief get the peak RAM used, in megabytes
    //! @return the peak RAM used, in megabytes
    size_t MemoryUsage::GetPeakRAMUsed() const
    {
      return m_PeakRAM;
    }

    //! @brief get the RAM in use, in megabytes
    //! @return the RAM in use, in megabytes
    size_t MemoryUsage::GetRAMInUse() const
    {
      return m_RAM;
    }

    //! @brief write peak memory usage, if available, to an output stream
    //! @param STREAM stream to write the output to
    //! @return stream reference
    std::ostream &MemoryUsage::WritePeakMemoryUsageInfo( std::ostream &STREAM)
    {
      MemoryUsage usage;
      if( IsDefined( usage.m_PeakVirtual))
      {
        STREAM << "Peak virtual memory used: " << usage.m_PeakVirtual << " MB\n";
      }
      if( IsDefined( usage.m_PeakRAM))
      {
        STREAM << "Peak RAM used: " << usage.m_PeakRAM << " MB\n";
      }
      return STREAM;
    }

    //! @brief write peak memory usage, if available, to an output stream
    //! @param STREAM stream to write the output to
    //! @return stream reference
    std::ostream &MemoryUsage::WriteCurrentMemoryUsageInfo( std::ostream &STREAM)
    {
      MemoryUsage usage;
      if( IsDefined( usage.m_Virtual))
      {
        STREAM << "Virtual memory used: " << usage.m_Virtual << " MB\n";
      }
      if( IsDefined( usage.m_RAM))
      {
        STREAM << "RAM used: " << usage.m_RAM << " MB\n";
      }
      return STREAM;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MemoryUsage::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_PeakVirtual, ISTREAM);
      io::Serialize::Read( m_Virtual, ISTREAM);
      io::Serialize::Read( m_PeakRAM, ISTREAM);
      io::Serialize::Read( m_RAM, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MemoryUsage::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_PeakVirtual, OSTREAM, INDENT);
      io::Serialize::Write( m_Virtual, OSTREAM, INDENT);
      io::Serialize::Write( m_PeakRAM, OSTREAM, INDENT);
      io::Serialize::Write( m_RAM, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace util
} // namespace bcl
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
#include "util/bcl_util_message.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

    //! @brief access to the only instance of the message class
    //! @return reference to MessageObject
    Message &GetMessenger()
    {
      static Message *s_message( new Message());
      return *s_message;
    }

  ///////////
  // enums //
  ///////////

    //! @brief MessageLevel as string
    //! @param MESSAGE_LEVEL the message level
    //! @return the MessageLevel as string
    const std::string &Message::GetLevelString( const Message::MessageLevel &MESSAGE_LEVEL)
    {
      static const std::string s_message_level_strings[] =
      {
        "Error", "Silent", "Critical", "Standard", "Verbose", "Debug",
        GetStaticClassName< MessageLevel>()
      };

      return s_message_level_strings[ MESSAGE_LEVEL];
    }

    //! @brief MessageVerbosity as string
    //! @param MESSAGE_VERBOSITY the message verbosity
    //! @return the MessageVerbosity as string
    const std::string &Message::GetVerbosityString( const Message::MessageVerbosity &MESSAGE_VERBOSITY)
    {
      static const std::string s_message_verbosity_strings[] =
      {
        "Summary", "Detail",
        GetStaticClassName< MessageVerbosity>()
      };

      return s_message_verbosity_strings[ MESSAGE_VERBOSITY];
    }

  //////////
  // data //
  //////////

    //! @brief command line flag to be used to set MessageLevel over the command line
    //! @return ShPtr to a FlagInterface which is used to set MessageLevel
    const ShPtr< command::FlagInterface> &Message::GetMessageLevelFlag()
    {
      static ShPtr< command::FlagInterface> s_flag_message_level;

      // first time function is called
      if( !s_flag_message_level.IsDefined())
      {
        s_flag_message_level = ShPtr< command::FlagInterface>
        (
          new command::FlagStatic
          (
            "message_level",
            "adjust the MessageLevel",
            ShPtrVector< command::ParameterInterface>::Create
            (
              GetMessageLevelParam(),
              GetMessageVerbosityParam()
            ),
            &Message::UpdateMessengerFromCommandLine
          )
        );
      }

      // end
      return s_flag_message_level;
    }

    //! @brief command line parameter to be used to set MessageLevel over the command line
    //! @return ShPtr to a parameter which is used to set MessageLevel
    const ShPtr< command::ParameterInterface> &Message::GetMessageLevelParam()
    {
      static const ShPtr< command::ParameterInterface> s_message_level_param
      (
        new command::Parameter
        (
          "level",
          "minimum level of messages that will be printed",
          command::ParameterCheckSerializable( MessageLevelEnum()),
          GetLevelString( e_Standard)
        )
      );

      // return
      return s_message_level_param;
    }

    //! @brief command line parameter to be used to set MessageVerbosity over the command line
    //! @return ShPtr to a parameter which is used to set MessageVerbosity
    const ShPtr< command::ParameterInterface> &Message::GetMessageVerbosityParam()
    {
      static const ShPtr< command::ParameterInterface> s_message_verbosity_param
      (
        new command::Parameter
        (
          "verbosity",
          "set to Detail to print the source file and line of origination for each message",
          command::ParameterCheckSerializable( MessageVerbosityEnum()),
          GetVerbosityString( e_Summary)
        )
      );

      // return
      return s_message_verbosity_param;
    }

  //////////
  // data //
  //////////

    //! @brief 3-char strings for the message level to be used when outputting messages to indicate level
    //! @return array with string for each MessageLevel enum
    const std::string *Message::GetMessageLevelTags()
    {
      //! @brief message level codes for "tagging" output messages
      static std::string *s_message_level_tags( NULL);
      if( s_message_level_tags == NULL)
      {
        s_message_level_tags = new std::string[ 6];
        s_message_level_tags[ 0] = "err";
        s_message_level_tags[ 1] = "slt";
        s_message_level_tags[ 2] = "crt";
        s_message_level_tags[ 3] = "std";
        s_message_level_tags[ 4] = "vrb";
        s_message_level_tags[ 5] = "dbg";
      }

      return s_message_level_tags;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief empty constructor - setting defaults
    Message::Message() :
      m_CurrentMessageLevel( e_Standard),
      m_CurrentMessageVerbosity( e_Summary)
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief store the current MessageLevel
    //! initialize the default message level to e_Standard
    Message::MessageLevel Message::GetCurrentMessageLevel() const
    {
      return m_CurrentMessageLevel;
    }

    //! @brief stores the MessageVerbosity flag
    //! initialize the default message verbosity to e_Summary
    Message::MessageVerbosity Message::GetMessageVerbosity() const
    {
      return m_CurrentMessageVerbosity;
    }

    //! @brief set MessageLevel to a new MESSAGE_LEVEL
    //! @param MESSAGE_LEVEL level that should be the current one
    void Message::SetMessageLevel( const Message::MessageLevel MESSAGE_LEVEL)
    {
      m_CurrentMessageLevel = MESSAGE_LEVEL;
    }

    //! @brief set MessageVerbosity to a new MESSAGE_VERBOSITY
    //! @param MESSAGE_VERBOSITY level of message output: summary or detail
    void Message::SetMessageVerbosity( const Message::MessageVerbosity MESSAGE_VERBOSITY)
    {
      m_CurrentMessageVerbosity = MESSAGE_VERBOSITY;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief sets the MessageLevel to what was given in the command line
    void Message::SetMessageLevelFromCommandLineFlag()
    {
      // level from command line
      const MessageLevelEnum new_level( GetMessageLevelParam()->GetValue());

      // check if it was changed
      if( m_CurrentMessageLevel != new_level)
      {
        SetMessageLevel( new_level);
      }
    }

    //! @brief sets the MessageLevel to what was given in the command line
    void Message::SetMessageVerbosityFromCommandLineFlag()
    {
      // verbosity from command line
      const MessageVerbosityEnum new_verbosity( GetMessageVerbosityParam()->GetValue());

      // check if it was changed
      if( m_CurrentMessageVerbosity != new_verbosity)
      {
        SetMessageVerbosity( new_verbosity);
      }
    }

    //! @brief is argument MessageLevel smaller or equal current MessageLevel
    //! @param MESSAGE_LEVEL the messagelevel that should be compared
    //! @return true if the MESSAGE_LEVEL is smaller or equal than the current MessageLevel
    bool Message::IsSmallerEqualCurrentMessageLevel( const MessageLevel MESSAGE_LEVEL) const
    {
      return MESSAGE_LEVEL <= m_CurrentMessageLevel;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief updates the message level and verbosity from the command line flags for the global messenger
    void Message::UpdateMessengerFromCommandLine()
    {
      GetMessenger().SetMessageLevelFromCommandLineFlag();
      GetMessenger().SetMessageVerbosityFromCommandLineFlag();
    }

    //! @brief send a message from a namespace to the user
    //! @param NAMESPACE the name of the current namespace
    //! @param MESSAGE the user message
    //! @return the stream where the message was written to
    std::ostream &Message::SendToUser( const std::string &NAMESPACE, const std::string &MESSAGE)
    {
      // write namespace, "> " and the message
      GetLogger().LogMessage( NAMESPACE + "> " + MESSAGE);
      return GetLogger();
    }

    //! @brief a message is outputted
    //! @param MESSAGE_LEVEL is the level of severity of the message
    //! @param NAMESPACE the namespace
    //! @param MESSAGE the user Message
    //! @return the stream where the message was written to
    std::ostream &Message::SendToUser
    (
      const MessageLevel &MESSAGE_LEVEL,
      const std::string &NAMESPACE,
      const std::string &MESSAGE
    )
    {
      // write arrow and message
      GetLogger().LogMessage( AssembleMessage( MESSAGE_LEVEL, NAMESPACE, MESSAGE));

      return GetLogger();
    }

    //! @brief a message is output with the current filename and functionname and line
    //! @param MESSAGE_LEVEL is the level of severity of the message
    //! @param NAMESPACE the namespace
    //! @param MESSAGE the user Message
    //! @param FILE_NAME the name of the file, where the Message came from
    //! @param LINE_NUMBER the line number the message came from
    //! @param FUNCTION_NAME the name of the function where the message came from
    //! @return the stream where the message was written to
    std::ostream &Message::SendToUser
    (
      const MessageLevel &MESSAGE_LEVEL,
      const std::string &NAMESPACE,
      const std::string &MESSAGE,
      const char *FILE_NAME,
      const int LINE_NUMBER,
      const char *FUNCTION_NAME
    )
    {
      if( GetMessenger().GetMessageVerbosity() == Message::e_Detail)
      {
        // remove path from filename
        std::string file_name( FILE_NAME);
        file_name = file_name.substr( file_name.find_last_of( PATH_SEPARATOR) + 1);
        // output message
        GetLogger().LogMessage
        (
          AssembleMessage( MESSAGE_LEVEL, NAMESPACE, MESSAGE) +
          ", function: " + FUNCTION_NAME + " | file: " + file_name + " | line: " + Format()( LINE_NUMBER)
        );
      }
      else
      {
        GetLogger().LogMessage( AssembleMessage( MESSAGE_LEVEL, NAMESPACE, MESSAGE));
      }

      return GetLogger();
    }

    //! @brief assemble a message from Namespace and Message
    //! @param NAMESPACE the namespace the message comes from
    //! @param MESSAGE the actual message
    //! @return string of the for {namespace}=>{message}
    std::string Message::AssembleMessage
    (
      const std::string &NAMESPACE,
      const std::string &MESSAGE
    )
    {
      return std::string( NAMESPACE + "=> " + MESSAGE);
    }

    //! @brief assemble a message from MessageLevel and Message
    //! @param MESSAGE_LEVEL the level for this message
    //! @param NAMESPACE the namespace the message comes from
    //! @param MESSAGE the actual message
    //! @return string of the for {level}={namespace}=>{message}
    std::string Message::AssembleMessage
    (
      const MessageLevel &MESSAGE_LEVEL,
      const std::string &NAMESPACE,
      const std::string &MESSAGE
    )
    {
      return std::string( GetMessageLevelTags()[ MESSAGE_LEVEL] + "=" + AssembleMessage( NAMESPACE, MESSAGE));
    }

  } // namespace util
} // namespace bcl
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

// include forward header of this class
#include "util/bcl_util_object_data_label.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_object_data_label_tokenizer.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically
#include <list>

namespace bcl
{
  namespace util
  {

  //////////
  // data //
  //////////

    //! @brief compares two object data labels by name
    //! @param A, B object data labels whose names should be compared
    //! @return true if the name of A is less than the name of B
    bool NameIsLessThan( const ObjectDataLabel &A, const ObjectDataLabel &B)
    {
      return A.GetName() < B.GetName();
    }

    //! @brief compares two object data label pointers by name
    //! @param A, B object data labels whose names should be compared
    //! @return true if the name of A is less than the name of B
    bool PtrNameIsLessThan( const ObjectDataLabel *const &A, const ObjectDataLabel *const &B)
    {
      return A->GetName() < B->GetName();
    }

    //! @brief default constructor
    ObjectDataLabel::ObjectDataLabel()
    {
    }

    //! @brief constructor from a string representing a label that is to be parsed
    //! @param LABEL the label (what does the name describe)
    ObjectDataLabel::ObjectDataLabel( const std::string &STR)
    {
      ObjectDataLabel::operator =( STR);
    }

    //! @brief constructor from a stream containing a data label
    //! @param STREAM stream to pass to get data label
    ObjectDataLabel::ObjectDataLabel( std::istream &STREAM)
    {
      long parenthesis_depth( 0);
      std::string line, new_line;

      while( ( line.size() == 0 || parenthesis_depth > 0) && STREAM.good())
      {
        STREAM >> std::ws;
        std::getline( STREAM, new_line);

        // this if-statement is a hack to allow reading in files that use bcl class identifiers
        // remove once these files have been updated
        if( StartsWith( new_line, "bcl::"))
        {
          continue;
        }

        // Skip comments
        if( StartsWith( new_line, "#"))
        {
          continue;
        }

        line += new_line;
        parenthesis_depth += GetParenthesisDepthChange( new_line);
      }
      ObjectDataLabel::operator =( line);
    }

    //! @brief constructor from arguments
    //! @param ARGS arguments for the data label
    ObjectDataLabel::ObjectDataLabel( const storage::Vector< ObjectDataLabel> &ARGS) :
      m_Name(),
      m_Value(),
      m_Arguments( ARGS.Begin(), ARGS.End())
    {
      SortArgumentsByName();
    }

    //! @brief constructor from a string representing a name and the arguments
    //! @param NAME the name of the object
    //! @param ARGS arguments for the data label
    ObjectDataLabel::ObjectDataLabel( const std::string &NAME, const storage::Vector< ObjectDataLabel> &ARGS) :
      m_Value( NAME),
      m_Arguments( ARGS.Begin(), ARGS.End())
    {
      SortArgumentsByName();
    }

    //! @brief constructor from a string representing a name and the arguments
    //! @param NAME the name of the object
    //! @param ARGS arguments for the data label
    ObjectDataLabel::ObjectDataLabel( const std::string &NAME, const std::vector< ObjectDataLabel> &ARGS) :
      m_Value( NAME),
      m_Arguments( ARGS)
    {
      SortArgumentsByName();
    }

    //! @brief constructor from a string representing a name, unparsed value, and type
    //! @param NAME identifies what the value represents
    //! @param VALUE the value of the string, which should be basic
    //! @param TYPE the type
    ObjectDataLabel::ObjectDataLabel
    (
      const std::string &NAME,
      const std::string &VALUE
    ) :
      m_Name( NAME),
      m_Value( VALUE)
    {
    }

    //! @brief constructor from a string representing a label, name and the arguments (implies type is object)
    //! @param NAME identifies what the value represents
    //! @param OBJECT_NAME the name or alias of the object's type
    //! @param ARGS arguments for the object
    ObjectDataLabel::ObjectDataLabel
    (
      const std::string &NAME,
      const std::string &OBJECT_NAME,
      const storage::Vector< ObjectDataLabel> &ARGS
    ) :
      m_Name( NAME),
      m_Value( OBJECT_NAME),
      m_Arguments( ARGS.Begin(), ARGS.End())
    {
      SortArgumentsByName();
    }

    //! @brief constructor from a string representing a label, name and a single argument (implies type is object)
    //! @param NAME identifies what the value represents
    //! @param OBJECT_NAME the name or alias of the object's type
    //! @param ARG argument for the object
    ObjectDataLabel::ObjectDataLabel
    (
      const std::string &NAME,
      const std::string &OBJECT_NAME,
      const ObjectDataLabel &ARG
    ) :
      m_Name( NAME),
      m_Value( OBJECT_NAME),
      m_Arguments( size_t( 1), ARG)
    {
    }

    //! @brief constructor from a string with a name and underlying label (implies type is object)
    //! @param NAME name given to the object
    //! @param OBJECT_LABEL the data label from the object
    ObjectDataLabel::ObjectDataLabel( const std::string &NAME, const ObjectDataLabel &OBJECT_LABEL) :
      m_Name(),
      m_Value( OBJECT_LABEL.m_Value),
      m_Arguments( OBJECT_LABEL.m_Arguments)
    {
      SetName( NAME);
      CopyArgumentSorting( OBJECT_LABEL);
    }

    //! @brief copy constructor
    //! @param PARENT the original data label
    ObjectDataLabel::ObjectDataLabel( const ObjectDataLabel &PARENT) :
      m_Name( PARENT.m_Name),
      m_Value( PARENT.m_Value),
      m_Arguments( PARENT.m_Arguments)
    {
      CopyArgumentSorting( PARENT);
    }

    //! @brief Clone function
    //! @return pointer to new ObjectDataLabel
    ObjectDataLabel *ObjectDataLabel::Clone() const
    {
      return new ObjectDataLabel( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ObjectDataLabel::GetClassIdentifier() const
    {
      return GetStaticClassName< ObjectDataLabel>();
    }

    //! @brief get the character that goes between each argument
    const char &ObjectDataLabel::GetArgumentDelimiter()
    {
      static const char s_Delimiter( ',');
      return s_Delimiter;
    }

    //! @brief get the character that goes between a label and its name / value
    //! @return the delimiter character
    const char &ObjectDataLabel::GetNameValueDelimiter()
    {
      static const char s_Delimiter( '=');
      return s_Delimiter;
    }

    //! @brief get the character that indicates to inline another file in this label
    //! @return the inline file character
    const char &ObjectDataLabel::GetInlineFileDelimiter()
    {
      static const char s_Delimiter( '@');
      return s_Delimiter;
    }

    //! @brief get all the delimiters
    //! @return a string containing all the delimiters
    const std::string &ObjectDataLabel::GetAllDelimiters()
    {
      static const std::string s_delimiters
      (
        "(\")" + std::string( 1, GetArgumentDelimiter()) + std::string( 1, GetNameValueDelimiter())
         + std::string( 1, GetInlineFileDelimiter())
      );
      return s_delimiters;
    }

    //! @return the primary identifier of this data label (either the property name or parameter value)
    const std::vector< ObjectDataLabel> &ObjectDataLabel::GetArguments() const
    {
      return m_Arguments;
    }

    //! @brief test whether this label is scalar (e.g. does not have any arguments)
    //! @return true if this label is scalar (e.g. does not have any arguments)
    bool ObjectDataLabel::IsScalar() const
    {
      return m_Arguments.empty();
    }

    //! @brief test whether this label is empty (e.g. does not have any arguments, name, or value)
    //! @return true if this label is empty (e.g. does not have any arguments, name, or value)
    bool ObjectDataLabel::IsEmpty() const
    {
      return m_Arguments.empty() && m_Name.empty() && m_Value.empty();
    }

    //! @brief get the number of arguments for this label
    //! @return the number of arguments for this label
    size_t ObjectDataLabel::GetNumberArguments() const
    {
      return m_Arguments.size();
    }

    //! @return the arguments of this data label (either the property name or parameter value)
    const ObjectDataLabel &ObjectDataLabel::GetArgument( const size_t &NUMBER) const
    {
      return m_Arguments[ NUMBER];
    }

    //! @brief set the tag / identifier of this data label (e.g. lhs=rhs -> lhs is the name)
    //! @param NEW_NAME new name
    //! @param ALLOW_CHANGE true to allow setting the name if this label was already given a non-empty name
    void ObjectDataLabel::SetName( const std::string &NEW_NAME, const bool &ALLOW_CHANGE)
    {
      std::string &target( m_Value.empty() && !IsScalar() ? m_Value : m_Name);
      if( target != NEW_NAME)
      {
        BCL_Assert
        (
          ALLOW_CHANGE || target.empty(),
          "Cannot change name (to: " + NEW_NAME + ") on label that was already named: " + ToNamedString()
        );
        target = NEW_NAME;
      }
    }

    //! @brief set the value of this data label (e.g. lhs=rhs -> rhs is the value)
    //! @param VALUE new value
    //! @param ALLOW_CHANGE true to allow setting the name if this label was already given a non-empty name
    void ObjectDataLabel::SetValue( const std::string &VALUE, const bool &ALLOW_CHANGE)
    {
      if( m_Value != VALUE)
      {
        BCL_Assert
        (
          ALLOW_CHANGE || m_Value.empty(),
          "Cannot change value (to: " + VALUE + ") on label that already had a value (from: " + m_Value + ")"
        );
        m_Value = VALUE;
      }
    }

    //! @brief get the beginning of the arguments
    //! @return the beginning of the arguments
    ObjectDataLabel::const_iterator ObjectDataLabel::Begin() const
    {
      return m_Arguments.begin();
    }

    //! @brief get the end of the arguments
    //! @return the end of the arguments
    ObjectDataLabel::const_iterator ObjectDataLabel::End() const
    {
      return m_Arguments.end();
    }

    //! @brief find an argument with the specified name.  If none is found, return End()
    //! @param NAME the value to find
    //! @param RECURSIVE whether to check sub-arguments (after checking arguments)
    ObjectDataLabel::const_iterator ObjectDataLabel::FindName( const std::string &NAME, const bool &RECURSIVE) const
    {
      // first, check arguments directly
      for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        if( itr->GetName() == NAME)
        {
          return itr;
        }
      }

      // handle recursion
      if( RECURSIVE)
      {
        for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
        {
          const_iterator itr_find( itr->FindName( NAME, true));
          if( itr_find != itr->End())
          {
            return itr_find;
          }
        }
      }

      return End();
    }

    //! @brief find an argument with the specified name.  If none is found, return End()
    //! @param NAME the value to find
    //! @param RECURSIVE whether to check sub-arguments (after checking arguments)
    ObjectDataLabel::const_iterator ObjectDataLabel::FindValue( const std::string &VALUE, const bool &RECURSIVE) const
    {
      // first, check arguments directly
      for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        if( itr->GetValue() == VALUE)
        {
          return itr;
        }
      }

      // handle recursion
      if( RECURSIVE)
      {
        for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
        {
          const_iterator itr_find( itr->FindValue( VALUE, true));
          if( itr_find != itr->End())
          {
            return itr_find;
          }
        }
      }

      return End();
    }

    //! @brief replace values == VALUE with REPLACEMENT, up to MAX_DEPTH below this argument
    //! @param VALUE the value to find
    //! @param REPLACEMENT what to replace the value with
    //! @param MAX_DEPTH maximum depth to perform replacement at(0 = just this object data label)
    //! @return the number of replacements performed
    size_t ObjectDataLabel::ReplaceValue
    (
      const std::string &VALUE,
      const std::string &REPLACEMENT,
      const size_t &MAX_DEPTH
    )
    {
      size_t n_replaced( 0);
      if( m_Value == VALUE)
      {
        if( IsScalar())
        {
          *this = ObjectDataLabel( m_Name, ObjectDataLabel( REPLACEMENT));
        }
        else
        {
          m_Value = REPLACEMENT;
        }
        ++n_replaced;
      }
      if( !MAX_DEPTH)
      {
        return n_replaced;
      }
      for( iterator itr( m_Arguments.begin()), itr_end( m_Arguments.end()); itr != itr_end; ++itr)
      {
        n_replaced += itr->ReplaceValue( VALUE, REPLACEMENT, MAX_DEPTH - 1);
      }
      return n_replaced;
    }

    //! @brief find an argument with the specified object data-label.  If none is found, return End()
    //! @param LABEL the label to find
    //! @param RECURSIVE whether to check sub-arguments (after checking arguments)
    //! @param CONSIDER_NAME whether to force a match on the name/tag given in LABEL; otherwise, any matching value
    //!        will be found
    ObjectDataLabel::const_iterator ObjectDataLabel::Find
    (
      const ObjectDataLabel &LABEL,
      const bool &RECURSIVE,
      const bool &CONSIDER_NAME
    ) const
    {
      if( CONSIDER_NAME)
      {
        // first, check arguments directly.
        for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
        {
          if( itr->TernaryCompare( LABEL) == 0)
          {
            return itr;
          }
        }
      }
      else
      {
        for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
        {
          if( itr->TernaryCompareWithoutName( LABEL) == 0)
          {
            return itr;
          }
        }
      }

      // handle recursion; name will always be considered on sub-arguments as well
      if( RECURSIVE)
      {
        for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
        {
          const_iterator itr_find( itr->Find( LABEL, true, CONSIDER_NAME));
          if( itr_find != itr->End())
          {
            return itr_find;
          }
        }
      }

      return End();
    }

    //! @brief compute the difference between this object data label and another
    //! @param OTHER the other object data label
    ObjectDataLabel ObjectDataLabel::Difference( const ObjectDataLabel &OTHER) const
    {
      if( !this->TernaryCompareWithoutName( OTHER))
      {
        return ObjectDataLabel();
      }
      if( m_Value != OTHER.m_Value)
      {
        return ObjectDataLabel( "", m_Value, m_Arguments);
      }
      storage::Vector< ObjectDataLabel> differing_arguments;
      std::list< ObjectDataLabel> other( OTHER.Begin(), OTHER.End());
      std::vector< const_iterator> nonidentical_arguments;
      // first, remove exact matches
      for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        std::list< ObjectDataLabel>::iterator itr_other_found( std::find( other.begin(), other.end(), *itr));
        if( itr_other_found != other.end())
        {
          other.erase( itr_other_found);
          continue;
        }
        // save the argument because it is not identical
        nonidentical_arguments.push_back( itr);
      }

      // compute difference for each argument using a greedy algorithm
      for
      (
        std::vector< const_iterator>::const_iterator
          itr_itr( nonidentical_arguments.begin()), itr_itr_end( nonidentical_arguments.end());
        itr_itr != itr_itr_end;
        ++itr_itr
      )
      {
        const ObjectDataLabel &label( **itr_itr);
        ObjectDataLabel best( label);
        size_t best_size( best.GetTreeSize());
        std::list< ObjectDataLabel>::iterator itr_best( other.end());
        for
        (
          std::list< ObjectDataLabel>::iterator itr_list( other.begin()), itr_list_end( other.end());
          itr_list != itr_list_end;
          ++itr_list
        )
        {
          // skip differently named arguments
          if( label.GetName() != itr_list->GetName())
          {
            ++itr_list;
            continue;
          }
          ObjectDataLabel difference_label( label.Difference( *itr_list));
          size_t current_size( difference_label.GetTreeSize());
          if( current_size && itr_list->GetValue() == label.GetValue())
          {
            current_size -= 1;
          }
          if( current_size < best_size)
          {
            best = difference_label;
            itr_best = itr_list;
            best_size = current_size;
          }
        }
        if( itr_best != other.end())
        {
          other.erase( itr_best);
        }
        differing_arguments.PushBack( best);
      }
      return ObjectDataLabel( m_Name, m_Value, differing_arguments);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief assign from a data-label string
    ObjectDataLabel &ObjectDataLabel::operator =( const std::string &LABEL)
    {
      // make sure that the parenthesis are balanced
      BCL_Assert( TryAssign( LABEL, GetLogger()), "invalid string");

      return *this;
    }

    //! @brief assignment from another object data label
    //! @param LABEL the other label
    //! @return a reference to this
    ObjectDataLabel &ObjectDataLabel::operator =( const ObjectDataLabel &LABEL)
    {
      m_Name = LABEL.m_Name;
      m_Value = LABEL.m_Value;
      m_Arguments = LABEL.m_Arguments;
      CopyArgumentSorting( LABEL);
      return *this;
    }

    //! @brief try to assign the label from a string
    //! @param STR string to attempt assignment from
    //! @param ERR error stream, stream to output errors to
    bool ObjectDataLabel::TryAssign( const std::string &STR, std::ostream &ERR)
    {
      m_Arguments.clear();       // remove any existing arguments
      m_Value.erase();           // and any existing value
      m_Name.erase();            // and any existing name
      m_SortedArguments.clear(); // and sorting order

      // parse the label
      ObjectDataLabelTokenizer tokenizer( STR);
      if( tokenizer.GetNextTokenType() != ObjectDataLabelTokenizer::e_End)
      {
        if( !ReadSubLabel( tokenizer, ERR))
        {
          return false;
        }
        if( tokenizer.GetNextTokenType() != ObjectDataLabelTokenizer::e_End)
        {
          switch( tokenizer.GetNextTokenType())
          {
            case ObjectDataLabelTokenizer::e_ArgDelimiter:
              ERR << "Arguments are not allowed without an associated object, that is, ',' cannot appear outside\n"
                  << "parenthetical scope for object data labels, but did in given string:\n";
              break;
            case ObjectDataLabelTokenizer::e_ScopeClose:
              ERR << "Excessive end-scope parenthesis in given string:\n";
              break;
            case ObjectDataLabelTokenizer::e_TagDelimiter:
              ERR << "Map object keys must be quoted but were not in given string:\n";
              break;
            default:
              ERR << "Multiple labels cannot appear on the same line but did in given string:\n";
              break;
          }
          ERR << STR << std::endl;
          return false;
        }
      }
      return true;
    }

    //! @brief ternary comparison operator (like strcmp), returns -1 if *this < A, 0 if *this == A, 1 if *this > A
    //! @param A the object data label to compare with
    //! @return -1 if *this < A, 0 if *this == A, 1 if *this > A
    int ObjectDataLabel::TernaryCompare( const ObjectDataLabel &LABEL) const
    {
      // compare names
      if( GetName() != LABEL.GetName())
      {
        return GetName() < LABEL.GetName() ? -1 : 1;
      }
      return TernaryCompareWithoutName( LABEL);
    }

    //! @brief Compare this label to another without considering the name/tag
    //! @param LABEL the object data label to compare with
    //! @return -1 if *this < LABEL, 0 if *this == LABEL, 1 if *this > LABEL
    int ObjectDataLabel::TernaryCompareWithoutName( const ObjectDataLabel &LABEL) const
    {
      // compare values
      if( GetValue() != LABEL.GetValue())
      {
        return GetValue() < LABEL.GetValue() ? -1 : 1;
      }
      // compare # arguments
      else if( GetNumberArguments() != LABEL.GetNumberArguments())
      {
        return GetNumberArguments() < LABEL.GetNumberArguments() ? -1 : 1;
      }
      else if( IsScalar())
      {
        return 0;
      }

      // consider all four possible cases:
      // 1. Neither this nor LABEL has an allocated m_SortedArguments vector
      // 2. this has an allocated m_SortedArguments vector but LABEL does not
      // 3. LABEL has an allocated m_SortedArguments vector but this does not
      // 4. both this and LABEL have an allocated m_SortedArgumentsVector
      if( m_SortedArguments.empty() == LABEL.m_SortedArguments.empty())
      {
        // case 1 & 4
        if( m_SortedArguments.empty())
        {
          // case 1
          for
          (
            const_iterator
              itr_a( m_Arguments.begin()), itr_b( LABEL.m_Arguments.begin()), itr_a_end( m_Arguments.end());
            itr_a != itr_a_end;
            ++itr_a, ++itr_b
          )
          {
            // compare argument a with argument b
            if( int compare_value = itr_a->TernaryCompare( *itr_b)) // values not equal
            {
              return compare_value;
            }
          }
        }
        else
        {
          // case 4, both sorted vectors are non-empty
          for
          (
            std::vector< ObjectDataLabel *>::const_iterator
              itr_a( m_SortedArguments.begin()),
              itr_b( LABEL.m_SortedArguments.begin()),
              itr_a_end( m_SortedArguments.end());
            itr_a != itr_a_end;
            ++itr_a, ++itr_b
          )
          {
            // compare argument a with argument b
            if( int compare_value = ( *itr_a)->TernaryCompare( **itr_b)) // values not equal
            {
              return compare_value;
            }
          }
        }
      }
      else
      {
        // disparate vector types, create iterators to the pointer and normal vector and compare them
        int multipler( 0);
        std::vector< ObjectDataLabel *>::const_iterator itr_ptr_arg, itr_ptr_arg_end;
        const_iterator itr_arg;
        if( m_SortedArguments.empty())
        {
          multipler = 1;
          itr_ptr_arg = LABEL.m_SortedArguments.begin();
          itr_ptr_arg_end = LABEL.m_SortedArguments.end();
          itr_arg = m_Arguments.begin();
        }
        else
        {
          multipler = -1;
          itr_ptr_arg = m_SortedArguments.begin();
          itr_ptr_arg_end = m_SortedArguments.end();
          itr_arg = LABEL.m_Arguments.begin();
        }
        while( itr_ptr_arg != itr_ptr_arg_end)
        {
          if( int compare_value = itr_arg->TernaryCompare( **itr_ptr_arg)) // values not equal
          {
            return multipler * compare_value;
          }
          ++itr_ptr_arg;
          ++itr_arg;
        }
      }

      // equal, so return 0
      return 0;
    }

    //! @brief get the string representation of this data label, without the label
    //! @return the property label converted to a string (always succeeds)
    std::string ObjectDataLabel::ToString() const
    {
      if( IsScalar())
      {
        // handle the case where there are no arguments, in which case m_Value must be present
        // write out the name or value, along with the arguments, if applicable, also write quotes if this was a string
        return QuoteStringIfDelimitersPresent( m_Value);
      }
      else if( m_Value.empty())
      {
        // no value, just a sequence, so write it out
        return ArgumentsToString();
      }
      // non-empty value and arguments
      // write out the value along with the arguments
      return QuoteStringIfDelimitersPresent( m_Value) + ArgumentsToString();
    }

    //! @return the property label converted to a string (always succeeds)
    //! @param MAX_DEPTH_FOR_SPLIT maximum depth of parameters to split onto separate lines
    std::string ObjectDataLabel::ToString
    (
      const size_t &LINE_LENGTH,
      const size_t &INDENT,
      const size_t &MAX_DEPTH_FOR_SPLIT
    ) const
    {
      std::ostringstream oss;

      // write out the name or value of this label with the correct indentation
      oss << std::string( INDENT, ' ');
      if( IsScalar())
      {
        oss << QuoteStringIfDelimitersPresent( m_Value);
      }
      else
      {
        // write the value string, but do not put empty quotes before a container
        if( !m_Value.empty())
        {
          oss << QuoteStringIfDelimitersPresent( m_Value);
        }

        // split the lines; ignore the indent because it can cause deeper labels to be split even though the
        // parent labels were not, which impairs readability
        if( MAX_DEPTH_FOR_SPLIT && GetLimitedLength( LINE_LENGTH, false, INDENT) > LINE_LENGTH)
        {
          // add the open parenthesis argument to the line to indicate that arguments follow
          oss << "(\n";

          // for all arguments except the last:
          for
          (
            size_t arg_number( 0), number_internal_arguments( m_Arguments.size() - 1);
            arg_number < number_internal_arguments;
            ++arg_number
          )
          {
            // indent by the proper string, end with the delimiter and a new line
            oss << m_Arguments[ arg_number].ToNamedString( LINE_LENGTH, INDENT + 2, MAX_DEPTH_FOR_SPLIT - 1)
                << GetArgumentDelimiter()
                << '\n';
          }

          // write the last argument, with no argument delimiter afterwards, since it is the last argument
          oss << m_Arguments.back().ToNamedString( LINE_LENGTH, INDENT + 2, MAX_DEPTH_FOR_SPLIT - 1)
              << '\n';

          // write the closing parenthesis
          oss << std::string( INDENT, ' ') << ")";
        }
        else
        {
          oss << ArgumentsToString();
        }
      }

      // return the string from the string stream
      return oss.str();
    }

    //! @brief get the string representation of this data label formatted for the logger
    //! This is equivalent to ToString( GetLogger().GetMaxLineWidth(), INDENT);
    //! @param INDENT amount to indent separate lines by
    //! @param MAX_DEPTH_FOR_SPLIT maximum depth of parameters to split onto separate lines
    //! @return the property label converted to a string (always succeeds)
    std::string ObjectDataLabel::ToStringForLogger( const size_t &INDENT, const size_t &MAX_DEPTH_FOR_SPLIT) const
    {
      return ToString( GetLogger().GetMaxLineWidth(), INDENT, MAX_DEPTH_FOR_SPLIT);
    }

    //! @brief get the string representation of this data label formatted with LoggerInterface::GetDefaultMaxLineWidth
    //! This is equivalent to ToString( LoggerInterface::GetDefaultMaxLineWidth(), INDENT);
    //! @param INDENT amount to indent separate lines by
    //! @param MAX_DEPTH_FOR_SPLIT maximum depth of parameters to split onto separate lines
    //! @return the property label converted to a string (always succeeds)
    std::string ObjectDataLabel::ToStringDefaultWidth( const size_t &INDENT, const size_t &MAX_DEPTH_FOR_SPLIT) const
    {
      return ToString( LoggerInterface::GetDefaultMaxLineWidth(), INDENT, MAX_DEPTH_FOR_SPLIT);
    }

    //! @brief get the string representation of this data label
    //! @return the property label converted to a string (always succeeds)
    std::string ObjectDataLabel::ToNamedString() const
    {
      return m_Name.empty() ? ToString() : QuoteStringIfDelimitersPresent( m_Name) + GetNameValueDelimiter() + ToString();
    }

    //! @brief get the string representation of this data label
    //! @param MAX_DEPTH_FOR_SPLIT maximum depth of parameters to split onto separate lines
    //! @return the property label converted to a string (always succeeds)
    std::string ObjectDataLabel::ToNamedString
    (
      const size_t &LINE_LENGTH,
      const size_t &INDENT,
      const size_t &MAX_DEPTH_FOR_SPLIT
    ) const
    {
      // If there is no name, then use the to-string function
      if( m_Name.empty())
      {
        return ToString( LINE_LENGTH, INDENT, MAX_DEPTH_FOR_SPLIT);
      }
      std::ostringstream oss;

      // write out the name or value of this label with the correct indentation
      oss << std::string( INDENT, ' ');

      // handle names that need quotes due to internal formatting
      oss << QuoteStringIfDelimitersPresent( m_Name) << GetNameValueDelimiter();

      if( IsScalar())
      {
        oss << QuoteStringIfDelimitersPresent( m_Value);
      }
      else
      {
        // write the value string, but do not put empty quotes before a container
        if( !m_Value.empty())
        {
          oss << QuoteStringIfDelimitersPresent( m_Value);
        }

        // split the lines
        if( MAX_DEPTH_FOR_SPLIT && GetLimitedLength( LINE_LENGTH, true, INDENT) > LINE_LENGTH)
        {
          // add the open parenthesis argument to the line to indicate that arguments follow
          oss << "(\n";

          // for all arguments except the last:
          for
          (
            size_t arg_number( 0), number_internal_arguments( m_Arguments.size() - 1);
            arg_number < number_internal_arguments;
            ++arg_number
          )
          {
            // indent by the proper string, end with the delimiter and a new line
            oss << m_Arguments[ arg_number].ToNamedString( LINE_LENGTH, INDENT + 2, MAX_DEPTH_FOR_SPLIT - 1)
                << GetArgumentDelimiter()
                << '\n';
          }

          // write the last argument, with no argument delimiter afterwards, since it is the last argument
          oss << m_Arguments.back().ToNamedString( LINE_LENGTH, INDENT + 2, MAX_DEPTH_FOR_SPLIT - 1)
              << '\n';

          // write the closing parenthesis
          oss << std::string( INDENT, ' ') << ")";
        }
        else // keep the arguments together on a single line
        {
          oss << ArgumentsToString();
        }
      }

      // return string stream
      return oss.str();
    }

    //! @brief get the string representation of this data label formatted for the logger
    //! This is equivalent to ToNamedString( GetLogger().GetMaxLineWidth(), INDENT);
    //! @param INDENT amount to indent separate lines by
    //! @param MAX_DEPTH_FOR_SPLIT maximum depth of parameters to split onto separate lines
    //! @return the property label converted to a string (always succeeds)
    std::string ObjectDataLabel::ToNamedStringForLogger
    (
      const size_t &INDENT,
      const size_t &MAX_DEPTH_FOR_SPLIT
    ) const
    {
      return ToNamedString( GetLogger().GetMaxLineWidth(), INDENT);
    }

    //! @brief get the string representation of this data label formatted with a default max line width
    //! This is equivalent to ToNamedString( LoggerInterface::GetDefaultMaxLineWidth(), INDENT);
    //! @param INDENT amount to indent separate lines by
    //! @param MAX_DEPTH_FOR_SPLIT maximum depth of parameters to split onto separate lines
    //! @return the property label converted to a string (always succeeds)
    std::string ObjectDataLabel::ToNamedStringDefaultWidth
    (
      const size_t &INDENT,
      const size_t &MAX_DEPTH_FOR_SPLIT
    ) const
    {
      return ToNamedString( LoggerInterface::GetDefaultMaxLineWidth(), INDENT, MAX_DEPTH_FOR_SPLIT);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief get the arguments of the data label printed as a string
    //! @param INCLUDE_PARENS whether to surround the string with () (default=true)
    //! @param ARGUMENT_DELIMITER what to delimit the arguments of the top level data label by
    //! @return the arguments, printed as a single-line string, surrounded by (), if desired
    std::string ObjectDataLabel::ArgumentsToString
    (
      const bool &INCLUDE_PARENS,
      const char &ARGUMENT_DELIMITER
    ) const
    {
      std::ostringstream oss;
      if( INCLUDE_PARENS)
      {
        oss << '(';
      }
      if( m_Arguments.size() > 0)
      {
        const_iterator itr( m_Arguments.begin());

        // write out the first argument
        oss << itr->ToNamedString();
        ++itr;

        // then prefix all later arguments with a delimiter
        for
        (
          const_iterator itr_end( m_Arguments.end());
          itr != itr_end;
          ++itr
        )
        {
          oss << ARGUMENT_DELIMITER << itr->ToNamedString();
        }
      }
      if( INCLUDE_PARENS)
      {
        oss << ')';
      }
      return oss.str();
    }

    //! @brief helper function to determine the length of the string returned by ToString
    //! @details used to determine whether to split labels out on separate lines
    //! @param LENGTH the maximum length of interest
    //! @param CONSIDER_NAME whether to consider the length of the name of this label
    //! @param INDENT the indent to subtract from the maximum length
    //! @return the size of the string returned by ToString (or ToNamedString if consider_name is true)
    //!         or any size > length, if the length would exceed LENGTH
    size_t ObjectDataLabel::GetLimitedLength
    (
      const size_t &LENGTH,
      const bool &CONSIDER_NAME,
      const size_t &INDENT
    ) const
    {
      size_t length( INDENT + m_Value.size() + 2);

      // return if length is already too large
      if( length > LENGTH)
      {
        return length;
      }

      // add the length of the name, if desired; also consider possible quotes and = sign
      if( CONSIDER_NAME && !m_Name.empty())
      {
        length += m_Name.size() + 3;
      }

      if( !IsScalar()) // there were arguments
      {
        // add space for the ()
        length += 2;

        // add the length of all the argument delimiters
        // which is 1 less than the # of arguments, since the last argument does not receive a delimiter
        length += m_Arguments.size() - 1;

        // add the lengths of all the arguments, stopping if length gets too long
        for
        (
          size_t arg_number( 0), number_args( m_Arguments.size());
          arg_number < number_args && length <= LENGTH;
          ++arg_number
        )
        {
          // add the length of the length of the next argument, automatically returning if the returned length would
          // be too large.  Don't add an indent to the arguments because the assumption is that they are still on the
          // same line
          length += m_Arguments[ arg_number].GetLimitedLength( LENGTH - length, true, 0);
        }
      }

      return length;
    }

    //! @brief add double quotes around a string, and escape any internal values
    //! @param STR the string of interest
    //! @return the string, quoted only if necessary
    std::string ObjectDataLabel::QuoteStringIfDelimitersPresent( const std::string &STR)
    {
      // check for all common delimiters, include quotes
      if( STR.find( '"') != std::string::npos)
      {
        // must quote the string and escape all explicit quote characters
        std::string new_str( 1, '"');
        for( size_t i( 0), sz( STR.size()); i < sz; ++i)
        {
          if( STR[ i] == '"')
          {
            new_str += '\\';
          }
          else if( STR[ i] == '\\')
          {
            new_str += '\\';
          }
          new_str += STR[ i];
        }
        new_str += '"';
        return new_str;
      }
      else if( STR.find_first_of( GetAllDelimiters()) != std::string::npos)
      {
        // must quote the string
        std::string new_str( 1, '"');
        new_str += STR;
        new_str += '"';
        return new_str;
      }
      else if( STR.empty())
      {
        return std::string( 2, '"');
      }
      // nothing to escape
      return STR;
    }

    //! @brief get the total size of the tree
    //! @return number of nodes in the tree
    size_t ObjectDataLabel::GetTreeSize() const
    {
      if( IsScalar())
      {
        return IsEmpty() ? 0 : 1;
      }
      size_t tree_size( 1);
      // add the tree sizes of all the arguments
      for
      (
        size_t arg_number( 0), number_args( m_Arguments.size());
        arg_number < number_args;
        ++arg_number
      )
      {
        tree_size += m_Arguments[ arg_number].GetTreeSize();
      }
      return tree_size;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ObjectDataLabel::Read( std::istream &ISTREAM)
    {
      *this = ObjectDataLabel( ISTREAM);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return output stream which was written to
    std::ostream &ObjectDataLabel::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      OSTREAM << ToNamedStringDefaultWidth( io::Serialize::s_Number_spaces_per_indent * INDENT);
      return OSTREAM;
    }

    //! @return the number of ( minus the number of ) in STR, ignoring quoted sections
    long ObjectDataLabel::GetParenthesisDepthChange( const std::string &STR)
    {
      long number_unclosed_parenthesis( 0); // number of unclosed parenthesis at the current position

      for
      (
        size_t i( 0), size( STR.size()); // position in the string
        i < size;
        ++i
      )
      {
        if( STR[ i] == '(')
        {
          ++number_unclosed_parenthesis;
        }
        else if( STR[ i] == ')')
        {
          --number_unclosed_parenthesis;
        }
        else if( STR[ i] == '"') // if we find a quote, then jump to the end quote
        {
          i = ObjectDataLabelTokenizer::FindEndQuote( STR, i);
          // handle unterminated strings
          if( i == std::string::npos)
          {
            break;
          }
        }
      }

      return number_unclosed_parenthesis;
    }

    //! @brief read a string using the tokenizer
    //! @param TOKENIZER tokenizer created by the parent
    //! @param ERR stream that errors should be output to
    //! @return true on success
    bool ObjectDataLabel::ReadSubLabel( ObjectDataLabelTokenizer &TOKENIZER, std::ostream &ERR)
    {
      // reset
      m_Arguments.clear();   // remove any existing arguments
      m_Name.erase();        // and any existing value
      m_Value.erase();       // and the existing tag

      // read the first scalar into the value
      if( TOKENIZER.GetNextTokenType() == ObjectDataLabelTokenizer::e_Scalar)
      {
        m_Value = TOKENIZER.Pop();
      }

      // the previously read value was a tag
      if( TOKENIZER.GetNextTokenType() == ObjectDataLabelTokenizer::e_TagDelimiter)
      {
        // pop off the delimiter
        TOKENIZER.Pop();

        // have a name value pair, any scalar previously read in is thus the name, any scalar that follows is the value
        std::swap( m_Name, m_Value);

        if( TOKENIZER.GetNextTokenType() == ObjectDataLabelTokenizer::e_Scalar)
        {
          m_Value = TOKENIZER.Pop();
          if( TOKENIZER.GetNextTokenType() == ObjectDataLabelTokenizer::e_TagDelimiter)
          {
            ERR << "Multiple parameter assignment is forbidden but was found in given string:\n"
                << TOKENIZER.GetString() << std::endl;
            return false;
          }
        }
      }

      if( TOKENIZER.GetNextTokenType() == ObjectDataLabelTokenizer::e_Scalar)
      {
        ERR << "Variable next to quoted string is forbidden, but was found in given string:\n"
            << TOKENIZER.GetString() << std::endl;
        return false;
      }

      // check whether there are arguments following the scalar
      if( TOKENIZER.GetNextTokenType() == ObjectDataLabelTokenizer::e_ScopeOpen)
      {
        // yep, so pop off the open scope
        TOKENIZER.Pop();

        std::list< ObjectDataLabel> arg_list;

        bool had_name( false);
        size_t number_args( 0); // track # of arguments to avoid O(N) call to std::list::size()

        // read in arguments until the end is reached or a scope close is found
        while
        (
          TOKENIZER.GetLastTokenType() != ObjectDataLabelTokenizer::e_ScopeClose
          && TOKENIZER.GetNextTokenType() != ObjectDataLabelTokenizer::e_End
        )
        {
          arg_list.push_back( ObjectDataLabel());
          if( !arg_list.back().ReadSubLabel( TOKENIZER, ERR))
          {
            return false;
          }
          had_name = had_name || arg_list.back().GetName().size();
          ++number_args;

          // pop the next token, which is either ) or ,
          TOKENIZER.Pop();
        }

        if( TOKENIZER.GetLastTokenType() == ObjectDataLabelTokenizer::e_ScopeClose)
        {
          if
          (
            TOKENIZER.GetScopeDepth()
            && TOKENIZER.GetNextTokenType() != ObjectDataLabelTokenizer::e_ArgDelimiter
            && TOKENIZER.GetNextTokenType() != ObjectDataLabelTokenizer::e_ScopeClose
            && TOKENIZER.GetNextTokenType() != ObjectDataLabelTokenizer::e_End
          )
          {
            ERR << "Missing comma near " << TOKENIZER.Pop() << " in given string:\n"
                << TOKENIZER.GetString() << std::endl;
            return false;
          }
        }
        else if( TOKENIZER.GetNextTokenType() == ObjectDataLabelTokenizer::e_End)
        {
          ERR << TOKENIZER.GetScopeDepth() << " unclosed parenthesis in given string:\n"
              << TOKENIZER.GetString() << std::endl;
          return false;
        }

        m_Arguments.resize( number_args);
        std::copy( arg_list.begin(), arg_list.end(), m_Arguments.begin());
        SortArgumentsByName();
      }
      return true;
    }

    //! @brief sort all arguments by name
    void ObjectDataLabel::SortArgumentsByName()
    {
      // short circuit if there is clearly no need to sort the arguments
      if( m_Arguments.size() <= size_t( 1))
      {
        return;
      }

      // keep track of whether any arguments have a name parameter
      bool have_name( false);

      for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        if( !itr->GetName().empty())
        {
          have_name = true;
        }
      }

      if( !have_name)
      {
        // sorting unnecessary if none of the arguments had a name
        return;
      }

      const size_t number_args( m_Arguments.size());
      m_SortedArguments.resize( 0);

      // test whether the labels are already sorted by name to avoid costly calls to Sort
      bool was_sorted( true);
      for
      (
        std::vector< ObjectDataLabel>::iterator itr( m_Arguments.begin()), itr_next( itr + 1), itr_end( m_Arguments.end());
        itr_next != itr_end;
        itr = itr_next, ++itr_next
      )
      {
        if( NameIsLessThan( *itr_next, *itr))
        {
          was_sorted = false;
          break;
        }
      }

      if( was_sorted)
      {
        return;
      }

      m_SortedArguments.reserve( number_args);

      // get the beginning of this object's arguments
      ObjectDataLabel *this_arguments_begin( &m_Arguments[ 0]);
      for( size_t i( 0); i < number_args; ++i)
      {
        m_SortedArguments.push_back( this_arguments_begin + i);
      }

      // sort argument pointers
      std::stable_sort( m_SortedArguments.begin(), m_SortedArguments.end(), &PtrNameIsLessThan);
    }

    //! @brief copy a sorted arguments vector from another object data label
    //! @param LABEL the label to copy it from
    void ObjectDataLabel::CopyArgumentSorting( const ObjectDataLabel &LABEL)
    {
      m_SortedArguments.clear();

      if( IsScalar() || LABEL.m_SortedArguments.empty())
      {
        return;
      }
      m_SortedArguments.reserve( LABEL.m_SortedArguments.size());

      // get the beginning of this object's arguments
      ObjectDataLabel *this_arguments_begin( &m_Arguments[ 0]);
      const ObjectDataLabel *that_arguments_begin( &LABEL.m_Arguments[ 0]);
      for
      (
        std::vector< ObjectDataLabel *>::const_iterator
          itr( LABEL.m_SortedArguments.begin()), itr_end( LABEL.m_SortedArguments.end());
        itr != itr_end;
        ++itr
      )
      {
        // add a pointer to the corresponding label for this object
        m_SortedArguments.push_back( this_arguments_begin + std::ptrdiff_t( *itr - that_arguments_begin));
      }
    }

    //! operator == (Comparison)
    bool operator ==( const ObjectDataLabel &LABEL_A, const ObjectDataLabel &LABEL_B)
    {
      return LABEL_A.TernaryCompare( LABEL_B) == 0;
    }

    //! operator != (Comparison)
    bool operator !=( const ObjectDataLabel &LABEL_A, const ObjectDataLabel &LABEL_B)
    {
      return LABEL_A.TernaryCompare( LABEL_B) != 0;
    }

    //! operator < (Comparison)
    bool operator <( const ObjectDataLabel &LABEL_A, const ObjectDataLabel &LABEL_B)
    {
      return LABEL_A.TernaryCompare( LABEL_B) < 0;
    }

    //! operator > (Comparison)
    bool operator >( const ObjectDataLabel &LABEL_A, const ObjectDataLabel &LABEL_B)
    {
      return LABEL_A.TernaryCompare( LABEL_B) > 0;
    }

  } // namespace util
} // namespace bcl
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
#include "util/bcl_util_object_data_label_tokenizer.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "io/bcl_io_ifstream.h"
#include "util/bcl_util_object_data_label.h"
#include "util/bcl_util_object_instances.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const SiPtr< const ObjectInterface> ObjectDataLabelTokenizer::s_Instance
    (
      GetObjectInstances().AddInstance( new ObjectDataLabelTokenizer( ""))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from string
    ObjectDataLabelTokenizer::ObjectDataLabelTokenizer( const std::string &STR, const size_t &POSITION) :
      m_String( STR),
      m_Size( STR.size()),
      m_Position( POSITION),
      m_ScopeDepth( 0),
      m_NextType( e_Start)
    {
      DetermineNextTokenType();
    }

    //! @brief Clone function
    //! @return pointer to new ObjectDataLabelTokenizer
    ObjectDataLabelTokenizer *ObjectDataLabelTokenizer::Clone() const
    {
      return new ObjectDataLabelTokenizer( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ObjectDataLabelTokenizer::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get token type
    //! @return the current token type
    const ObjectDataLabelTokenizer::Type &ObjectDataLabelTokenizer::GetLastTokenType() const
    {
      return m_LastType;
    }

    //! @brief get type of token that will next be returned by Pop
    //! @return the type of token that will next be returned by Pop
    const ObjectDataLabelTokenizer::Type &ObjectDataLabelTokenizer::GetNextTokenType() const
    {
      return m_NextType;
    }

    //! @brief get token type
    //! @return the current token type
    const size_t &ObjectDataLabelTokenizer::GetScopeDepth() const
    {
      return m_ScopeDepth;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get the string for the next token
    std::string ObjectDataLabelTokenizer::Pop()
    {
      m_LastType = m_NextType;
      if( m_LastType == e_End)
      {
        return std::string();
      }

      // handle scalars
      if( m_LastType == e_Scalar)
      {
        size_t start_position( m_Position);
        const char first_char( m_String[ m_Position]);

        // handle quoted strings
        if( first_char == '"')
        {
          ++start_position;
          m_Position = ObjectDataLabelTokenizer::FindEndQuote( m_String, m_Position);
          if( m_Position == std::string::npos)
          {
            // no end quote, return complete string
            return m_String.substr( start_position);
          }

          const size_t end_position( m_Position++);

          // determine next type
          DetermineNextTokenType();

          return m_String.substr( start_position, end_position - start_position);
        }

        // handle normal scalars, keep track of the last position that was not a space

        // find the next delimiter
        m_Position = m_String.find_first_of( ObjectDataLabel::GetAllDelimiters(), m_Position);

        // if there are no remaining delimiters, reset to the end of the string
        if( m_Position == std::string::npos)
        {
          m_Position = m_Size;
        }

        // determine the next token type
        DetermineNextTokenType();

        // find the last non-space character, since the scalar does not include leading/trailing spaces
        size_t last_non_space_pos( m_Position - 1);
        while( isspace( m_String[ last_non_space_pos]))
        {
          --last_non_space_pos;
        }

        // return the substring containing only the scalar
        return m_String.substr( start_position, last_non_space_pos + 1 - start_position);
      }

      // handle delimiters
      // get the current character and advance the position
      const char character( m_String[ m_Position++]);

      // determine the next token type
      DetermineNextTokenType();

      // handle scope depth changes
      if( m_LastType == e_ScopeOpen)
      {
        ++m_ScopeDepth;
      }
      else if( m_LastType == e_ScopeClose)
      {
        --m_ScopeDepth;
      }

      // return the string with the character
      return std::string( 1, character);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief find the end of the quoted string that starts at a given position, ignoring escaped quotes
    //! @param STR string of interest
    //! @param START_QUOTE_POS the position of the starting quote
    //! @return the end quote position, or std::string::npos if no such quote exists
    size_t ObjectDataLabelTokenizer::FindEndQuote( const std::string &STR, const size_t &START_QUOTE_POS)
    {
      const size_t last_char_pos( STR.size() - 1);
      const char quote( STR[ START_QUOTE_POS]);

      // check for quotes at the end of the string
      if( START_QUOTE_POS == last_char_pos)
      {
        return std::string::npos;
      }

      // find the next quote character
      size_t pos( STR.find( quote, START_QUOTE_POS + 1));

      // continue while the quote character is preceded by an escape character
      while( pos < last_char_pos && STR[ pos - 1] == '\\')
      {
        pos = STR.find( quote, pos + 1);
      }

      // if pos is exactly at the last character, make sure that it was not preceded by an escape
      if( pos == last_char_pos && STR[ pos - 1] == '\\')
      {
        return std::string::npos;
      }

      // return position
      return pos;
    }

    //! @brief  validate a string, write out errors to a given stream
    //! @param  STR the string to validate
    //! @param  ERR the stream to write errors to
    //! @return true if the string was valid
    bool ObjectDataLabelTokenizer::Validate( const std::string &STR, std::ostream &ERR)
    {
      ObjectDataLabelTokenizer validator( STR);
      Type previous_type( e_Start);
      while( validator.GetNextTokenType() != e_End)
      {
        // check for invalid combinations of last/next tokens
        if( validator.GetScopeDepth() == 0)
        {
          if( validator.GetNextTokenType() == e_ArgDelimiter)
          {
            ERR << "Arguments are not allowed without an associated object, that is, ',' cannot appear outside\n"
                << "parenthetical scope for object data labels, but did in given string:\n"
                << STR << std::endl;
            return false;
          }
          else if( validator.GetNextTokenType() == e_ScopeClose)
          {
            ERR << "Excessive end-scope parenthesis in given string:\n"
                << STR << std::endl;
            return false;
          }
        }
        if( validator.GetLastTokenType() == e_Scalar && validator.GetNextTokenType() == e_Scalar)
        {
          ERR << "Variable next to quoted string is forbidden, but was found in given string:\n" << STR << std::endl;
          return false;
        }
        if( validator.GetLastTokenType() == e_ScopeClose)
        {
          if
          (
            validator.GetNextTokenType() != e_ArgDelimiter
            && validator.GetNextTokenType() != e_ScopeClose
            && validator.GetNextTokenType() != e_End
          )
          {
            ERR << "Missing comma near " << validator.Pop() << " in given string:\n"
                << STR << std::endl;
            return false;
          }
        }

        // check that the string =X= never appears
        // people with programming might try to use this to set multiple variables equal to some value,
        // but this currently is not supported
        if
        (
          previous_type == e_TagDelimiter
          && validator.GetLastTokenType() == e_Scalar
          && validator.GetNextTokenType() == e_TagDelimiter
        )
        {
          ERR << "Multiple parameter assignment is forbidden but was found in given string:\n"
              << STR << std::endl;
          return false;
        }

        previous_type = validator.GetLastTokenType();
        validator.Pop();
      }
      if( validator.GetScopeDepth())
      {
        ERR << validator.GetScopeDepth() << " unclosed parenthesis in given string:\n"
            << STR << std::endl;
        return false;
      }
      return true;
    }

    //! @brief update the next token type, after skipping spaces
    void ObjectDataLabelTokenizer::DetermineNextTokenType()
    {
      while( m_Position < m_Size && isspace( m_String[ m_Position]))
      {
        ++m_Position;
      }

      m_LastType = m_NextType;
      if( m_Position >= m_Size)
      {
        m_NextType = e_End;
      }
      else
      {
        // switch to find delimiter-based enums
        bool is_file( false);
        switch( m_String[ m_Position])
        {
          case '(':
            m_NextType = e_ScopeOpen;
            break;
          case ')':
            m_NextType = e_ScopeClose;
            break;
          case ',':
            m_NextType = e_ArgDelimiter;
            break;
          case '=':
            m_NextType = e_TagDelimiter;
            break;
          case '@':
            is_file = true;
            break;
          default:
            // all other letters indicate scalars
            m_NextType = e_Scalar;
            break;
        };
        if( is_file)
        {
          const size_t end_filename_pos( m_String.find_first_of( ObjectDataLabel::GetAllDelimiters(), m_Position + 1));
          const std::string filename( RStrip( m_String.substr( m_Position + 1, end_filename_pos - m_Position - 1), " \t\n\r"));

          io::IFStream ist;
          io::File::MustOpenIFStream( ist, filename);
          const std::string complete_file
          (
            Join( "\n", util::StringListFromIStream( ist))
          );
          m_String = m_String.substr( 0, m_Position) + complete_file + m_String.substr( end_filename_pos);
          m_Size = m_String.size();
          DetermineNextTokenType();
        }
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ObjectDataLabelTokenizer::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_String, ISTREAM);
      m_Size = m_String.size();

      size_t new_position( 0);
      io::Serialize::Read( new_position, ISTREAM);

      // get back to the previous state
      m_Position = 0;
      m_LastType = m_NextType = e_Start;
      while( m_Position < new_position)
      {
        Pop();
      }

      DetermineNextTokenType();
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ObjectDataLabelTokenizer::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_String, OSTREAM, INDENT) << '\t' << m_Position;
      return OSTREAM;
    }

  } // namespace util
} // namespace bcl
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
#include "util/bcl_util_object_instances.h"

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

  //////////
  // data //
  //////////

    //! @brief Singleton set of known objects, used to avoid costly list searches
    storage::Map< std::string, SiPtr< const ObjectInterface> > &GetKnownObjects()
    {
      static storage::Map< std::string, SiPtr< const ObjectInterface> > s_known_objects;

      return s_known_objects;
    }

    //! @brief vector of names of all known objects
    //! @return a vector of names of all known bcl objects
    storage::Vector< std::string> &GetKnownObjectsNamesNonConst()
    {
      static storage::Vector< std::string> s_objects;
      return s_objects;
    }

    //! @brief vector of names of all known objects
    //! @return a vector of names of all known bcl objects
    const storage::Vector< std::string> &ObjectInstances::GetKnownObjectNames() const
    {
      return GetKnownObjectsNamesNonConst();
    }

    //! @brief returns a pointer to object for the DESCRIPTOR / a NULL ptr if there is no such DESCRIPTOR in the map
    //! @param DESCRIPTOR the bcl descriptor as returned by GetStaticClassname
    //! @return ptr to ObjectInterface
    const SiPtr< const ObjectInterface> &ObjectInstances::GetPtrToObjectFromIdentifier
    (
      const std::string &DESCRIPTOR
    ) const
    {
      // find the corresponding enum
      storage::Map< std::string, SiPtr< const ObjectInterface> >::const_iterator
        itr( GetKnownObjects().Find( DESCRIPTOR));

      //check that entry with DESCRIPTOR was found
      BCL_Assert
      (
        itr != GetKnownObjects().End(),
        "entry with descriptor \"" + DESCRIPTOR + "\" was not found"
      );

      //if for some reason the map contains an empty pointer
      if( !itr->second.IsDefined())
      {
        static const SiPtr< const ObjectInterface> s_empty;
        BCL_MessageCrt( "Object instance for " + DESCRIPTOR + " was NULL!");
        return s_empty;
      }

      return itr->second;
    }

    //! @brief get the address used for writing in files when a shared pointer is not intended for combining with
    //!        other shared pointers upon reading
    const size_t &ObjectInstances::GetUncombinedPointerAddress()
    {
      static const size_t s_uncombined( 12345678);
      return s_uncombined;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief returns a pointer to object for the DESCRIPTOR / an empty ptr if there is no such DESCRIPTOR in the map
    //! @param DESCRIPTOR the bcl descriptor as returned by GetStaticClassname
    //! @return ptr to ObjectInterface (must be a base of the actual object, since a dynamic cast is used) to a copy of
    //!         the class instance on every bcl object
    ObjectInterface *ObjectInstances::GetNewPtrToObjectFromIdentifier( const std::string &DESCRIPTOR) const
    {
      return GetPtrToObjectFromIdentifier( DESCRIPTOR)->Clone();
    }

    //! @brief write out an object behind a shared pointer
    //! @param PTR the shared pointer to write out
    //! @param OSTREAM the stream to write the pointer out to
    //! @param INDENT the amount of indent to use
    void ObjectInstances::WriteShPtr
    (
      const ObjectInterface *PTR,
      std::ostream &OSTREAM,
      const size_t INDENT
    ) const
    {
      // write the address; completely handles null pointers
      WritePtrAddress( size_t( PTR), OSTREAM, INDENT);
      if( !PTR)
      {
        return;
      }

      // check that the descriptor map contains the classname to guarantee a safe reading
      if( !GetKnownObjects().Has( PTR->GetClassIdentifier()))
      {
        BCL_MessageCrt
        (
          "entry with name: \"" + PTR->GetClassIdentifier() +
          "\" does not exist => reading from ShPtr will be impossible"
        );
      }

      //end
      PTR->WriteObject( OSTREAM, INDENT);
    }

    //! @brief write out a pointer address
    //! @param ADDRESS addressed memory location
    //! @param STREAM the stream to write the pointer out to
    //! @param INDENT the amount of indent to use
    void ObjectInstances::WritePtrAddress( const size_t &ADDRESS, std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write address of pointer
      io::Serialize::Write( ADDRESS ? GetUncombinedPointerAddress() : 0, OSTREAM, INDENT) << '\n';
    }

  /////////////////
  // data access //
  /////////////////

  ////////////////
  // operations //
  ////////////////

    //! @brief add pointer of class instance to enumeration of ObjectInstances
    //! @param INTERFACE pointer to the instance passed derived from ObjectInterface
    //! @return INTERFACE
    SiPtr< const ObjectInterface> ObjectInstances::AddInstance( const ObjectInterface *INTERFACE)
    {
      if( INTERFACE == NULL)
      {
        return SiPtr< const ObjectInterface>();
      }
      return AddInstanceWithName( INTERFACE, INTERFACE->GetClassIdentifier());
    }

    //! @brief add pointer of class instance to enumeration of ObjectInstances
    //! @param INTERFACE pointer to the instance passed derived from ObjectInterface
    //! @return SiPtr to the interface
    SiPtr< const ObjectInterface> ObjectInstances::AddInstanceWithName
    (
      const ObjectInterface *INTERFACE,
      const std::string &NAME
    )
    {
      if( GetKnownObjects().Insert( std::make_pair( NAME, SiPtr< const ObjectInterface>( INTERFACE))).second)
      {
        GetKnownObjectsNamesNonConst().PushBack( NAME);
      }
      else
      {
        BCL_MessageCrt( NAME + " was tried to insert twice into object instances, skipped the second time!");
      }

      return SiPtr< const ObjectInterface>( INTERFACE);
    }

    //! @brief add pointer of class instance to enumeration of ObjectInstances, if it is not already there
    //! @param INTERFACE pointer to the instance passed derived from ObjectInterface
    //! @return true if the instance was added
    bool ObjectInstances::TryAddInstance( ObjectInterface *INTERFACE)
    {
      if( INTERFACE == NULL)
      {
        return false;
      }

      if( !GetKnownObjects().Has( INTERFACE->GetClassIdentifier()))
      {
        // add the object to object instances, if it is not already there
        AddInstance( INTERFACE);
        return true;
      }
      return false;
    }

  } // namespace util

  //! @brief get the only instance of ObjectInstances
  util::ObjectInstances &GetObjectInstances()
  {
    static util::ObjectInstances s_instances;
    return s_instances;
  }

} // namespace bcl
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
#include "util/bcl_util_object_interface.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_assert.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

    //! @brief virtual destructor
    ObjectInterface::~ObjectInterface()
    {
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read the identifier
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ObjectInterface::ReadIdentifier( std::istream &ISTREAM) const
    {
      //extract the class name form the stream
      std::string identifier( ExtractIdentifier( ISTREAM));

       // assert that this is a stream for this object
      BCL_Assert
      (
        GetClassIdentifier() == identifier,
        "wrong stream; expected \"" + GetClassIdentifier() + "\" but encountered \"" + identifier + "\""
      );

      //end
      return ISTREAM;
    }

    //! @brief write the identifier
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return output stream which was written to
    std::ostream &ObjectInterface::WriteIdentifier( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write class name as identifier of the string
      OSTREAM << std::string( INDENT * 2, ' ') << GetClassIdentifier() << '\n';

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ObjectInterface::ReadObject( std::istream &ISTREAM)
    {
      //read and check stream identifier
      ReadIdentifier( ISTREAM);

      //call the Read function of the class
      Read( ISTREAM);

      //end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return output stream which was written to
    std::ostream &ObjectInterface::WriteObject( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write identifier
      WriteIdentifier( OSTREAM, INDENT);

      // call the write function of the class
      Write( OSTREAM, INDENT + 1);

      // end
      return OSTREAM;
    }

    // extracts the identifier string for the given ISTREAM
    std::string ObjectInterface::ExtractIdentifier( std::istream &ISTREAM)
    {
      //instantiate identifier string that will be extracted
      std::string identifier;

      // repeat getline until eof is reached or the extracted stream is not empty
      while( !ISTREAM.eof() && ISTREAM.good())
      {
        ISTREAM >> std::ws;
        //get identifier string
        ISTREAM >> identifier;

        if( !identifier.empty())
        {
          ISTREAM >> std::ws;
          break;
        }
      }

      // return identifier
      return identifier;
    }

    //! @brief extract namespace name without the "bcl::" part from a given object identifier
    //! @param IDENTIFIER object identifier such as "bcl::util::ClassName"
    //! @return namespace name without the "bcl::" part such as "util" or "assemble"
    const std::string ObjectInterface::ExtractNamespaceName( const std::string &IDENTIFIER)
    {
      static const std::string s_name( bcl::GetNamespaceIdentifier() + "::");
      BCL_Assert
      (
        IDENTIFIER.size() >= s_name.size() + 1
        && std::equal( IDENTIFIER.begin(), IDENTIFIER.begin() + s_name.size(), s_name.begin()),
        "The provided identifier does not have at least one namespace name after BCL : " + IDENTIFIER
      );

      const size_t namespace_end_pos( IDENTIFIER.find( "::", s_name.size() + 1));

      // return the second one after "bcl"
      return IDENTIFIER.substr( s_name.size(), namespace_end_pos - s_name.size());
    }

    //! @brief operator >> to read ReadWrite derived Object from std::istream
    //! @param ISTREAM inpustream the object is read from
    //! @param OBJECT the object derived from ReadWrite
    //! @return the inputstream where the object was read from
    std::istream &operator >>
    (
      std::istream &ISTREAM,
      ObjectInterface &OBJECT
    )
    {
      return OBJECT.ReadObject( ISTREAM);
    }

    //! @brief operator >> to read const ReadWrite derived Object from std::istream which will assert, since this is not possible
    //! @param ISTREAM inpustream the object is read from
    //! @param OBJECT the object derived from ReadWrite
    //! @return the inputstream where the object was read from
    std::istream &operator >>
    (
      std::istream &ISTREAM,
      const ObjectInterface &OBJECT
    )
    {
      BCL_Exit( "do not call this function with a const OBJECT", -1);
      return ISTREAM;
    }

    //! @brief write ReadWrite derived Object to std::ostream
    //! @param OSTREAM outputstream the object is written to
    //! @param OBJECT object derived from ReadWrite
    //! @return outputstream the object was written to
    std::ostream &operator <<
    (
      std::ostream &OSTREAM,
      const ObjectInterface &OBJECT
    )
    {
      return OBJECT.WriteObject( OSTREAM, 0);
    }

  } // namespace util
} // namespace bcl

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
#include "util/bcl_util_ptr_interface.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_call_stack.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    //! @brief Writes out a warning message that a pointer cast failed between the two types
    //! @param PTR_TYPE the pointer type that was received
    //! @param CAST_TYPE the type the pointer was cast to
    //! @note this function provides a convenient hook for debuggers, which need only set a single breakpoint in the cpp
    void NotifyUserBadPointerCast( const std::string &PTR_TYPE, const std::string &CAST_TYPE)
    {
      BCL_MessageCrt
      (
        "was not able to cast pointer from " + PTR_TYPE + " to " + CAST_TYPE
        + " Callstack: " + CallStack().String()
      );
    }

    //! @brief get the string for an empty pointer, e.g. "NULL"
    //! @return the string for an empty pointer, e.g. "NULL"
    const std::string &GetNullDescriptor()
    {
      static const std::string s_null_descriptor( "NULL");
      return s_null_descriptor;
    }

  } // namespace util
} // namespace bcl
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
#include "util/bcl_util_runtime_environment_default.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RuntimeEnvironmentDefault::RuntimeEnvironmentDefault()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @return the class name if this function is overwritten
    const std::string &RuntimeEnvironmentDefault::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief access to the process number or slot number
    //! @return the process or slot, this app is running in
    int RuntimeEnvironmentDefault::GetProcessNumber() const
    {
      return GetUndefined< int>();
    }

  ////////////////
  // operations //
  ////////////////

    //! @copydoc RuntimeEnvironmentInterface::ResolveFileName()
    std::string RuntimeEnvironmentDefault::ResolveFileName( const std::string &FILE_NAME) const
    {
      return FILE_NAME;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RuntimeEnvironmentDefault::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &RuntimeEnvironmentDefault::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Initializes the environment
    //! @return A boolean representing whether the initialization was a success or failure.
    bool RuntimeEnvironmentDefault::Initialize() const
    {
      return ( true);
    }

    //! @brief Finalizes the environment
    //! @return A boolean representing whether finalize was a success or failure.
    //! @param STATUS indicates if an error occured
    bool RuntimeEnvironmentDefault::Finalize( const int STATUS) const
    {
      return STATUS == 0;
    }

  } // namespace util
} // namespace bcl
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
#include "util/bcl_util_runtime_environment_interface.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_runtime_environments.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

  ////////////////
  // operations //
  ////////////////

    //! @brief Returns the operating system specific path separator
    //! @return the operating-system specific file path separator
    const std::string &RuntimeEnvironmentInterface::GetPathSeperator() const
    {
      static const std::string s_path_seperator( "/");

      // return
      return s_path_seperator;
    }

    //! @brief report progress
    //! @param FRACTION_PROGRESS progress between [0,1]
    //! @return true if successful reported
    bool RuntimeEnvironmentInterface::ReportProgress( const double FRACTION_PROGRESS) const
    {
      BCL_MessageStd( Format().W( 3).FFP( 0)( 100 * FRACTION_PROGRESS) + "% finished");
      return true;
    }

    //! @brief add name to the process number
    //! @param NAME name to be added to process number
    //! @return NAME + processnumber as string
    std::string RuntimeEnvironmentInterface::AddProcessNumberToName( const std::string &NAME) const
    {
      return std::string( NAME + Format()( GetProcessNumber()));
    }

    //! @brief returns active RuntimeEnvironment
    //! @return active RuntimeEnvironment
    const RuntimeEnvironmentInterface &GetRuntimeEnvironment()
    {
      return GetRuntimeEnvironments().GetCurrentEnvironment();
    }

  } // namespace util
} // namespace bcl
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
#include "util/bcl_util_runtime_environments.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "util/bcl_util_enumerate.hpp"
#include "util/bcl_util_runtime_environment_default.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

  //////////
  // data //
  //////////

    //! command line flag to be used to set RuntimeEnvironment over the command line
    ShPtr< command::FlagInterface> &RuntimeEnvironments::GetFlagRuntimeEnvironment()
    {
      static ShPtr< command::FlagInterface> s_runtime_environment_flag
      (
        new command::FlagStatic
        (
          "runtime_environment",
          "change the runtime environment this executable runs in",
          command::Parameter
          (
            "environment",
            "choice of environment",
            command::ParameterCheckEnumerate< RuntimeEnvironments>(),
            "Default"
          ),
          &RuntimeEnvironments::UpdateCurrentRuntimeEnvironmentFromCommandLineFlag
        )
      );

      return s_runtime_environment_flag;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RuntimeEnvironments::RuntimeEnvironments() :
      Enumerate< ShPtr< RuntimeEnvironmentInterface>, RuntimeEnvironments>( false),
      e_Default( AddEnum( "Default", ShPtr< RuntimeEnvironmentInterface>( new RuntimeEnvironmentDefault()))),
      m_CurrentRuntimeEnvironment( *e_Default)
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RuntimeEnvironments::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief Initialize the runtime environment from the command line flag
    void RuntimeEnvironments::UpdateCurrentRuntimeEnvironmentFromCommandLineFlag()
    {
      // name of runtime environment from commandline
      const std::string &re_name( GetFlagRuntimeEnvironment()->GetFirstParameter()->GetValue());

      // construct the environment enum from the name given in the commandline
      ShPtr< RuntimeEnvironmentInterface> new_runtime_environment( *RuntimeEnvironment( re_name));

      // the current environment is correct
      if( new_runtime_environment == GetEnums().m_CurrentRuntimeEnvironment)
      {
        return;
      }

      // update the current runtime environment to the commandline selected one
      GetEnums().m_CurrentRuntimeEnvironment = new_runtime_environment;

      // message to user
      BCL_MessageTop( "updating runtime environment to " + re_name);

      // initialize environment
      BCL_Assert
      (
        GetEnums().m_CurrentRuntimeEnvironment->Initialize(),
        "unable to initialize the runtime environment " + re_name
      );
    }

    //! @brief Initialize and return the current environment
    //! @return one and only reference to one of the environments
    const RuntimeEnvironmentInterface &RuntimeEnvironments::GetCurrentEnvironment() const
    {
      // return the current runtime environment
      return *m_CurrentRuntimeEnvironment;
    }

    //! @brief returns RuntimeEnvironments
    //! @return RuntimeEnvironments
    RuntimeEnvironments &GetRuntimeEnvironments()
    {
      return RuntimeEnvironments::GetEnums();
    }

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< RuntimeEnvironmentInterface>, RuntimeEnvironments>;

  } // namespace util
} // namespace bcl
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
#include "util/bcl_util_serializable_interface.h"

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief virtual destructor
    SerializableInterface::~SerializableInterface()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief return the label for this object
    //! @param INCLUDE_DATA whether to include members defined as data, otherwise, include only initialization info
    //! @return a label for this object
    ObjectDataLabel SerializableInterface::GetLabel( const bool &INCLUDE_DATA) const
    {
      return GetCompleteSerializer().GetLabel( INCLUDE_DATA);
    }

    //! @brief return the string for this object
    //! @param INCLUDE_DATA whether to include members defined as data, otherwise, include only initialization info
    //! @return the initialization string for this object
    std::string SerializableInterface::GetString( const bool &INCLUDE_DATA) const
    {
      return GetLabel( INCLUDE_DATA).ToString();
    }

    //! @brief determine the type of value that can be parsed
    DataType::Type SerializableInterface::GetSerializedType() const
    {
      return DataType::e_StaticObject;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief Try to read the members of this object from the given label
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    bool SerializableInterface::TryRead( const ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      return ValidateRead( LABEL, ERROR_STREAM).IsAllowed();
    }

    //! @brief Read the members of this object from the given label, assert on failure, error messages sent to logger
    //! @param LABEL the label containing members that should be read of this class
    void SerializableInterface::AssertRead( const ObjectDataLabel &LABEL)
    {
      BCL_MessageDbg( "Trying to read " + LABEL.ToString() + " into " + this->GetClassIdentifier());
      std::stringstream err_stream;
      io::ValidationResult result( ValidateRead( LABEL, err_stream));
      if( !result.IsAllowed())
      {
        if( result.IsHelp())
        {
          BCL_ExitWithoutCallstack( "Requested help: " + err_stream.str(), 0);
        }
        else if( !result.IsComplete())
        {
          BCL_Exit( "Failed to read " + this->GetClassIdentifier() + ", error was: " + err_stream.str() + ", exiting", -1);
        }
      }
    }

    //! @brief set the value of the corresponding member based on the label
    //! @param LABEL label that is used to set the string
    //! @param ERROR_STREAM stream to write errors to
    //! @return result of validation
    io::ValidationResult SerializableInterface::ValidateRead
    (
      const ObjectDataLabel &LABEL,
      std::ostream &ERROR_STREAM
    )
    {
      if( LABEL.GetName() == io::ValidationResult::GetHelpString())
      {
        WriteHelp( ERROR_STREAM);
        return io::ValidationResult( io::e_Help);
      }

      // call PreReadHook: overridable function to reset any necessary members
      io::ValidationResult result( PreReadHook( LABEL, ERROR_STREAM));
      if( result.IsInvalid())
      {
        ERROR_STREAM << "Failed at pre-read hook";
        return result;
      }
      else if( result.IsComplete())
      {
        return io::ValidationResult( io::e_Allowed);
      }

      // get the serializer for this object
      io::Serializer serializer( GetCompleteSerializer());

      // read in arguments
      result = serializer.ReadArguments( LABEL, ERROR_STREAM);
      if( !result)
      {
        ERROR_STREAM << "Failed at reading arguments";
        // failed to read object correctly
        return result;
      }
      else if( serializer.GetWasDataRead())
      {
        // call function to finalize object if data members were read
        result = io::ValidationResult( ReadDataSuccessHook( LABEL, ERROR_STREAM));
        if( result.IsInvalid())
        {
          ERROR_STREAM << "Failed at ReadDataSuccessHook for " << this->GetClassIdentifier();
        }
        return result;
      }

      // only read initializers
      result = io::ValidationResult( ReadInitializerSuccessHook( LABEL, ERROR_STREAM));
      if( result.IsInvalid())
      {
        ERROR_STREAM << "Failed at ReadInitializerSuccessHook for " << this->GetClassIdentifier();
      }
      return result;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SerializableInterface::Read( std::istream &ISTREAM)
    {
      AssertRead( ObjectDataLabel( ISTREAM));
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &SerializableInterface::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return GetCompleteSerializer().Write( OSTREAM, INDENT);
    }

    //! @brief write this help for this object
    //! @param STREAM stream to write to
    //! @param INDENT indentation to give each line
    std::ostream &SerializableInterface::WriteHelp( std::ostream &STREAM, const size_t INDENT) const
    {
      return GetCompleteSerializer().WriteHelp( STREAM, INDENT);
    }

    //! @brief write this help for this object
    //! @param WRITER a fixed line width writer
    io::FixedLineWidthWriter &SerializableInterface::WriteHelp( io::FixedLineWidthWriter &WRITER) const
    {
      return GetCompleteSerializer().WriteHelp( WRITER);
    }

    //! @brief get the complete serializer, including command line identifier and type
    //! @return the complete serializer, including command line identifier and type
    io::Serializer SerializableInterface::GetCompleteSerializer() const
    {
      return GetSerializer().SetCommandLineIdentifier( GetAlias()).SetTypeAndFinalize( GetSerializedType());
    }

    //! These functions are sometimes overridden to initialize data members from initialization variables setup during read

    //! @brief a function that derived classes can override to perform some action on a class before any data members
    //!        are read, e.g. resetting certain data members so that a post-read command can update parameters that
    //!        were not set
    //! @param SERIALIZER The serializer that will be used
    //! @param ERR_STREAM stream to write any errors encountered to
    //! @return result of any validation performed internally
    io::ValidationResult SerializableInterface::PreReadHook( const ObjectDataLabel &SERIALIZER, std::ostream &ERR_STREAM)
    {
      return io::ValidationResult( io::e_Allowed);
    }

    //! @brief a function that derived classes can override to set up additional data members whenever Read is called successfully
    //! @param SERIALIZER the serializer object with initialization information
    //! @param ERR_STREAM stream to write out errors to
    //! @return true, unless there were new errors
    bool SerializableInterface::ReadInitializerSuccessHook( const ObjectDataLabel &SERIALIZER, std::ostream &ERR_STREAM)
    {
      return true;
    }

    //! @brief a function that derived classes can override to take additional actions whenever Read is called successfully
    //!        AND any data members are specified for the given class
    //! @param SERIALIZER the serializer object with initialization information
    //! @param ERR_STREAM stream to write out errors to
    //! @return true, unless there were new errors
    bool SerializableInterface::ReadDataSuccessHook( const ObjectDataLabel &SERIALIZER, std::ostream &ERR_STREAM)
    {
      return true;
    }

    //! @brief Get a set of all class names used by the serializer. Useful for introspection
    //! @param TYPES set to insert type names into
    //! @param INCLUDE_OPTIONAL true to also count optional members
    //! @param INCLUDE_DATA true to also include data-containing members
    void SerializableInterface::InsertDataTypes
    (
      storage::Map< std::string, size_t> &TYPES,
      const bool &INCLUDE_OPTIONAL,
      const bool &INCLUDE_DATA,
      const size_t &MAX_DEPTH
    ) const
    {
      ++TYPES[ this->GetClassIdentifier()];
      this->GetCompleteSerializer().InsertDataTypes( TYPES, INCLUDE_OPTIONAL, INCLUDE_DATA, MAX_DEPTH);
    }

    //! @brief Get addresses of all objects serialized as part of this
    //! @notes used to ensure uniqueness of serialized objects
    std::vector< size_t> SerializableInterface::GetSerializationAddresses() const
    {
      if( this->GetSerializedType() == DataType::e_DynamicObject)
      {
        return io::SerializationInterface::GetSerializationAddresses();
      }
      return this->GetCompleteSerializer().GetSerializationAddresses();
    }

    //! @brief infer functional type of the serializable interface, relative to a specific interface
    //! @param INTERFACE_STR the type name string of the interface of interest (by default, any interface)
    //! @return the functional type of the serializable interface, relative to a specific interface
    FunctionalType::Type SerializableInterface::InferFunctionalType( const std::string &INTERFACE_STR) const
    {
      storage::Map< std::string, size_t> type_counts;
      this->InsertDataTypes( type_counts, true, false, 2);
      size_t count_match_iface( 0), count_unmatched_iface( 0);
      for
      (
        storage::Map< std::string, size_t>::const_iterator itr( type_counts.Begin()), itr_end( type_counts.End());
        itr != itr_end;
        ++itr
      )
      {
        if( itr->first != this->GetClassIdentifier())
        {
          ( StartsWith( itr->first, INTERFACE_STR) ? count_match_iface : count_unmatched_iface) += itr->second;
        }
      }
      return FunctionalType::GetType( count_match_iface, count_unmatched_iface);
    }

    //! @brief operator >> to read ReadWrite derived Object from std::istream
    //! @param ISTREAM input stream the object is read from
    //! @param OBJECT the object derived from ReadWrite
    //! @return the input stream where the object was read from
    std::istream &operator >>( std::istream &ISTREAM, SerializableInterface &OBJECT)
    {
      return OBJECT.Read( ISTREAM);
    }

    //! @brief operator >> to read const ReadWrite derived Object from std::istream which will assert, since this is not possible
    //! @param ISTREAM input stream the object is read from
    //! @param OBJECT the object derived from ReadWrite
    //! @return the input stream where the object was read from
    std::istream &operator >>( std::istream &ISTREAM, const SerializableInterface &OBJECT)
    {
      BCL_Exit( "cannot read a const OBJECT", -1);
      return ISTREAM;
    }

    //! @brief write ReadWrite derived Object to std::ostream
    //! @param OSTREAM output stream the object is written to
    //! @param OBJECT object derived from ReadWrite
    //! @return output stream the object was written to
    std::ostream &operator <<( std::ostream &OSTREAM, const SerializableInterface &OBJECT)
    {
      return OBJECT.Write( OSTREAM, 0);
    }

  } // namespace util
} // namespace bcl
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
#include "util/bcl_util_sh_ptr.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    template class BCL_API ShPtr< ObjectInterface>;
  } // namespace util
} // namespace bcl
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

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically
#include <iostream>

namespace bcl
{
  namespace util
  {

    bool WriteFiasco( const std::string &FILE_NAME)
    {
      static size_t counter( 0);
      ++counter;

      std::cerr << "Starting to initialize file - number: [" << counter
                << "] filename: [" << FILE_NAME << "]" << '\n';

      return true;
    }

  } // namespace util
} // namespace bcl
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
#include "util/bcl_util_stopwatch.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! set the Timer Start to the current time
    //! @param PRINT_ON_DESTRUCTION true if the current time of this stopwatch should be printed when destroyed
    Stopwatch::Stopwatch( const bool PRINT_ON_DESTRUCTION) :
      m_Description(),
      m_TimerStart(), // set to zero
      m_TotalTime(),  // set to zero
      m_LastNotificationTime(),
      m_NotificationInterval( std::numeric_limits< size_t>::max(), 0),
      m_MessageLevel( Message::e_Verbose),
      m_PrintOnDestruction( PRINT_ON_DESTRUCTION)
    {
      Start();
    }

    //! @brief constructor with description and notification options
    //! set the Timer Start to the current time
    //! @param DESCRIPTION string that describes what the stopwatch is timing
    //! @param NOTIFICATION_INTERVAL the minimum interval between messages
    //! @param MESSAGE_LEVEL the message level to use when printing messages
    //! @param PRINT_ON_DESTRUCTION true if the time on the stopwatch should be printed when destroyed
    //! @param AUTO_START whether to start the stopwatch upon construction
    Stopwatch::Stopwatch
    (
      const std::string &DESCRIPTION,
      const Time &NOTIFICATION_INTERVAL,
      const Message::MessageLevel &MESSAGE_LEVEL,
      const bool PRINT_ON_DESTRUCTION,
      const bool AUTO_START
    ) :
      m_Description( DESCRIPTION),
      m_TimerStart(), // set to zero
      m_TotalTime(),  // set to zero
      m_LastNotificationTime(),
      m_NotificationInterval( NOTIFICATION_INTERVAL),
      m_MessageLevel( MESSAGE_LEVEL),
      m_PrintOnDestruction( PRINT_ON_DESTRUCTION)
    {
      if( AUTO_START)
      {
        Start();
      }
    }

    //! @brief constructor with description and message level
    //! set the Timer Start to the current time
    //! @param DESCRIPTION string that describes what the stopwatch is timing
    //! @param MESSAGE_LEVEL the message level to use when printing messages
    //! @param PRINT_ON_DESTRUCTION true if the time on the stopwatch should be printed when destroyed
    Stopwatch::Stopwatch
    (
      const std::string &DESCRIPTION,
      const Message::MessageLevel &MESSAGE_LEVEL,
      const bool PRINT_ON_DESTRUCTION,
      const bool AUTO_START
    ) :
      m_Description( DESCRIPTION),
      m_TimerStart(), // set to zero
      m_TotalTime(),  // set to zero
      m_LastNotificationTime(),
      m_NotificationInterval( std::numeric_limits< size_t>::max(), 0),
      m_MessageLevel( MESSAGE_LEVEL),
      m_PrintOnDestruction( PRINT_ON_DESTRUCTION)
    {
      if( AUTO_START)
      {
        Start();
      }
    }

    //! virtual copy constructor
    Stopwatch *Stopwatch::Clone() const
    {
      return new Stopwatch( *this);
    }

    //! virtual destructor - will output the lifetime of this object and though the duration of the scope in which it was initilized
    Stopwatch::~Stopwatch()
    {
      if( m_PrintOnDestruction)
      {
        Stop();
        // only write the message if the stopwatch was actually used.
        if( m_TotalTime.GetSecondsFractional())
        {
          WriteMessage();
        }
      }
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Stopwatch::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return description for that Stopwatch
    //! @return description
    const std::string &Stopwatch::GetDescription() const
    {
      return m_Description;
    }

    //! @brief return the starting time
    //! @return Time object of the current time, when this Stopwatch was constructed
    const Time &Stopwatch::GetLastStartTime() const
    {
      return m_TimerStart;
    }

    //! @brief return the total accumulative Time
    Time Stopwatch::GetTotalTime() const
    {
      // if running,
      return IsRunning() ? m_TotalTime + GetProcessDuration() : m_TotalTime;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief is stopwatch running
    //! @return true, is the stopwatch is currently running
    bool Stopwatch::IsRunning() const
    {
      return !m_TimerStart.IsZero();
    }

    //! @brief start this stopwatch again.
    void Stopwatch::Start()
    {
      // check if not already running
      if( IsRunning())
      {
        return;
      }

      m_TimerStart = Time::GetCurrent();
    }

    //! @brief stop this stopwatch and print out a message if enough time has passed
    void Stopwatch::Stop()
    {
      // check if already stopped
      if( !IsRunning())
      {
        return;
      }

      // add difference between timer start and current time to total time
      m_TotalTime += Time::GetCurrent() - m_TimerStart;
      m_TimerStart.SetZero(); // set start time to zero, to indicate, that the stopwatch is not running

      // check for how long no notification was issued
      if( ( m_TotalTime - m_LastNotificationTime) > m_NotificationInterval)
      {
        m_LastNotificationTime = m_TotalTime;
        WriteMessage();
      }
    }

    //! @brief stops the timer, and sets total time to zero
    void Stopwatch::Reset()
    {
      Stop();
      m_TotalTime.SetZero();
      m_LastNotificationTime.SetZero();
    }

    //! @brief return the time since the stop watch was last started
    Time Stopwatch::GetProcessDuration() const
    {
      // calculate the current duration if it is running
      if( IsRunning())
      {
        return Time::GetCurrent() - m_TimerStart;
      }
      // return a zero Time object
      else
      {
        return Time();
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write a message that contains the total time that this stopwatch has accumulated
    void Stopwatch::WriteMessage() const
    {
      Time total_time( GetTotalTime());
      BCL_Message
      (
        m_MessageLevel,
        m_Description + " has run for " + util::Format()( total_time.GetSecondsFractional()) + " seconds"
      );
    }

    //! @brief Write StopWatch to std::ostream
    std::ostream &Stopwatch::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Description,          OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_TimerStart,           OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_TotalTime,            OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_LastNotificationTime, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NotificationInterval, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MessageLevel,         OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PrintOnDestruction,   OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief Read StopWatch fom std::istream
    std::istream &Stopwatch::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Description,          ISTREAM);
      io::Serialize::Read( m_TimerStart,           ISTREAM);
      io::Serialize::Read( m_TotalTime,            ISTREAM);
      io::Serialize::Read( m_LastNotificationTime, ISTREAM);
      io::Serialize::Read( m_NotificationInterval, ISTREAM);
      io::Serialize::Read( m_MessageLevel,         ISTREAM);
      io::Serialize::Read( m_PrintOnDestruction,   ISTREAM);

      // end
      return ISTREAM;
    }

  } // namespace util
} // namespace bcl
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
#include "util/bcl_util_string_functions.h"

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    //! @brief test whether a given string begins with another string
    //! @param STRING the full string
    //! @param TEST_PREFIX the prefix
    //! @return true if STRING begins with TEST_PREFIX
    bool StartsWith( const std::string &STRING, const std::string &TEST_PREFIX)
    {
      return STRING.size() >= TEST_PREFIX.size()
             && std::equal( STRING.begin(), STRING.begin() + TEST_PREFIX.size(), TEST_PREFIX.begin());
    }

    //! @brief test whether a given string ends with another string
    //! @param STRING the full string
    //! @param TEST_SUFFIX the prefix
    //! @return true if STRING ends with TEST_SUFFIX
    bool EndsWith( const std::string &STRING, const std::string &TEST_SUFFIX)
    {
      return STRING.size() >= TEST_SUFFIX.size()
             && std::equal( STRING.end() - TEST_SUFFIX.size(), STRING.end(), TEST_SUFFIX.begin());
    }

    //! @brief Convert characters ]A-Z] in a string converted to lower case
    //! @param STRING the full string
    //! @return the string with all upper-case characters converted to lower case
    std::string ToLower( const std::string &STRING)
    {
      std::string lowered( STRING);
      for( std::string::iterator itr( lowered.begin()), itr_end( lowered.end()); itr != itr_end; ++itr)
      {
        *itr = tolower( *itr);
      }
      return lowered;
    }

    //! @brief Convert characters [a-z] in a string converted to upper case
    //! @param STRING the full string
    //! @return the string with all lower-case characters converted to upper case
    std::string ToUpper( const std::string &STRING)
    {
      std::string raised( STRING);
      for( std::string::iterator itr( raised.begin()), itr_end( raised.end()); itr != itr_end; ++itr)
      {
        *itr = toupper( *itr);
      }
      return raised;
    }

    //! @brief Remove a specified set of starting and trailing characters from the string
    //! @param STRING the full string
    //! @param CHARS_TO_STRIP the set of characters to strip from the beginning and end of the string
    //! @return the string stripped of starting and ending characters found in CHARS_TO_STRIP
    std::string Strip( const std::string &STRING, const std::string &CHARS_TO_STRIP)
    {
      //searches for first character not to be stripped
      const std::string::size_type pos1( STRING.find_first_not_of( CHARS_TO_STRIP));

      if( pos1 == std::string::npos)
      {
        // all characters stripped
        return std::string();
      }

      //searches for the last character that will not be stripped
      const std::string::size_type pos2( STRING.find_last_not_of( CHARS_TO_STRIP));

      //returns substring from pos1 of length pos2 - pos1 + 1
      return STRING.substr( pos1, pos2 - pos1 + 1);
    }

    //! @brief Remove a specified set of starting characters from the string
    //! @param STRING the full string
    //! @param CHARS_TO_STRIP the set of characters to strip from the beginning of the string
    //! @return the string stripped of starting characters found in CHARS_TO_STRIP
    std::string LStrip( const std::string &STRING, const std::string &CHARS_TO_STRIP)
    {
      //searches for first character not to be stripped
      const std::string::size_type pos( STRING.find_first_not_of( CHARS_TO_STRIP));

      if( pos == std::string::npos)
      {
        // all characters stripped
        return std::string();
      }
      return STRING.substr( pos);
    }

    //! @brief Remove a specified set of trailing characters from the string
    //! @param STRING the full string
    //! @param CHARS_TO_STRIP the set of characters to strip from the end of the string
    //! @return the string stripped of ending characters found in CHARS_TO_STRIP
    std::string RStrip( const std::string &STRING, const std::string &CHARS_TO_STRIP)
    {
      //searches for the last character that will not be stripped
      const std::string::size_type pos( STRING.find_last_not_of( CHARS_TO_STRIP));

      if( pos == std::string::npos)
      {
        // all characters stripped
        return std::string();
      }

      //returns substring from first character to the last non-stripped character
      return STRING.substr( 0, pos + 1);
    }

    //! @brief repeat a string a particular number of times
    //! @param STRING the string to repeat
    //! @param REPETITIONS the number of times to repeat the string
    //! @return the string, repeated REPETITIONS times
    //! @note complexity is O(STRING.size() * log2(REPETITIONS))
    std::string Repeat( const std::string &STRING, const size_t &REPETITIONS)
    {
      size_t repetitions_binary( REPETITIONS);
      std::string repeated( REPETITIONS & 1 ? STRING : std::string());
      repeated.reserve( STRING.size() * REPETITIONS);

      // at each iteration in the loop below, string_doubler doubles in size
      std::string string_doubler( STRING);
      string_doubler.reserve( STRING.size() * REPETITIONS);
      while( repetitions_binary >>= 1)
      {
        string_doubler += string_doubler;

        // if the remaining # of repeats is odd, add the doubled string to the repeated string
        if( repetitions_binary & 1)
        {
          repeated += string_doubler;
        }
      }
      return repeated;
    }

    //! @brief join a vector of strings with a different internal string
    //! @param JOINER the string to insert between consecutive elements of the vector of strings
    //! @param STRINGS array of strings to join with the JOINER
    //! @return STRINGS(0){JOINER}{STRINGS(1)}{JOINER}...{STRINGS(N)}
    //! @note complexity is O(STRING.size() * log2(REPETITIONS))
    std::string Join( const std::string &JOINER, const storage::Vector< std::string> &STRINGS)
    {
      std::string joined;
      if( STRINGS.IsEmpty())
      {
        return joined;
      }

      // initialize with the 1st element
      storage::Vector< std::string>::const_iterator itr( STRINGS.Begin()), itr_end( STRINGS.End());
      joined = *itr;

      // join the remaining elements
      for( ++itr; itr != itr_end; ++itr)
      {
        joined += JOINER;
        joined += *itr;
      }
      return joined;
    }

    //! returns a vector of std::string from an array of char *
    storage::Vector< std::string> StringListFromCharacterArray( const int NUMBER_ARGUMENTS, const char *ARGUMENTS[])
    {
      //list of strings
      storage::Vector< std::string> string_list;
      string_list.AllocateMemory( NUMBER_ARGUMENTS);

      // copy arguments
      for( int i( 0); i < NUMBER_ARGUMENTS; ++i)
      {
        string_list.PushBack( std::string( ARGUMENTS[ i]));
      }

      //end
      return string_list;
    }

    //! reads all strings from an istream and stores them consecutively in a storage vector of strings
    storage::Vector< std::string> StringListFromIStream( std::istream &ISTREAM)
    {
      //initialize string_list
      storage::Vector< std::string> string_list;

      //read strings from the ISTREAM till the end
      while( !ISTREAM.eof())
      {
        // get the next character
        char c( ISTREAM.get());
        if( !isspace( c) && !ISTREAM.eof())
        {
          // add a new string at the end of the vector and get a reference on it
          string_list.PushBack();
          std::string &current_argument( string_list.LastElement());

          // keep track of whether we are in quotes, escapes, etc. quotes and escapes
          bool in_quotes( false);
          bool in_escape( false);
          char quotes_char( '"'); // if in_quotes is true, this is the type of quotes that we are presently in

          do
          {
            if( in_escape)
            {
              in_escape = false;
              current_argument += c;
            }
            else if( c == '\\')
            {
              in_escape = true;
            }
            else if( in_quotes)
            {
              if( c == quotes_char)
              {
                in_quotes = false;
              }
              else
              {
                current_argument += c;
              }
            }
            else if( c == '\'' || c == '"')
            {
              in_quotes = true;
              quotes_char = c;
            }
            else
            {
              current_argument += c;
            }

            c = char( ISTREAM.get());
          } while( !ISTREAM.eof() && ( !isspace( c) || in_quotes));
        }
      }

      //end
      return string_list;
    }

    //! @brief SplittedStringLineListFromIStream reads all lines from an istream, splits them, and then stores them
    //! @param ISTREAM the stream from which the lines will be read
    //! @param SPLITTER the character to be used to split the line strings
    //! @return Vector of strings which are each consecutive line coming from ISTREAM
    storage::Vector< storage::Vector< std::string> > SplittedStringLineListFromIStream
    (
      std::istream &ISTREAM, const std::string &SPLITTER
    )
    {
      // create const storage::Vector "lines" and initialize with all the unsplitted lines in "ISTREAM"
      const storage::Vector< std::string> lines( StringLineListFromIStream( ISTREAM));

      // create storage::Vector of storage::Vector of strings "splitted_lines" which will hold all the splitted lines
      // of "lines"
      storage::Vector< storage::Vector< std::string> > splitted_lines;

      // iterate through "lines" in order to split the lines into the individual strings that make up each line
      for
      (
        storage::Vector< std::string>::const_iterator iter( lines.Begin()), iter_end( lines.End());
        iter != iter_end;
        ++iter
      )
      {
        // add the strings of the line currently denoted by "itr" to "splitted_lines"
        if( !iter->empty())
        {
          splitted_lines.PushBack( SplitString( *iter, SPLITTER));
        }
      }

      // return the lines of "ISTREAM" which have been split into the strings that make them up
      return splitted_lines;
    }

    //! returns a vector of substrings of STRING, that don't contain SPLITTER
    storage::Vector< std::string> SplitString( const std::string &STRING, const std::string &SPLITTER)
    {
      // initialize results
      storage::Vector< std::string> result;
      // Skip delimiters at beginning.
      std::string::size_type last_pos( STRING.find_first_not_of( SPLITTER, 0));
      // Find first "non-delimiter".
      std::string::size_type pos( STRING.find_first_of( SPLITTER, last_pos));

      // as long as something was found
      while( std::string::npos != pos || std::string::npos != last_pos)
      {
        // Found a token, add it to the vector.
        result.PushBack( STRING.substr( last_pos, pos - last_pos));
        // Skip delimiters.  Note the "not_of"
        last_pos = STRING.find_first_not_of( SPLITTER, pos);
        // Find next "non-delimiter"
        pos = STRING.find_first_of( SPLITTER, last_pos);
      }

      // end
      return result;
    }

    //! returns a string which has all REPLACE replaced with REPLACEMENT in ORIGINAL
    std::string ReplaceString( const std::string &ORIGINAL, const std::string &REPLACE, const std::string &REPLACEMENT)
    {
      // split class identifier at REPLACE
      const storage::Vector< std::string> split_original( SplitString( ORIGINAL, REPLACE));

      // is there anything to replace
      if( split_original.GetSize() <= 1)
      {
        return ORIGINAL;
      }

      // initialize filename with the first part without REPLACEMENT in front
      std::string result( split_original.FirstElement());
      // loop over remaining parts and add them to filename after REPLACEMENT
      for
      (
        storage::Vector< std::string>::const_iterator name_itr( split_original.Begin() + 1),
          name_itr_end( split_original.End());
        name_itr != name_itr_end;
        ++name_itr
      )
      {
        result += REPLACEMENT + *name_itr;
      }

      return result;
    }

    //! trims spaces from beginning and end of a copied string
    std::string TrimString( const std::string &STRING)
    {
      return Strip( STRING, " \n\t\r\0");
    }

    //! test whether string is numerical (double, float, size_t, int) with '.', '-', leading and tailing spaces are allowed
    bool IsNumerical( const std::string &STRING)
    {
      return LengthOfFloatingPointType( STRING) == STRING.size();
    }

    //! @brief Find the number of characters that belong to a character string that can be converted to a numerical type
    //! @brief (e.g. float, double, int, size_t, etc.), including whitespace
    //! @param STRING the string to search
    //! @param START is where number (or whitespace) should begin in the string
    //! @return the number of characters in the number starting at STRING[ START]
    //! @note returns 0 if STRING[ START]...STRING[START+X]
    //! @note this function can handle precisely those numerical formats that can be loaded with stream input operators,
    size_t LengthOfFloatingPointType( const std::string &STRING, const size_t &START)
    {
      size_t i( START);

      // Operating on the string directly is inefficient because it requires that every
      // loop check that i < STRING.size().  With the character array, the last character
      // is guaranteed to be a '\0' and because '\0' is an invalid character for a double,
      // none of the loops would continue if they reached it.
      const char *my_string = STRING.c_str();

      while( my_string[ i] == ' ' || my_string[ i] == '\t') // Go past optional spaces
      {
        ++i;
      }

      if( my_string[ i] == '-' || my_string[ i] == '+') // Optional sign after spaces
      {
        ++i;
      }

      if( isdigit( my_string[ i]))
      {
        ++i; // MYSTRING[ i] is a digit; no need to check it again

        while( isdigit( my_string[ i])) // go through the digits
        {
          ++i;
        }

        if( my_string[ i] == '.') // allow a decimal
        {
          ++i; // step over the decimal
          while( isdigit( my_string[ i])) // and then any digits that follow
          {
            ++i;
          }
        }
      }
      else if( my_string[ i] == '.' && isdigit( my_string[ i + 1])) // number that starts with decimal, e.g. .9
      {
        ++i; // step over the decimal
        while( isdigit( my_string[ i])) // go through the numbers after the decimal point
        {
          ++i;
        }
      }
      else // something other than a number
      {
        return size_t( 0);
      }

      // lastly, account for the optional exponent
      // 'i' will move forward iff there is a sign or number after the e/E
      if( my_string[ i] == 'e' || my_string[ i] == 'E') //
      {
        if( isdigit( my_string[ i + 1])) // implied positive sign after e or E, e.g. 1e9 == 1e+9
        {
          i += 2;
        }
        else if( ( my_string[ i + 1] == '-' || my_string[ i + 1] == '+') && isdigit( my_string[ i + 2])) // [eE][+-][:digit] after the number
        { // [eE]?[+-][[:digits:]+] after the mantissa
          i += 3;
        }

        while( isdigit( my_string[ i])) // walk past the remainder of the exponent, if any
        {
          ++i;
        }
      }

      while( my_string[ i] == ' ' || my_string[ i] == '\t') // Go past optional spaces at end
      {
        ++i;
      }

      // i is now 1 past the last character that could be considered part of the number, so
      // returning i-START gives the number of characters in the array.
      return ( i - START);
    }

    //! Find the number of consecutive characters in a string that can be converted to a single integer type
    //! (e.g. long, int, short), including whitespace
    //! @param STRING the string to search
    //! @param START is where the number (or whitespace) should begin in the string
    //! @return an integer X such that STRING[ START] to STRING[ START+X-1] is convertable to a long/integer/short
    //! @note returns 0 if there is no such X.
    //! @note this function can handle precisely those numerical formats that can be loaded with stream input operators
    size_t LengthOfIntegerType( const std::string &STRING, const size_t &START)
    {
      size_t i = START;

      // Operating on the string directly is inefficient because it requires that every
      // loop check that i < STRING.size().  With the character array, the last character
      // is guaranteed to be a '\0' and because '\0' is an invalid character for a numeric type,
      // none of the loops would continue if they reached it.
      const char *const &my_string( STRING.c_str());

      while( my_string[ i] == ' ' || my_string[ i] == '\t') // Optional spaces
      {
        i++;
      }

      if( my_string[ i] == '+' || my_string[ i] == '-') // Optional sign after spaces
      {
        i++;
      }

      if( isdigit( my_string[ i]))
      {
        ++i;
        while( isdigit( my_string[ i]))
        {
          ++i;
        }
        while( my_string[ i] == ' ' || my_string[ i] == '\t') // Go past optional spaces at end
        {
          i++;
        }
        return ( i - START);
      }
      else // no digits found; go back to start
      {
        return size_t( 0);
      }
    }

    //! Find the number of consecutive characters in a string that can be converted to a single unsigned integer type
    //! (e.g. size_t or unsigned long/int/short), including whitespace
    //! @param STRING the string to search
    //! @param START is where the number (or whitespace) should begin in the string
    //! @return an integer X such that STRING[ START] to STRING[ START+X-1] is convertable to a long/integer/short
    //! @note returns 0 if there is no such X.
    //! @note this function can handle precisely those numerical formats that can be loaded with stream input operators
    size_t LengthOfUnsignedIntegerType( const std::string &STRING, const size_t &START)
    {
      size_t i = START;

      // Operating on the string directly is inefficient because it requires that every
      // loop check that i < STRING.size().  With the character array, the last character
      // is guaranteed to be a '\0' and because '\0' is an invalid character for a numeric type,
      // none of the loops would continue if they reached it.
      const char *const &my_string = STRING.c_str();

      while( my_string[ i] == ' ' || my_string[ i] == '\t') // Optional spaces
      {
        ++i;
      }

      if( my_string[ i] == '+') // Optional sign after spaces
      {
        ++i;
      }

      if( isdigit( my_string[ i])) // now at the first digit.
      {
        ++i;
        while( isdigit( my_string[ i]))
        {
          ++i;
        }
        while( my_string[ i] == ' ' || my_string[ i] == '\t') // Go past optional spaces at end
        {
          i++;
        }
        return ( i - START);
      }
      else // no digits found
      {
        return size_t( 0);
      }
    }

    //! searches for spaces and removes them from the string
    std::string RemoveSpacesFromString( const std::string &STRING)
    {
      std::string cleaned_string( STRING);

      // erase every space
      cleaned_string.erase
      (
        std::remove( cleaned_string.begin(), cleaned_string.end(), ' '),
        cleaned_string.end()
      );

      return cleaned_string;
    }

    //! @brief ConvertStringToBoolean takes a string and converts it to a boolean value
    //! @param STRING the string that will be converted to a boolean, should be either "true" or "false"
    //! @return boolean either true of false depending on STRING
    bool
    ConvertStringToBoolean( const std::string &STRING)
    {
      // remove tabs, spaces, new lines, etc.
      const std::string string( TrimString( STRING));

      // make sure the string is either true or false
      BCL_Assert
      (
        string == "true" || string == "false",
        "string to convert to boolean is neither \"true\" nor \"false\" but is \"" + string + "\""
      );

      // return whether the string was true
      return string == "true";
    }

    //! @brief wrap a string for a given line length into multiple lines
    //! @details like split string, only that it will also split, if the resulting string is longer that LINE_LENGTH
    //! @param STRING the string to be wrapped
    //! @param LINE_LENGTH length of lines
    //! @param WRAPPER_CHARS all characters to wrap at
    //! @return vector of lines
    storage::Vector< std::string> WrapString( const std::string &STRING, const size_t LINE_LENGTH, const std::string &WRAPPER_CHARS)
    {
      // stores the wrapped lines
      storage::Vector< std::string> wrapped_lines;

      const std::string::size_type string_length( STRING.length());
      std::string::size_type current_start( 0);
      std::string::size_type current_end( 0);
      std::string::size_type skip( 0); // if wrap is at space, remove that space

      // as long as end was not hit
      while( current_start != string_length)
      {
        // potentially go till end of line
        current_end = current_start + LINE_LENGTH;

        // current end beyond the string
        if( current_end > string_length)
        {
          current_end = string_length;
        }
        // not at the end yet, try to wrap at space
        else
        {
          // find the last space before the max line length
          const std::string::size_type current_space_pos( STRING.rfind( WRAPPER_CHARS, current_end));

          // try to wrap the string at the space position
          if( current_space_pos != std::string::npos && current_space_pos > current_start)
          {
            current_end = current_space_pos;
            skip = 1;
          }
        }

        // print the substring from start with length
        wrapped_lines.PushBack( STRING.substr( current_start, current_end - current_start));

        // update the start to the current end
        current_start = current_end + skip;
        skip = 0;
      }

      // end
      return wrapped_lines;
    }

    //! ignore lines that start with '#' and empty lines
    void ChopHeader( std::istream &ISTREAM)
    {
      std::string line;
      while( ISTREAM.peek() == '#' || ISTREAM.peek() == '\n')
      {
        std::getline( ISTREAM, line);
      }
    }

    //! returns the maximum matching string in terms of pair< size_t, size_t>( index, length) ~ substr( index, length) within first and second string
    storage::Triplet< std::string, storage::Pair< size_t, size_t>, storage::Pair< size_t, size_t> >
    MaximumMatchingSubstring( const std::string &STRING_A, const std::string &STRING_B, const size_t &MINIMUM_SIZE)
    {
      // search second within first string
      // start with maximum string and reduce size iteratively
      size_t size( std::min( STRING_B.size(), STRING_A.size()));
      while( size >= MINIMUM_SIZE)
      {
        // shift window through first sequence
        for( size_t i = 0; i < STRING_B.size() - size + 1; ++i)
        {
          size_t index( STRING_A.find( STRING_B.substr( i, size)));
          if( index < std::string::npos)
          {
            return storage::Triplet< std::string, storage::Pair< size_t, size_t>, storage::Pair< size_t, size_t> >
            (
              STRING_B.substr( i, size),
              storage::Pair< size_t, size_t>( index, size),
              storage::Pair< size_t, size_t>( i, size)
            );
          }
        }
        --size;
      }
      return storage::Triplet< std::string, storage::Pair< size_t, size_t>, storage::Pair< size_t, size_t> >
      (
        std::string( ""),
        storage::Pair< size_t, size_t>( GetUndefined< size_t>(), GetUndefined< size_t>()),
        storage::Pair< size_t, size_t>( GetUndefined< size_t>(), GetUndefined< size_t>())
      );
    }

    //! @brief find the matching brackets revers
    //! @param STRING the string of interest
    //! @param POSITION the position of the closing bracket
    //! @param BRACKET_OPEN the character for opening bracket
    //! @param BRACKET_CLOSE the character for closing bracket
    //! @return the Position of the matching opening bracket
    size_t
    RFindMatchingBracket
    (
      const std::string &STRING,
      const size_t POSITION,
      const char BRACKET_OPEN,
      const char BRACKET_CLOSE
    )
    {
      if( STRING[ POSITION] != BRACKET_CLOSE)
      {
        return POSITION;
      }

      //number of brackets that need to be closed and the current pos in the string
      size_t number_brackets_to_close( 1), current_pos( POSITION);
      //iterate as long as there is a a bracket to close or the rend is reached
      while( number_brackets_to_close > 0 && current_pos != std::string::npos)
      {
        //step back
        --current_pos;
        //increase or decrease number_brackets_to_close
        if( STRING[ current_pos] == BRACKET_CLOSE)
        {
          ++number_brackets_to_close;
        }
        else if( STRING[ current_pos] == BRACKET_OPEN)
        {
          --number_brackets_to_close;
        }
      }

      //check that every pair of brackets are closed
      BCL_Assert
      (
        number_brackets_to_close == 0,
        std::string( "could not find every matching pair of \'")
          + BRACKET_OPEN + "\' and \'" + BRACKET_CLOSE + "\' in \""
          + STRING + "\" starting at " + util::Format()( POSITION)
      );

      //return position pointing to the matching openeing bracket
      return current_pos;
    }

    //! @brief StringLineListFromIStream reads all lines from an istream and stores them consecutively as strings
    //! @param ISTREAM the stream from which the lines will be read
    //! @return Vector of strings which are each consecutive line coming from ISTREAM
    storage::Vector< std::string> StringLineListFromIStream( std::istream &ISTREAM)
    {
      // create storage::Vector "string_lines" to hold the consecutive lines being read in from "ISTREAM"
      storage::Vector< std::string> string_lines;

      // create std::string "current_line" to hold the current line of the accessibility file
      std::string current_line;

      // read in all the line of "ISTREAM"
      while( std::getline( ISTREAM, current_line))
      {
        // add "current_line" to "string_lines"
        string_lines.PushBack( TrimString( current_line));
      }

      // return "string_lines" which has all the lines which were in "ISTREAM"
      return string_lines;
    }

  } // namespace util
} // namespace bcl
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
#include "util/bcl_util_string_numeric_conversion.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically
#include <iostream>
#include <map>
#include <errno.h>
#include <stdlib.h>
#include <string.h>

namespace bcl
{
  namespace util
  {
    // floating point types

    //! @brief a helper function to test for a string that is nan, inf, +inf, or -inf (all case-insensitive)
    //! @param TYPE variable to try to convert STRING into
    //! @param STR string of interest
    //! @return true on success
    template< typename t_DataType>
    bool TryConvertNanInfString( t_DataType &TYPE, const std::string &STR)
    {
      if( STR.size() == size_t( 3))
      {
        if( !strcasecmp( STR.c_str(), "nan"))
        {
          TYPE = GetUndefined< t_DataType>();
          return true;
        }
        else if( !strcasecmp( STR.c_str(), "inf"))
        {
          TYPE = std::numeric_limits< t_DataType>::infinity();
          return true;
        }
      }
      else if( STR.size() == size_t( 4) && !strcasecmp( STR.c_str() + 1, "inf"))
      {
        // check for +/- inf as well
        if( STR[ 0] == '+')
        {
          TYPE = std::numeric_limits< t_DataType>::infinity();
          return true;
        }
        else if( STR[ 0] == '-')
        {
          TYPE = -std::numeric_limits< t_DataType>::infinity();
          return true;
        }
      }
      return false;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString( double &TYPE, const std::string &STR, std::ostream &ERR_STREAM)
    {
      // test that the string is non-empty
      if( STR.empty())
      {
        ERR_STREAM << "Tried to convert an empty string into a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // create an int to hold any existing status of errno
      const int current_err( errno);
      errno = 0;

      // declare a variable to hold the end position for string, which strtod will set
      char *number_end( NULL);

      // get a pointer to the position where a number is expected
      const char *start_ptr( STR.c_str());

      // convert to value
      TYPE = strtod( start_ptr, &number_end);

      // check / reset errno, which is set by strtod in case the value is too large to represent
      if( errno == ERANGE)
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        errno = current_err;
        return false;
      }
      errno = current_err;
      const size_t length( number_end - start_ptr);

      // check for length shorter than STR -> not all the string represents the proper data type
      if( length != STR.size())
      {
        // check for nans and infs (possibly signed) specifically, since the current mingw cross-compiler does not
        // handle them
        if( length == size_t( 0) && TryConvertNanInfString( TYPE, STR))
        {
          return true;
        }
        ERR_STREAM << STR << " could not be converted to a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }
      return true;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString( float &TYPE, const std::string &STR, std::ostream &ERR_STREAM)
    {
      // test that the string is non-empty
      if( STR.empty())
      {
        ERR_STREAM << "Tried to convert an empty string into a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // create an int to hold any existing status of errno
      const int current_err( errno);
      errno = 0;

      // declare a variable to hold the end position for string, which strtof will set
      char *number_end( NULL);

      // get a pointer to the position where a number is expected
      const char *start_ptr( STR.c_str());

      // convert to value
      TYPE = strtof( start_ptr, &number_end);

      // check / reset errno, which is set by strtod in case the value is too large to represent
      if( errno == ERANGE)
      {
        // handle the underflow condition; not really a problem, just set to 0
        if( STR.find( "e-") != std::string::npos)
        {
          ERR_STREAM << "Warning; interpreting " << STR << " as 0 because it is too small to fit in a float\n";
          errno = current_err;
          return true;
        }
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        errno = current_err;
        return false;
      }
      errno = current_err;
      const size_t length( number_end - start_ptr);

      // check for length shorter than STR -> not all the string represents the proper data type
      if( length != STR.size())
      {
        // check for nan's specifically, since the current mingw cross-compiler does not handle them
        // check for nans and infs (possibly signed) specifically, since the current mingw cross-compiler does not
        // handle them
        // check for nans and infs (possibly signed) specifically, since the current mingw cross-compiler does not
        // handle them
        if( length == size_t( 0) && TryConvertNanInfString( TYPE, STR))
        {
          return true;
        }
        ERR_STREAM << STR << " could not be converted to a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }
      return true;
    }

    // character types

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString( unsigned char &TYPE, const std::string &STR, std::ostream &ERR_STREAM)
    {
      if( STR.size() != size_t( 1))
      {
        ERR_STREAM << "Expected a single character, but got \"" << STR << "\"";
        return false;
      }
      memcpy( &TYPE, &STR[ 0], 1);
      return true;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString( char &TYPE, const std::string &STR, std::ostream &ERR_STREAM)
    {
      if( STR.size() != size_t( 1))
      {
        ERR_STREAM << "Expected a single character, but got \"" << STR << "\"";
        return false;
      }
      TYPE = STR[ 0];
      return true;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString
    (
      char &TYPE,
      const char &MIN,
      const char &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    )
    {
      if( STR.size() != size_t( 1))
      {
        ERR_STREAM << "Expected a single character, but got \"" << STR << "\"";
        return false;
      }
      TYPE = STR[ 0];
      return true;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString( signed char &TYPE, const std::string &STR, std::ostream &ERR_STREAM)
    {
      if( STR.size() != size_t( 1))
      {
        ERR_STREAM << "Expected a single character, but got \"" << STR << "\"";
        return false;
      }
      memcpy( &TYPE, &STR[ 0], 1);
      return true;
    }

    // integer types

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString( short &TYPE, const std::string &STR, std::ostream &ERR_STREAM)
    {
      // read into a long
      long new_value;

      // attempt the conversion
      if
      (
        TryConvertFromString
        (
          new_value,
          long( std::numeric_limits< short>::min()),
          long( std::numeric_limits< short>::max()),
          STR,
          ERR_STREAM
        )
      )
      {
        TYPE = short( new_value);
        return true;
      }
      return false;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString( int &TYPE, const std::string &STR, std::ostream &ERR_STREAM)
    {
      // read into a long
      long new_value;

      // attempt the conversion
      if
      (
        TryConvertFromString
        (
          new_value,
          long( std::numeric_limits< int>::min()),
          long( std::numeric_limits< int>::max()),
          STR,
          ERR_STREAM
        )
      )
      {
        TYPE = int( new_value);
        return true;
      }
      return false;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString( long &TYPE, const std::string &STR, std::ostream &ERR_STREAM)
    {
      // test that the string is non-empty
      if( STR.empty())
      {
        ERR_STREAM << "Tried to convert an empty string into a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // declare a variable to hold the end position for string, which strtol will set
      char *number_end( NULL);

      // get a pointer to the position where a number is expected
      const char *start_ptr( STR.c_str());

      // create an int to hold any existing status of errno
      const int current_err( errno);
      errno = 0;

      // convert to value
      TYPE = strtol( start_ptr, &number_end, 10);

      // check / reset errno, which is set by strtod in case the value is too large to represent
      if( errno == ERANGE)
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        errno = current_err;
        return false;
      }
      errno = current_err;
      const size_t length( number_end - start_ptr);

      // check for length shorter than STR -> not all the string represents the proper data type
      if( length != STR.size())
      {
        ERR_STREAM << STR << " could not be converted to a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // in some library implementations strtol does not change errno if the number is representable as an unsigned
      // long, even though it may not be represented as a long.  Likewise, it is prudent to check that the signs of the
      // string and number agree
      if( STR[ 0] == '-')
      {
        if( TYPE > 0)
        {
          ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
          return false;
        }
      }
      else if( TYPE < 0)
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      return true;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString( long long &TYPE, const std::string &STR, std::ostream &ERR_STREAM)
    {
      // test that the string is non-empty
      if( STR.empty())
      {
        ERR_STREAM << "Tried to convert an empty string into a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // declare a variable to hold the end position for string, which strtol will set
      char *number_end( NULL);

      // get a pointer to the position where a number is expected
      const char *start_ptr( STR.c_str());

      // create an int to hold any existing status of errno
      const int current_err( errno);
      errno = 0;

      // convert to value
      TYPE = strtoll( start_ptr, &number_end, 10);

      // check / reset errno, which is set by strtod in case the value is too large to represent
      if( errno == ERANGE)
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        errno = current_err;
        return false;
      }
      errno = current_err;
      const size_t length( number_end - start_ptr);

      // check for length shorter than STR -> not all the string represents the proper data type
      if( length != STR.size())
      {
        ERR_STREAM << STR << " could not be converted to a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // in some library implementations strtol does not change errno if the number is representable as an unsigned
      // long, even though it may not be represented as a long.  Likewise, it is prudent to check that the signs of the
      // string and number agree
      if( STR[ 0] == '-')
      {
        if( TYPE > 0)
        {
          ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
          return false;
        }
      }
      else if( TYPE < 0)
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      return true;
    }

    // unsigned types

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString( unsigned short &TYPE, const std::string &STR, std::ostream &ERR_STREAM)
    {
      // read into a long
      unsigned long new_value;

      // attempt the conversion
      if
      (
        TryConvertFromString
        (
          new_value,
          (unsigned long)( std::numeric_limits< unsigned short>::min()),
          (unsigned long)( std::numeric_limits< unsigned short>::max()),
          STR,
          ERR_STREAM
        )
      )
      {
        TYPE = ( unsigned short)( new_value);
        return true;
      }
      return false;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString( unsigned int &TYPE, const std::string &STR, std::ostream &ERR_STREAM)
    {
      // read into a long
      unsigned long new_value;

      // attempt the conversion
      if
      (
        TryConvertFromString
        (
          new_value,
          (unsigned long)( std::numeric_limits< unsigned int>::min()),
          (unsigned long)( std::numeric_limits< unsigned int>::max()),
          STR,
          ERR_STREAM
        )
      )
      {
        TYPE = ( unsigned int)( new_value);
        return true;
      }
      return false;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString( unsigned long &TYPE, const std::string &STR, std::ostream &ERR_STREAM)
    {
      // test that the string is non-empty
      if( STR.empty())
      {
        ERR_STREAM << "Tried to convert an empty string into a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // declare a variable to hold the end position for string, which strtol will set
      char *number_end( NULL);

      // get a pointer to the position where a number is expected
      const char *start_ptr( STR.c_str());

      // create an int to hold any existing status of errno
      const int current_err( errno);
      errno = 0;

      // convert to value
      TYPE = strtoul( start_ptr, &number_end, 10);

      // check / reset errno, which is set by strtod in case the value is too large to represent
      if( errno == ERANGE)
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        errno = current_err;
        return false;
      }
      errno = current_err;
      const size_t length( number_end - start_ptr);

      // check for length shorter than STR -> not all the string represents the proper data type
      if( length != STR.size())
      {
        ERR_STREAM << STR << " could not be converted to a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // in some library implementations strtol does not change errno if the number is representable as an unsigned
      // long, even though it may not be represented as a long.  Likewise, it is prudent to check that the signs of the
      // string and number agree
      if( STR[ 0] == '-')
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      return true;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString( unsigned long long &TYPE, const std::string &STR, std::ostream &ERR_STREAM)
    {
      // test that the string is non-empty
      if( STR.empty())
      {
        ERR_STREAM << "Tried to convert an empty string into a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // declare a variable to hold the end position for string, which strtol will set
      char *number_end( NULL);

      // get a pointer to the position where a number is expected
      const char *start_ptr( STR.c_str());

      // create an int to hold any existing status of errno
      const int current_err( errno);
      errno = 0;

      // convert to value
      TYPE = strtoull( start_ptr, &number_end, 10);

      // check / reset errno, which is set by strtod in case the value is too large to represent
      if( errno == ERANGE)
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        errno = current_err;
        return false;
      }
      errno = current_err;
      const size_t length( number_end - start_ptr);

      // check for length shorter than STR -> not all the string represents the proper data type
      if( length != STR.size())
      {
        ERR_STREAM << STR << " could not be converted to a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // in some library implementations strtol does not change errno if the number is representable as an unsigned
      // long, even though it may not be represented as a long.  Likewise, it is prudent to check that the signs of the
      // string and number agree
      if( STR[ 0] == '-')
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      return true;
    }

    // bool

    //! @brief create a static map containing all valid boolean strings
    std::map< std::string, bool> GetBoolStrMap()
    {
      // these are from http://yaml.org/type/bool.html, except 0/1, which are bcl extensions
      std::map< std::string, bool> bool_map;
      bool_map[ "0"]     = bool_map[ "n"]     = bool_map[ "N"]     = false;
      bool_map[ "no"]    = bool_map[ "No"]    = bool_map[ "NO"]    = false;
      bool_map[ "false"] = bool_map[ "False"] = bool_map[ "FALSE"] = false;
      bool_map[ "off"]   = bool_map[ "Off"]   = bool_map[ "OFF"]   = false;
      bool_map[ "1"]     = bool_map[ "y"]     = bool_map[ "Y"]     = true;
      bool_map[ "yes"]   = bool_map[ "Yes"]   = bool_map[ "YES"]   = true;
      bool_map[ "true"]  = bool_map[ "True"]  = bool_map[ "TRUE"]  = true;
      bool_map[ "on"]    = bool_map[ "On"]    = bool_map[ "ON"]    = true;
      return bool_map;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString( bool &TYPE, const std::string &STR, std::ostream &ERR_STREAM)
    {
      static const std::map< std::string, bool> s_bool_map( GetBoolStrMap());
      std::map< std::string, bool>::const_iterator itr( s_bool_map.find( STR));
      if( itr == s_bool_map.end())
      {
        ERR_STREAM << "Expected a bool value (0/n/N/no/No/No/false/False/FALSE/off/Off/OFF/1/y/Y/yes/Yes/YES/true/True/TRUE/on/On/ON) but received: " << STR;
        return false;
      }
      TYPE = itr->second;
      return true;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString
    (
      double &TYPE,
      const double &MIN,
      const double &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    )
    {
      // test that the string is non-empty
      if( STR.empty())
      {
        ERR_STREAM << "Tried to convert an empty string into a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // declare a variable to hold the end position for string, which strtod will set
      char *number_end( NULL);

      // get a pointer to the position where a number is expected
      const char *start_ptr( STR.c_str());

      // create an int to hold any existing status of errno
      const int current_err( errno);
      errno = 0;

      // convert to value
      TYPE = strtod( start_ptr, &number_end);

      // check / reset errno, which is set by strtod in case the value is too large to represent
      if( errno == ERANGE)
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        errno = current_err;
        return false;
      }
      errno = current_err;
      const size_t length( number_end - start_ptr);

      // check for length shorter than STR -> not all the string represents the proper data type
      if( length != STR.size())
      {
        // check for nans and infs (possibly signed) specifically, since the current mingw cross-compiler does not
        // handle them
        if( length != size_t( 0) || !TryConvertNanInfString( TYPE, STR))
        {
          ERR_STREAM << STR << " could not be converted to a " << GetStaticClassName( TYPE) << '\n';
          return false;
        }
      }

      if( IsNaN( TYPE) || TYPE < MIN || TYPE > MAX)
      {
        ERR_STREAM << STR << " is outside the allowed range: [" << MIN << ',' << MAX << "]\n";
        return false;
      }
      return true;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString
    (
      float &TYPE,
      const float &MIN,
      const float &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    )
    {
      // test that the string is non-empty
      if( STR.empty())
      {
        ERR_STREAM << "Tried to convert an empty string into a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // declare a variable to hold the end position for string, which strtod will set
      char *number_end( NULL);

      // get a pointer to the position where a number is expected
      const char *start_ptr( STR.c_str());

      // create an int to hold any existing status of errno
      const int current_err( errno);
      errno = 0;

      // convert to value
      TYPE = strtof( start_ptr, &number_end);

      // check / reset errno, which is set by strtod in case the value is too large to represent
      if( errno == ERANGE)
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        errno = current_err;
        return false;
      }
      errno = current_err;
      const size_t length( number_end - start_ptr);

      // check for length shorter than STR -> not all the string represents the proper data type
      if( length != STR.size())
      {
        // check for nans and infs (possibly signed) specifically, since the current mingw cross-compiler does not
        // handle them
        if( length != size_t( 0) || !TryConvertNanInfString( TYPE, STR))
        {
          ERR_STREAM << STR << " could not be converted to a " << GetStaticClassName( TYPE) << '\n';
          return false;
        }
      }

      if( IsNaN( TYPE) || TYPE < MIN || TYPE > MAX)
      {
        ERR_STREAM << STR << " is outside the allowed range: [" << MIN << ',' << MAX << "]\n";
        return false;
      }
      return true;
    }

    // character types

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString
    (
      unsigned char &TYPE,
      const unsigned char &MIN,
      const unsigned char &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    )
    {
      if( STR.size() != size_t( 1))
      {
        ERR_STREAM << "Expected a single character, but got \"" << STR << "\"";
        return false;
      }
      memcpy( &TYPE, &STR[ 0], 1);
      if( TYPE < MIN || TYPE > MAX)
      {
        ERR_STREAM << "Expected a character in the range " << MIN << '-' << MAX << ", but got \"" << STR << "\"";
        return false;
      }
      return true;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString
    (
      signed char &TYPE,
      const signed char &MIN,
      const signed char &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    )
    {
      if( STR.size() != size_t( 1))
      {
        ERR_STREAM << "Expected a single character, but got \"" << STR << "\"";
        return false;
      }
      memcpy( &TYPE, &STR[ 0], 1);
      if( TYPE < MIN || TYPE > MAX)
      {
        ERR_STREAM << "Expected a character in the range " << MIN << '-' << MAX << ", but got \"" << STR << "\"";
        return false;
      }
      return true;
    }

    // integer types

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString
    (
      short &TYPE,
      const short &MIN,
      const short &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    )
    {
      // read into a long
      long new_value;

      // attempt the conversion
      if( TryConvertFromString( new_value, long( MIN), long( MAX), STR, ERR_STREAM))
      {
        TYPE = short( new_value);
        return true;
      }
      return false;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString
    (
      int &TYPE,
      const int &MIN,
      const int &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    )
    {
      // read into a long
      long new_value;

      // attempt the conversion
      if( TryConvertFromString( new_value, long( MIN), long( MAX), STR, ERR_STREAM))
      {
        TYPE = int( new_value);
        return true;
      }
      return false;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString
    (
      long &TYPE,
      const long &MIN,
      const long &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    )
    {
      // test that the string is non-empty
      if( STR.empty())
      {
        ERR_STREAM << "Tried to convert an empty string into a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // declare a variable to hold the end position for string, which strtol will set
      char *number_end( NULL);

      // get a pointer to the position where a number is expected
      const char *start_ptr( STR.c_str());

      // create an int to hold any existing status of errno
      const int current_err( errno);
      errno = 0;

      // convert to value
      TYPE = strtol( start_ptr, &number_end, 10);

      // check / reset errno, which is set by strtod in case the value is too large to represent
      if( errno == ERANGE)
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        errno = current_err;
        return false;
      }
      errno = current_err;
      const size_t length( number_end - start_ptr);

      // check for length shorter than STR -> not all the string represents the proper data type
      if( length != STR.size())
      {
        ERR_STREAM << STR << " could not be converted to a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // in some library implementations strtol does not change errno if the number is representable as an unsigned
      // long, even though it may not be represented as a long.  Likewise, it is prudent to check that the signs of the
      // string and number agree
      if( STR[ 0] == '-')
      {
        if( TYPE > 0)
        {
          ERR_STREAM << STR << " is outside the allowed range: [" << MIN << ',' << MAX << "]\n";
          return false;
        }
      }
      else if( TYPE < 0 || TYPE < MIN || TYPE > MAX)
      {
        ERR_STREAM << STR << " is outside the allowed range: [" << MIN << ',' << MAX << "]\n";
        return false;
      }

      return true;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString
    (
      long long &TYPE,
      const long long &MIN,
      const long long &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    )
    {
      // test that the string is non-empty
      if( STR.empty())
      {
        ERR_STREAM << "Tried to convert an empty string into a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // declare a variable to hold the end position for string, which strtol will set
      char *number_end( NULL);

      // get a pointer to the position where a number is expected
      const char *start_ptr( STR.c_str());

      // create an int to hold any existing status of errno
      const int current_err( errno);
      errno = 0;

      // convert to value
      TYPE = strtoll( start_ptr, &number_end, 10);

      // check / reset errno, which is set by strtod in case the value is too large to represent
      if( errno == ERANGE)
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        errno = current_err;
        return false;
      }
      errno = current_err;
      const size_t length( number_end - start_ptr);

      // check for length shorter than STR -> not all the string represents the proper data type
      if( length != STR.size())
      {
        ERR_STREAM << STR << " could not be converted to a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // in some library implementations strtol does not change errno if the number is representable as an unsigned
      // long, even though it may not be represented as a long.  Likewise, it is prudent to check that the signs of the
      // string and number agree
      if( STR[ 0] == '-')
      {
        if( TYPE > 0)
        {
          ERR_STREAM << STR << " is outside the allowed range: [" << MIN << ',' << MAX << "]\n";
          return false;
        }
      }
      else if( TYPE < 0 || TYPE < MIN || TYPE > MAX)
      {
        ERR_STREAM << STR << " is outside the allowed range: [" << MIN << ',' << MAX << "]\n";
        return false;
      }

      return true;
    }

    // unsigned types

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString
    (
      unsigned short &TYPE,
      const unsigned short &MIN,
      const unsigned short &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    )
    {
      // read into a long
      unsigned long new_value;

      // attempt the conversion
      if( TryConvertFromString( new_value, ( unsigned long)( MIN), ( unsigned long)( MAX), STR, ERR_STREAM))
      {
        TYPE = ( unsigned short)( new_value);
        return true;
      }
      return false;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString
    (
      unsigned int &TYPE,
      const unsigned int &MIN,
      const unsigned int &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    )
    {
      // read into a long
      unsigned long new_value;

      // attempt the conversion
      if( TryConvertFromString( new_value, ( unsigned long)( MIN), ( unsigned long)( MAX), STR, ERR_STREAM))
      {
        TYPE = ( unsigned int)( new_value);
        return true;
      }
      return false;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString
    (
      unsigned long &TYPE,
      const unsigned long &MIN,
      const unsigned long &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    )
    {
      // test that the string is non-empty
      if( STR.empty())
      {
        ERR_STREAM << "Tried to convert an empty string into a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // declare a variable to hold the end position for string, which strtol will set
      char *number_end( NULL);

      // get a pointer to the position where a number is expected
      const char *start_ptr( STR.c_str());

      // create an int to hold any existing status of errno
      const int current_err( errno);
      errno = 0;

      // convert to value
      TYPE = strtoul( start_ptr, &number_end, 10);

      // check / reset errno, which is set by strtod in case the value is too large to represent
      if( errno == ERANGE)
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        errno = current_err;
        return false;
      }
      errno = current_err;
      const size_t length( number_end - start_ptr);

      // check for length shorter than STR -> not all the string represents the proper data type
      if( length != STR.size())
      {
        ERR_STREAM << STR << " could not be converted to a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // in some library implementations strtol does not change errno if the number is representable as an unsigned
      // long, even though it may not be represented as a long.  Likewise, it is prudent to check that the signs of the
      // string and number agree
      if( STR[ 0] == '-' || TYPE < MIN || TYPE > MAX)
      {
        ERR_STREAM << STR << " is outside the allowed range: [" << MIN << ',' << MAX << "]\n";
        return false;
      }

      return true;
    }

    //! @brief overloaded function to attempt conversion of a string into a specific data type
    //! @param TYPE variable to try to convert STRING into
    //! @param MIN the minimum allowed value
    //! @param MAX the maximum allowed value
    //! @param STRING string of interest
    //! @param ERR_STREAM stream to write
    //! @return true if the conversion was successful
    bool TryConvertFromString
    (
      unsigned long long &TYPE,
      const unsigned long long &MIN,
      const unsigned long long &MAX,
      const std::string &STR,
      std::ostream &ERR_STREAM
    )
    {
      // test that the string is non-empty
      if( STR.empty())
      {
        ERR_STREAM << "Tried to convert an empty string into a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // declare a variable to hold the end position for string, which strtol will set
      char *number_end( NULL);

      // get a pointer to the position where a number is expected
      const char *start_ptr( STR.c_str());

      // create an int to hold any existing status of errno
      const int current_err( errno);
      errno = 0;

      // convert to value
      TYPE = strtoull( start_ptr, &number_end, 10);

      // check / reset errno, which is set by strtod in case the value is too large to represent
      if( errno == ERANGE)
      {
        ERR_STREAM << STR << " is outside the range of " << GetStaticClassName( TYPE) << '\n';
        errno = current_err;
        return false;
      }
      errno = current_err;
      const size_t length( number_end - start_ptr);

      // check for length shorter than STR -> not all the string represents the proper data type
      if( length != STR.size())
      {
        ERR_STREAM << STR << " could not be converted to a " << GetStaticClassName( TYPE) << '\n';
        return false;
      }

      // in some library implementations strtol does not change errno if the number is representable as an unsigned
      // long, even though it may not be represented as a long.  Likewise, it is prudent to check that the signs of the
      // string and number agree
      if( STR[ 0] == '-' || TYPE < MIN || TYPE > MAX)
      {
        ERR_STREAM << STR << " is outside the allowed range: [" << MIN << ',' << MAX << "]\n";
        return false;
      }

      return true;
    }

  } // namespace util
} // namespace bcl
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
#include "util/bcl_util_string_replacement.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically
#include <cstring>

namespace bcl
{
  namespace util
  {

    // instantiate s_Instance
    const SiPtr< const ObjectInterface> StringReplacement::s_Instance
    (
      GetObjectInstances().AddInstance( new StringReplacement())
    );

    //! @brief GetMatchContextDescription provides the name of ENUM
    //! @param ENUM - the context for which a name is desired
    const std::string &StringReplacement::GetMatchContextDescription( const StringReplacement::MatchContext &ENUM)
    {
      static const std::string s_descriptors[] =
      {
        "Any",
        "Prefix",
        "Suffix",
        "Word",
        "Exact",
        GetStaticClassName< MatchContext>()
      };
      return s_descriptors[ size_t( ENUM)];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! Construct from data members
    StringReplacement::StringReplacement
    (
      const MatchContext &CONTEXT,
      const std::string  &MATCH,
      const std::string  &REPLACE
    ) :
      m_Context( CONTEXT),
      m_Match( MATCH),
      m_Replacement( REPLACE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new StringReplacer
    StringReplacement *StringReplacement::Clone() const
    {
      return new StringReplacement( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &StringReplacement::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief ReplaceEachIn makes the replacement in a given string
    //! @param STRING the string to perform the replacements in
    //! @return the number of replacements made
    //! Non-recursive, so
    //! StringReplacement x( e_Any, "ab", "abab");
    //! string y( "abcd");
    //! x.ReplaceEachIn( y); // returns 1; now y == "ababcd"
    size_t StringReplacement::ReplaceEachIn( std::string &STRING) const
    {
      size_t count( 0);

      // for each letter in the (non-constant) STRING
      // string size is non-constant, so don't move the .size() out of the for loop
      for( size_t i( 0); i < STRING.size(); i++)
      {
        if( Matches( STRING, i))
        {
          ++count;
          Replace( STRING, i);
        }
      }

      return count;
    }

    //! @brief ReplaceAllIn recursively makes the replacement
    //! @param STRING the string to perform the replacements in
    //! Recursive, so
    //! StringReplacement x( e_Any, "abc", "ab");
    //! string y( "abcc");
    //! ReplaceAllIn( y); // y == "ab"
    void StringReplacement::ReplaceAllIn( std::string &STRING) const
    {
      // must ensure that replacement string does not match the match string or else this would go on forever
      BCL_Assert
      (
        FindNextMatch( m_Replacement, 0) == std::string::npos,
        "Call to ReplaceAllIn( " + STRING + ") is recursive " + GetClassIdentifier()
        + " = " + util::Format()( *this)
      );

      // keep calling the non-recursive version until no more replacements are found
      while( ReplaceEachIn( STRING));
    }

    //! @brief ReplaceEachWithExclusions makes the replacement in a given string, excluding things that are matched
    //!        by other StringReplacements
    //! @param STRING the string to perform the replacements in
    //! @param EXCLUSIONS a list of string replacements.  If the string at a given position matches any exclusion, the
    //!        replacement will not happen
    //! @return the number of replacements made
    //! Non-recursive, so
    //! StringReplacement x( e_Any, "ab", "abab");
    //! string y( "abcd");
    //! x.ReplaceEachIn( y); // returns 1; now y == "ababcd"
    size_t StringReplacement::ReplaceEachWithExclusions
    (
      std::string &STRING,
      const storage::List< StringReplacement> &EXCLUSIONS
    ) const
    {
      size_t count( 0);

      storage::List< StringReplacement>::const_iterator itr_end( EXCLUSIONS.End());

      // for each letter in the (non-constant) STRING
      // string size is non-constant, so don't move the .size() out of the for loop
      for( size_t i( 0); i < STRING.size(); i++)
      {
        if( Matches( STRING, i)) // move to the next position because this position was excluded
        {
          bool matched_exclusion( false);
          for( storage::List< StringReplacement>::const_iterator itr( EXCLUSIONS.Begin()); itr != itr_end; ++itr)
          {
            if( itr->Matches( STRING, i))
            {
              i += itr->GetMatch().size() - 1;
              matched_exclusion = true;
              break;
            }
          }
          if( !matched_exclusion)
          {
            ++count;
            Replace( STRING, i);
          }
        }
      }

      return count;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief finds all matching positions for a given search string and context
    //! @param STRING_TO_SEARCH the string we are searching inside
    //! @return position for all occurrences
    std::vector< size_t> StringReplacement::FindAllMatches
    (
      const std::string &STRING_TO_SEARCH
    ) const
    {
      std::vector< size_t> matching_positions;
      for( size_t pos( 0); pos < STRING_TO_SEARCH.size(); ++pos)
      {
        if( Matches( STRING_TO_SEARCH, pos))
        {
          matching_positions.push_back( pos);
        }
      }

      return matching_positions;
    }

    //! @brief finds all matching positions for a given search string and context
    //! @param STRING_TO_SEARCH the string we are searching inside
    //! @param START_POSITION the position to test for the match
    //! @return position of the next match
    size_t StringReplacement::FindNextMatch
    (
      const std::string &STRING_TO_SEARCH,
      const size_t &START_POSITION
    ) const
    {
      for
      (
        size_t pos( START_POSITION), last_pos( STRING_TO_SEARCH.size());
        pos < last_pos;
        pos++
      )
      {
        if( Matches( STRING_TO_SEARCH, pos))
        {
          return pos;
        }
      }

      return std::string::npos;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &StringReplacement::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_Context, ISTREAM);
      io::Serialize::Read( m_Match, ISTREAM);
      io::Serialize::Read( m_Replacement, ISTREAM);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &StringReplacement::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_Context, OSTREAM, INDENT);
      io::Serialize::Write( m_Match, OSTREAM, INDENT);
      io::Serialize::Write( m_Replacement, OSTREAM, INDENT);
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Matches check whether a given string matches this string replacement at a given position
    //! @param STRING the string we are interested in
    //! @param POSITION the position to test for the match
    //! @return true if STRING matches this string replacement at POSITION
    bool StringReplacement::Matches( const std::string &STRING, const size_t &POSITION) const
    {
      if( m_Context == e_Exact)
      {
        // only match an exact context at the first position if the strings match exactly
        return ( POSITION == 0 && STRING == m_Match);
      }
      else if( POSITION + m_Match.size() <= STRING.size())
      {
        // STRING could contain the string to be matched
        if( m_Context == e_Word || m_Context == e_Prefix) // check the start context
        {
          // Check the letter at POSITION-1, it if exists, make sure it is a valid variable character
          if( POSITION > 0 && !IsNonVariableCharacter( STRING[ POSITION - 1]))
          {
            return false; // start context did not match
          }
        }

        if( m_Context == e_Word || m_Context == e_Suffix)
        {
          // Check the letter just after where the match would end, it if exists
          if
          (
            POSITION + m_Match.size() < STRING.size()
            && !IsNonVariableCharacter( STRING[ POSITION + m_Match.size()])
          )
          {
            return false; // end context did not match
          }
        }

        return !strncmp( STRING.c_str() + POSITION, m_Match.c_str(), m_Match.size());
      }

      return false; // string was not big enough
    }

    //! @brief Replace makes the replacement in a given string at a given position
    //! @param STRING the string we are interested in
    //! @param POSITION the position to test for the match.  This is moved forward to the end of the replaced string
    //!        if a replacement was performed
    void StringReplacement::Replace( std::string &STRING, size_t &POSITION) const
    {
      if( m_Replacement.size() == m_Match.size())
      {
        // same size, no need to take substrings, just overwrite what is already there
        strncpy( &STRING[ POSITION], m_Replacement.c_str(), m_Replacement.size());
      }
      else
      {
        if( POSITION == 0)
        {
          STRING = m_Replacement + SafeSubstr( STRING, POSITION + m_Match.size());
        }
        else
        {
          // different sizes, use substrings
          STRING = STRING.substr( 0, POSITION) + m_Replacement + SafeSubstr( STRING, POSITION + m_Match.size());
        }
      }

      POSITION += m_Replacement.size() - 1;
    }

  ////////////////
  // comparison //
  ////////////////

    //! @brief operator < test StringReplacements for specificity
    //! @param A a string replacement
    //! @param B another string replacement
    //! @return true iff A should be performed before B
    //! A string replacement has priority over another if it either has a more specific context, or an
    //! equivalent context but matches a larger string
    bool operator <( const StringReplacement &A, const StringReplacement &B)
    {
      return size_t( A.m_Context) > size_t( B.m_Context)
             ||
             (
               A.m_Match.size() > B.m_Match.size()
               &&
               (
                 A.m_Context == B.m_Context
                 ||
                 ( A.m_Context == StringReplacement::e_Prefix && B.m_Context == StringReplacement::e_Suffix)
               )
             );
    }

  } // namespace util
} // namespace bcl
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
#include "util/bcl_util_template_instantiations.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_geometry_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

    BCL_EXPIMP_TEMPLATE template class BCL_API FunctionInterface< storage::VectorND< 2, SiPtr< const coord::GeometryInterface> >, bool>;

    BCL_EXPIMP_TEMPLATE template class BCL_API FunctionInterface< storage::VectorND< 2, linal::Vector< double> >, storage::VectorND< 2, linal::Vector< double> > >;

  } // namespace util
} // namespace bcl
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
#include "util/bcl_util_time.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

// http://gethelp.devx.com/techtips/cpp_pro/10min/10min0500.asp
// http://www.linuxjournal.com/node/5574/print
// http://synesis.com.au/software/unixem.html

#if defined (__GNUC__) && !defined(__MINGW32__)
  #include <sys/time.h>
  #include <unistd.h>
  #define msecond_sleep( MSECONDS) usleep( MSECONDS * 1000)
#elif defined (_MSC_VER)
  #include <sys/timeb.h>
  #include <sys/types.h>
  #include <time.h>
  #include <windows.h>
  #include <winsock.h>
  void gettimeofday( struct timeval *t, void *timezone)
  {       struct _timeb timebuffer;
          _ftime64_s( &timebuffer);
          t->tv_sec = long( timebuffer.time); // explicit cast to long, otherwise C4244 compiler warning in VS
          t->tv_usec = 1000 * timebuffer.millitm;
  };
  //in windows sleep is written with capital Sleep
  #define msecond_sleep( MSECONDS) Sleep( MSECONDS)
#elif defined(__MINGW32__)
  #include <sys/time.h>
  #include <windows.h>
  //in windows sleep is written with capital Sleep
  #define msecond_sleep( MSECONDS) Sleep( MSECONDS)
#endif

namespace bcl
{
  namespace util
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    Time::Time() :
      m_Seconds( 0),
      m_Microseconds( 0)
    {
    }

    //! construct Time from Seconds and Microseconds
    Time::Time( const size_t &SECONDS, const size_t &MICROSECONDS) :
      m_Seconds( SECONDS + MICROSECONDS / s_MicroSecondsPerSecond),
      m_Microseconds( MICROSECONDS % s_MicroSecondsPerSecond)
    {
    }

    //! virtual copy constructor
    Time *Time::Clone() const
    {
      return new Time( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Time::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! get seconds
    const size_t &Time::GetSeconds() const
    {
      return m_Seconds;
    }

    //! get microseconds
    const size_t &Time::GetMicroSeconds() const
    {
      return m_Microseconds;
    }

    //! set this to current time
    Time &Time::SetToCurrentTime()
    {
      return *this = GetCurrent();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief is the time zero
    //! @return true is time is zero
    bool Time::IsZero() const
    {
      return m_Seconds == 0 && m_Microseconds == 0;
    }

    //! @brief set the time to zero time
    void Time::SetZero()
    {
      m_Seconds = 0;
      m_Microseconds = 0;
    }

    //! return the current system time
    Time Time::GetCurrent()
    {
      struct timeval tv;
      gettimeofday( &tv, NULL);
      return Time( tv.tv_sec, tv.tv_usec);
    }

    //! return time in the format like Tue Nov 18 10:54:53 2008
    std::string Time::GetTimeAsDate() const
    {
      const time_t currenttime( m_Seconds);
      return TrimString( std::string( ctime( &currenttime)));
    }

    //! returns the time in the format : hour:minute:second
    std::string Time::GetTimeAsHourMinuteSecond() const
    {
      static const Format s_format( Format().W( 2).R().Fill( '0'));
      std::string time_string;
      time_string = s_format( long( m_Seconds) / long( 3600)) + ":"
                    + s_format( ( long( m_Seconds) % long( 3600)) / long( 60)) + ":" + s_format( long( m_Seconds) % long( 60));

      return time_string;
    }

    //! return time as days hour:minute:second
    std::string Time::GetTimeAsDayHourMinuteSecond() const
    {
      static const Format s_format( Format().W( 2).R().Fill( '0'));
      std::string time_string;
      time_string = s_format(   long( m_Seconds) / long( 3600 * 24)) + "d " +
                    s_format( ( long( m_Seconds) % long( 3600 * 24)) / long( 3600)) + ":" +
                    s_format( ( long( m_Seconds) % long( 3600)) / long( 60)) + ":" +
                    s_format(   long( m_Seconds) % long( 60));

      return time_string;
    }

    //! returns the time in the format hour:minute:second.milliseconds
    std::string Time::GetTimeAsHourMinuteSecondMilliSeconds() const
    {
      const std::string time_string
      (
        Format().W( 2).R().Fill( '0')( int( m_Seconds) / 3600) +
        ":" + Format().W( 2).R().Fill( '0')( ( int( m_Seconds) % 3600) / 60) +
        ":" + Format().W( 2).R().Fill( '0')( int( m_Seconds) % 60) +
        "." + Format().W( 3).R().Fill( '0')( m_Microseconds / 1000)
      );

      return time_string;
    }

    //! @brief get seconds as fractional seconds.fractional
    //! @return double as seconds, and microseconds as fraction of a second
    double Time::GetSecondsFractional() const
    {
      return double( m_Seconds) + double( m_Microseconds) / double( s_MicroSecondsPerSecond);
    }

    //! @brief get the time in total milliseconds
    //! @return sum of seconds an microsenconds converted to milliseconds
    size_t Time::GetTotalMilliseconds() const
    {
      return m_Seconds * s_MilliSecondsPerSecond + m_Microseconds / s_MicroSecondsPerMiliSecond;
    }

    //! @brief get the time in total microseconds
    //! @return sum of seconds converted to microseconds and microseconds
    size_t Time::GetTotalMicroseconds() const
    {
      return m_Seconds * s_MicroSecondsPerSecond + m_Microseconds;
    }

    //! @brief delays the process
    //! @param DELAY the time to delay by
    void Time::Delay( const Time &DELAY)
    {
      msecond_sleep( DELAY.GetTotalMilliseconds());
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief subtract time from this time
    Time &Time::operator -=( const Time &TIME)
    {
      //if TIME is larger then this, this is set to zero and returned
      if( TIME.m_Seconds > m_Seconds || ( TIME.m_Seconds == m_Seconds && TIME.m_Microseconds > m_Microseconds))
      {
        return *this = Time();
      }

      //Seconds are larger but microseonds are smaller
      if( TIME.m_Microseconds > m_Microseconds)
      {
        return *this = Time( m_Seconds - TIME.m_Seconds - 1, s_MicroSecondsPerSecond - ( TIME.m_Microseconds - m_Microseconds));
      }

      //Seconds and Microseconds are larger
      return *this = Time( m_Seconds - TIME.m_Seconds, m_Microseconds - TIME.m_Microseconds);
    }

    //! @brief add time to this time
    Time &Time::operator +=( const Time &TIME)
    {
      // add microseconds
      m_Microseconds += TIME.m_Microseconds;

      // add seconds
      m_Seconds += TIME.m_Seconds;

      // if Microseconds exceed
      if( m_Microseconds >= s_MicroSecondsPerSecond)
      {
        m_Seconds += m_Microseconds / s_MicroSecondsPerSecond;
        m_Microseconds %= s_MicroSecondsPerSecond;
      }

      // end
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! Write Time to std::ostream
    std::ostream &Time::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write time as seconds.xxxxxx
      io::Serialize::Write( m_Seconds, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_Microseconds, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! Read Time from std::istream
    std::istream &Time::Read( std::istream &ISTREAM)
    {
      //read member
      io::Serialize::Read( m_Seconds, ISTREAM);
      io::Serialize::Read( m_Microseconds, ISTREAM);

      // end
      return ISTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Create time from a fractional minutes string, e.g. 35.75 means 35 minutes and 45 seconds
    //! @param TIME the time object to initialize
    //! @param INIT the string to use to initialize the time object
    //! @param ERR_STREAM stream to output errors to
    //! @return true on success
    bool Time::CreateTimeFromMinutesString( Time &TIME, const std::string &INIT, std::ostream &ERR_STREAM)
    {
      double minutes;
      if( !TryConvertFromString( minutes, INIT, ERR_STREAM))
      {
        return false;
      }
      // convert minute to seconds
      TIME.m_Seconds = size_t( minutes * s_SecondsPerMinute);
      TIME.m_Microseconds = size_t( minutes * s_SecondsPerMinute * s_MicroSecondsPerSecond) % s_MicroSecondsPerSecond;
      return true;
    }

    //! @brief Convert the given time object to a minutes string
    //! @return the converted string
    std::string Time::ConvertTimeToMinutesString( const Time &TIME)
    {
      return
        Format()
        (
          ( double( TIME.m_Seconds) + double( TIME.m_Microseconds) / double( s_MicroSecondsPerSecond))
          / double( s_SecondsPerMinute)
        );
    }

    //! @brief Time object from number of hours
    //! @param HOURS the number of hours
    //! @return a Time object with seconds for that many hours
    Time Time::CreateTimeFromHours( const size_t HOURS)
    {
      return Time( HOURS * s_SecondsPerHour, 0);
    }

    //! @brief create time object from compiler __DATE__ and __TIME__ macro
    //! @param DATE format has to be: Mmm dd yyyy as it comes from the __DATE__ macro
    //! @param TIME format has to be: hh::mm::ss as is come from the __TIME__ macro
    Time Time::CreateTimeFromCompilerMacro( const std::string &DATE, const std::string &TIME)
    {
      // split data and time into its components
      const storage::Vector< std::string> splitted_date( SplitString( DATE));
      const storage::Vector< std::string> splitted_time( SplitString( TIME, ":"));

      // array containing all three letter months
      static const std::string s_months[] =
      {
        "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
      };

      // get current timeinfo and modify it to the user's choice
      time_t resulting_time;
      struct tm *timeinfo;
      time( &resulting_time);
      timeinfo = localtime( &resulting_time);
      timeinfo->tm_year = ConvertStringToNumericalValue< size_t>( splitted_date( 2)) - 1900;
      timeinfo->tm_mon = size_t( std::find( s_months, s_months + 12, splitted_date( 0)) - s_months);
      timeinfo->tm_mday = ConvertStringToNumericalValue< size_t>( splitted_date( 1));
      timeinfo->tm_hour = ConvertStringToNumericalValue< size_t>( splitted_time( 0));
      timeinfo->tm_min  = ConvertStringToNumericalValue< size_t>( splitted_time( 1));
      timeinfo->tm_sec  = ConvertStringToNumericalValue< size_t>( splitted_time( 2));
      BCL_Assert( timeinfo->tm_mon < 12, "unknown month supplied: " + splitted_date( 0));

      // call mktime to convert to real time
      resulting_time = mktime( timeinfo);

      // create time object from seconds
      return Time( size_t( resulting_time), 0);
    }

    //! @brief operator for subtracting times.
    //! if time 1 smaller time 2 difference will be zero
    //! @param TIME_LHS time that is usually larger than the time that is subtracted
    //! @param TIME_RHS smaller time that is subtracted
    //! @return difference of times - if left was smaller than right, it will be time zero
    Time operator -( const Time &TIME_LHS, const Time &TIME_RHS)
    {
      return Time( TIME_LHS).operator -=( TIME_RHS);
    }

    //! @brief operator for adding times.
    //! @param TIME_LHS left hand side time
    //! @param TIME_RHS right hand side time
    //! @return summation of times
    Time operator +( const Time &TIME_LHS, const Time &TIME_RHS)
    {
      return Time( TIME_LHS).operator +=( TIME_RHS);
    }

    //! @brief equal two time objects
    //! @param TIME_LHS left hand side time
    //! @param TIME_RHS right hand side time
    //! @return true if LHS time is equal to RHS time, false otherwise
    bool operator ==( const Time &TIME_LHS, const Time &TIME_RHS)
    {
      return TIME_LHS.GetSeconds() == TIME_RHS.GetSeconds() && TIME_LHS.GetMicroSeconds() == TIME_RHS.GetMicroSeconds();
    }

    //! @brief not equal than two time objects
    //! @param TIME_LHS left hand side time
    //! @param TIME_RHS right hand side time
    //! @return true if LHS time not equal than RHS time, false otherwise
    bool operator !=( const Time &TIME_LHS, const Time &TIME_RHS)
    {
      return TIME_LHS.GetSeconds() != TIME_RHS.GetSeconds() || TIME_LHS.GetMicroSeconds() != TIME_RHS.GetMicroSeconds();
    }

    //! @brief smaller than two time objects
    //! @param TIME_LHS left hand side time
    //! @param TIME_RHS right hand side time
    //! @return true if LHS time is smaller than RHS time, false otherwise
    bool operator <( const Time &TIME_LHS, const Time &TIME_RHS)
    {
      if( TIME_LHS.GetSeconds() < TIME_RHS.GetSeconds())
      {
        return true;
      }

      if( TIME_LHS.GetSeconds() == TIME_RHS.GetSeconds())
      {
        return TIME_LHS.GetMicroSeconds() < TIME_RHS.GetMicroSeconds();
      }

      // end
      return false;
    }

    //! @brief smaller equal than two time objects
    //! @param TIME_LHS left hand side time
    //! @param TIME_RHS right hand side time
    //! @return true if LHS time is smaller or equal than RHS time, false otherwise
    bool operator <=( const Time &TIME_LHS, const Time &TIME_RHS)
    {
      if( TIME_LHS.GetSeconds() < TIME_RHS.GetSeconds())
      {
        return true;
      }

      if( TIME_LHS.GetSeconds() == TIME_RHS.GetSeconds())
      {
        return TIME_LHS.GetMicroSeconds() <= TIME_RHS.GetMicroSeconds();
      }

      // end
      return false;
    }

    //! @brief larger than two time objects
    //! @param TIME_LHS left hand side time
    //! @param TIME_RHS right hand side time
    //! @return true if LHS time is larger than RHS time, false otherwise
    bool operator >( const Time &TIME_LHS, const Time &TIME_RHS)
    {
      if( TIME_LHS.GetSeconds() > TIME_RHS.GetSeconds())
      {
        return true;
      }

      if( TIME_LHS.GetSeconds() == TIME_RHS.GetSeconds())
      {
        return TIME_LHS.GetMicroSeconds() > TIME_RHS.GetMicroSeconds();
      }

      // end
      return false;
    }

    //! @brief larger equal than two time objects
    //! @param TIME_LHS left hand side time
    //! @param TIME_RHS right hand side time
    //! @return true if LHS time is larger or equal than RHS time, false otherwise
    bool operator >=( const Time &TIME_LHS, const Time &TIME_RHS)
    {
      if( TIME_LHS.GetSeconds() > TIME_RHS.GetSeconds())
      {
        return true;
      }

      if( TIME_LHS.GetSeconds() == TIME_RHS.GetSeconds())
      {
        return TIME_LHS.GetMicroSeconds() >= TIME_RHS.GetMicroSeconds();
      }

      // end
      return false;
    }

  } // namespace util
} // namespace bcl
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
#include "util/bcl_util_undefined.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically
#if defined (__APPLE__)
  #include <math.h>
#elif defined (__GNUC__)
  #include <math.h>     // _isnan function in Linux located in math.h
#elif defined (_MSC_VER) //this block exist to avoid unneccesarry warnings with Visual C++
  #include <float.h>    // isnan function in Windows located in float.h
  #define isnan _isnan  // in float.h isnan is called _isnan
  #define isinf !_finite
#endif
//this block is necessary since for mingw and cygwin, isnan and isinf is defined in namespace std
#if defined (__MINGW32__) || defined(__CYGWIN__)
  #include <cmath> // either include <math.h> and use isnan or include <cmath> and use std::isnan
  using std::isnan;
  using std::isinf;
#endif

namespace bcl
{
  namespace util
  {

    //! @brief specialization of template GetUndefined for float
    template<>
    const float &GetUndefined< float>()
    {
      static const float s_undefined( std::numeric_limits< float>::quiet_NaN());
      return s_undefined;
    }

    //! @brief specialization of template GetUndefined for double
    template<>
    const double &GetUndefined< double>()
    {
      static const double s_undefined( std::numeric_limits< double>::quiet_NaN());
      return s_undefined;
    }

    //! @brief specialization of template GetUndefined for short
    template<>
    const short &GetUndefined< short>()
    {
      static const short s_undefined( std::numeric_limits< short>::max());
      return s_undefined;
    }

    //! @brief specialization of template GetUndefined for unsigned short
    template<>
    const unsigned short &GetUndefined< unsigned short>()
    {
      static const unsigned short s_undefined( std::numeric_limits< unsigned short>::max());
      return s_undefined;
    }

    //! @brief specialization of template GetUndefined for int
    template<>
    const int &GetUndefined< int>()
    {
      static const int s_undefined( std::numeric_limits< int>::max());
      return s_undefined;
    }

    //! @brief specialization of template GetUndefined for unsigned int
    template<>
    const unsigned int &GetUndefined< unsigned int>()
    {
      static const unsigned int s_undefined( std::numeric_limits< unsigned int>::max());
      return s_undefined;
    }

    //! @brief specialization of template GetUndefined for long
    template<>
    const long &GetUndefined< long>()
    {
      static const long s_undefined( std::numeric_limits< long>::max());
      return s_undefined;
    }

    //! @brief specialization of template GetUndefined for unsigned long
    template<>
    const unsigned long &GetUndefined< unsigned long>()
    {
      static const unsigned long s_undefined( std::numeric_limits< unsigned long>::max());
      return s_undefined;
    }

    //! @brief specialization of template GetUndefined for long long
    template<>
    const long long &GetUndefined< long long>()
    {
      static const long long s_undefined( std::numeric_limits< long>::max());
      return s_undefined;
    }

    //! @brief specialization of template GetUndefined for unsigned long
    template<>
    const unsigned long long &GetUndefined< unsigned long long>()
    {
      static const unsigned long long s_undefined( std::numeric_limits< unsigned long long>::max());
      return s_undefined;
    }

    //! @brief specialization of template GetUndefined for char
    template<>
    const char &GetUndefined< char>()
    {
      static const char s_undefined( std::numeric_limits< char>::max());
      return s_undefined;
    }

    const double &GetUndefinedDouble()
    {
      static const double s_undefined( std::numeric_limits< double>::quiet_NaN());
      return s_undefined;
    }

    const size_t &GetUndefinedSize_t()
    {
      static const size_t s_undefined( std::numeric_limits< size_t>::max());
      return s_undefined;
    }

    //! @brief specialization of template GetUndefined for bool
    template<>
    BCL_API
    const bool &GetUndefined< bool>()
    {
      // it is not really possible to define an undefined bool; but in most contexts, false will do
      static const bool s_undefined( false);
      return s_undefined;
    }

    //! function to return whether the supplied UndefinedObject is defined or not
    bool IsDefined( const UndefinedObject &OBJECT)
    {
      return false;
    }

    //! function to return whether the supplied UndefinedObject is defined or not for double
    bool IsDefined( const double DOUBLE)
    {
      return isnan( DOUBLE) == 0 && isinf( DOUBLE) == 0;
    }

    //! function to return whether the supplied UndefinedObject is defined or not for float
    bool IsDefined( const float FLOAT)
    {
      return isnan( FLOAT) == 0 && isinf( FLOAT) == 0;
    }

    //! function to return whether the supplied UndefinedObject is defined or not for unsigned int
    bool IsDefined( const unsigned int INT)
    {
      return INT != GetUndefined< unsigned int>();
    }

    //! function to return whether the supplied UndefinedObject is defined or not for int
    bool IsDefined( const int INT)
    {
      return INT != GetUndefined< int>();
    }

    //! function to return whether the supplied UndefinedObject is defined or not for unsigned long
    bool IsDefined( const unsigned long LONG)
    {
      return LONG != GetUndefined< unsigned long>();
    }

    //! function to return whether the supplied UndefinedObject is defined or not for long
    bool IsDefined( const long LONG)
    {
      return LONG != GetUndefined< long>();
    }

    //! function to return whether the supplied UndefinedObject is defined or not for unsigned short
    bool IsDefined( const unsigned short SHORT)
    {
      return SHORT != GetUndefined< unsigned short>();
    }

    //! function to return whether the supplied UndefinedObject is defined or not for short
    bool IsDefined( const short SHORT)
    {
      return SHORT != GetUndefined< short>();
    }

    //! function to return whether the supplied UndefinedObject is defined or not for unsigned long long
    bool IsDefined( const unsigned long long LONG)
    {
      return LONG != GetUndefined< unsigned long long>();
    }

    //! function to return whether the supplied UndefinedObject is defined or not for long long
    bool IsDefined( const long long LONG)
    {
      return LONG != GetUndefined< long long>();
    }

    //! function to return whether the supplied value is not a number
    bool IsNaN( const double &DOUBLE)
    {
      return DOUBLE != DOUBLE;
    }

    //! function to return whether the supplied value is not a number
    bool IsNaN( const float &FLOAT)
    {
      return FLOAT != FLOAT;
    }

  } // namespace util
} // namespace bcl
