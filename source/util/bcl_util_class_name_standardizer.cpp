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
