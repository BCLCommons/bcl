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
