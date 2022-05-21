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

#ifndef BCL_UTIL_OBJECT_DATA_LABEL_H_
#define BCL_UTIL_OBJECT_DATA_LABEL_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_util_object_interface.h"

// external includes - sorted alphabetically
#include <string>
#include <vector>

namespace bcl
{
  namespace util
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ObjectDataLabel
    //! @brief a structured label, commonly used for I/O of classes via the LabeledObjectInterface
    //! @details labels have a simple function-call-like syntax similar to OGDL (http://ogdl.sourceforge.net/), such as
    //! RDF( steps = 5, step size = 0.5, temperature = 0.1, atom property = Atom_VdwSurfaceArea)
    //!
    //! This class parses data labels into a tree of names and basic values, and can write out property trees
    //! given such an interface
    //!
    //! Example: When parsed, DoStuff( coefficient = 1.0, DoMoreStuff( 5.0, 6.0, 100)) will become
    //! m_NameOrValue == DoStuff
    //! m_Arguments[ 0].m_Label == 1.0
    //! m_Arguments[ 0].m_NameOrValue == 1.0
    //! m_Arugments[ 1].m_NameOrValue == DoMoreStuff
    //! m_Arugments[ 1].m_Arguments[ 0].m_NameOrValue == 5.0
    //! m_Arugments[ 1].m_Arguments[ 1].m_NameOrValue == 6.0
    //! m_Arugments[ 1].m_Arguments[ 2].m_NameOrValue == 100
    //! m_Arugments[ 1].m_Arguments[ {0,1,2} ].m_Arguments.GetSize() == 0
    //!
    //! @see @link example_util_object_data_label.cpp @endlink
    //! @author mendenjl
    //! @date 06/21/10
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ObjectDataLabel :
      public ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      std::string m_Name;  //!< name AKA tag, the part to the left hand side of the = sign, if there was one
      std::string m_Value; //!< value, everything between the = sign and the opening parenthesis

      //! if m_Value refers to a parameterized object, this vector will hold its arguments
      std::vector< ObjectDataLabel> m_Arguments;

      //! Vector of si-ptrs to arguments, sorted by name
      //! used when comparing labels, but only if at least 1 of m_Arguments is a map-like node
      std::vector< ObjectDataLabel *> m_SortedArguments;

      typedef std::vector< ObjectDataLabel>::iterator iterator;

    public:

    //////////
    // data //
    //////////

      typedef std::vector< ObjectDataLabel>::const_iterator const_iterator;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ObjectDataLabel();

      //! @brief constructor from a string representing a label that is to be parsed
      //! @param STR string to parse to make the object data label
      ObjectDataLabel( const std::string &STR);

      //! @brief constructor from a stream containing a data label
      //! @param STREAM stream to pass to get data label
      explicit ObjectDataLabel( std::istream &STREAM);

      //! @brief constructor from arguments
      //! @param ARGS arguments for the data label
      explicit ObjectDataLabel( const storage::Vector< ObjectDataLabel> &ARGS);

      //! @brief constructor from a string representing a name and the arguments (implies type is object)
      //! @param OBJECT_NAME the name of the object
      //! @param ARGS arguments for the data label
      ObjectDataLabel( const std::string &OBJECT_NAME, const storage::Vector< ObjectDataLabel> &ARGS);

      //! @brief constructor from a string representing a name and the arguments (implies type is object)
      //! @param OBJECT_NAME the name of the object
      //! @param ARGS arguments for the data label
      ObjectDataLabel( const std::string &OBJECT_NAME, const std::vector< ObjectDataLabel> &ARGS);

      //! @brief constructor from a string representing a name, unparsed value, and type
      //! @param NAME identifies what the value represents
      //! @param VALUE the value of the string
      ObjectDataLabel( const std::string &NAME, const std::string &VALUE);

      //! @brief constructor from a string with a name and underlying label (implies type is object)
      //! @param NAME identifies what the label represents
      //! @param OBJECT_LABEL the data label from the object
      ObjectDataLabel( const std::string &NAME, const ObjectDataLabel &OBJECT_LABEL);

      //! @brief constructor from a string representing a label, name and the arguments (implies type is object)
      //! @param NAME identifies what the value represents
      //! @param OBJECT_NAME the name or alias of the object's type
      //! @param ARGS arguments for the object
      ObjectDataLabel
      (
        const std::string &NAME,
        const std::string &OBJECT_NAME,
        const storage::Vector< ObjectDataLabel> &ARGS
      );

      //! @brief constructor from a string representing a label, name and a single argument (implies type is object)
      //! @param NAME identifies what the value represents
      //! @param OBJECT_NAME the name or alias of the object's type
      //! @param ARG argument for the object
      ObjectDataLabel
      (
        const std::string &NAME,
        const std::string &OBJECT_NAME,
        const ObjectDataLabel &ARG
      );

      //! @brief copy constructor
      //! @param PARENT the original data label
      ObjectDataLabel( const ObjectDataLabel &PARENT);

      //! @brief Clone function
      //! @return pointer to new ObjectDataLabel
      ObjectDataLabel *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get the character that goes between each argument
      //! @return the delimiter character
      static const char &GetArgumentDelimiter();

      //! @brief get the character that goes between a label and its name / value
      //! @return the delimiter character
      static const char &GetNameValueDelimiter();

      //! @brief get the character that indicates to inline another file in this label
      //! @return the inline file character
      static const char &GetInlineFileDelimiter();

      //! @brief get all the delimiters
      //! @return a string containing all the delimiters
      static const std::string &GetAllDelimiters();

      //! @return the primary identifier of this data label (either the property name or parameter value)
      const std::string &GetName() const
      {
        return m_Name;
      }

      //! @return the value (e.g. lhs=rhs -> rhs is the value)
      const std::string &GetValue() const
      {
        return m_Value;
      }

      //! @brief set the tag / identifier of this data label (e.g. lhs=rhs -> lhs is the name)
      //! @param NEW_NAME new name
      //! @param ALLOW_CHANGE true to allow setting the name if this label was already given a non-empty name
      void SetName( const std::string &NEW_NAME, const bool &ALLOW_CHANGE = false);

      //! @brief set the value of this data label (e.g. lhs=rhs -> rhs is the value)
      //! @param VALUE new value
      //! @param ALLOW_CHANGE true to allow setting the value if this label was already given a non-empty value
      void SetValue( const std::string &VALUE, const bool &ALLOW_CHANGE = false);

      //! @brief get the beginning of the arguments
      //! @return the beginning of the arguments
      const_iterator Begin() const;

      //! @brief get the beginning of the arguments
      //! @return the beginning of the arguments
      const_iterator End() const;

      //! @brief find an argument with the specified name.  If none is found, return End()
      //! @param NAME the value to find
      //! @param RECURSIVE whether to check sub-arguments (after checking arguments)
      const_iterator FindName( const std::string &NAME, const bool &RECURSIVE = true) const;

      //! @brief find an argument with the specified value.  If none is found, return End()
      //! @param VALUE the value to find
      //! @param RECURSIVE whether to check sub-arguments (after checking arguments)
      const_iterator FindValue( const std::string &VALUE, const bool &RECURSIVE = true) const;

      //! @brief replace values == VALUE with REPLACEMENT, up to MAX_DEPTH below this argument
      //! @param VALUE the value to find
      //! @param REPLACEMENT what to replace the value with
      //! @param MAX_DEPTH maximum depth to perform replacement at(0 = just this object data label)
      //! @return the number of replacements performed
      size_t ReplaceValue( const std::string &VALUE, const std::string &REPLACEMENT, const size_t &MAX_DEPTH = size_t( -1));

      //! @brief find an argument with the specified object data-label.  If none is found, return End()
      //! @param LABEL the label to find
      //! @param RECURSIVE whether to check sub-arguments (after checking arguments)
      //! @param CONSIDER_NAME whether to force a match on the name/tag given in LABEL; otherwise, any matching value
      //!        will be found
      const_iterator Find
      (
        const ObjectDataLabel &LABEL,
        const bool &RECURSIVE = true,
        const bool &CONSIDER_NAME = true
      ) const;

      //! @brief compute the difference between this object data label and another
      //! @param OTHER the other object data label
      ObjectDataLabel Difference( const ObjectDataLabel &OTHER) const;

      //! @return the primary identifier of this data label (either the property name or parameter value)
      const std::vector< ObjectDataLabel> &GetArguments() const;

      //! @brief test whether this label is scalar (e.g. does not have any arguments)
      //! @return true if this label is scalar (e.g. does not have any arguments)
      bool IsScalar() const;

      //! @brief test whether this label is empty (e.g. does not have any arguments, name, or value)
      //! @return true if this label is empty (e.g. does not have any arguments, name, or value)
      bool IsEmpty() const;

      //! @brief get the number of arguments for this label
      //! @return the number of arguments for this label
      size_t GetNumberArguments() const;

      //! @return the primary identifier of this data label (either the property name or parameter value)
      const ObjectDataLabel &GetArgument( const size_t &POSITION) const;

    ////////////////
    // operations //
    ////////////////

      //! @brief assignment from another object data label
      //! @param LABEL the other label
      //! @return a reference to this
      ObjectDataLabel &operator =( const ObjectDataLabel &LABEL);

      //! @brief assign from a data-label string
      ObjectDataLabel &operator =( const std::string &LABEL);

      //! @brief try to assign the label from a string
      //! @param STR string to attempt assignment from
      //! @param ERR error stream, stream to output errors to
      bool TryAssign( const std::string &STR, std::ostream &ERR_STREAM);

      //! @brief ternary comparison operator (like strcmp), returns -1 if *this < A, 0 if *this == A, 1 if *this > A
      //! @param A the object data label to compare with
      //! @return -1 if *this < A, 0 if *this == A, 1 if *this > A
      int TernaryCompare( const ObjectDataLabel &LABEL) const;

      //! @brief Compare this label to another without considering the name/tag
      //! @param LABEL the object data label to compare with
      //! @return -1 if *this < LABEL, 0 if *this == LABEL, 1 if *this > LABEL
      int TernaryCompareWithoutName( const ObjectDataLabel &LABEL) const;

      //! @brief get the string representation of this data label, without the label
      //! @param LINE_LENGTH maximum desired length of each line
      //! @param INDENT amount to indent separate lines by
      //! @param MAX_DEPTH_FOR_SPLIT maximum depth of parameters to split onto separate lines
      //! @return the property label converted to a string (always succeeds)
      std::string ToString
      (
        const size_t &LINE_LENGTH,
        const size_t &INDENT = 0,
        const size_t &MAX_DEPTH_FOR_SPLIT = size_t( -1)
      ) const;

      //! @brief operator for string conversion
      operator std::string () const
      {
        return ToString( size_t( -1));
      }

      //! @brief get the string representation of this data label formatted for the logger
      //! This is equivalent to ToString( GetLogger().GetMaxLineWidth(), INDENT);
      //! @param INDENT amount to indent separate lines by
      //! @param MAX_DEPTH_FOR_SPLIT maximum depth of parameters to split onto separate lines
      //! @return the property label converted to a string (always succeeds)
      std::string ToStringForLogger
      (
        const size_t &INDENT = 0,
        const size_t &MAX_DEPTH_FOR_SPLIT = size_t( -1)
      ) const;

      //! @brief get the string representation of this data label formatted with LoggerInterface::GetDefaultMaxLineWidth
      //! This is equivalent to ToString( LoggerInterface::GetDefaultMaxLineWidth(), INDENT);
      //! @param INDENT amount to indent separate lines by
      //! @param MAX_DEPTH_FOR_SPLIT maximum depth of parameters to split onto separate lines
      //! @return the property label converted to a string (always succeeds)
      std::string ToStringDefaultWidth
      (
        const size_t &INDENT = 0,
        const size_t &MAX_DEPTH_FOR_SPLIT = size_t( -1)
      ) const;

      //! @brief get the string representation of this data label, without the label
      //! @return the property label converted to a string (always succeeds)
      std::string ToString() const;

      //! @brief get the string representation of this data label
      //! @return the property label converted to a string (always succeeds)
      std::string ToNamedString() const;

      //! @brief get the string representation of this data label
      //! @param LINE_LENGTH maximum desired length of each line
      //! @param INDENT amount to indent separate lines by
      //! @param MAX_DEPTH_FOR_SPLIT maximum depth of parameters to split onto separate lines
      //! @return the property label converted to a string (always succeeds)
      std::string ToNamedString
      (
        const size_t &LINE_LENGTH,
        const size_t &INDENT = 0,
        const size_t &MAX_DEPTH_FOR_SPLIT = size_t( -1)
      ) const;

      //! @brief get the string representation of this data label formatted for the logger
      //! This is equivalent to ToNamedString( GetLogger().GetMaxLineWidth(), INDENT);
      //! @param INDENT amount to indent separate lines by
      //! @param MAX_DEPTH_FOR_SPLIT maximum depth of parameters to split onto separate lines
      //! @return the property label converted to a string (always succeeds)
      std::string ToNamedStringForLogger
      (
        const size_t &INDENT = 0,
        const size_t &MAX_DEPTH_FOR_SPLIT = size_t( -1)
      ) const;

      //! @brief get the string representation of this data label formatted with a default max line width
      //! This is equivalent to ToNamedString( LoggerInterface::GetDefaultMaxLineWidth(), INDENT);
      //! @param INDENT amount to indent separate lines by
      //! @param MAX_DEPTH_FOR_SPLIT maximum depth of parameters to split onto separate lines
      //! @return the property label converted to a string (always succeeds)
      std::string ToNamedStringDefaultWidth
      (
        const size_t &INDENT = 0,
        const size_t &MAX_DEPTH_FOR_SPLIT = size_t( -1)
      ) const;

      //! @brief get the arguments of the data label printed as a string
      //! @param INCLUDE_PARENS whether to surround the string with () (default=true)
      //! @param ARGUMENT_DELIMITER what to delimit the arguments of the top level data label by
      //! @return the arguments, printed as a single-line string, surrounded by (), if desired
      std::string ArgumentsToString
      (
        const bool &INCLUDE_PARENS = true,
        const char &ARGUMENT_DELIMITER = GetArgumentDelimiter()
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! @return the number of ( minus the number of ) in STR, ignoring quoted sections
      static long GetParenthesisDepthChange( const std::string &STR);

      //! @brief add double quotes around a string, and escape any internal values
      //! @param STR the string of interest
      //! @return the string, quoted only if necessary
      static std::string QuoteStringIfDelimitersPresent( const std::string &STR);

      //! @brief helper function to determine the length of the string returned by ToString
      //! @details used to determine whether to split labels out on separate lines
      //! @param LENGTH the maximum length of interest
      //! @param CONSIDER_NAME whether to consider the length of the name of this label
      //! @param INDENT the indent to subtract from the maximum length
      //! @return the size of the string returned by ToString (or ToNamedString if consider_name is true)
      //!         or any size > length, if the length would exceed LENGTH
      size_t GetLimitedLength( const size_t &LENGTH, const bool &CONSIDER_NAME, const size_t &INDENT) const;

      //! @brief get the total size of the tree
      //! @return number of nodes in the tree
      size_t GetTreeSize() const;

      //! @brief read a string using the tokenizer
      //! @param TOKENIZER tokenizer created by the parent
      //! @param ERR_STREAM stream that errors should be output to
      //! @return true on success
      bool ReadSubLabel( ObjectDataLabelTokenizer &TOKENIZER, std::ostream &ERR_STREAM);

      //! @brief sort all arguments by name
      void SortArgumentsByName();

      //! @brief copy a sorted arguments vector from another object data label
      //! @param LABEL the label to copy it from
      void CopyArgumentSorting( const ObjectDataLabel &LABEL);

    }; // class ObjectDataLabel

    //! operator == (Comparison)
    bool operator ==( const ObjectDataLabel &LABEL_A, const ObjectDataLabel &LABEL_B);

    //! operator != (Comparison)
    bool operator !=( const ObjectDataLabel &LABEL_A, const ObjectDataLabel &LABEL_B);

    //! operator < (Comparison)
    bool operator <( const ObjectDataLabel &LABEL_A, const ObjectDataLabel &LABEL_B);

    //! operator > (Comparison)
    bool operator >( const ObjectDataLabel &LABEL_A, const ObjectDataLabel &LABEL_B);

  } // namespace util
} // namespace bcl

#endif // BCL_UTIL_OBJECT_DATA_LABEL_H_

