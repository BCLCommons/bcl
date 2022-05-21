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

#ifndef BCL_UTIL_OBJECT_DATA_LABEL_TOKENIZER_H_
#define BCL_UTIL_OBJECT_DATA_LABEL_TOKENIZER_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_object_interface.h"

// external includes - sorted alphabetically
#include <iostream>
#include <string>

namespace bcl
{
  namespace util
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ObjectDataLabelTokenizer
    //! @brief Tokenizer used to efficiently parse object data labels using a subset of OGDL
    //!
    //! @see @link example_util_object_data_label_tokenizer.cpp @endlink
    //! @author mendenjl
    //! @date Sep 28, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API ObjectDataLabelTokenizer :
      public ObjectInterface
    {
    public:
    //////////
    // enum //
    //////////

      //! Token type
      enum Type
      {
        e_Start,        //!< Start of label
        e_Scalar,       //!< A scalar value or string
        e_ScopeOpen,    //!< Opening of an inner scope
        e_ScopeClose,   //!< Closing of an inner scope
        e_TagDelimiter, //!< Delimiter between tag and value
        e_ArgDelimiter, //!< Delimiter between arguments
        e_End           //!< End of label
      };

    //////////
    // data //
    //////////

    private:

      std::string m_String;     //!< The string being parsed
      size_t      m_Size;       //!< The size of the string being parsed
      size_t      m_Position;   //!< Current position in the string
      size_t      m_ScopeDepth; //!< Depth of scope (e.g. # of unclosed scope opens)
      Type        m_LastType;   //!< Type of token last returned by Pop()
      Type        m_NextType;   //!< Type of token that will next be returned by Pop()

    public:

      //! single instance of that class
      static const SiPtr< const ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from string
      explicit ObjectDataLabelTokenizer( const std::string &STR, const size_t &POSITION = 0);

      //! @brief Clone function
      //! @return pointer to new ObjectDataLabelTokenizer
      ObjectDataLabelTokenizer *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get the complete string
      //! @return the complete string
      const std::string &GetString() const
      {
        return m_String;
      }

      //! @brief get type of token last returned by Pop
      //! @return the type of token last returned by Pop
      const Type &GetLastTokenType() const;

      //! @brief get type of token that will next be returned by Pop
      //! @return the type of token that will next be returned by Pop
      const Type &GetNextTokenType() const;

      //! @brief get token type
      //! @return the current token type
      const size_t &GetScopeDepth() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief get the string for the next token
      //! @return the string for the next token
      std::string Pop();

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief find the end of the quoted string that starts at a given position, ignoring escaped quotes
      //! @param STR string of interest
      //! @param START_QUOTE_POS the position of the starting quote
      //! @return the end quote position, or std::string::npos if no such quote exists
      static size_t FindEndQuote( const std::string &STR, const size_t &START_QUOTE_POS);

      //! @brief  validate a string, write out errors to a given stream
      //! @param  STR the string to validate
      //! @param  ERR the stream to write errors to
      //! @return true if the string was valid
      static bool Validate( const std::string &STR, std::ostream &ERR);

    private:

      //! @brief update the next token type, after skipping spaces
      void DetermineNextTokenType();

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    };

  } // namespace util
} // namespace bcl

#endif // BCL_UTIL_OBJECT_DATA_LABEL_TOKENIZER_H_

