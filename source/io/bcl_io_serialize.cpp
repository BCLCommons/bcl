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
#include "io/bcl_io_serialize.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_serializable_interface.h"
#include "util/bcl_util_string_functions.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace io
  {

    //! @brief helper function used for writing any numeric type, using nan if the value is undefined
    //! @param VALUE the numeric value to write to the stream
    //! @param OSTREAM the stream to write the value to
    //! @param INDENT - number of indentations
    //! @return reference to the stream
    template< typename t_NumericType>
    std::ostream &WriteNumericType
    (
      const t_NumericType &VALUE,
      std::ostream &OSTREAM,
      const size_t INDENT
    )
    {
      // insert indents
      Serialize::InsertIndent( OSTREAM, INDENT);

      if( util::IsDefined( VALUE))
      {
        // write the value
        OSTREAM << VALUE;
      }
      else
      {
        // write nan, since the actual undefined value is machine-dependent
        OSTREAM << "nan";
      }

      // end
      return OSTREAM;
    }

  //////////////////////
  // write operations //
  //////////////////////

    //! @brief write OBJECT to ostream
    //! @param OBJECT ObjectInterface derived class
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return reference to OSTREAM
    std::ostream &Serialize::Write
    (
      const util::ObjectInterface &OBJECT,
      std::ostream &OSTREAM,
      const size_t INDENT
    )
    {
      // write the object with the given indents
      OBJECT.WriteObject( OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief write OBJECT to ostream
    //! use this function, if there is ambiguity when writing std::string or other bcl objects derived from std class
    //! @param OBJECT ObjectInterface derived class
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return reference to OSTREAM
    std::ostream &Serialize::WriteObject
    (
      const util::ObjectInterface &OBJECT,
      std::ostream &OSTREAM,
      const size_t INDENT
    )
    {
      // write the object with the given indents
      OBJECT.WriteObject( OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief write OBJECT to ostream
    //! @param OBJECT ObjectInterface derived class
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return reference to OSTREAM
    std::ostream &Serialize::Write
    (
      const util::SerializableInterface &OBJECT,
      std::ostream &OSTREAM,
      const size_t INDENT
    )
    {
      return OBJECT.Write( OSTREAM, INDENT);
    }

    //! @brief write OBJECT to ostream
    //! use this function, if there is ambiguity when writing std::string or other bcl objects derived from std class
    //! @param OBJECT util::ObjectInterface derived class
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return reference to OSTREAM
    std::ostream &Serialize::WriteSerializable
    (
      const util::SerializableInterface &OBJECT,
      std::ostream &OSTREAM,
      const size_t INDENT
    )
    {
      return OBJECT.Write( OSTREAM, INDENT);
    }

    //! @brief write double to ostream
    //! @param DOUBLE the double to be written
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return reference to OSTREAM
    std::ostream &Serialize::Write
    (
      const double &DOUBLE,
      std::ostream &OSTREAM,
      const size_t INDENT
    )
    {
      // insert indents
      InsertIndent( OSTREAM, INDENT);

      if( !util::IsDefined( DOUBLE))
      {
        // write nan, since the actual undefined value is machine-dependent
        OSTREAM << "nan";
      }
      else
      {
        // write DOUBLE
        OSTREAM << DOUBLE;
      }

      // end
      return OSTREAM;
    }

    //! @brief write float to ostream
    //! @param FLOAT the float to be written
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return reference to OSTREAM
    std::ostream &Serialize::Write
    (
      const float &FLOAT,
      std::ostream &OSTREAM,
      const size_t INDENT
    )
    {
      // insert indents
      InsertIndent( OSTREAM, INDENT);

      if( !util::IsDefined( FLOAT))
      {
        // write nan, since the actual undefined value is machine-dependent
        OSTREAM << "nan";
      }
      else
      {
        // write FLOAT
        OSTREAM << FLOAT;
      }

      // end
      return OSTREAM;
    }

    //! @brief write boolean to ostream
    //! @param BOOLEAN the boolean to be written
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return reference to OSTREAM
    std::ostream &Serialize::Write
    (
      const bool &BOOLEAN,
      std::ostream &OSTREAM,
      const size_t INDENT
    )
    {
      // insert indents
      InsertIndent( OSTREAM, INDENT);

      // write BOOLEAN
      OSTREAM << int( BOOLEAN);

      // end
      return OSTREAM;
    }

    //! @brief write short to ostream
    //! @param SHORT the short to be written
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return reference to OSTREAM
    std::ostream &Serialize::Write
    (
      const short &SHORT,
      std::ostream &OSTREAM,
      const size_t INDENT
    )
    {
      return WriteNumericType( SHORT, OSTREAM, INDENT);
    }

    //! @brief write int to ostream
    //! @param INT the int to be written
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return reference to OSTREAM
    std::ostream &Serialize::Write
    (
      const int &INT,
      std::ostream &OSTREAM,
      const size_t INDENT
    )
    {
      return WriteNumericType( INT, OSTREAM, INDENT);
    }

    //! @brief write size_t to ostream
    //! @param SIZE_T the size_t to be written
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return reference to OSTREAM
    std::ostream &Serialize::Write
    (
      const unsigned long &SIZE_T,
      std::ostream &OSTREAM,
      const size_t INDENT
    )
    {
      return WriteNumericType( SIZE_T, OSTREAM, INDENT);
    }

    //! @brief write size_t to ostream
    //! @param SIZE_T the size_t to be written
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return reference to OSTREAM
    std::ostream &Serialize::Write
    (
      const unsigned int &SIZE_T,
      std::ostream &OSTREAM,
      const size_t INDENT
    )
    {
      return WriteNumericType( SIZE_T, OSTREAM, INDENT);
    }

    //! @brief write size_t to ostream
    //! @param SIZE_T the size_t to be written
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return reference to OSTREAM
    std::ostream &Serialize::Write
    (
      const unsigned short &SIZE_T,
      std::ostream &OSTREAM,
      const size_t INDENT
    )
    {
      return WriteNumericType( SIZE_T, OSTREAM, INDENT);
    }

    //! @brief write long to ostream
    //! @param LONG the long to be written
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return reference to OSTREAM
    std::ostream &Serialize::Write
    (
      const long &LONG,
      std::ostream &OSTREAM,
      const size_t INDENT
    )
    {
      return WriteNumericType( LONG, OSTREAM, INDENT);
    }

    //! @brief write long long int to ostream
    //! @param LONG the long long int to be written
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return reference to OSTREAM
    std::ostream &Serialize::Write
    (
      const long long int &LONG,
      std::ostream &OSTREAM,
      const size_t INDENT
    )
    {
      return WriteNumericType( LONG, OSTREAM, INDENT);
    }

    //! @brief write long long unsigned int to ostream
    //! @param LONG the long long unsigned int to be written
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return reference to OSTREAM
    std::ostream &Serialize::Write
    (
      const long long unsigned int &LONG,
      std::ostream &OSTREAM,
      const size_t INDENT
    )
    {
      return WriteNumericType( LONG, OSTREAM, INDENT);
    }

    //! @brief write char to ostream
    //! @param CHAR the char to be written
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return reference to OSTREAM
    std::ostream &Serialize::Write
    (
      const char &CHAR,
      std::ostream &OSTREAM,
      const size_t INDENT
    )
    {
      // insert indents
      InsertIndent( OSTREAM, INDENT);

      // space characters must be escaped
      if( isspace( CHAR))
      {
        OSTREAM << '\'' << CHAR << '\'';
      }
      else if( CHAR == '\'' || CHAR == '\\') // single quotes and escapes must also be escaped
      {
        OSTREAM << '\\' << CHAR;
      }
      else
      {
        // write char directly
        OSTREAM << CHAR;
      }

      // end
      return OSTREAM;
    }

    //! @brief write char to ostream
    //! @param CHAR the char to be written
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return reference to OSTREAM
    std::ostream &Serialize::Write
    (
      const unsigned char &CHAR,
      std::ostream &OSTREAM,
      const size_t INDENT
    )
    {
      // insert indents
      InsertIndent( OSTREAM, INDENT);

      // put CHAR between ' '
      OSTREAM << ( unsigned int)( CHAR);

      // end
      return OSTREAM;
    }

    //! @brief write string to ostream
    //! @param STRING the string to be written
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return reference to OSTREAM
    std::ostream &Serialize::Write
    (
      const std::string &STRING,
      std::ostream &OSTREAM,
      const size_t INDENT
    )
    {
      // insert indents
      InsertIndent( OSTREAM, INDENT);

      // write the string
      OSTREAM << "\"";

      // iterate over string and escape all '"' and '\'
      for( std::string::const_iterator itr( STRING.begin()), itr_end( STRING.end()); itr != itr_end; ++itr)
      {
        if( *itr == '"' || *itr == '\\')
        {
          OSTREAM.put( '\\');
        }
        OSTREAM.put( *itr);
      }

      OSTREAM << "\"";

      // end
      return OSTREAM;
    }

    //! @brief write Wrapper< std::string> to ostream
    //! @param WRAPPER_STRING the data to be written
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return reference to OSTREAM
    std::ostream &Serialize::Write
    (
      const util::Wrapper< std::string> &WRAPPER_STRING,
      std::ostream &OSTREAM,
      const size_t INDENT
    )
    {
      // write the object
      WRAPPER_STRING.WriteObject( OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  /////////////////////
  // read operations //
  /////////////////////

    //! @brief read OBJECT from istream
    //! @param OBJECT util::ObjectInterface derived class
    //! @param ISTREAM input stream to read from
    //! @return reference to ISTREAM
    std::istream &Serialize::Read
    (
      util::ObjectInterface &OBJECT,
      std::istream &ISTREAM
    )
    {
      // read OBJECT
      return OBJECT.ReadObject( ISTREAM);
    }

    //! @brief read OBJECT from istream
    //! @param OBJECT util::ObjectInterface derived class
    //! @param ISTREAM input stream to read from
    //! @return reference to ISTREAM
    std::istream &Serialize::Read
    (
      util::SerializableInterface &OBJECT,
      std::istream &ISTREAM
    )
    {
      // read OBJECT
      return OBJECT.Read( ISTREAM);
    }

    //! @brief read double from istream
    //! @param DOUBLE the double to be read
    //! @param ISTREAM input stream to read from
    //! @return reference to ISTREAM
    std::istream &Serialize::Read
    (
      double &DOUBLE,
      std::istream &ISTREAM
    )
    {
      DOUBLE = util::ReadNumberFromStream< double>( ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief read float from istream
    //! @param FLOAT the float to be read
    //! @param ISTREAM input stream to read from
    //! @return reference to ISTREAM
    std::istream &Serialize::Read
    (
      float &FLOAT,
      std::istream &ISTREAM
    )
    {
      FLOAT = util::ReadNumberFromStream< float>( ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief read boolean from istream
    //! @param BOOLEAN the boolean to be read
    //! @param ISTREAM input stream to read from
    //! @return reference to ISTREAM
    std::istream &Serialize::Read
    (
      bool &BOOLEAN,
      std::istream &ISTREAM
    )
    {
      // read bool
      ISTREAM >> BOOLEAN;

      // end
      return ISTREAM;
    }

    //! @brief read short from istream
    //! @param SHORT the short to be read
    //! @param ISTREAM input stream to read from
    //! @return reference to ISTREAM
    std::istream &Serialize::Read
    (
      short &SHORT,
      std::istream &ISTREAM
    )
    {
      SHORT = util::ReadNumberFromStream< short>( ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief read int from istream
    //! @param INT the int to be read
    //! @param ISTREAM input stream to read from
    //! @return reference to ISTREAM
    std::istream &Serialize::Read
    (
      int &INT,
      std::istream &ISTREAM
    )
    {
      INT = util::ReadNumberFromStream< int>( ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief read long int from istream
    //! @param LONG_INT the long int to be read
    //! @param ISTREAM input stream to read from
    //! @return reference to ISTREAM
    std::istream &Serialize::Read
    (
      long int &LONG_INT,
      std::istream &ISTREAM
    )
    {
      LONG_INT = util::ReadNumberFromStream< long int>( ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief read long long int from istream
    //! @param LONG_LONG_INT the long long int to be read
    //! @param ISTREAM input stream to read from
    //! @return reference to ISTREAM
    std::istream &Serialize::Read
    (
      long long int &LONG_LONG_INT,
      std::istream &ISTREAM
    )
    {
      LONG_LONG_INT = util::ReadNumberFromStream< long long int>( ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief read unsigned short from istream
    //! @param U_SHORT the short to be read
    //! @param ISTREAM input stream to read from
    //! @return reference to ISTREAM
    std::istream &Serialize::Read
    (
      unsigned short &U_SHORT,
      std::istream &ISTREAM
    )
    {
      U_SHORT = util::ReadNumberFromStream< unsigned short>( ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief read unsigned int from istream
    //! @param U_INT the unsigned int to be read
    //! @param ISTREAM input stream to read from
    //! @return reference to ISTREAM
    std::istream &Serialize::Read
    (
      unsigned int &U_INT,
      std::istream &ISTREAM
    )
    {
      U_INT = util::ReadNumberFromStream< unsigned int>( ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief read long unsigned int from istream
    //! @param U_LONG_INT the long unsigned int to be read
    //! @param ISTREAM input stream to read from
    //! @return reference to ISTREAM
    std::istream &Serialize::Read
    (
      long unsigned int &U_LONG_INT,
      std::istream &ISTREAM
    )
    {
      U_LONG_INT = util::ReadNumberFromStream< long unsigned int>( ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief read size_t from istream
    //! @param U_LONG_LONG_INT the long long unsigned int to be read
    //! @param ISTREAM input stream to read from
    //! @return reference to ISTREAM
    std::istream &Serialize::Read
    (
      long long unsigned int &U_LONG_LONG_INT,
      std::istream &ISTREAM
    )
    {
      U_LONG_LONG_INT = util::ReadNumberFromStream< long long unsigned int>( ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief read char from istream
    //! @param CHAR the char to be read
    //! @param ISTREAM input stream to read from
    //! @return reference to ISTREAM
    std::istream &Serialize::Read
    (
      char &CHAR,
      std::istream &ISTREAM
    )
    {
      // read the next character, skipping WS
      ISTREAM >> CHAR;

      // if the character is a single quote, read the character afterward
      if( CHAR == '\'')
      {
        // and take the next character directly
        ISTREAM.get( CHAR);
        // also get the last character and discard it
        BCL_Assert
        (
          ISTREAM.good() && ISTREAM.get() == '\'',
          "expected end single quote after character"
        );
      }
      else if( CHAR == '\\') // disregard an escape character
      {
        // get the next character directly
        ISTREAM.get( CHAR);
      }

      // end
      return ISTREAM;
    }

    //! @brief read char from istream
    //! @param CHAR the char to be read
    //! @param ISTREAM input stream to read from
    //! @return reference to ISTREAM
    std::istream &Serialize::Read
    (
      unsigned char &CHAR,
      std::istream &ISTREAM
    )
    {
      // temp character
      unsigned int tmp;

      ISTREAM >> tmp;

      CHAR = ( unsigned char)( tmp);

      // end
      return ISTREAM;
    }

    //! @brief read string from istream
    //! @param STRING the string to be read
    //! @param ISTREAM input stream to read from
    //! @return reference to ISTREAM
    std::istream &Serialize::Read
    (
      std::string &STRING,
      std::istream &ISTREAM
    )
    {
      // find leading "
      // temp character
      char tmp;

      BCL_Assert
      (
        ISTREAM >> tmp && tmp == '"',
        "could not find \" indicating a beginning of a string in given stream up to the end but found: \'" +
        util::Format()( tmp) + "\'"
      );

      // read until find " without \ in front of it
      STRING.clear();
      while( ISTREAM.get( tmp))
      {
        // end of string in ISTREAM
        if( tmp == '"')
        {
          break;
        }

        // escape character
        if( tmp == '\\')
        {
          ISTREAM.get( tmp);
        }

        // put char in string
        STRING.push_back( tmp);
      }

      // end
      return ISTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief insert INDENT to stream
    //! @param OSTREAM output stream
    //! @param INDENT number of indentation to be inserted
    //! @return reference to output stream
    std::ostream &Serialize::InsertIndent( std::ostream &OSTREAM, const size_t INDENT)
    {
      // insert indents
      OSTREAM << std::string( INDENT * s_Number_spaces_per_indent, ' ');

      // end
      return OSTREAM;
    }

  } // namespace io
} // namespace bcl
