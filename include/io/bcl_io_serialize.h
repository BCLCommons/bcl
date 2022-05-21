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

#ifndef BCL_IO_SERIALIZE_H_
#define BCL_IO_SERIALIZE_H_

// include the namespace header
#include "bcl_io.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_assert.h"
#include "util/bcl_util_class_descriptor.h"
#include "util/bcl_util_format.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace io
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Serialize
    //! @brief TODO: add a brief comment
    //!
    //! @see @link example_io_serialize.cpp @endlink
    //! @author woetzen
    //! @date Aug 15, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Serialize
    {

    public:

    //////////
    // data //
    //////////

      //! number of space is standard indentation
      static const size_t s_Number_spaces_per_indent = 2;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    private:

      //! @brief default constructor
      Serialize()
      {
      }

    public:

    //////////////////////
    // write operations //
    //////////////////////

      //! @brief write OBJECT to ostream
      //! @param OBJECT ObjectInterface derived class
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return reference to OSTREAM
      static std::ostream &Write
      (
        const util::ObjectInterface &OBJECT,
        std::ostream &OSTREAM,
        const size_t INDENT
      );

      //! @brief write OBJECT to ostream
      //! use this function, if there is ambiguity when writing std::string or other bcl objects derived from std class
      //! @param OBJECT util::ObjectInterface derived class
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return reference to OSTREAM
      static std::ostream &WriteObject
      (
        const util::ObjectInterface &OBJECT,
        std::ostream &OSTREAM,
        const size_t INDENT
      );

      //! @brief write OBJECT to ostream
      //! @param OBJECT ObjectInterface derived class
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return reference to OSTREAM
      static std::ostream &Write
      (
        const util::SerializableInterface &OBJECT,
        std::ostream &OSTREAM,
        const size_t INDENT
      );

      //! @brief write OBJECT to ostream
      //! use this function, if there is ambiguity when writing std::string or other bcl objects derived from std class
      //! @param OBJECT util::ObjectInterface derived class
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return reference to OSTREAM
      static std::ostream &WriteSerializable
      (
        const util::SerializableInterface &OBJECT,
        std::ostream &OSTREAM,
        const size_t INDENT
      );

      //! @brief write double to ostream
      //! @param DOUBLE the double to be written
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return reference to OSTREAM
      static std::ostream &Write
      (
        const double &DOUBLE,
        std::ostream &OSTREAM,
        const size_t INDENT
      );

      //! @brief write float to ostream
      //! @param FLOAT the float to be written
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return reference to OSTREAM
      static std::ostream &Write
      (
        const float &FLOAT,
        std::ostream &OSTREAM,
        const size_t INDENT
      );

      //! @brief write boolean to ostream
      //! @param BOOLEAN the boolean to be written
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return reference to OSTREAM
      static std::ostream &Write
      (
        const bool &BOOLEAN,
        std::ostream &OSTREAM,
        const size_t INDENT
      );

      //! @brief write short to ostream
      //! @param SHORT the short to be written
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return reference to OSTREAM
      static std::ostream &Write
      (
        const short &SHORT,
        std::ostream &OSTREAM,
        const size_t INDENT
      );

      //! @brief write int to ostream
      //! @param INT the int to be written
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return reference to OSTREAM
      static std::ostream &Write
      (
        const int &INT,
        std::ostream &OSTREAM,
        const size_t INDENT
      );

      //! @brief write size_t to ostream
      //! @param SIZE_T the size_t to be written
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return reference to OSTREAM
      static std::ostream &Write
      (
        const unsigned long &SIZE_T,
        std::ostream &OSTREAM,
        const size_t INDENT
      );

      //! @brief write size_t to ostream
      //! @param SIZE_T the size_t to be written
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return reference to OSTREAM
      static std::ostream &Write
      (
        const unsigned int &SIZE_T,
        std::ostream &OSTREAM,
        const size_t INDENT
      );

      //! @brief write size_t to ostream
      //! @param SIZE_T the size_t to be written
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return reference to OSTREAM
      static std::ostream &Write
      (
        const unsigned short &SIZE_T,
        std::ostream &OSTREAM,
        const size_t INDENT
      );

      //! @brief write long to ostream
      //! @param LONG the long to be written
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return reference to OSTREAM
      static std::ostream &Write
      (
        const long &LONG,
        std::ostream &OSTREAM,
        const size_t INDENT
      );

      //! @brief write long long int to ostream
      //! @param LONG the long long int to be written
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return reference to OSTREAM
      static std::ostream &Write
      (
        const long long int &LONG,
        std::ostream &OSTREAM,
        const size_t INDENT
      );

      //! @brief write long long unsigned int to ostream
      //! @param LONG the long long unsigned int to be written
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return reference to OSTREAM
      static std::ostream &Write
      (
        const long long unsigned int &LONG,
        std::ostream &OSTREAM,
        const size_t INDENT
      );

      //! @brief write char to ostream
      //! @param CHAR the char to be written
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return reference to OSTREAM
      static std::ostream &Write
      (
        const char &CHAR,
        std::ostream &OSTREAM,
        const size_t INDENT
      );

      //! @brief write char to ostream
      //! @param CHAR the char to be written
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return reference to OSTREAM
      static std::ostream &Write
      (
        const unsigned char &CHAR,
        std::ostream &OSTREAM,
        const size_t INDENT
      );

      //! @brief write string to ostream
      //! @param STRING the string to be written
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return reference to OSTREAM
      static std::ostream &Write
      (
        const std::string &STRING,
        std::ostream &OSTREAM,
        const size_t INDENT
      );

      //! @brief write any data without indent to ostream
      //! @param DATA the data to be written
      //! @param OSTREAM output stream to write to
      //! @return reference to OSTREAM
      template< typename t_DataType>
      static std::ostream &Write
      (
        const t_DataType &DATA,
        std::ostream &OSTREAM
      )
      {
        // call according write function with indent 0
        return Write( DATA, OSTREAM, 0);
      }

      //! @brief write data to ostream
      //! @param PAIR the std::pair to be written
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return reference to OSTREAM
      template< typename t_First, typename t_Second>
      static std::ostream &Write
      (
        const std::pair< t_First, t_Second> &PAIR,
        std::ostream &OSTREAM,
        const size_t INDENT
      )
      {
        // insert indents
        InsertIndent( OSTREAM, INDENT) << GetStaticClassName( PAIR) << '\n';

        // write data
        Serialize::Write( PAIR.first, OSTREAM, INDENT + 1) << '\n';
        Serialize::Write( PAIR.second, OSTREAM, INDENT + 1);

        // end
        return OSTREAM;
      }

      //! @brief write data to ostream using a formatting
      //! @param DATA the data to be written
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @param FORMAT use given format
      //! @return reference to OSTREAM
      template< typename t_DataType>
      static std::ostream &Write
      (
        const t_DataType &DATA,
        std::ostream &OSTREAM,
        const size_t INDENT,
        const util::Format &FORMAT
      )
      {
        // format given OSTREAM
        FORMAT.SetFlags( OSTREAM);

        // write data
        Serialize::Write( DATA, OSTREAM, INDENT);

        // remove formatting on OSTREAM
        FORMAT.UnsetFlags( OSTREAM);

        // end
        return OSTREAM;
      }

      //! @brief write Wrapper< std::string> to ostream
      //! @param WRAPPER_STRING the data to be written
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return reference to OSTREAM
      static std::ostream &Write
      (
        const util::Wrapper< std::string> &WRAPPER_STRING,
        std::ostream &OSTREAM,
        const size_t INDENT
      );

    /////////////////////
    // read operations //
    /////////////////////

      //! @brief read OBJECT from istream
      //! @param OBJECT util::ObjectInterface derived class
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      static std::istream &Read
      (
        util::ObjectInterface &OBJECT,
        std::istream &ISTREAM
      );

      //! @brief read OBJECT from istream
      //! @param OBJECT util::ObjectInterface derived class
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      static std::istream &Read
      (
        util::SerializableInterface &OBJECT,
        std::istream &ISTREAM
      );

      //! @brief read double from istream
      //! @param DOUBLE the double to be read
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      static std::istream &Read
      (
        double &DOUBLE,
        std::istream &ISTREAM
      );

      //! @brief read float from istream
      //! @param FLOAT the float to be read
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      static std::istream &Read
      (
        float &FLOAT,
        std::istream &ISTREAM
      );

      //! @brief read boolean from istream
      //! @param BOOLEAN the boolean to be read
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      static std::istream &Read
      (
        bool &BOOLEAN,
        std::istream &ISTREAM
      );

      //! @brief read short from istream
      //! @param SHORT the short to be read
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      static std::istream &Read
      (
        short &SHORT,
        std::istream &ISTREAM
      );

      //! @brief read int from istream
      //! @param INT the int to be read
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      static std::istream &Read
      (
        int &INT,
        std::istream &ISTREAM
      );

      //! @brief read long int from istream
      //! @param LONG_INT the long int to be read
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      static std::istream &Read
      (
        long int &LONG_INT,
        std::istream &ISTREAM
      );

      //! @brief read long long int from istream
      //! @param LONG_LONG_INT the long long int to be read
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      static std::istream &Read
      (
        long long int &LONG_LONG_INT,
        std::istream &ISTREAM
      );

      //! @brief read unsigned short from istream
      //! @param U_SHORT the short to be read
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      static std::istream &Read
      (
        unsigned short &U_SHORT,
        std::istream &ISTREAM
      );

      //! @brief read unsigned int from istream
      //! @param U_INT the unsigned int to be read
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      static std::istream &Read
      (
        unsigned int &U_INT,
        std::istream &ISTREAM
      );

      //! @brief read long unsigned int from istream
      //! @param U_LONG_INT the long unsigned int to be read
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      static std::istream &Read
      (
        long unsigned int &U_LONG_INT,
        std::istream &ISTREAM
      );

      //! @brief read size_t from istream
      //! @param U_LONG_LONG_INT the long long unsigned int to be read
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      static std::istream &Read
      (
        long long unsigned int &U_LONG_LONG_INT,
        std::istream &ISTREAM
      );

      //! @brief read char from istream
      //! @param CHAR the char to be read
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      static std::istream &Read
      (
        char &CHAR,
        std::istream &ISTREAM
      );

      //! @brief read char from istream
      //! @param CHAR the char to be read
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      static std::istream &Read
      (
        unsigned char &CHAR,
        std::istream &ISTREAM
      );

      //! @brief read string from istream
      //! @param STRING the string to be read
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      static std::istream &Read
      (
        std::string &STRING,
        std::istream &ISTREAM
      );

      //! @brief read pair from istream
      //! @param PAIR the std::pair to be read
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      template< typename t_First, typename t_Second>
      static std::istream &Read
      (
        std::pair< t_First, t_Second> &PAIR,
        std::istream &ISTREAM
      )
      {
        // read identifier
        std::string identifier;
        ISTREAM >> identifier;

        BCL_Assert
        (
          identifier == GetStaticClassName( PAIR),
          "identifier was not " + GetStaticClassName( PAIR) + " but " + identifier
        );

        // read data
        Serialize::Read( PAIR.first, ISTREAM);
        Serialize::Read( PAIR.second, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief read other data from istream
      //! @param DATA the data to be read
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      template< typename t_DataType>
      static std::istream &Read
      (
        t_DataType &DATA,
        std::istream &ISTREAM
      )
      {
        // read data
        ISTREAM >> DATA;

        // end
        return ISTREAM;
      }

      //! @brief read other data from istream
      //! @param DATA the data to be read
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      template< typename t_DataType>
      static std::istream &Read
      (
        const t_DataType &DATA,
        std::istream &ISTREAM
      )
      {
        BCL_Exit( "reading of const DATA is not possible!", -1);

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
      static std::ostream &InsertIndent
      (
        std::ostream &OSTREAM,
        const size_t INDENT
      );

    }; // class Serialize

  } // namespace io
} // namespace bcl

#endif // BCL_IO_SERIALIZE_H_
