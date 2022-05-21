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

#ifndef BCL_IO_BINARY_SERIALIZE_H_
#define BCL_IO_BINARY_SERIALIZE_H_

// include the namespace header
#include "bcl_io.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_assert.h"
#include "util/bcl_util_format.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace io
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BinarySerialize
    //! @brief functions to read and write basic data types in binary format, always as little endian
    //!
    //! @see @link example_io_binary_serialize.cpp @endlink
    //! @author mendenjl
    //! @date Feb 02, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API BinarySerialize
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    private:

      //! @brief default constructor
      BinarySerialize();

    public:

    //////////////////////
    // write operations //
    //////////////////////

      //! @brief write double to ostream
      //! @param DOUBLE the double to be written
      //! @param OSTREAM output stream to write to
      //! @return reference to OSTREAM
      static std::ostream &Write( const double &DOUBLE, std::ostream &OSTREAM);

      //! @brief write float to ostream
      //! @param FLOAT the float to be written
      //! @param OSTREAM output stream to write to
      //! @return reference to OSTREAM
      static std::ostream &Write( const float &FLOAT, std::ostream &OSTREAM);

      //! @brief write boolean to ostream
      //! @param BOOLEAN the boolean to be written
      //! @param OSTREAM output stream to write to
      //! @return reference to OSTREAM
      static std::ostream &Write( const bool &BOOLEAN, std::ostream &OSTREAM);

      //! @brief write short to ostream
      //! @param SHORT the short to be written
      //! @param OSTREAM output stream to write to
      //! @return reference to OSTREAM
      static std::ostream &Write( const short &SHORT, std::ostream &OSTREAM);

      //! @brief write int to ostream
      //! @param INT the int to be written
      //! @param OSTREAM output stream to write to
      //! @return reference to OSTREAM
      static std::ostream &Write( const int &INT, std::ostream &OSTREAM);

      //! @brief write size_t to ostream
      //! @param SIZE_T the size_t to be written
      //! @param OSTREAM output stream to write to
      //! @return reference to OSTREAM
      static std::ostream &Write( const unsigned long &SIZE_T, std::ostream &OSTREAM);

      //! @brief write size_t to ostream
      //! @param SIZE_T the size_t to be written
      //! @param OSTREAM output stream to write to
      //! @return reference to OSTREAM
      static std::ostream &Write( const unsigned int &SIZE_T, std::ostream &OSTREAM);

      //! @brief write size_t to ostream
      //! @param SIZE_T the size_t to be written
      //! @param OSTREAM output stream to write to
      //! @return reference to OSTREAM
      static std::ostream &Write( const unsigned short &SIZE_T, std::ostream &OSTREAM);

      //! @brief write long to ostream
      //! @param LONG the long to be written
      //! @param OSTREAM output stream to write to
      //! @return reference to OSTREAM
      static std::ostream &Write( const long &LONG, std::ostream &OSTREAM);

      //! @brief write char to ostream
      //! @param CHAR the char to be written
      //! @param OSTREAM output stream to write to
      //! @return reference to OSTREAM
      static std::ostream &Write( const char &CHAR, std::ostream &OSTREAM);

      //! @brief write char to ostream
      //! @param CHAR the char to be written
      //! @param OSTREAM output stream to write to
      //! @return reference to OSTREAM
      static std::ostream &Write( const unsigned char &CHAR, std::ostream &OSTREAM);

      //! @brief write string to ostream
      //! @param STRING the string to be written
      //! @param OSTREAM output stream to write to
      //! @return reference to OSTREAM
      static std::ostream &Write( const std::string &STRING, std::ostream &OSTREAM);

    /////////////////////
    // read operations //
    /////////////////////

      //! @brief read double from istream
      //! @param DOUBLE the double to be read
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      static std::istream &Read( double &DOUBLE, std::istream &ISTREAM);

      //! @brief read float from istream
      //! @param FLOAT the float to be read
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      static std::istream &Read( float &FLOAT, std::istream &ISTREAM);

      //! @brief read boolean from istream
      //! @param BOOLEAN the boolean to be read
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      static std::istream &Read( bool &BOOLEAN, std::istream &ISTREAM);

      //! @brief read short from istream
      //! @param SHORT the short to be read
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      static std::istream &Read( short &SHORT, std::istream &ISTREAM);

      //! @brief read int from istream
      //! @param INT the int to be read
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      static std::istream &Read( int &INT, std::istream &ISTREAM);

      //! @brief read long int from istream
      //! @param LONG_INT the long int to be read
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      static std::istream &Read( long int &LONG_INT, std::istream &ISTREAM);

      //! @brief read unsigned short from istream
      //! @param U_SHORT the short to be read
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      static std::istream &Read( unsigned short &U_SHORT, std::istream &ISTREAM);

      //! @brief read unsigned int from istream
      //! @param U_INT the unsigned int to be read
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      static std::istream &Read( unsigned int &U_INT, std::istream &ISTREAM);

      //! @brief read long unsigned int from istream
      //! @param U_LONG_INT the long unsigned int to be read
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      static std::istream &Read( long unsigned int &U_LONG_INT, std::istream &ISTREAM);

      //! @brief read char from istream
      //! @param CHAR the char to be read
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      static std::istream &Read( char &CHAR, std::istream &ISTREAM);

      //! @brief read char from istream
      //! @param CHAR the char to be read
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      static std::istream &Read( unsigned char &CHAR, std::istream &ISTREAM);

      //! @brief read string from istream
      //! @param STRING the string to be read
      //! @param ISTREAM input stream to read from
      //! @return reference to ISTREAM
      static std::istream &Read( std::string &STRING, std::istream &ISTREAM);

      //! @brief read from istream into character buffer (invokes fewer calls to read than calling read on each item individually)
      //! @param BUFFER the vector to be read into, should already be allocated to the correct size
      //! @param ISTREAM input stream to read from
      static std::istream &ReadVector( linal::Vector< char> &BUFFER, std::istream &ISTREAM);

      //! @brief read from istream into vector (invokes fewer calls to read than calling read on each item individually)
      //! @param BUFFER the vector to be read into, should already be allocated to the correct size
      //! @param ISTREAM input stream to read from
      static std::istream &ReadVector( linal::Vector< float> &BUFFER, std::istream &ISTREAM);

    }; // class BinarySerialize

  } // namespace io
} // namespace bcl

#endif // BCL_IO_BINARY_SERIALIZE_H_
