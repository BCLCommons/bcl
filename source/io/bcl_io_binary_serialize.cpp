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
#include "io/bcl_io_binary_serialize.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace io
  {
    //! @brief helper class used by ConvertToLittleEndian (below)
    //! @tparam N the sizeof the data structure in bytes
    //! @param DATA the data, converted to a character pointer
    template< int N>
    void MakeLittleEndianHelper( unsigned char *DATA);

    //! specializations for base cases of MakeLittleEndianHelper, 0 bytes
    template<>
    void MakeLittleEndianHelper< 0>( unsigned char *DATA);

    //! specializations for base cases of MakeLittleEndianHelper, 1 byte
    template<>
    void MakeLittleEndianHelper< 1>( unsigned char *DATA);

    //! @brief reverse all the bits in a character
    //! @param ORIGINAL_CHAR the original character, whose bits will be reversed
    void ReverseBitsInCharacter( unsigned char &ORIGINAL_CHAR)
    {
      // create an array of chars such that s_bit_reverse_table[ char_1] = char_1 with bits reversed
      static const unsigned char s_bit_reverse_table[ 256] =
      {
         0, 128, 64, 192, 32, 160,  96, 224, 16, 144, 80, 208, 48, 176, 112, 240,
         8, 136, 72, 200, 40, 168, 104, 232, 24, 152, 88, 216, 56, 184, 120, 248,
         4, 132, 68, 196, 36, 164, 100, 228, 20, 148, 84, 212, 52, 180, 116, 244,
        12, 140, 76, 204, 44, 172, 108, 236, 28, 156, 92, 220, 60, 188, 124, 252,
         1, 129, 65, 193, 33, 161,  97, 225, 17, 145, 81, 209, 49, 177, 113, 241,
         9, 137, 73, 201, 41, 169, 105, 233, 25, 153, 89, 217, 57, 185, 121, 249,
         5, 133, 69, 197, 37, 165, 101, 229, 21, 149, 85, 213, 53, 181, 117, 245,
        13, 141, 77, 205, 45, 173, 109, 237, 29, 157, 93, 221, 61, 189, 125, 253,
         2, 130, 66, 194, 34, 162,  98, 226, 18, 146, 82, 210, 50, 178, 114, 242,
        10, 138, 74, 202, 42, 170, 106, 234, 26, 154, 90, 218, 58, 186, 122, 250,
         6, 134, 70, 198, 38, 166, 102, 230, 22, 150, 86, 214, 54, 182, 118, 246,
        14, 142, 78, 206, 46, 174, 110, 238, 30, 158, 94, 222, 62, 190, 126, 254,
         3, 131, 67, 195, 35, 163,  99, 227, 19, 147, 83, 211, 51, 179, 115, 243,
        11, 139, 75, 203, 43, 171, 107, 235, 27, 155, 91, 219, 59, 187, 123, 251,
         7, 135, 71, 199, 39, 167, 103, 231, 23, 151, 87, 215, 55, 183, 119, 247,
        15, 143, 79, 207, 47, 175, 111, 239, 31, 159, 95, 223, 63, 191, 127, 255
      };

      ORIGINAL_CHAR = s_bit_reverse_table[ int( ORIGINAL_CHAR)];
    }

    //! @brief helper class used by ConvertToLittleEndian (below)
    //! @tparam N the sizeof the data structure in bytes
    //! @param DATA the data, converted to a character pointer
    template< int N>
    void MakeLittleEndianHelper( unsigned char *DATA)
    {
      // swap the highest and lowest bit
      std::swap( DATA[ 0], DATA[ N - 1]);

      // reverse the order of the bits in those characters
      ReverseBitsInCharacter( DATA[ 0]);
      ReverseBitsInCharacter( DATA[ N - 1]);

      // recursively call make little endian helper; the compiler can unroll this to avoid any function call overhead
      MakeLittleEndianHelper< N - 2>( DATA + 1);
    }

    template<>
    void MakeLittleEndianHelper< 0>( unsigned char *DATA)
    {
    }

    template<>
    void MakeLittleEndianHelper< 1>( unsigned char *DATA)
    {
      ReverseBitsInCharacter( DATA[ 0]);
    }

    // determine whether the system is big-endian
    // big endian systems start off with the most significant bit, which is never 1, so
    // if the first byte is 0, the system is not big-endian; the size of an address never
    // occupies the highest byte of a size_t, so if the first byte is zero, the system is
    // a big endian
    const bool g_IsBigEndian( reinterpret_cast< const char *>( &g_SizeOfAddress)[ 0] == char( 0));

    //! @brief writes to a binary file; always converts DATA to little endian type
    template< typename t_DataType>
    std::ostream &WriteToBinaryFile( std::ostream &WRITE, const t_DataType &DATA)
    {
      if( g_IsBigEndian) // convert to little-endian format before writing
      {
        t_DataType data( DATA);
        MakeLittleEndianHelper< sizeof( DATA)>( reinterpret_cast< unsigned char *>( &data));
        return WRITE.write( reinterpret_cast< char *>( &data), sizeof( data));
      }

      // little endian format; write it out directly
      return WRITE.write( reinterpret_cast< const char *>( &DATA), sizeof( DATA));
    }

    //! @brief reads from a binary file; always converts DATA to little endian type
    template< typename t_DataType>
    std::istream &ReadFromBinaryFile( std::istream &READ, t_DataType &DATA)
    {
      READ.read( reinterpret_cast< char *>( &DATA), sizeof( DATA));

      if( g_IsBigEndian) // for big-endian systems, reverse the bits
      {
        MakeLittleEndianHelper< sizeof( DATA)>( reinterpret_cast< unsigned char *>( &DATA));
      }
      return READ;
    }

  //////////////////////
  // write operations //
  //////////////////////

    //! @brief write double to ostream
    //! @param DOUBLE the double to be written
    //! @param OSTREAM output stream to write to
    //! @return reference to OSTREAM
    std::ostream &BinarySerialize::Write
    (
      const double &DOUBLE,
      std::ostream &OSTREAM
    )
    {
      return WriteToBinaryFile( OSTREAM, DOUBLE);
    }

    //! @brief write float to ostream
    //! @param FLOAT the float to be written
    //! @param OSTREAM output stream to write to
    //! @return reference to OSTREAM
    std::ostream &BinarySerialize::Write
    (
      const float &FLOAT,
      std::ostream &OSTREAM
    )
    {
      return WriteToBinaryFile( OSTREAM, FLOAT);
    }

    //! @brief write boolean to ostream
    //! @param BOOLEAN the boolean to be written
    //! @param OSTREAM output stream to write to
    //! @return reference to OSTREAM
    std::ostream &BinarySerialize::Write
    (
      const bool &BOOLEAN,
      std::ostream &OSTREAM
    )
    {
      return WriteToBinaryFile( OSTREAM, ( unsigned char)( BOOLEAN));
    }

    //! @brief write short to ostream
    //! @param SHORT the short to be written
    //! @param OSTREAM output stream to write to
    //! @return reference to OSTREAM
    std::ostream &BinarySerialize::Write
    (
      const short &SHORT,
      std::ostream &OSTREAM
    )
    {
      // to ensure cross-platform compatibility of binary files, write shorts out with 32 bits
      return WriteToBinaryFile( OSTREAM, int32_t( SHORT));
    }

    //! @brief write int to ostream
    //! @param INT the int to be written
    //! @param OSTREAM output stream to write to
    //! @return reference to OSTREAM
    std::ostream &BinarySerialize::Write
    (
      const int &INT,
      std::ostream &OSTREAM
    )
    {
      // to ensure cross-platform compatibility of binary files, write ints out with 64 bits
      return WriteToBinaryFile( OSTREAM, int64_t( INT));
    }

    //! @brief write size_t to ostream
    //! @param SIZE_T the size_t to be written
    //! @param OSTREAM output stream to write to
    //! @return reference to OSTREAM
    std::ostream &BinarySerialize::Write
    (
      const unsigned long &SIZE_T,
      std::ostream &OSTREAM
    )
    {
      // to ensure cross-platform compatibility of binary files, write longs out with 64 bits
      return WriteToBinaryFile( OSTREAM, uint64_t( SIZE_T));
    }

    //! @brief write size_t to ostream
    //! @param SIZE_T the size_t to be written
    //! @param OSTREAM output stream to write to
    //! @return reference to OSTREAM
    std::ostream &BinarySerialize::Write
    (
      const unsigned int &SIZE_T,
      std::ostream &OSTREAM
    )
    {
      // to ensure cross-platform compatibility of binary files, write ints out with 64 bits
      return WriteToBinaryFile( OSTREAM, uint64_t( SIZE_T));
    }

    //! @brief write size_t to ostream
    //! @param SIZE_T the size_t to be written
    //! @param OSTREAM output stream to write to
    //! @return reference to OSTREAM
    std::ostream &BinarySerialize::Write
    (
      const unsigned short &SIZE_T,
      std::ostream &OSTREAM
    )
    {
      // to ensure cross-platform compatibility of binary files, write shorts out with 32 bits
      return WriteToBinaryFile( OSTREAM, uint32_t( SIZE_T));
    }

    //! @brief write long to ostream
    //! @param LONG the long to be written
    //! @param OSTREAM output stream to write to
    //! @return reference to OSTREAM
    std::ostream &BinarySerialize::Write
    (
      const long &LONG,
      std::ostream &OSTREAM
    )
    {
      return WriteToBinaryFile( OSTREAM, uint64_t( LONG));
    }

    //! @brief write char to ostream
    //! @param CHAR the char to be written
    //! @param OSTREAM output stream to write to
    //! @return reference to OSTREAM
    std::ostream &BinarySerialize::Write
    (
      const char &CHAR,
      std::ostream &OSTREAM
    )
    {
      return WriteToBinaryFile( OSTREAM, CHAR);
    }

    //! @brief write char to ostream
    //! @param CHAR the char to be written
    //! @param OSTREAM output stream to write to
    //! @return reference to OSTREAM
    std::ostream &BinarySerialize::Write
    (
      const unsigned char &CHAR,
      std::ostream &OSTREAM
    )
    {
      return WriteToBinaryFile( OSTREAM, CHAR);
    }

    //! @brief write string to ostream
    //! @param STRING the string to be written
    //! @param OSTREAM output stream to write to
    //! @return reference to OSTREAM
    std::ostream &BinarySerialize::Write
    (
      const std::string &STRING,
      std::ostream &OSTREAM
    )
    {
      // write the size
      WriteToBinaryFile( OSTREAM, uint64_t( STRING.size()));
      // write the characters
      for( std::string::const_iterator itr( STRING.begin()), itr_end( STRING.end()); itr != itr_end; ++itr)
      {
        WriteToBinaryFile( OSTREAM, *itr);
      }

      // end
      return OSTREAM;
    }

  /////////////////////
  // read operations //
  /////////////////////

    //! @brief read double from istream
    //! @param DOUBLE the double to be read
    //! @param ISTREAM input stream to read from
    //! @return reference to ISTREAM
    std::istream &BinarySerialize::Read
    (
      double &DOUBLE,
      std::istream &ISTREAM
    )
    {
      return ReadFromBinaryFile( ISTREAM, DOUBLE);
    }

    //! @brief read float from istream
    //! @param FLOAT the float to be read
    //! @param ISTREAM input stream to read from
    //! @return reference to ISTREAM
    std::istream &BinarySerialize::Read
    (
      float &FLOAT,
      std::istream &ISTREAM
    )
    {
      return ReadFromBinaryFile( ISTREAM, FLOAT);
    }

    //! @brief read boolean from istream
    //! @param BOOLEAN the boolean to be read
    //! @param ISTREAM input stream to read from
    //! @return reference to ISTREAM
    std::istream &BinarySerialize::Read
    (
      bool &BOOLEAN,
      std::istream &ISTREAM
    )
    {
      // read bool back into a single byte;
      unsigned char tmp;
      ReadFromBinaryFile( ISTREAM, tmp);
      BOOLEAN = bool( tmp);
      return ISTREAM;
    }

    //! @brief read short from istream
    //! @param SHORT the short to be read
    //! @param ISTREAM input stream to read from
    //! @return reference to ISTREAM
    std::istream &BinarySerialize::Read
    (
      short &SHORT,
      std::istream &ISTREAM
    )
    {
      // read back into a 32 bit value
      int32_t tmp;
      ReadFromBinaryFile( ISTREAM, tmp);
      SHORT = short( tmp);
      // end
      return ISTREAM;
    }

    //! @brief read int from istream
    //! @param INT the int to be read
    //! @param ISTREAM input stream to read from
    //! @return reference to ISTREAM
    std::istream &BinarySerialize::Read
    (
      int &INT,
      std::istream &ISTREAM
    )
    {
      int64_t tmp;
      ReadFromBinaryFile( ISTREAM, tmp);
      INT = int( tmp);

      // end
      return ISTREAM;
    }

    //! @brief read long int from istream
    //! @param LONG_INT the long int to be read
    //! @param ISTREAM input stream to read from
    //! @return reference to ISTREAM
    std::istream &BinarySerialize::Read
    (
      long int &LONG_INT,
      std::istream &ISTREAM
    )
    {
      int64_t tmp;
      ReadFromBinaryFile( ISTREAM, tmp);
      LONG_INT = long( tmp);

      // end
      return ISTREAM;
    }

    //! @brief read unsigned short from istream
    //! @param U_SHORT the short to be read
    //! @param ISTREAM input stream to read from
    //! @return reference to ISTREAM
    std::istream &BinarySerialize::Read
    (
      unsigned short &U_SHORT,
      std::istream &ISTREAM
    )
    {
      // read back into a 32 bit value
      uint32_t tmp;
      ReadFromBinaryFile( ISTREAM, tmp);
      U_SHORT = ( unsigned short)( tmp);
      // end
      return ISTREAM;
    }

    //! @brief read unsigned int from istream
    //! @param U_INT the unsigned int to be read
    //! @param ISTREAM input stream to read from
    //! @return reference to ISTREAM
    std::istream &BinarySerialize::Read
    (
      unsigned int &U_INT,
      std::istream &ISTREAM
    )
    {
      uint64_t tmp;
      ReadFromBinaryFile( ISTREAM, tmp);
      U_INT = ( unsigned int)( tmp);

      // end
      return ISTREAM;
    }

    //! @brief read long unsigned int from istream
    //! @param U_LONG_INT the long unsigned int to be read
    //! @param ISTREAM input stream to read from
    //! @return reference to ISTREAM
    std::istream &BinarySerialize::Read
    (
      long unsigned int &U_LONG_INT,
      std::istream &ISTREAM
    )
    {
      uint64_t tmp;
      ReadFromBinaryFile( ISTREAM, tmp);
      U_LONG_INT = ( unsigned long)( tmp);

      // end
      return ISTREAM;
    }

    //! @brief read char from istream
    //! @param CHAR the char to be read
    //! @param ISTREAM input stream to read from
    //! @return reference to ISTREAM
    std::istream &BinarySerialize::Read
    (
      char &CHAR,
      std::istream &ISTREAM
    )
    {
      return ReadFromBinaryFile( ISTREAM, CHAR);
    }

    //! @brief read char from istream
    //! @param CHAR the char to be read
    //! @param ISTREAM input stream to read from
    //! @return reference to ISTREAM
    std::istream &BinarySerialize::Read
    (
      unsigned char &CHAR,
      std::istream &ISTREAM
    )
    {
      return ReadFromBinaryFile( ISTREAM, CHAR);
    }

    //! @brief read string from istream
    //! @param STRING the string to be read
    //! @param ISTREAM input stream to read from
    //! @return reference to ISTREAM
    std::istream &BinarySerialize::Read
    (
      std::string &STRING,
      std::istream &ISTREAM
    )
    {
      // read the size
      uint64_t size;
      ReadFromBinaryFile( ISTREAM, size);
      STRING.clear();
      STRING.reserve( size);
      for( size_t i( 0); i < size; ++i)
      {
        // read the character
        char tmp;
        ReadFromBinaryFile( ISTREAM, tmp);

        // append it to the string
        STRING += tmp;
      }

      // end
      return ISTREAM;
    }

    //! @brief read from istream into character buffer (invokes fewer calls to read than calling read on each item individually)
    //! @param BUFFER the vector to be read into, should already be allocated to the correct size
    //! @param ISTREAM input stream to read from
    std::istream &BinarySerialize::ReadVector( linal::Vector< char> &BUFFER, std::istream &ISTREAM)
    {
      if( BUFFER.GetSize() == size_t( 0))
      {
        // nothing to read, return
        return ISTREAM;
      }
      const size_t size( BUFFER.GetSize());
      ISTREAM.read( BUFFER.Begin(), size);
      if( g_IsBigEndian)
      {
        unsigned char *buffer( reinterpret_cast< unsigned char *>( BUFFER.Begin()));
        for( unsigned char *const buffer_end( buffer + size); buffer < buffer_end; ++buffer)
        {
          ReverseBitsInCharacter( *buffer);
        }
      }
      return ISTREAM;
    }

    //! @brief read from istream into vector (invokes fewer calls to read than calling read on each item individually)
    //! @param BUFFER the vector to be read into, should already be allocated to the correct size
    //! @param ISTREAM input stream to read from
    std::istream &BinarySerialize::ReadVector( linal::Vector< float> &BUFFER, std::istream &ISTREAM)
    {
      if( BUFFER.GetSize() == size_t( 0))
      {
        // nothing to read, return
        return ISTREAM;
      }
      const size_t size( BUFFER.GetSize() * sizeof( float));
      ISTREAM.read( reinterpret_cast< char *const>( BUFFER.Begin()), size);
      if( g_IsBigEndian)
      {
        unsigned char *buffer( reinterpret_cast< unsigned char *const>( BUFFER.Begin()));
        for( unsigned char *const buffer_end( buffer + size); buffer < buffer_end; buffer += sizeof( float))
        {
          MakeLittleEndianHelper< sizeof( float)>( buffer);
        }
      }
      return ISTREAM;
    }

  } // namespace io
} // namespace bcl
