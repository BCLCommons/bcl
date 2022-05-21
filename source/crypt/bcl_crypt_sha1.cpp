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
#include "crypt/bcl_crypt_sha1.h"
#include "util/bcl_util_object_instances.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace crypt
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> Sha1::s_Instance
    (
      GetObjectInstances().AddInstance( new Sha1())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new Sha1
    Sha1 *Sha1::Clone() const
    {
      return new Sha1( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Sha1::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief default constructor
    Sha1::Sha1State::Sha1State( const std::string &STRING)
    {
      // initial values for the message digest, which is the result
      m_MessageDigest[ 0] = 0x67452301;
      m_MessageDigest[ 1] = 0xEFCDAB89;
      m_MessageDigest[ 2] = 0x98BADCFE;
      m_MessageDigest[ 3] = 0x10325476;
      m_MessageDigest[ 4] = 0xC3D2E1F0;

      // create an array to hold the block
      unsigned char message_block[ e_BlockSize];

      const size_t string_size( STRING.size());

      // To allow for binary-message hashing, SHA-1 specifies that the length of the string must be stored in bits,
      // and that this length must fit in 8 bytes, thus, the maximum number is 2 ^ 61 ( e.g. 1 << 61)
      BCL_Assert
      (
        uint64_t( string_size) < ( UINT64_C( 1) << UINT64_C( 61)),
        "String was too long to be hashed by SHA1!"
      );

      // get the string size and its pointer, converted to unsigned char, which is required by SHA1
      const unsigned char *char_ptr( ( const unsigned char *)STRING.c_str());

      uint32_t w[ e_WordBlockSize]; // Temporary word block used in processing

      size_t message_block_index( 0);
      for( size_t i( 0); i < string_size; ++i)
      {
        message_block[ message_block_index] = char_ptr[ i];
        if( ++message_block_index == e_BlockSize)
        {
          ProcessMessageBlock( message_block, w);
          message_block_index = 0;
        }
      }

      // pad the string
      const size_t prepadded_message_size( e_BlockSize - e_MessagePadding);

      // Check to see if the current message block is too small to hold
      // the initial padding bits and length.  If so, we will pad the
      // block, process it, and then continue padding into a second block.
      message_block[ message_block_index++] = 0x80;
      if( message_block_index > prepadded_message_size)
      {
        for( ; message_block_index < e_BlockSize; ++message_block_index)
        {
          message_block[ message_block_index] = 0;
        }

        ProcessMessageBlock( message_block, w);
        message_block_index = 0;
      }
      for( ; message_block_index < prepadded_message_size; ++message_block_index)
      {
        message_block[ message_block_index] = 0;
      }

      // add the length (in bits) to the end of the sha1 message as a 64 bit integer
      const uint64_t length( uint64_t( string_size) * 8);
      message_block[ 56] = (length >> 56)       ;
      message_block[ 57] = (length >> 48) & 0xFF;
      message_block[ 58] = (length >> 40) & 0xFF;
      message_block[ 59] = (length >> 32) & 0xFF;
      message_block[ 60] = (length >> 24) & 0xFF;
      message_block[ 61] = (length >> 16) & 0xFF;
      message_block[ 62] = (length >>  8) & 0xFF;
      message_block[ 63] = (length      ) & 0xFF;

      ProcessMessageBlock( message_block, w);
    }

    //! @brief convert SHA1 hashed string to hex string
    std::string Sha1::Sha1State::ToHexString() const
    {
      static const char s_hex[ 16] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F'};
      std::string characters( e_OutputSize * 2, ' ');
      for( size_t i( 0), c( 0); i < e_OutputSizeInts; ++i, ++c)
      {
        characters[ c]   = s_hex[ (  m_MessageDigest[ i] >> 28)      ];
        characters[ ++c] = s_hex[ (( m_MessageDigest[ i] >> 24) & 15)];
        characters[ ++c] = s_hex[ (( m_MessageDigest[ i] >> 20) & 15)];
        characters[ ++c] = s_hex[ (( m_MessageDigest[ i] >> 16) & 15)];
        characters[ ++c] = s_hex[ (( m_MessageDigest[ i] >> 12) & 15)];
        characters[ ++c] = s_hex[ (( m_MessageDigest[ i] >>  8) & 15)];
        characters[ ++c] = s_hex[ (( m_MessageDigest[ i] >>  4) & 15)];
        characters[ ++c] = s_hex[ (  m_MessageDigest[ i]        & 15)];
      }
      return characters;
    }

    //! @brief  Performs a circular left shift operation
    inline uint32_t CircularShift( const int &BITS, const uint32_t &WORD)
    {
      return ( WORD << BITS) | ( WORD >> ( 32 - BITS));
    }

    //! @brief Process the next 512 bits of the message
    //! @param MESSAGE_BLOCK characters that make up the message
    //! @param STORAGE a temporary array used during processing
    void Sha1::Sha1State::ProcessMessageBlock( unsigned char *MESSAGE_BLOCK, uint32_t *STORAGE)
    {
      for( size_t c( 0), i( 0); c < e_BlockSize; ++i, ++c)
      {
        // This will swap endian on big endian and keep endian on little endian.
        uint32_t value( ( uint32_t)( MESSAGE_BLOCK[ c]) << 24);
        value |= uint32_t( MESSAGE_BLOCK[ ++c]) << 16;
        value |= uint32_t( MESSAGE_BLOCK[ ++c]) << 8;
        value |= uint32_t( MESSAGE_BLOCK[ ++c]);
        STORAGE[ i] = value;
      }

      for( size_t i( 16); i < 80; ++i)
      {
        STORAGE[ i] =
          CircularShift( 1, STORAGE[ i - 3] ^ STORAGE[ i - 8] ^ STORAGE[ i - 14] ^ STORAGE[ i - 16]);
      }

      uint32_t a( m_MessageDigest[ 0]);
      uint32_t b( m_MessageDigest[ 1]);
      uint32_t c( m_MessageDigest[ 2]);
      uint32_t d( m_MessageDigest[ 3]);
      uint32_t e( m_MessageDigest[ 4]);

      for( size_t t( 0); t < 20; ++t)
      {
        uint32_t temp( CircularShift( 5, a) + ( ( b & c) | ( ~b & d)) + e + STORAGE[ t] + 0x5A827999);
        e = d;
        d = c;
        c = CircularShift( 30, b);
        b = a;
        a = temp;
      }

      for( size_t t( 20); t < 40; ++t)
      {
        uint32_t temp( CircularShift( 5, a) + ( b ^ c ^ d) + e + STORAGE[ t] + 0x6ED9EBA1);
        e = d;
        d = c;
        c = CircularShift( 30, b);
        b = a;
        a = temp;
      }

      for( size_t t( 40); t < 60; ++t)
      {
        uint32_t temp( CircularShift( 5, a) + ( ( b & c) | ( b & d) | ( c & d)) + e + STORAGE[ t] + 0x8F1BBCDC);
        e = d;
        d = c;
        c = CircularShift( 30, b);
        b = a;
        a = temp;
      }

      for( size_t t( 60); t < 80; ++t)
      {
        uint32_t temp( CircularShift( 5, a) + ( b ^ c ^ d) + e + STORAGE[ t] + 0xCA62C1D6);
        e = d;
        d = c;
        c = CircularShift( 30, b);
        b = a;
        a = temp;
      }

      m_MessageDigest[ 0] += a;
      m_MessageDigest[ 1] += b;
      m_MessageDigest[ 2] += c;
      m_MessageDigest[ 3] += d;
      m_MessageDigest[ 4] += e;
    }

    //! @brief encrypt a plain text message and return an encrypted message
    //! @param ORIGINAL_MESSAGE plain text message as a string
    //! @return hashed message (in hex)
    std::string Sha1::Hash( const std::string &ORIGINAL_MESSAGE) const
    {
      return Sha1State( ORIGINAL_MESSAGE).ToHexString();
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Sha1::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Sha1::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  } // namespace crypt
} // namespace bcl
