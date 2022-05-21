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
#include "bzip2/bcl_bzip2.h"
#include "io/bcl_io_stream_buffer_interface.h"

// external includes - sorted alphabetically
#include "bzlib.h"
#include <cstring>
#include <fstream>

namespace bcl
{
  namespace bzip2
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class StreamBuffer
    //! @brief Class implements a bz2 compressions stream
    //! @details Class implements a bz2 compressions stream for the stream buffer interface
    //!
    //! @author alexanns, woetzen, mendenjl
    //! @date 09/20/2008
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API StreamBuffer :
      public io::StreamBufferInterface
    {

    private:

    //////////
    // data //
    //////////

      static const unsigned int s_InBufferSize = 1024;    //!< size of data buffer
      static const unsigned int s_OutBufferSize = 2048;   //!< size of data buffer

      std::filebuf       m_FileBuffer;                //!< the actual filebuffer
      std::vector< char> m_Buffer;                    //!< buffer
      bz_stream          m_BZStream;                  //!< the bz stream

      std::ios_base::openmode m_Mode;                 //!< I/O mode
      bool                    m_Opened;               //!< open/close state of stream

      // for outbuffer
      std::vector< char> m_OutBuffer;

      // for input buffer
      std::vector< char> m_InBuffer;
      char* m_PutbackEnd;
      char* m_InBegin;
      char* m_InEnd;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief default constructor
      StreamBuffer();

      //! @brief copy constructor
      StreamBuffer( const StreamBuffer &BUFFER);

      //! @brief Clone constructor
      StreamBuffer *Clone() const;

      //! @brief destructor
      ~StreamBuffer();

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the default filename extension
      //! @return extension of a filename of that stream type as string
      const std::string &GetDefaultFileExtension() const;

      //! @brief return the StreamBufferClass
      //! @return StreamBufferClass representing the type of compression
      const io::StreamBufferClass &GetStreamBufferClass() const;

      //! @brief is_open indicates whether a StreamBuffer is currently associated with a file
      //! @return boolean true if StreamBuffer is open and false if StreamBuffer is not open
      bool IsOpen();

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      //! @brief goto a particular location in the stream (either absolute or relative to the current position)
      //! @param OFF offset to jump to
      //! @param WAY (std::ios_base:ssekdir::(beg/cur/end)) denotes whether the offset is relative to the stream
      //!            beginning, current position, or end (respectively)
      //! @param WHICH openmode what the current openmode is
      //! Overrides std::streambuf::seekoff, so BCL function naming guidelines cannot be used
      std::streampos seekoff
      (
        std::streamoff OFF,
        std::ios_base::seekdir WAY,
        std::ios_base::openmode WHICH = std::ios_base::in | std::ios_base::out
      );

      //! @brief open opens a StreamBuffer from filename and open_mode
      //! @param NAME name of the file which this StreamBuffer will be associated with
      //! @param OPEN_MODE the manner with which this StreamBuffer should be opened
      //!        for explanation on the types and use of open modes please see
      //!        <http://www.cplusplus.com/reference/iostream/ios_base/openmode.html>
      //! @return returns a pointer to this StreamBuffer if successful; otherwise returns a null pointer
      io::StreamBufferInterface *Open( const char *NAME, const std::ios_base::openmode OPEN_MODE);

      //! @brief close dissociates this StreamBuffer from the file it is currently bound to
      //! @return returns a pointer to this StreamBuffer if successful; otherwise returns a null pointer
      io::StreamBufferInterface *Close();

      //! @brief check if given open mode is valid for that stream buffer class
      //! @param OPEN_MODE desired open mode
      //! @return flase is open mode is impossible for that stream class - true otherwise
      bool IsValidOpenMode( const std::ios_base::openmode OPEN_MODE);

    protected:

    /////////////////////////////////////////////////////
    // Functions related to reading compressed streams //
    /////////////////////////////////////////////////////

      //! @brief underflow is for handling the instance where the buffer has been emptied but more data is to come
      //!        i.e. data is being read from the buffer more quickly than the buffer is being filled
      //! @return returns int signaling the state of the buffer
      int underflow();

    /////////////////////////////////////////////////////
    // Functions related to writing compressed streams //
    /////////////////////////////////////////////////////

      //! @brief sync the stream; essentially a virtual version of Flush inherited from streambuff
      //!        i.e. outputs any remaining characters in the buffer to the file ()
      //! @return 0 on success, -1 otherwise
      //! Overrides std::streambuf::sync, so BCL function naming guidelines cannot be used
      int sync();

      //! @brief writes up to NUM characters from BUFFER to the output stream
      //! @return the number of characters written
      //! Overrides std::streambuf::xsputn, so BCL function naming guidelines cannot be used
      std::streamsize xsputn( const char *BUFFER, std::streamsize NUM);

      //! @brief overflow (used for output buffer only) writes given integer to current position of the put pointer
      //!        see http://www.cplusplus.com/reference/iostream/filebuf/overflow.html for more information
      //! @param c integer to be written at the current position of the put pointer
      //! @return returns an int signaling the state of the buffer
      //! Overrides std::streambuf::overflow, so BCL function naming guidelines cannot be used
      int overflow( int c);

      //! @brief Flush the stream to the output; this should ordinarily only be called at closing, otherwise it will
      //!        degrade compression.
      //! @return true on success
      bool Flush();

      //! @brief compress the current block
      //! @return true on success
      bool CompressCurrentBlock();

      //! @brief finalize the output stream by completing the current Burrows-Wheeler block
      //! @return true on success
      bool FinalizeOutputStream();

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read restraint from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write restraint to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief initialize; calls InitInput or InitOutput and catches errors
      void Initialize();

      //! @brief initialize if used as input buffer
      //! @return int return value from BZ2_bzDecompressInit
      int InitInput();

      //! @brief initialize if used as output buffer
      //! @return int return value from BZ2_bzCompressInit
      int InitOutput();

      //! @brief set up pointers to the output buffer
      void SetOutputBufferPointers()
      {
        m_BZStream.next_out = &m_OutBuffer[ 0];
        m_BZStream.avail_out = m_OutBuffer.size();
      }

      //! @brief instance of this stream buffer class in stream_buffer_classes enum
      static const io::StreamBufferClass e_BZ2;

    }; // class StreamBuffer

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    StreamBuffer::StreamBuffer() :
      m_Mode(),
      m_Opened( false)
    {
    }

    //! @brief copy constructor
    StreamBuffer::StreamBuffer( const StreamBuffer &BUFFER) :
      m_Mode(),
      m_Opened( false)
    {
    }

    //! @brief Clone constructor
    StreamBuffer *StreamBuffer::Clone() const
    {
      return new StreamBuffer( *this);
    }

    //! @brief destructor
    StreamBuffer::~StreamBuffer()
    {
      Close();
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &StreamBuffer::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the default filename extension
    //! @return extension of a filename of that stream type as string
    const std::string &StreamBuffer::GetDefaultFileExtension() const
    {
      // std::string s_file_extension is the file extension associated with gzip
      static const std::string s_file_extension( "bz2");
      return s_file_extension;
    }

    //! @brief return the StreamBufferClass
    //! @return StreamBufferClass representing the type of compression
    const io::StreamBufferClass &StreamBuffer::GetStreamBufferClass() const
    {
      return e_BZ2;
    }

  ///////////////
  // operators //
  ///////////////

  ////////////////
  // operations //
  ////////////////

    //! @brief open opens a StreamBuffer from filename and open_mode
    //! @param NAME name of the file which this StreamBuffer will be associated with
    //! @param OPEN_MODE the manner with which this StreamBuffer should be opened
    //!        for explanation on the types and use of open modes please see
    //!        <http://www.cplusplus.com/reference/iostream/ios_base/openmode.html>
    //! @return returns a pointer to this StreamBuffer if successful; otherwise returns a null pointer
    io::StreamBufferInterface *StreamBuffer::Open
    (
      const char *NAME, const std::ios_base::openmode OPEN_MODE
    )
    {
      // not valid for already opened buffer
      if( IsOpen())
      {
        return ( StreamBuffer *) 0;
      }

      m_Mode = OPEN_MODE;

      // no append nor read/write mode
      if( !IsValidOpenMode( OPEN_MODE))
      {
        return ( StreamBuffer *) 0;
      }

      // open the target file stream.  Binary mode is necessary because compressed data is binary and will contain ASCII
      // control characters
      m_FileBuffer.open( NAME, OPEN_MODE | std::ios::binary);

      if( !m_FileBuffer.is_open())
      {
        return ( StreamBuffer *) 0;
      }

      Initialize();

      // successfully opened
      m_Opened = true;

      // end
      return this;
    }

    //! @brief check if stream is open
    bool StreamBuffer::IsOpen()
    {
      return m_Opened;
    }

    //! @brief close dissociates this StreamBuffer from the file it is currently bound to
    //! @return returns a pointer to this StreamBuffer if successful; otherwise returns a null pointer
    io::StreamBufferInterface *StreamBuffer::Close()
    {
      if( IsOpen())
      {
        // when used as output
        if( m_Mode & std::ios::out)
        {
          // Compress remaining characters in the stream
          CompressCurrentBlock();

          // Finalize the output stream
          FinalizeOutputStream();

          // delete the compressor stream
          BZ2_bzCompressEnd( &m_BZStream);
        }

        // when used as input
        if( m_Mode & std::ios::in)
        {
          // uninit the bz2 stream
          BZ2_bzDecompressEnd( &m_BZStream);
        }

        m_FileBuffer.close();
        m_Opened = false;
      }

      return ( StreamBuffer *) 0;
    }

    //! @brief compress the current block
    //! @return true on success
    bool StreamBuffer::CompressCurrentBlock()
    {
      // Precondition for calling BZ2_bzCompress is that next_in be set to the data to be compressed, avail in needs
      // to be set to the number of characters available in the stream
      const int num( pptr() - pbase());
      m_BZStream.next_in = pbase();

      // loop until all characters have been compressed
      for( m_BZStream.avail_in = num; m_BZStream.avail_in;)
      {
        // call bzCompress to compress the data in the stream
        SetOutputBufferPointers();
        BZ2_bzCompress( &m_BZStream, BZ_RUN);

        // write the data to the file stream buffer
        const int n_characters_compressed( m_OutBuffer.size() - m_BZStream.avail_out);

        // test that all data could be written; this could fail if the disk is full or there is a write error
        if( m_FileBuffer.sputn( &m_OutBuffer[ 0], n_characters_compressed) != n_characters_compressed)
        {
          return false;
        }
      }

      // bump the underlying file pointer now that the additional bytes have been written to the stream
      pbump( -num);
      return true;
    }

    //! @brief Flush the stream to the output; this should ordinarily only be called at closing, otherwise it will
    //!        degrade compression.
    //! @return true on success
    //! Overrides std::streambuf::overflow, so BCL function naming guidelines cannot be used
    bool StreamBuffer::Flush()
    {
      m_BZStream.avail_in = 0;
      m_BZStream.next_in = NULL;

      for( bool done( false); !done;)
      {
        // call bzCompress to Flush the stream by stopping current compression block.
        SetOutputBufferPointers();
        done = BZ2_bzCompress( &m_BZStream, BZ_FLUSH) == BZ_RUN_OK;

        // write the data to the (file/string) stream buffer
        const int n_characters_compressed( m_OutBuffer.size() - m_BZStream.avail_out);

        // test that all data could be written; this could fail if the disk is full or there is a write error
        if( m_FileBuffer.sputn( &m_OutBuffer[ 0], n_characters_compressed) != n_characters_compressed)
        {
          return false;
        }
      }

      return true;
    }

    //! @brief finalize the output stream by completing the current Burrows-Wheeler block
    //! @return true on success
    bool StreamBuffer::FinalizeOutputStream()
    {
      m_BZStream.avail_in = 0;
      m_BZStream.next_in = NULL;

      for( bool done( false); !done;)
      {
        // call bzCompress to finish the compression
        SetOutputBufferPointers();
        done = BZ2_bzCompress( &m_BZStream, BZ_FINISH) == BZ_STREAM_END;

        // write the data to the (file/string) stream buffer
        const int n_characters_compressed( m_OutBuffer.size() - m_BZStream.avail_out);

        // test that all data could be written; this could fail if the disk is full or there is a write error
        if( m_FileBuffer.sputn( &m_OutBuffer[ 0], n_characters_compressed) != n_characters_compressed)
        {
          return false;
        }
      }

      return true;
    }

    //! @brief underflow is for handling the instance where the buffer has been emptied but more data is to come
    //!        i.e. data is being read from the buffer more quickly than the buffer is being filled
    //! @return returns int signaling the state of the buffer
    int StreamBuffer::underflow()
    {
      // function does nothing for output streams, since they are not reading from the buffer
      if( m_Mode & std::ios::out)
      {
        return std::streambuf::underflow();
      }

      // calculate the number of characters to put back == smaller of # in the read-in buffer, or the area reserved
      // for decompression
      char *buffer_begin( m_Buffer.begin().base()), *buffer_end( m_Buffer.end().base());
      const size_t putback_num( std::min( gptr() - eback(), m_PutbackEnd - buffer_begin));

      // copy the putback data into the putback area
      std::memcpy( m_PutbackEnd - putback_num, gptr() - putback_num, putback_num);

      // Read into the bzstream's buffers; keep calling BZ2_Decompress until something is actually decompressed
      do
      {
        // read from the file buffer if the input buffer is empty
        if( m_InBegin == m_InEnd)
        {
          std::streamsize num_read( m_FileBuffer.sgetn( m_InBuffer.begin().base(), m_InBuffer.size()));
          if( num_read == 0)
          {
            // no characters read; presumably at end of file (maybe should do some error checking here to ensure that is
            // the case though)
            return traits_type::eof();
          }
          m_InBegin = m_InBuffer.begin().base();
          m_InEnd = m_InBegin + num_read;
        }

        // setup the bz stream's pointers for decompression
        m_BZStream.next_in = m_InBegin;
        m_BZStream.avail_in = m_InEnd - m_InBegin;
        m_BZStream.next_out = m_PutbackEnd;
        m_BZStream.avail_out = buffer_end - m_PutbackEnd;

        // attempt decompression
        const int ret( BZ2_bzDecompress( &m_BZStream));
        if( ret == BZ_STREAM_END)
        {
          if( buffer_end - m_PutbackEnd == ( int) m_BZStream.avail_out)
          {
            return traits_type::eof();
          }
        }
        else if( ret != BZ_OK)
        {
          BCL_Exit( "Error reading BZ2 file (corrupted or non-bz2 file?)", -1);
        }

        // update the input buffer pointer
        m_InBegin = m_InEnd - m_BZStream.avail_in;
      } while( m_BZStream.avail_out + m_PutbackEnd == buffer_end);

      // update base streambuf classes' stream buffer pointers
      setg( m_PutbackEnd - putback_num, m_PutbackEnd, buffer_end - m_BZStream.avail_out);

      // return next char
      return traits_type::to_int_type( *gptr());
    }

    //! @brief overflow (used for output buffer only) writes given integer to current position of the put pointer
    //!        see http://www.cplusplus.com/reference/iostream/filebuf/overflow.html for more information
    //! @param c integer to be written at the current position of the put pointer
    //! @return returns an int signaling the state of the buffer
    //! Overrides std::streambuf::overflow, so BCL function naming guidelines cannot be used
    int StreamBuffer::overflow( int c)
    {
      if( m_Mode & std::ios::in)
      {
        return std::streambuf::overflow( c);
      }

      if( !traits_type::eq_int_type( c, traits_type::eof()))
      {
        // put this character in the last position, which will be pptr whenever this function is called, which points
        // to a valid position due to the overallocation
        *pptr() = c;
        pbump( 1);
      }

      return CompressCurrentBlock() ? traits_type::not_eof( c) : traits_type::eof();
    }

    //! @brief sync synchronizes what has been written to the file and what has been put into the buffer
    //!        i.e. outputs any remaining characters in the buffer to the file ()
    //! @return returns integer signaling the state of the buffer
    int StreamBuffer::sync()
    {
      if( m_Mode & std::ios::in)
      {
        return std::streambuf::sync();
      }

      return CompressCurrentBlock() && Flush() ? 0 : -1;
    }

    //! @brief goto a particular location in the stream (either absolute or relative to the current position)
    //! @param OFF offset to jump to
    //! @param WAY (std::ios_base:ssekdir::(beg/cur/end)) denotes whether the offset is relative to the stream
    //!            beginning, current position, or end (respectively)
    //! @param WHICH openmode what the current openmode is
    //! Overrides std::streambuf::seekoff, so BCL function naming guidelines cannot be used
    std::streampos StreamBuffer::seekoff
    (
      std::streamoff OFF,
      std::ios_base::seekdir WAY,
      std::ios_base::openmode WHICH
    )
    {
      // TODO: align offset based on input size and output size
      if( WAY == std::ios_base::cur && OFF == 0)
      {
        return std::streampos( 0);
      }
      m_FileBuffer.pubseekoff( OFF, WAY, WHICH);

      // ensure that this is an offset that can be directly jumped to; currently only the beginning is supported
      BCL_Assert
      (
        OFF == 0 && ( WAY == std::ios_base::beg || WAY == std::ios_base::cur),
        this->GetClassIdentifier() + " can only seek the beginning of the file, tried to go to position " +
        util::Format()( OFF) + " relative to " + util::Format()( int( WAY)) + " compare to " + util::Format()( int( std::ios_base::beg))
      );

      Initialize();
      return std::streampos( OFF);
    }

    //! @brief writes up to NUM characters from BUFFER to the output stream
    //! @return the number of characters written
    //! Overrides std::streambuf::xsputn, so BCL function naming guidelines cannot be used
    std::streamsize StreamBuffer::xsputn( const char *BUFFER, std::streamsize NUM)
    {
      // nothing to be done for input streams
      if( m_Mode & std::ios::in)
      {
        return std::streambuf::xsputn( BUFFER, NUM);
      }

      std::streamsize chars_copied( 0);

      // loop while we haven't reached the end of BUFFER
      while( chars_copied < NUM)
      {
        // get the number of chars between epptr and pptr, inclusive
        const std::streamsize block_len( std::min( NUM - chars_copied, std::streamsize( epptr() - pptr() + 1)));

        // write them
        memcpy( pptr(), BUFFER + chars_copied, size_t( block_len));

        // update write pointer
        pbump( int( block_len));

        // compress the current block if necessary
        if( pptr() >= epptr() && !CompressCurrentBlock())
        {
          break;
        }

        // update the count of total characters copied
        chars_copied += block_len;
      }

      return chars_copied;
    }

    //! @brief check if given open mode is valid for that stream buffer class
    //! @param OPEN_MODE desired open mode
    //! @return false is open mode is impossible for that stream class - true otherwise
    bool StreamBuffer::IsValidOpenMode( const std::ios_base::openmode OPEN_MODE)
    {
      if //< true if "OPEN_MODE" is not valid
      (
           ( OPEN_MODE & std::ios::ate) //< unsupported mode
        || ( OPEN_MODE & std::ios::app) //< unsupported mode
        || ( ( OPEN_MODE & std::ios::in) && ( OPEN_MODE & std::ios::out)) //< unsupported combination of modes
      )
      {
        return false;
      }

      // open mode is valid so return true
      return true;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read restraint from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &StreamBuffer::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write restraint to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &StreamBuffer::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief initialize; calls InitInput or InitOutput and catches errors
    void StreamBuffer::Initialize()
    {
      // set modeptr according to OPEN_MODE
      int initialization_return_value( 0);
      if( m_Mode & std::ios::in)
      {
        initialization_return_value = InitInput();
      }
      else if( m_Mode & std::ios::out)
      {
        initialization_return_value = InitOutput();
      }
      if( initialization_return_value == BZ_CONFIG_ERROR)
      {
        BCL_Exit( "libbz2 was not compiled correctly.", -1);
      }
      else if( initialization_return_value == BZ_MEM_ERROR)
      {
        BCL_Exit( "libbz2 Allocation problem", -1);
      }
      else if( initialization_return_value != BZ_OK)
      {
        BCL_Exit( "Unknown error creating bz2 decompressor stream buffer.", -1);
      }
    }

    //! @brief initialize if used as input buffer
    int StreamBuffer::InitInput()
    {
      // allocate the buffers
      m_Buffer.resize( s_InBufferSize);
      m_InBuffer.resize( s_InBufferSize);
      m_InBegin = m_InEnd = m_InBuffer.begin().base();

      // set the get pointers so as to force an underflow on the first read
      m_PutbackEnd = &m_Buffer[ 64];
      setg( m_PutbackEnd, m_PutbackEnd, m_PutbackEnd);

      // zero out the decompressor stream; BZ2_bzCompressInit should work regardless, but this safeguards against poorly
      // implemented BZ2_bzDecompressInit functions that may depend on the default bz_stream object constructor having
      // just been called
      memset( &m_BZStream, 0, sizeof( m_BZStream));

      // use standard memory allocation routines
      m_BZStream.bzalloc = NULL;
      m_BZStream.bzfree = NULL;
      m_BZStream.opaque = NULL;

      // init a decompressor stream with default values; no verbosity, and prefer speed over memory conservation
      return BZ2_bzDecompressInit( &m_BZStream, 0, 0);
    }

    //! @brief initialize if used as output buffer
    int StreamBuffer::InitOutput()
    {
      // allocate the working and output buffers
      m_Buffer.resize( s_OutBufferSize);
      m_OutBuffer.resize( s_OutBufferSize);

      // set the buffer pointers; reserve an extra character at the end of the buffer
      setp( m_Buffer.begin().base(), &*--m_Buffer.end());

      // zero out the compressor stream; BZ2_bzCompressInit should work regardless, but this safeguards against poorly
      // implemented BZ2_bzCompressInit functions that may depend on the default bz_stream object constructor having
      // just been called
      memset( &m_BZStream, 0, sizeof( m_BZStream));

      // use standard memory allocation routines
      m_BZStream.opaque = NULL;
      m_BZStream.bzfree = NULL;
      m_BZStream.bzalloc = NULL;

      // Default parameters (constant for now);
      // Verbosity == 0 -> no verbose output
      // WorkFactor == 0 -> use default work factor to determine when to use the backup algorithm; default works well
      // BlockSize == 9 -> Maximum compression at some minor expense in memory usage
      const unsigned int verbosity( 0), work_factor( 0), block_size( 9);

      // initialize the stream
      return BZ2_bzCompressInit( &m_BZStream, block_size, verbosity, work_factor);
    }

    //! @brief instance of this stream buffer class in stream_buffer_classes enum
    const io::StreamBufferClass StreamBuffer::e_BZ2
    (
      io::StreamBufferClasses::GetEnums().AddEnum( "BZ2", util::ShPtr< io::StreamBufferInterface>( new StreamBuffer()))
    );

  } // namespace bzip2
} // namespace bcl
