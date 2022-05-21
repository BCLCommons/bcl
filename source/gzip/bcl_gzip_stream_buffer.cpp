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
#include "gzip/bcl_gzip.h"
#include "io/bcl_io_stream_buffer_interface.h"

// external includes - sorted alphabetically
#include "zlib.h"
#include <cstdio>
#include <cstring>

namespace bcl
{
  namespace gzip
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class StreamBuffer
    //! @brief Class implements a GZ compression stream
    //! @details Class implements a GZ compression stream
    //!
    //! @author alexanns, woetzen
    //! @date 09/20/2008
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API StreamBuffer :
      public io::StreamBufferInterface
    {

    //////////
    // data //
    //////////

    private:

      static const int s_BufferSize = 1024;    //!< size of data buffer
      gzFile                  m_File;                 //!< file handle for compressed file
      char                    m_Buffer[s_BufferSize]; //!< data buffer
      bool                    m_Opened;               //!< open/close state of stream
      std::ios_base::openmode m_Mode;                 //!< I/O mode

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

      //! @brief return the io::StreamBufferClass
      //! @return io::StreamBufferClass representing the type of compression
      const io::StreamBufferClass &GetStreamBufferClass() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief open opens a StreamBuffer from filename and open_mode
      //! @param NAME name of the file which this StreamBuffer will be associated with
      //! @param OPEN_MODE the manner with which this StreamBuffer should be opened
      //!        for explanation on the types and use of open modes please see
      //!        <http://www.cplusplus.com/reference/iostream/ios_base/openmode.html>
      //! @return returns a pointer to this StreamBuffer if successful; otherwise returns a null pointer
      io::StreamBufferInterface *Open
      (
        const char *NAME, const std::ios_base::openmode OPEN_MODE
      );

      //! @brief is_open indicates whether a StreamBuffer is currently associated with a file
      //! @return boolean true if StreamBuffer is open and false if StreamBuffer is not open
      bool IsOpen();

      //! @brief close dissociates this StreamBuffer from the file it is currently bound to
      //! @return returns a pointer to this StreamBuffer if successful; otherwise returns a null pointer
      io::StreamBufferInterface *Close();

      //! @brief Seek a given offset in the file
      std::streampos seekoff
      (
        std::streamoff OFF,
        std::ios_base::seekdir WAY,
        std::ios_base::openmode WHICH = std::ios_base::in | std::ios_base::out
      );

      //! @brief overflow (used for output buffer only) writes given integer to current position of the put pointer
      //!        see http://www.cplusplus.com/reference/iostream/filebuf/overflow.html for more information
      //! @param CHARACTER integer to be written at the current position of the put pointer
      //! @return returns an int signaling the state of the buffer
      int overflow( int CHARACTER = EOF);

      //! @brief underflow is for handling the instance where the buffer has been emptied but more data is to come
      //!        i.e. data is being read from the buffer more quickly than the buffer is being filled
      //! @return returns int signalling the state of the buffer
      int underflow();

      //! @brief sync synchronizes what has been written to the file and what has been put into the buffer
      //!        i.e. outputs any remaining characters in the buffer to the file ()
      //! @return returns integer signalling the state of the buffer
      int sync();

      //! @brief check if given open mode is valid for that stream buffer class
      //! @param OPEN_MODE desired open mode
      //! @return flase is open mode is impossible for that stream class - true otherwise
      bool IsValidOpenMode( const std::ios_base::openmode OPEN_MODE);

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

      //! @brief flush_buffer outputs the data remaining in the buffer
      //! @return returns an int signaling the state of the buffer
      int Flush();

      //! @brief instance of this stream buffer class in stream_buffer_classes enum
      static const io::StreamBufferClass e_GZ;

    }; // class StreamBuffer

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    StreamBuffer::StreamBuffer() :
      m_Opened( false)
    {
      setp( m_Buffer, m_Buffer + ( s_BufferSize - 1));
      setg( m_Buffer + 4, m_Buffer + 4, m_Buffer + 4);
    }

    //! @brief copy constructor
    StreamBuffer::StreamBuffer( const StreamBuffer &BUFFER) :
      m_Opened( false)
    {
      setp( m_Buffer, m_Buffer + ( s_BufferSize - 1));
      setg( m_Buffer + 4, m_Buffer + 4, m_Buffer + 4);
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
      static const std::string s_file_extension( "gz");
      return s_file_extension;
    }

    //! @brief return the io::StreamBufferClass
    //! @return io::StreamBufferClass representing the type of compression
    const io::StreamBufferClass &StreamBuffer::GetStreamBufferClass() const
    {
      return e_GZ;
    }

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

      // open the file in read or write mode according to Mode
      m_File = gzopen( NAME, ( m_Mode & std::ios::in ? "rb" : "wb"));

      // check if opening was sucessful
      if( m_File == 0)
      {
        return ( StreamBuffer *) 0;
      }

      // sucessfully opened
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
      // if the file is open and sync (flush) is successful
      if( m_Opened && !sync())
      {
        m_Opened = false;
        // close the file
        if( gzclose( m_File) == Z_OK)
        {
          return this;
        }
      }

      return ( StreamBuffer *) 0;
    }

    //! @brief underflow is for handling the instance where the buffer has been emptied but more data is to come
    //!        i.e. data is being read from the buffer more quickly than the buffer is being filled
    //! @return returns int signalling the state of the buffer
    int StreamBuffer::underflow()
    {
      // This is only used when reading; not in output
      if( gptr() && ( gptr() < egptr()))
      {
        return *reinterpret_cast< unsigned char *>( gptr());
      }

      // test that the file is open for input
      if( !( m_Mode & std::ios::in) || !m_Opened)
      {
        return EOF;
      }

      // copy the last bytes (up to 4) into the first 4 bytes of the buffer
      const int n_putback( std::min( long( 4), long( gptr() - eback())));
      memcpy( m_Buffer + ( 4 - n_putback), gptr() - n_putback, n_putback);

      // Fill m_Buffer with the decompressed stream from m_File
      int characters_read( gzread( m_File, m_Buffer + 4, s_BufferSize - 4));
      if( characters_read <= 0) // ERROR or EOF
      {
        return EOF;
      }

      // update buffer pointers to the read in segment
      setg( m_Buffer + 4 - n_putback, m_Buffer + 4, m_Buffer + 4 + characters_read);

      // return next character
      return *reinterpret_cast< unsigned char *>( gptr());
    }

    //! @brief overflow (used for output buffer only) writes given integer to current position of the put pointer
    //!        see http://www.cplusplus.com/reference/iostream/filebuf/overflow.html for more information
    //! @param CHARACTER integer to be written at the current position of the put pointer
    //! @return returns an int signaling the state of the buffer
    int StreamBuffer::overflow( int CHARACTER)
    {
      if( !( m_Mode & std::ios::out) || !m_Opened)
      {
        return EOF;
      }

      if( CHARACTER != EOF)
      {
        // add the character to the end of the stream
        *pptr() = CHARACTER;
        pbump( 1);
      }

      // attempt flushing to output; on failure, signal EOF
      return Flush() == EOF ? EOF : CHARACTER;
    }

    //! @brief sync synchronizes what has been written to the file and what has been put into the buffer
    //!        i.e. outputs any remaining characters in the buffer to the file ()
    //! @return returns integer signaling the state of the buffer
    int StreamBuffer::sync()
    {
      return pptr() && pptr() > pbase() && Flush() == EOF ? -1 : 0;
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

    //! @brief Seek a given offset in the file
    std::streampos StreamBuffer::seekoff
    (
      std::streamoff OFF,
      std::ios_base::seekdir WAY,
      std::ios_base::openmode WHICH
    )
    {
      // seek the offset in the underlying file buffer
      if( WHICH & std::ios_base::in)
      {
        // input files
        setp( m_Buffer, m_Buffer + ( s_BufferSize - 1));
      }
      else if( WHICH & std::ios_base::out)
      {
        // output files
        setg( m_Buffer + 4, m_Buffer + 4, m_Buffer + 4);
      }
      // call gzseek, only trick is to translate std::ios_base::seekdir enum into the lseek variety
      if( WAY == std::ios_base::beg)
      {
        return gzseek( m_File, OFF, SEEK_SET);
      }
      else if( WAY == std::ios_base::cur)
      {
        return gzseek( m_File, OFF, SEEK_CUR);
      }
      return gzseek( m_File, OFF, SEEK_END);
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

    //! @brief flush_buffer outputs the data remaining in the buffer
    //! @return returns an int signaling the state of the buffer
    int StreamBuffer::Flush()
    {
      // Write everything in the compressed buffer (base class std::streambuf) to m_File
      const int bytes( pptr() - pbase());
      if( gzwrite( m_File, pbase(), bytes) != bytes)
      {
        // could not complete write; signal EOF
        return EOF;
      }

      // move back to the beginning of the write buffer
      pbump( -bytes);

      // return # of bytes written
      return bytes;
    }

    //! @brief instance of this stream buffer class in stream_buffer_classes enum
    const io::StreamBufferClass StreamBuffer::e_GZ
    (
      io::StreamBufferClasses::GetEnums().AddEnum( "GZ", util::ShPtr< io::StreamBufferInterface>( new StreamBuffer()))
    );

  } // namespace gzip
} // namespace bcl
