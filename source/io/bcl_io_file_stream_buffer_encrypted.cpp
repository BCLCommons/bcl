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
#include "io/bcl_io_stream_buffer_classes.h"
#include "io/bcl_io_stream_buffer_interface.h"

// external includes - sorted alphabetically
#include <fstream>

namespace bcl
{
  namespace io
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FileStreamBufferEncrypted
    //! @brief Class implements a encrypted file stream
    //! @details Class implements an xor encrypted stream for the stream buffer interface
    //! http://www.ifj.edu.pl/~rmaj/c++/ex/streambuf_2_encryption/x.cpp.html
    //!
    //! @author woetzen
    //! @date 04/08/2011
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FileStreamBufferEncrypted :
      public StreamBufferInterface
    {

    private:

    //////////
    // data //
    //////////

      std::filebuf          m_FileBuffer; //!< the actual filebuffer
      traits_type::int_type m_Key;        //!< key for encryption and decryption
      traits_type::pos_type m_Pos;        //!< this filter (decryption) specific

      std::ios_base::openmode m_Mode;   //!< I/O mode
      bool                    m_Opened; //!< open/close state of stream

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief default constructor
      FileStreamBufferEncrypted();

      //! @brief copy constructor
      FileStreamBufferEncrypted( const FileStreamBufferEncrypted &BUFFER);

      //! @brief Clone constructor
      FileStreamBufferEncrypted *Clone() const;

      //! @brief destructor
      ~FileStreamBufferEncrypted();

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
      const StreamBufferClass &GetStreamBufferClass() const;

    ////////////////
    // operations //
    ////////////////

      traits_type::int_type underflow();

      traits_type::int_type uflow();

      traits_type::int_type overflow( traits_type::int_type c);

      //! @brief open opens a StreamBuffer from filename and open_mode
      //! @param NAME name of the file which this StreamBuffer will be associated with
      //! @param OPEN_MODE the manner with which this StreamBuffer should be opened
      //!        for explanation on the types and use of open modes please see
      //!        <http://www.cplusplus.com/reference/iostream/ios_base/openmode.html>
      //! @return returns a pointer to this StreamBuffer if successful; otherwise returns a null pointer
      StreamBufferInterface *Open
      (
        const char *NAME, const std::ios_base::openmode OPEN_MODE
      );

      //! @brief is_open indicates whether a StreamBuffer is currently associated with a file
      //! @return boolean true if StreamBuffer is open and false if StreamBuffer is not open
      bool IsOpen();

      //! @brief close dissociates this StreamBuffer from the file it is currently bound to
      //! @return returns a pointer to this StreamBuffer if successful; otherwise returns a null pointer
      StreamBufferInterface *Close();

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

      traits_type::char_type decode( traits_type::char_type CHAR, const traits_type::pos_type POS)
      {
        const traits_type::char_type decoded( CHAR ^ ( POS ^ m_Key));
        return traits_type::to_int_type( decoded) == traits_type::eof() ? CHAR : decoded;
      }

      traits_type::char_type encode( traits_type::char_type CHAR, const traits_type::pos_type POS)
      {
        const traits_type::char_type encoded( CHAR ^ ( POS ^ m_Key));
        return traits_type::to_int_type( encoded) == traits_type::eof() ? CHAR : encoded;
      }

      //! @brief instance of this stream buffer class in stream_buffer_classes enum
      static const StreamBufferClass e_Crypt;

    }; // class FileStreamBufferEncrypted

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FileStreamBufferEncrypted::FileStreamBufferEncrypted() :
      m_FileBuffer(),
      m_Key( 1234),
      m_Pos( 0),
      m_Mode(),
      m_Opened( false)
    {
    }

    //! @brief copy constructor
    FileStreamBufferEncrypted::FileStreamBufferEncrypted( const FileStreamBufferEncrypted &BUFFER) :
      m_FileBuffer(),
      m_Key( BUFFER.m_Key),
      m_Pos( 0),
      m_Mode(),
      m_Opened( false)
    {
    }

    //! @brief Clone constructor
    FileStreamBufferEncrypted *FileStreamBufferEncrypted::Clone() const
    {
      return new FileStreamBufferEncrypted( *this);
    }

    //! @brief destructor
    FileStreamBufferEncrypted::~FileStreamBufferEncrypted()
    {
      Close();
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FileStreamBufferEncrypted::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the default filename extension
    //! @return extension of a filename of that stream type as string
    const std::string &FileStreamBufferEncrypted::GetDefaultFileExtension() const
    {
      // std::string s_file_extension is the file extension associated with gzip
      static const std::string s_file_extension( "crypt");
      return s_file_extension;
    }

    //! @brief return the StreamBufferClass
    //! @return StreamBufferClass representing the type of compression
    const StreamBufferClass &FileStreamBufferEncrypted::GetStreamBufferClass() const
    {
      return e_Crypt;
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
    StreamBufferInterface *FileStreamBufferEncrypted::Open
    (
      const char *NAME, const std::ios_base::openmode OPEN_MODE
    )
    {
      // not valid for already opened buffer
      if( IsOpen())
      {
        return ( FileStreamBufferEncrypted *) 0;
      }

      m_Mode = OPEN_MODE;

      // no append nor read/write mode
      if( !IsValidOpenMode( OPEN_MODE))
      {
        return ( FileStreamBufferEncrypted *) 0;
      }

      // open the target file stream
      m_FileBuffer.open( NAME, OPEN_MODE | std::ios::binary);

      if( !m_FileBuffer.is_open())
      {
        return ( FileStreamBufferEncrypted *) 0;
      }

      // sucessfully opened
      m_Opened = true;

      // end
      return this;
    }

    //! @brief check if stream is open
    bool FileStreamBufferEncrypted::IsOpen()
    {
      return m_Opened;
    }

    //! @brief close dissociates this StreamBuffer from the file it is currently bound to
    //! @return returns a pointer to this StreamBuffer if successful; otherwise returns a null pointer
    StreamBufferInterface *FileStreamBufferEncrypted::Close()
    {
      if( IsOpen())
      {
        m_FileBuffer.close();
        m_Opened = false;
      }

      return ( FileStreamBufferEncrypted *) 0;
    }

    //! @brief underflow is for handling the instance where the buffer has been emptied but more data is to come
    //!        i.e. data is being read from the buffer more quickly than the buffer is being filled
    //! @return returns int signaling the state of the buffer
    FileStreamBufferEncrypted::traits_type::int_type FileStreamBufferEncrypted::underflow()
    {
      const traits_type::int_type fromUnder( m_FileBuffer.sgetc()); // read from underlying buf
      if( fromUnder == traits_type::eof())
      {
        return traits_type::eof(); // EOF
      }
      const traits_type::int_type coded( traits_type::to_int_type( decode( traits_type::to_char_type( fromUnder), m_Pos)));

      return coded;
    }

    FileStreamBufferEncrypted::traits_type::int_type FileStreamBufferEncrypted::uflow()
    {
      const traits_type::int_type fromUnder( m_FileBuffer.sbumpc()); // read from underlying buf
      if( fromUnder == traits_type::eof())
      {
        return traits_type::eof(); // EOF
      }
      const traits_type::int_type coded( traits_type::to_int_type( decode( traits_type::to_char_type( fromUnder), m_Pos)));
      m_Pos += 1;

      return coded;
    }

    //! @brief overflow (used for output buffer only) writes given integer to current position of the put pointer
    //!        see http://www.cplusplus.com/reference/iostream/filebuf/overflow.html for more information
    //! @param c integer to be written at the current position of the put pointer
    //! @return returns an int signaling the state of the buffer
    FileStreamBufferEncrypted::traits_type::int_type FileStreamBufferEncrypted::overflow( traits_type::int_type c)
    {
      if( c == traits_type::eof())
      {
        return traits_type::eof(); // EOF
      }

      c = m_FileBuffer.sputc( encode( traits_type::to_char_type( c), m_Pos));
      m_Pos += 1;

      return c;
    }

    //! @brief check if given open mode is valid for that stream buffer class
    //! @param OPEN_MODE desired open mode
    //! @return false is open mode is impossible for that stream class - true otherwise
    bool FileStreamBufferEncrypted::IsValidOpenMode( const std::ios_base::openmode OPEN_MODE)
    {
      // open mode is valid so return true
      return true;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read restraint from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FileStreamBufferEncrypted::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write restraint to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &FileStreamBufferEncrypted::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief instance of this stream buffer class in stream_buffer_classes enum
    const StreamBufferClass FileStreamBufferEncrypted::e_Crypt
    (
      StreamBufferClasses::GetEnums().AddEnum( "Crypt", util::ShPtr< io::StreamBufferInterface>( new FileStreamBufferEncrypted()))
    );

  } // namespace io
} // namespace bcl
