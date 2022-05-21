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

#ifndef BCL_IO_FILE_STREAM_BUFFER_H_
#define BCL_IO_FILE_STREAM_BUFFER_H_

// include the namespace header
#include "bcl_io.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_io_stream_buffer_interface.h"

// external includes - sorted alphabetically
#include <fstream>

namespace bcl
{
  namespace io
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FileStreamBuffer
    //! @brief FileStreamBuffer is a StreamBufferInterface which is to be used for uncompressed input and output
    //!
    //! @see @link example_io_file_stream_buffer.cpp @endlink
    //! @author alexanns, woetzen
    //! @date 09/20/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FileStreamBuffer :
      public StreamBufferInterface,
      std::filebuf
    {

    //////////
    // data //
    //////////

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FileStreamBuffer();

      //! @brief copy constructor
      FileStreamBuffer( const FileStreamBuffer &BUFFER);

      //! @brief Clone constructor
      FileStreamBuffer *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the default filename extension
      //! @return extension of a filename of that stream type as string
      virtual const std::string &GetDefaultFileExtension() const;

      //! @brief return the StreamBufferClass
      //! @return StreamBufferClass representing the type of compression
      virtual const StreamBufferClass &GetStreamBufferClass() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief open opens a FileStreamBuffer from filename and open_mode
      //! @param NAME name of the file which this FileStreamBuffer will be associated with
      //! @param OPEN_MODE the manner with which this FileStreamBuffer should be opened
      //!        for explanation on the types and use of open modes please see
      //!        <http://www.cplusplus.com/reference/iostream/ios_base/openmode.html>
      //! @return returns a pointer to this FileStreamBuffer if successful; otherwise returns a null pointer
      virtual StreamBufferInterface *Open
      (
        const char *NAME, const std::ios_base::openmode OPEN_MODE
      );

      //! @brief is_open indicates whether a FileStreamBuffer is currently associated with a file
      //! @return boolean true if FileStreamBuffer is open and false if FileStreamBuffer is not open
      virtual bool IsOpen();

      //! @brief close dissociates this FileStreamBuffer from the file it is currently bound to
      //! @return returns a pointer to this FileStreamBuffer if successful; otherwise returns a null pointer
      virtual StreamBufferInterface *Close();

      //! @brief check if given open mode is valid for that stream buffer class
      //! @param OPEN_MODE desired open mode
      //! @return flase is open mode is impossible for that stream class - true otherwise
      virtual bool IsValidOpenMode( const std::ios_base::openmode OPEN_MODE);

    protected:

      //! @brief showmanyc
      std::streamsize showmanyc()
      {
        return std::filebuf::showmanyc();
      }

      //! @brief underflow
      int underflow()
      {
        return std::filebuf::underflow();
      }

      //! @brief uflow
      int uflow()
      {
        return std::filebuf::uflow();
      }

      //! @brief pbackfail
      int pbackfail( int C = EOF)
      {
        return std::filebuf::pbackfail( C);
      }

      //! @brief overflow
      int overflow( int C = EOF)
      {
        return std::filebuf::overflow( C);
      }

      //! @brief setbuf
      std::streambuf *setbuf( char *s, std::streamsize n)
      {
        return std::filebuf::setbuf( s, n);
      }

      //! @brief seekoff
      std::streampos seekoff
      (
        std::streamoff off,
        std::ios_base::seekdir way,
        std::ios_base::openmode which = std::ios_base::in | std::ios_base::out
      )
      {
        return std::filebuf::seekoff( off, way, which);
      }

      //! @brief seekpos
      std::streampos seekpos
      (
        std::streampos sp, std::ios_base::openmode which = std::ios_base::in | std::ios_base::out
      )
      {
        return std::filebuf::seekpos( sp, which);
      }

      //! @brief sync
      int sync()
      {
        return std::filebuf::sync();
      }

      //! @brief imbue
      void imbue( const std::locale &loc)
      {
        return std::filebuf::imbue( loc);
      }
      
      std::streamsize xsputn( const char *p, std::streamsize num)
      {
        return std::filebuf::xsputn( p, num);
      }

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

    }; // class FileStreamBuffer

  } // namespace io
} // namespace bcl
#endif // BCL_IO_FILE_STREAM_BUFFER_H_
