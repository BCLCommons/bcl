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
#include "io/bcl_io_file_stream_buffer.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace io
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FileStreamBuffer::FileStreamBuffer() :
      std::filebuf()
    {
    }

    //! @brief copy constructor
    FileStreamBuffer::FileStreamBuffer( const FileStreamBuffer &BUFFER)
    {
    }

    //! @brief Clone constructor
    FileStreamBuffer *FileStreamBuffer::Clone() const
    {
      return new FileStreamBuffer();
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FileStreamBuffer::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the default filename extension
    //! @return extension of a filename of that stream type as string
    const std::string &FileStreamBuffer::GetDefaultFileExtension() const
    {
      // std::string s_file_extension is the file extension associated with uncompressed files
      static const std::string s_file_extension( "");

      return s_file_extension;
    }

    //! @brief return the StreamBufferClass
    //! @return StreamBufferClass representing the type of compression
    const StreamBufferClass &FileStreamBuffer::GetStreamBufferClass() const
    {
      return GetStreamBufferClasses().e_Uncompressed;
    }

  ///////////////
  // operators //
  ///////////////

  ////////////////
  // operations //
  ////////////////

    //! @brief open opens a FileStreamBuffer from filename and open_mode
    //! @param NAME name of the file which this FileStreamBuffer will be associated with
    //! @param OPEN_MODE the manner with which this FileStreamBuffer should be opened
    //!        for explanation on the types and use of open modes please see
    //!        <http://www.cplusplus.com/reference/iostream/ios_base/openmode.html>
    //! @return returns a pointer to this FileStreamBuffer if successful; otherwise returns a null pointer
    StreamBufferInterface *FileStreamBuffer::Open
    (
      const char *NAME, const std::ios_base::openmode OPEN_MODE
    )
    {
      std::filebuf::open( NAME, OPEN_MODE);
      return this;
    }

    //! @brief check if stream is open
    bool FileStreamBuffer::IsOpen()
    {
      return std::filebuf::is_open();
    }

    //! @brief close dissociates this FileStreamBuffer from the file it is currently bound to
    //! @return returns a pointer to this FileStreamBuffer if successful; otherwise returns a null pointer
    StreamBufferInterface *FileStreamBuffer::Close()
    {
      std::filebuf::close();
      return this;
    }

    //! @brief check if given open mode is valid for that stream buffer class
    //! @param OPEN_MODE desired open mode
    //! @return false is open mode is impossible for that stream class - true otherwise
    bool FileStreamBuffer::IsValidOpenMode( const std::ios_base::openmode OPEN_MODE)
    {
      if //< true if "OPEN_MODE" is not valid
      (
        ( ( OPEN_MODE & std::ios::in) && ( OPEN_MODE & std::ios::out)) //< unsupported combination of modes
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
    std::istream &FileStreamBuffer::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write restraint to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &FileStreamBuffer::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace io
} // namespace bcl
