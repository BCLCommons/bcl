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

#ifndef BCL_IO_STREAM_BUFFER_INTERFACE_H_
#define BCL_IO_STREAM_BUFFER_INTERFACE_H_

// include the namespace header
#include "bcl_io.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_io_stream_buffer_classes.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically
#include <streambuf>

namespace bcl
{
  namespace io
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class StreamBufferInterface
    //! @brief class defines the basic functionality that is needed to be provided for any StreamBuffer class
    //!
    //! @remarks example unnecessary
    //! @author alexanns, woetzen
    //! @date 09/20/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API StreamBufferInterface :
      public std::streambuf,
      public util::ObjectInterface
    {
    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief destructor
      virtual ~StreamBufferInterface();

    /////////////////
    // data access //
    /////////////////

      //! @brief return the default filename extension
      //! @return extension of a filename of that stream type as string
      virtual const std::string &GetDefaultFileExtension() const = 0;

      //! @brief add the default file extension to a filename, separated by a "." if the extension is non-empty
      //! @param FILENAME name of the file to add extension to
      //! @return filename with proper compression extension appended
      std::string AddExtension( const std::string &FILENAME) const;

      //! @brief return the StreamBufferClass
      //! @return StreamBufferClass representing the type of compression
      virtual const StreamBufferClass &GetStreamBufferClass() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief open opens a StreamBufferInterface from filename and open_mode
      //! @param NAME name of the file which this StreamBufferInterface will be associated with
      //! @param OPEN_MODE the manner with which this StreamBufferInterface should be opened
      //!        for explanation on the types and use of open modes please see
      //!        <http://www.cplusplus.com/reference/iostream/ios_base/openmode.html>
      //! @return returns a pointer to this StreamBufferInterface if successful; otherwise returns a null pointer
      virtual StreamBufferInterface *Open
      (
        const char *NAME, const std::ios_base::openmode OPEN_MODE
      ) = 0;

      //! @brief is_open indicates whether a StreamBufferInterface is currently associated with a file
      //! @return boolean true if StreamBufferInterface is open and false if StreamBufferInterface is not open
      virtual bool IsOpen() = 0;

      //! @brief close dissociates this StreamBufferInterface from the file it is currently bound to
      //! @return returns a pointer to this StreamBufferInterface if successful; otherwise returns a null pointer
      virtual StreamBufferInterface *Close() = 0;

      //! @brief check if given open mode is valid for that stream buffer class
      //! @param OPEN_MODE desired open mode
      //! @return false is open mode is impossible for that stream class - true otherwise
      virtual bool IsValidOpenMode( const std::ios_base::openmode OPEN_MODE) = 0;

    }; // class StreamBufferInterface

  } // namespace io
} // namespace bcl

#endif // BCL_IO_STREAM_BUFFER_INTERFACE_H_
