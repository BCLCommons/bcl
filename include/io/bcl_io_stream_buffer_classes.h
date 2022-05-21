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

#ifndef BCL_IO_STREAM_BUFFER_CLASSES_H_
#define BCL_IO_STREAM_BUFFER_CLASSES_H_

// include the namespace header
#include "bcl_io.h"

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util.h"
#include "util/bcl_util_enumerate.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace io
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class StreamBufferClasses
    //! @brief StreamBufferClasses provides the enumeration of all StreamBufferInterface derived classes.
    //! @details  This class also has the "open" function which provides the functionality of switching compressions over
    //! the command line.
    //!
    //! @see @link example_io_stream_buffer_classes.cpp @endlink
    //! @author alexanns, woetzen
    //! @date 09/26/08
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API StreamBufferClasses :
      public util::Enumerate< util::ShPtr< StreamBufferInterface>, StreamBufferClasses>
    {
      friend class util::Enumerate< util::ShPtr< StreamBufferInterface>, StreamBufferClasses>;

    public:

    //////////
    // data //
    //////////

      // declare all stream buffer classes
      const StreamBufferClass e_Uncompressed; //!< FileStreamBuffer

      //! @brief GetFlagFileCompression gives commandline flag to adjust default file compression
      //! @return ShPtr to a FlagInterface which is used to get the desired compression from the command line
      util::ShPtr< command::FlagInterface> &GetFlagFileCompression() const;

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct all StreamBufferClasses
      StreamBufferClasses();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief open opens a StreamBufferInterface from filename and open_mode
      //! @param NAME name of the file which this StreamBufferInterface will be associated with
      //! @param OPEN_MODE the manner with which this StreamBufferInterface should be opened
      //!        for explanation on the types and use of open modes please see
      //!        <http://www.cplusplus.com/reference/iostream/ios_base/openmode.html>
      //! @return returns a pointer to this StreamBufferInterface if successful; otherwise returns a null pointer
      util::ShPtr< StreamBufferInterface> Open
      (
        const char *NAME, const std::ios_base::openmode OPEN_MODE
      ) const;

      //! @brief get the compression from the commandline
      //! @return StreamBufferClass that is to be used as given in commandline
      const StreamBufferClass &GetCompressionFromCommandLine() const;

      //! @brief get StreamBufferClass from file extension
      //! @param FILE_EXTENSION extension of filename like "gz"
      //! @return StreamBufferClass for given file extension
      const StreamBufferClass &GetCompressionFromExtension( const std::string &FILE_EXTENSION) const;

    }; // class StreamBufferClasses

    //! @brief construct on access function for all StreamBufferClasses
    //! @return reference to only instances of StreamBufferClasses
    BCL_API
    StreamBufferClasses &GetStreamBufferClasses();

  } // namespace io

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< io::StreamBufferInterface>, io::StreamBufferClasses>;

  } // namespace util
} // namespace bcl

#endif //BCL_IO_STREAM_BUFFER_CLASSES_H_
