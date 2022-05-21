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

#ifndef BCL_IO_STREAM_INTERFACE_H_
#define BCL_IO_STREAM_INTERFACE_H_

// include the namespace header
#include "bcl_io.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_io_stream_buffer_interface.h"

// external includes - sorted alphabetically
#include <ios>

namespace bcl
{
  namespace io
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class StreamInterface
    //! @brief class defines the basic functionalities that must be provided by any Stream class
    //!
    //! @remarks example unnecessary
    //! @author alexanns, woetzen
    //! @date 09/20/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API StreamInterface :
      public virtual std::ios
    {

    //////////
    // data //
    //////////

    protected:

      // ShPtr to a StreamBufferInterface which determines the compression
      util::ShPtr< StreamBufferInterface> m_Buffer;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      StreamInterface();

      //! @brief construct form filename and MODE
      StreamInterface( const char *name, const std::ios_base::openmode MODE);

      //! @brief destructor
      virtual ~StreamInterface();

      //! @brief return the compression for that stream
      const StreamBufferClass &GetCompression() const;

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      //! @brief open opens a StreamInterface from filename and open_mode
      //! @param NAME name of the file which this StreamInterface will be associated with
      //! @param OPEN_MODE the manner with which this StreamInterface should be opened
      //!        for explanation on the types and use of open modes please see
      //!        <http://www.cplusplus.com/reference/iostream/ios_base/openmode.html>
      //! @return returns a pointer to this StreamInterface if successful; otherwise returns a null pointer
      virtual void open
      (
        const char *NAME, const std::ios_base::openmode OPEN_MODE
      );

      //! @brief is_open indicates whether a StreamInterface is currently associated with a file
      //! @return boolean true if StreamInterface is open and false if StreamInterface is not open
      bool is_open();

      //! @brief close dissociates this StreamInterface from the file it is currently bound to
      //! @return returns a pointer to this StreamInterface if successful; otherwise returns a null pointer
      void close();

      //! @brief rdbuf returns a pointer to the buffer object bound to the StreamInterface
      //! @return returns a StreamBufferInterface which is associated with the StreamInterface
      StreamBufferInterface *rdbuf() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

    private:

    }; // class StreamInterface

  } // namespace io
} // namespace bcl

#endif // BCL_IO_STREAM_INTERFACE_H_
