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

#ifndef BCL_IO_IFSTREAM_H_
#define BCL_IO_IFSTREAM_H_

// include the namespace header
#include "bcl_io.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_io_stream_interface.h"

// external includes - sorted alphabetically
#include <istream>

namespace bcl
{
  namespace io
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class IFStream
    //! @brief class is a StreamInterface which can be used for reading in input from a file.
    //!
    //! @see @link example_io_ifstream.cpp @endlink
    //! @author alexanns, woetzen
    //! @date 09/15/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API IFStream :
      public StreamInterface,
      public std::istream
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      IFStream();

      //! @brief construct form filename and open_mode
      IFStream( const char *name, const std::ios_base::openmode open_mode = std::ios::in);

      //! @brief destructor
      ~IFStream();

      //! @brief acces stream buffer
      //__streambuf_type
      StreamBufferInterface *rdbuf() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief open is used to open the IFStream
      void open( const char *name, const std::ios::openmode open_mode = std::ios::in);

    }; // class IFStream

  } // namespace io
} // namespace bcl

#endif // BCL_IO_IFSTREAM_H_ 
