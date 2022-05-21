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

#ifndef BCL_IO_OFSTREAM_H_
#define BCL_IO_OFSTREAM_H_

// include the namespace header
#include "bcl_io.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_io_stream_interface.h"

// external includes - sorted alphabetically
#include <ostream>

namespace bcl
{
  namespace io
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class OFStream
    //! @brief class is a StreamInterface which can be used for writing output to a file.
    //!
    //! @see @link example_io_ofstream.cpp @endlink
    //! @author alexanns, woetzen
    //! @date 09/15/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API OFStream :
      public StreamInterface,
      public std::ostream
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      OFStream();

      //! @brief construct form filename and open_mode
      OFStream( const char *name, const std::ios_base::openmode open_mode = std::ios::out);

      //! @brief destructor
      ~OFStream();

      //! @brief acces stream buffer
      StreamBufferInterface *rdbuf() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief open is used to open the OFStream
      void open( const char *name, const std::ios::openmode open_mode = std::ios::out);

    }; // class OFStream

  } // namespace io
} // namespace bcl

#endif // BCL_IO_OFSTREAM_H_ 
