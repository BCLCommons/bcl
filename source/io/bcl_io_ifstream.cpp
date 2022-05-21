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
#include "io/bcl_io_ifstream.h"

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
    IFStream::IFStream() :
      StreamInterface(),
      std::istream( StreamInterface::m_Buffer.IsDefined() ? &( *StreamInterface::m_Buffer) : nullptr)
    {
    }

    //! @brief construct from filename and open_mode
    IFStream::IFStream( const char *name, const std::ios_base::openmode open_mode) :
      StreamInterface(),
      std::istream( StreamInterface::m_Buffer.IsDefined() ? &( *StreamInterface::m_Buffer) : nullptr)
    {
      open( name, open_mode | std::ios::in);
    }

    //! @brief destructor
    IFStream::~IFStream()
    {
      if( is_open())
      {
        close();
      }
    }

    //! @brief acces stream buffer
    //__streambuf_type
    StreamBufferInterface *IFStream::rdbuf() const
    {
      return StreamInterface::rdbuf();
    }

  ///////////////
  // operators //
  ///////////////

  ////////////////
  // operations //
  ////////////////

    //! @brief open is used to open the IFStream
    void IFStream::open( const char *name, const std::ios::openmode open_mode)
    {
      // call base class open
      StreamInterface::open( name, open_mode | std::ios::in);

      // set the stream buffer
      std::istream::rdbuf( &( *m_Buffer));
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace io
} // namespace bcl
