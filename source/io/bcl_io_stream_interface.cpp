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
#include "io/bcl_io_ofstream.h"

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
    StreamInterface::StreamInterface() :
      std::ios(),
      m_Buffer( GetStreamBufferClasses().e_Uncompressed->HardCopy())
    {
    }

    //! @brief construct form filename and open_mode
    StreamInterface::StreamInterface( const char *name, const std::ios_base::openmode open_mode) :
      std::ios(),
      m_Buffer( GetStreamBufferClasses().Open( name, open_mode))
    {
    }

    //! @brief destructor
    StreamInterface::~StreamInterface()
    {
    }

    //! @brief return the compression for that stream
    const StreamBufferClass &StreamInterface::GetCompression() const
    {
      return m_Buffer->GetStreamBufferClass();
    }

  ///////////////
  // operators //
  ///////////////

  ////////////////
  // operations //
  ////////////////

    //! @brief is_open indicates whether a StreamInterface is currently associated with a file
    //! @return boolean true if StreamInterface is open and false if StreamInterface is not open
    bool StreamInterface::is_open()
    {
      return m_Buffer->IsOpen();
    }

    //! @brief close dissociates this StreamInterface from the file it is currently bound to
    //! @return returns a pointer to this StreamInterface if successful; otherwise returns a null pointer
    void StreamInterface::close()
    {
      if( m_Buffer->IsOpen())
      {
        if( !m_Buffer->Close())
        {
          clear( rdstate() | std::ios::badbit);
        }
      }
    }

    //! @brief rdbuf returns a pointer to the buffer object bound to the StreamInterface
    //! @return returns a StreamBufferInterface which is associated with the StreamInterface
    StreamBufferInterface *StreamInterface::rdbuf() const
    {
      return const_cast< StreamBufferInterface *>( &( *m_Buffer));
    }

    //! @brief open opens a StreamInterface from filename and open_mode
    //! @param NAME name of the file which this StreamInterface will be associated with
    //! @param OPEN_MODE the manner with which this StreamInterface should be opened
    //!        for explanation on the types and use of open modes please see
    //!        <http://www.cplusplus.com/reference/iostream/ios_base/openmode.html>
    //! @return returns a pointer to this StreamInterface if successful; otherwise returns a null pointer
    void StreamInterface::open
    (
      const char *NAME, const std::ios_base::openmode OPEN_MODE
    )
    {
      // close and clear the stream, if open
      m_Buffer->Close();

      // set m_Buffer
      m_Buffer      = GetStreamBufferClasses().Open( NAME, OPEN_MODE);
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace io
} // namespace bcl
