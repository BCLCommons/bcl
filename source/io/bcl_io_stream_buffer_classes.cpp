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
#include "io/bcl_io_stream_buffer_classes.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_storage_file.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_file_stream_buffer.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace io
  {
  //////////
  // data //
  //////////

    //! @brief enum, that adds FileCompression flag to default app flags
    static const util::ShPtr< command::FlagInterface> e_FileCompressionFlag
    (
      command::GetAppDefaultFlags().AddDefaultFlag
      (
        GetStreamBufferClasses().GetFlagFileCompression(),
        command::e_Io
      )
    );

    //! @brief GetFlagFileCompression gives commandline flag to adjust default file compression
    //! @return ShPtr to a FlagInterface which is used to get the desired compression from the command line
    util::ShPtr< command::FlagInterface> &StreamBufferClasses::GetFlagFileCompression() const
    {
      static util::ShPtr< command::FlagInterface> s_flag_file_compression
      (
        new command::FlagStatic
        (
          "file_compression",
          "the type of file compression to be used for files",
          command::Parameter
          (
            "compression_type",
            "compression algorithm to be used",
            command::ParameterCheckEnumerate< StreamBufferClasses>(),
            e_Uncompressed.GetName()
          )
        )
      );

      return s_flag_file_compression;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct all StreamBufferClasses
    StreamBufferClasses::StreamBufferClasses() :
      util::Enumerate< util::ShPtr< StreamBufferInterface>, StreamBufferClasses>( false),
      e_Uncompressed( AddEnum( "Uncompressed", util::ShPtr< StreamBufferInterface>( new FileStreamBuffer())))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &StreamBufferClasses::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief open opens a StreamBufferInterface from filename and open_mode
    //! @param NAME name of the file which this StreamBufferInterface will be associated with
    //! @param OPEN_MODE the manner with which this StreamBufferInterface should be opened
    //!        for explanation on the types and use of open modes please see
    //!        <http://www.cplusplus.com/reference/iostream/ios_base/openmode.html>
    //! @return returns a pointer to this StreamBufferInterface if successful; otherwise returns a null pointer
    util::ShPtr< StreamBufferInterface> StreamBufferClasses::Open
    (
      const char *NAME, const std::ios_base::openmode OPEN_MODE
    ) const
    {
      // initialize filename and get extension
      std::string filename( NAME);
      const std::string extension( File::GetLastExtension( filename));

      util::ShPtr< StreamBufferInterface> buffer;

      // check for input or output
      if( OPEN_MODE & std::ios::in) // input - just open according to extension
      {
        buffer = GetCompressionFromExtension( extension)->HardCopy();
      }
      // output open according to extension, if uncompressed apply commandline compression and append file extension
      else if( OPEN_MODE & std::ios::out)
      {
        const StreamBufferClass compression_extension( GetCompressionFromExtension( extension));
        const StreamBufferClass commandline_compression( GetCompressionFromCommandLine());

        if( compression_extension == e_Uncompressed && commandline_compression != e_Uncompressed)
        {
          // append commandline compression file extension to filename
          filename += ".";
          filename += ( *commandline_compression)->GetDefaultFileExtension();
          buffer = commandline_compression->HardCopy();
        }
        else
        {
          buffer = compression_extension->HardCopy();
        }
      }
      else
      {
        BCL_Exit( "open mode has to be at least in or out", -1);
      }

      // check if the "OPEN_MODE" is allowed for "stream_buffer_interface"
      if( !buffer->IsValidOpenMode( OPEN_MODE))
      {
        BCL_MessageStd
        (
          buffer->GetStreamBufferClass().GetName() +
          " does not support the given open mode: " + util::Format()( OPEN_MODE)
        );
      }

      // open buffer
      buffer->Open( filename.c_str(), OPEN_MODE);

      // return the stream
      return buffer;
    }

    //! @brief get StreamBufferClass from file extension
    //! @param FILE_EXTENSION extension of filename like "gz"
    //! @return StreamBufferClass for given file extension
    const StreamBufferClass &StreamBufferClasses::GetCompressionFromExtension( const std::string &FILE_EXTENSION) const
    {
      for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        if( ( **itr)->GetDefaultFileExtension() == FILE_EXTENSION)
        {
          return *itr;
        }
      }

      // no matching compression was found, return uncompressed
      return e_Uncompressed;
    }

    //! @brief get the compression from the commandline
    //! @return StreamBufferClass that is to be used as given in commandline
    const StreamBufferClass &StreamBufferClasses::GetCompressionFromCommandLine() const
    {
      return GetEnumFromName( GetFlagFileCompression()->GetFirstParameter()->GetValue());
    }

    //! @brief construct on access function for all StreamBufferClasses
    //! @return reference to only instances of StreamBufferClasses
    StreamBufferClasses &GetStreamBufferClasses()
    {
      return StreamBufferClasses::GetEnums();
    }

  } // namespace io

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< io::StreamBufferInterface>, io::StreamBufferClasses>;

  } // namespace util
} // namespace bcl
