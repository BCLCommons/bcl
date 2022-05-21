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

#ifndef BCL_IO_FILE_H_
#define BCL_IO_FILE_H_

// include the namespace header
#include "bcl_io.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_io_ifstream.h"
#include "bcl_io_ofstream.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_runtime_environment_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace io
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class File
    //! @brief File is a collection of functions which are useful when working with files and streams
    //!
    //! @see @link example_io_file.cpp @endlink
    //! @author alexanns
    //! @date 09/10/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API File
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    private:

      //! @brief default constructor
      File();

    public:

      //! @brief static function to return the GetPostFixes of this class
      //! @return name of this class as std::string
      static const std::string &GetPostFixes();

      //! @brief static function to return the GetExtensionDelimiter of this class
      //! @return name of this class as std::string
      static const std::string &GetExtensionDelimiter();

      //! @brief commandline flag that enables to try to open alternative input files, if a file with given name cannot be opened for reading
      //! @return the flag
      static const util::ShPtr< command::FlagInterface> &GetFlagCheckCompressedAlternatives();

      //! @brief check for alternatives input files that are compressed as given in commandline
      //! @return true if alternative compressed files should be used
      static bool CheckCompressedAlternatives();

    ////////////////
    // operations //
    ////////////////

      //! @brief count the number of slashes in given PATH name
      //! @param PATH path that possibly containing slashes
      //! @return number of slashes in path
      static size_t CountSlashes( const std::string &PATH);

      //! @brief list all existing compressed alternatives for a given filename
      //! @param FILENAME the filename of interest
      //! @return list of filenames that exist and have different possible compression extensions appended
      static storage::Vector< std::string> ListCompressedAlternatives( const std::string &FILENAME);

      //! @brief TryOpenFstream tries to open a stream to a file; if it cannot it just returns false
      //! @param FSTREAM the file stream to be used
      //! @param FILENAME the file which will be attempted to be opened
      //! @param OPEN_MODE the manner in which FILENAME should be opened
      //!        for explanation on the types and use of open modes please see
      //!        <http://www.cplusplus.com/reference/iostream/ios_base/openmode.html>
      //! @return boolean: true if FILENAME can be and is bound to FSTREAM; false if not
      template< typename t_Fstream>
      static bool TryOpenFstream
      (
        t_Fstream &FSTREAM,
        const std::string &FILENAME,
        const std::ios_base::openmode OPEN_MODE
      )
      {
        // resolve filename from runtime environmenet
        const std::string resolved_filename
        (
          util::GetRuntimeEnvironment().ResolveFileName( FILENAME)
        );

        // if resolve unsuccessful
        if( resolved_filename.empty())
        {
          // report error
          BCL_MessageCrt( "unable to resolve filename in environment: " + FILENAME)

          // could not resolve FILENAME in the environment -> return failure
          return false;
        }

        // bind "FILENAME" to "FSTREAM"
        File::OpenFstream( FSTREAM, resolved_filename, OPEN_MODE);

        // check that "FILENAME" could be associated with "FSTREAM" and neither failbit or badbit is set
        // please see <http://www.cplusplus.com/reference/iostream/ios_base/iostate.html> for explanation about
        // failbit, badbit and other stream states
        if( FSTREAM && FSTREAM.is_open())
        {
          return true;
        }

        // "FILENAME" could not be bound to "FSTREAM" or stream state of "FSTREAM" is failed/bad
        return false;
      }

      //! @brief TryOpenOFStream tries to open a stream to a file; if it cannot it just returns false
      //! @param OFSTREAM the file stream to be used
      //! @param FILENAME the file which will be attempted to be opened
      //! @param OPEN_MODE the manner in which FILENAME should be opened
      //!        for explanation on the types and use of open modes please see
      //!        <http://www.cplusplus.com/reference/iostream/ios_base/openmode.html>
      //! @return boolean: true if FILENAME can be and is bound to FSTREAM; false if not
      template< typename t_Fstream>
      static bool TryOpenOFStream
      (
        t_Fstream &OFSTREAM,
        const std::string &FILENAME,
        const std::ios_base::openmode OPEN_MODE = std::ios::out
      )
      {
        return TryOpenFstream( OFSTREAM, FILENAME, OPEN_MODE | std::ios::out);
      }

      //! @brief TryOpenIFStream tries to open a stream to a file; if it cannot it just returns false
      //! @param IFSTREAM the file stream to be used
      //! @param FILENAME the file which will be attempted to be opened
      //! @param OPEN_MODE the manner in which FILENAME should be opened
      //!        for explanation on the types and use of open modes please see
      //!        <http://www.cplusplus.com/reference/iostream/ios_base/openmode.html>
      //! @return boolean: true if FILENAME can be and is bound to FSTREAM; false if not
      template< typename t_Fstream>
      static bool TryOpenIFStream
      (
        t_Fstream &IFSTREAM,
        const std::string &FILENAME,
        const std::ios_base::openmode OPEN_MODE = std::ios::in
      )
      {
        // convert all '/' into the system's path separator
        std::string filename( FILENAME);
        if( PATH_SEPARATOR != '/')
        {
          std::replace( filename.begin(), filename.end(), '/', PATH_SEPARATOR);
        }

        if( TryOpenFstream( IFSTREAM, filename, OPEN_MODE | std::ios::in))
        {
          return true;
        }

        // no compressed alternatives should be checked
        if( !CheckCompressedAlternatives())
        {
          return false;
        }

        BCL_MessageStd( "check if there are compressed files for the file: " + FILENAME);

        // list all filenames of compressed files for that filename
        const storage::Vector< std::string> different_compression_files( ListCompressedAlternatives( FILENAME));
        if( different_compression_files.IsEmpty())
        {
          return false;
        }

        // iterate over availabe compression types till opening was successfull
        for
        (
          storage::Vector< std::string>::const_iterator
            itr( different_compression_files.Begin()), itr_end( different_compression_files.End());
          itr != itr_end;
          ++itr
        )
        {
          // check if this itr is a filename that can be opened
          if( TryOpenFstream( IFSTREAM, *itr, OPEN_MODE | std::ios::in))
          {
            BCL_MessageStd( "compressed alternative was found and opened: " + *itr);
            return true;
          }
        }

        // no alternative was found
        return false;
      }

      //! @brief MustOpenFstream tries to open a stream to a file to write to.
      //!        If it cannot, then it gives an bcl::Assert to die.
      //! @param OFSTREAM the file stream to be used
      //! @param FILENAME the file which will be attempted to be opened
      //! @param OPEN_MODE the manner in which FILENAME should be opened
      //!        for explanation on the types and use of open modes please see
      //!        <http://www.cplusplus.com/reference/iostream/ios_base/openmode.html>
      template< typename t_Fstream>
      static void MustOpenOFStream
      (
        t_Fstream &OFSTREAM,
        const std::string &FILENAME,
        const std::ios_base::openmode OPEN_MODE = std::ios::out
      )
      {
        // make sure that "FILENAME" can be bound to "OFSTREAM"
        BCL_Assert( TryOpenOFStream( OFSTREAM, FILENAME, OPEN_MODE | std::ios::out), "Could not open " + FILENAME + " for writing");
      }

      //! @brief MustOpenFstream tries to open a stream to a file to read from.
      //!        If it cannot, then it gives an bcl::Assert to die.
      //! @param IFSTREAM the file stream to be used
      //! @param FILENAME the file which will be attempted to be opened
      //! @param OPEN_MODE the manner in which FILENAME should be opened
      //!        for explanation on the types and use of open modes please see
      //!        <http://www.cplusplus.com/reference/iostream/ios_base/openmode.html>
      template< typename t_Fstream>
      static void MustOpenIFStream
      (
        t_Fstream &IFSTREAM,
        const std::string &FILENAME,
        const std::ios_base::openmode OPEN_MODE = std::ios::in
      )
      {
        // make sure that "FILENAME" can be bound to "IFSTREAM"
        BCL_Assert( TryOpenIFStream( IFSTREAM, FILENAME, OPEN_MODE | std::ios::in), "Could not open " + FILENAME + " for reading");
      }

      //! @brief TryOpenStreamBuf tries to open a streambuf to a filename; if it cannot it just returns false
      //! @param STREAM_BUF the std::streambuf to be used
      //! @param FILENAME the file which will be attempted to be opened
      //! @param OPEN_MODE the manner in which FILENAME should be opened
      //!        for explanation on the types and use of open modes please see
      //!        <http://www.cplusplus.com/reference/iostream/ios_base/openmode.html>
      //! @return boolean: true if FILENAME can be and is bound to STREAM_BUF; false if not
      template< typename t_StreamBufType>
      static bool TryOpenStreamBuf
      (
        t_StreamBufType &STREAM_BUF, const std::string &FILENAME, const std::ios_base::openmode OPEN_MODE
      )
      {
        // open the stream buffer
        STREAM_BUF.open( FILENAME.c_str(), OPEN_MODE);

        // if the buffer is open, return true
        if( STREAM_BUF && STREAM_BUF.is_open())
        {
          return true;
        }

        // else false
        return false;
      }

      //! @brief Size gives the size of a file
      //! @param FILENAME is the name of the file whose size is desired
      //! @return returns a pair of boolean (true if file size could be gotten; false if not) and long the file size
      static storage::Pair< bool, long> Size( const std::string &FILENAME);

      //! @brief GetFullExtension gives all extensions of a filename
      //!        (i.e. everything after the first s_ExtensionDelimiter)
      //! @param FILENAME is the filename whose extension is desired
      //! @return returns a std::string which is the full extension of FILENAME
      static std::string GetFullExtension( const std::string &FILENAME);

      //! @brief GetLastExtension gives the last extension of a filename
      //!        (i.e. everything after the last s_ExtensionDelimiter)
      //! @param FILENAME is the filename whose last extension is desired
      //! @return returns a std::string which is the last extension of a filename
      static std::string GetLastExtension( const std::string &FILENAME);

      //! @brief Remove known compression type extensions from a filename
      //! @param FILENAME is the filename for which the compression extension, if any, should be removed
      //! @return returns a std::string which lacks the compression extension
      static std::string RemoveCompressionExtension( const std::string &FILENAME);

      //! @brief CleanFilename removes the path and any postfixes (e.x. '*' or '~') from a filename
      //! @param FILENAME is the filename which is desired to be cleaned
      //! @param POSTFIX the string of possible characters considered as a postfix
      //!        everything after the first occurrence of any of the characters in POSTFIX will be removed
      //! @return returns a string which is FILENAME without the path or any postfixes
      static std::string CleanFilename( const std::string &FILENAME, const std::string &POSTFIX);

      //! @brief RemoveLastExtension removes the last extension from a filename
      //! @param FILENAME the filename whose last extension is desired to be removed
      //! @return returns a string which is FILENAME with its last extension removed
      static std::string RemoveLastExtension( const std::string &FILENAME);

      //! @brief RemoveFullExtension removes the full extension from a filename
      //! @param FILENAME the filename whose full extension is desired to be removed
      //! @return returns a string which is FILENAME with its full extension removed
      static std::string RemoveFullExtension( const std::string &FILENAME);

      //! @brief RemoveCommonComponents removes common directories & filenames from QUERY that also exist in TEMPLATE
      //! @param QUERY the filename to find the non-common components in
      //! @param TEMPLATE another filename with a common prefix/suffix to be removed
      //! @return returns a string which is FILENAME with its full extension removed
      static std::string RemoveCommonComponents( const std::string &QUERY, const std::string &TEMPLATE);

      //! @brief static function to convert a Class Identifier into a filename
      //! this converts any string into a filename by replacing the namespace-"::", template-"<>", "/" and "\"
      //! @param CLASS_IDENTIFIER string to convert
      //! @return string for filename
      static const std::string ConvertClassIdentifierToFilename( const std::string &CLASS_IDENTIFIER);

      //! @brief FilesMatch checks if the text contents (strings) of each line of two files are the same
      //! @param FILE_NAME_A the first file which will be compared
      //! @param FILE_NAME_B the file which will be compared to "FILE_NAME_A"
      //! @return returns a bool if each line of "FILE_NAME_A" is the same as "CORRECT_FILE_NAME"
      static bool FilesMatch( const std::string &FILE_NAME_A, const std::string &FILE_NAME_B);

      //! @brief FilesMatchWithinAbsoluteTolerance is like FilesMatch, except that strings that are numerical in both
      //! files are considered equal if the numbers they represent are within an absolute tolerance
      //! @param FILE_NAME_A the first file which will be compared
      //! @param FILE_NAME_B the file which will be compared to "FILE_NAME_A"
      //! @param ABSOLUTE_TOLERANCE the absolute tolerance to use, should be ~= precision of #s in the file
      //! @return returns a bool if each line of "FILE_NAME_A" is the same as "CORRECT_FILE_NAME"
      static bool FilesMatchWithinAbsoluteTolerance
      (
        const std::string &FILE_NAME_A,
        const std::string &FILE_NAME_B,
        const double &ABSOLUTE_TOLERANCE
      );

      //! @brief FilesMatch checks if the binary contents (bits) of each of two files are the same
      //! @param FILE_NAME_A the first file which will be compared
      //! @param FILE_NAME_B the file which will be compared to "FILE_NAME_A"
      //! @return returns a bool if each bit of "FILE_NAME_A" is the same as "FILE_NAME_B"
      static bool BinaryFilesMatch( const std::string &FILE_NAME_A, const std::string &FILE_NAME_B);

      //! @brief StreamsMatch tests whether the text contents (strings) of each whitespace-delimited token of two streams are the same
      //! @param STREAM_A, STREAM_B the streams to compare
      //! @return returns a bool if each token of "STREAM_A" is the same as "STREAM_B"
      static bool StreamsMatch( std::istream &STREAM_A, std::istream &STREAM_B);

      //! @brief StreamsMatchWithinAbsoluteTolerance is like StreamsMatch, except that strings that are numerical in
      //! both streams are considered equal if the numbers they represent are within an absolute tolerance
      //! @param FILE_NAME_A, STREAM_B the streams to compare
      //! @param ABSOLUTE_TOLERANCE the absolute tolerance to use, should be ~= precision of #s in the file
      //! @return returns true if each string of "STREAM_A" is the same as "STREAM_B" or within the given numerical tolerance
      static bool StreamsMatchWithinAbsoluteTolerance
      (
        std::istream &STREAM_A,
        std::istream &STREAM_B,
        const double &ABSOLUTE_TOLERANCE
      );

      //! @brief BinaryStreamsMatch tests whether the contents (bits) of each file are identical
      //! @param STREAM_A, STREAM_B the streams to compare
      //! @return returns a bool if each token of "STREAM_A" is the same as "STREAM_B"
      static bool BinaryStreamsMatch( std::istream &STREAM_A, std::istream &STREAM_B);

      //! @brief ObjectsMatchWithinAbsoluteTolerance is like StreamsMatchWithin Absolute Tolerance, except
      //! objects are converted into streams and then StreamMathwithinAbolute tolerance is called
      //! @param OBJECT_A, OBJECT_B the objects to compare
      //! @param ABSOLUTE_TOLERANCE the absolute tolerance to use, should be ~= precision of #s in the file
      //! @return returns true if each string of "STREAM_A" is the same as "STREAM_B" or within the given numerical tolerance
      static bool ObjectsMatchWithinAbsoluteTolerance
      (
        const util::ObjectInterface &OBJECT_A,
        const util::ObjectInterface &OBJECT_B,
        const double &ABSOLUTE_TOLERANCE
      );

    //////////////////////
    // input and output //
    //////////////////////

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief OpenFstream tries to bind an file stream to a file; it makes no checks or guarantees about success
      //! @param FSTREAM the file steam to be used
      //! @param FILENAME the file which will be attempted to be opened
      //! @param OPEN_MODE the manner in which FILENAME should be opened
      //!        for explanation on the types and use of open modes please see
      //!        <http://www.cplusplus.com/reference/iostream/ios_base/openmode.html>
      template< typename t_Fstream>
      static void OpenFstream
      (
        t_Fstream &FSTREAM,
        const std::string &FILENAME,
        const std::ios_base::openmode OPEN_MODE
      )
      {
        // make sure "FSTREAM" is not bound to something already
        File::CloseClearFStream( FSTREAM);

        // attempt to bind "FILENAME" to "FSTREAM"
        FSTREAM.open( FILENAME.c_str(), OPEN_MODE);
      }

      //! @brief CloseClearFStream dissociates an file stream from its file and resets any set error flags
      //! @param FSTREAM the file stream to be dissociated from the file it is bound to
      template< typename t_Fstream>
      static void CloseClearFStream( t_Fstream &FSTREAM)
      {
        // close the file associated with "FSTREAM" so that it is unbound from "FSTREAM"
        FSTREAM.close();

        // reset all the error state flags of "FSTREAM"
        File::ClearStream( FSTREAM);
      }

      //! @brief ClearStream resets any error flags that have been set on a std::ios to goodbit
      //!        see <http://www.cplusplus.com/reference/iostream/ios_base/iostate.html> for info on stream states
      //! @param STREAM the std::ios which will have its error flags reset
      static void ClearStream( std::ios &STREAM);

      //! @brief RemovePath removes the path from a filename
      //! @param FILENAME is the filename whose path is desired to be removed
      //! @return string which is FILENAME without its path
      static std::string RemovePath( const std::string &FILENAME);

      //! @brief split a given filename into path and filename
      //! @param FILENAME the file name to be separated
      //! @return first is the path, second is the single filename
      static storage::VectorND< 2, std::string> SplitToPathAndFileName( const std::string &FILENAME);

      //! @brief RemoveAfterFirstOfAnyPostfix removes everything after any character in string is found in a filename
      //! @param FILENAME the filename which has a postfix to be removed
      //! @param POSTFIX the string of possible characters considered as a postfix
      //!        everything after the first occurrence of any of the characters in POSTFIX will be removed
      //! @return returns string which is FILENAME without its postfix
      static std::string RemoveAfterFirstOfAnyPostfix( const std::string &FILENAME, const std::string &POSTFIX);

      //! @brief checks if the string is an absolute path
      //! @param FILENAME the filename which is checked for being an absolute path
      //! @return true if it an absolute path ('/' as first char in linux, {drive_letter}:\ as first in windows)
      static bool IsAbsolutePath( const std::string &FILENAME);

      //! @brief makes the string into an absolute path
      //! @param FILENAME the filename which is checked for being an absolute path
      //! @return the file name as an absolute path
      static std::string MakeAbsolutePath( const std::string &FILENAME);

    }; // class File

  } // namespace io
} // namespace bcl

#endif // BCL_IO_FILE_H_
