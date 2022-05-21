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
#include "io/bcl_io_file.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_static.h"
#include "io/bcl_io_directory_entry.h"
#include "math/bcl_math.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically
#include <fstream>
#include <string.h>
#include <unistd.h>

namespace bcl
{
  namespace io
  {

  //////////
  // data //
  //////////

    //! @brief enum, that adds check compressed alternatives flag to default app flags
    static const util::ShPtr< command::FlagInterface> e_CheckCompressedAlternativesFlag
    (
      command::GetAppDefaultFlags().AddDefaultFlag
      (
        File::GetFlagCheckCompressedAlternatives(),
        command::e_Io
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    File::File()
    {
    }

    //! @brief static function to return the GetPostFixes of this class
    //! @return name of this class as std::string
    const std::string &File::GetPostFixes()
    {
      static const std::string s_post_fixes( "~*");
      return s_post_fixes;
    }

    //! @brief static function to return the GetExtensionDelimiter of this class
    //! @return name of this class as std::string
    const std::string &File::GetExtensionDelimiter()
    {
      static const std::string s_extension_delimiter( ".");
      return s_extension_delimiter;
    }

    //! @brief commandline flag that enables to try to open alternative input files, if a file with given name cannot be opened for reading
    //! @return the flag
    const util::ShPtr< command::FlagInterface> &File::GetFlagCheckCompressedAlternatives()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "file_compressed_alternatives",
          "set to enable the search for files that cannot be opened, but for which a file with a compression"
          " extension might exist as alternative"
        )
      );

      // end
      return s_flag;
    }

    //! @brief check for alternatives input files that are compressed as given in commandline
    //! @return true if alternative compressed files should be used
    bool File::CheckCompressedAlternatives()
    {
      return GetFlagCheckCompressedAlternatives()->GetFlag();
    }

  ///////////////
  // operators //
  ///////////////

  ////////////////
  // operations //
  ////////////////

    //! @brief count the number of slashes in given PATH name
    //! @param PATH path that possibly containing slashes
    //! @return number of slashes in path
    size_t File::CountSlashes( const std::string &PATH)
    {
      size_t slash_count( 0);
      for( std::string::const_iterator itr( PATH.begin()), itr_end( PATH.end()); itr != itr_end; ++itr)
      {
        if( *itr == '\\' || *itr == '/')
        {
          ++slash_count;
        }
      }

      // end
      return slash_count;
    }

    //! @brief list all existing compressed alternatives for a given filename
    //! @param FILENAME the filename of interest
    //! @return list of filenames that exist and have different possible compression extensions appended
    storage::Vector< std::string> File::ListCompressedAlternatives( const std::string &FILENAME)
    {
      storage::Vector< std::string> compressed_alternatives;

      // it already is a compressed file, so there should not be any alternatives
      if
      (
           GetStreamBufferClasses().GetCompressionFromExtension( GetLastExtension( FILENAME))
        != GetStreamBufferClasses().e_Uncompressed
      )
      {
        return compressed_alternatives;
      }

      // iterate over compression types
      for
      (
        StreamBufferClasses::const_iterator
          itr( GetStreamBufferClasses().Begin()), itr_end( GetStreamBufferClasses().End());
        itr != itr_end;
        ++itr
      )
      {
        const std::string compressed_filename( FILENAME + GetExtensionDelimiter() + ( **itr)->GetDefaultFileExtension());
        if( DirectoryEntry( compressed_filename).DoesExist())
        {
          compressed_alternatives.PushBack( compressed_filename);
        }
      }

      // end
      return compressed_alternatives;
    }

    //! @brief Size gives the size of a file
    //! @param FILENAME is the name of the file whose size is desired
    //! @return returns a pair of boolean (true if file size could be gotten; false if not) and size_t the file size
    storage::Pair< bool, long> File::Size( const std::string &FILENAME)
    {
      // true if file does not exist
      if( !DirectoryEntry( FILENAME).DoesExist())
      {
        // give warning
        BCL_MessageStd( FILENAME + " does not exist, so size can not be determined");

        // file does not exist so return pair of false boolean and size zero
        return storage::Pair< bool, long>( false, 0);
      }

      // create std::ifsteam "ifstream"
      // Do not use io::IFStream; it will not work if the file contains compressed data
      // because seekg does not work on the compressed streams
      std::ifstream ifstream;
      ifstream.open( FILENAME.c_str(), std::ios::binary | std::ios::in);

      // true if file is not opened
      if( !ifstream)
      {
        // give warning
        BCL_MessageStd( FILENAME + " could not be opened, so size can not be determined");

        // file cannot be opened so return pair of false boolean and size zero
        return storage::Pair< bool, long>( false, 0);
      }

      // move the get pointer to the beginning of the stream with zero offset
      // see http://www.cplusplus.com/reference/iostream/istream/seekg.html for info on seekg
      ifstream.seekg( 0, std::ios_base::beg);

      // true if seekg failed and set either the failbit or the badbit
      if( !ifstream)
      {
        // give warning
        BCL_MessageStd( "seekg could not find begin of \"ifstream\" . failbit or badbit set");

        // beginning of the stream with zero offset could not be gotten so return pair of false boolean and size zero
        return storage::Pair< bool, long>( false, 0);
      }

      // create long "begin" and initialize with result of tellg() ( tellg() gives the position of the get pointer)
      // see http://www.cplusplus.com/reference/iostream/istream/tellg.html for info on tellg
      const std::istream::pos_type begin( ifstream.tellg());

      // move the get pointer to the end of the stream with zero offset
      // see http://www.cplusplus.com/reference/iostream/istream/seekg.html for info on seekg
      ifstream.seekg( 0, std::ios_base::end);

      // true if seekg failed and set either the failbit or the badbit
      if( !ifstream)
      {
        // give warning
        BCL_MessageStd( "seekg could not find end of \"ifstream\" and set failbit or badbit");

        // end of the stream with zero offset could not be gotten so return pair of false boolean and size zero
        return storage::Pair< bool, long>( false, 0);
      }

      // create long "end" and initialize with result of tellg() ( tellg() gives the position of the get pointer)
      // see http://www.cplusplus.com/reference/iostream/istream/tellg.html for info on tellg
      const std::istream::pos_type end( ifstream.tellg());

      // close and clear "ifstream"
      CloseClearFStream( ifstream);

      // create long "file_size" and initialize with the difference between "end" and "begin"
      const long file_size( long( end - begin));

      // return boolean true indicating success and "file_size"
      return storage::Pair< bool, long>( true, file_size);
    }

    //! @brief GetFullExtension gives all extensions of a filename
    //!        (i.e. everything after the first s_ExtensionDelimiter)
    //! @param FILENAME is the filename whose extension is desired
    //! @return returns a std::string which is the full extension of FILENAME
    std::string File::GetFullExtension( const std::string &FILENAME)
    {
      // create string "cleaned_filename" and initialize with filename without the path or s_PostFixes
      const std::string cleaned_filename( File::CleanFilename( FILENAME, GetPostFixes()));

      // create size_t "dot_position"; initialize with first occurance of s_ExtensionDelimiter in "stripped_filename"
      // ignore first position in "stripped_filename" since hidden files start with s_ExtensionDelimiter
      const size_t dot_position( cleaned_filename.find_first_of( GetExtensionDelimiter(), 1));

      // true if s_ExtensionDelimiter does not occur in "FILENAME"
      if( dot_position == std::string::npos)
      {
        // there is no extension so return empty string
        return std::string();
      }

      // create size_t "extension_length" and initialize with the difference between "info_position" and "dot_position"
      const size_t extension_length( cleaned_filename.size() - dot_position);

      // create string "extension"; initialize with substring of length "extension_length"
      // use "dot_position" + 1 so that s_ExtensionDelimiter is not included in "extension"
      const std::string extension( cleaned_filename.substr( dot_position + 1, extension_length));

      // return "extension"
      return extension;
    }

    //! @brief GetLastExtension gives the last extension of a filename
    //!        (i.e. everything after the last s_ExtensionDelimiter)
    //! @param FILENAME is the filename whose last extension is desired
    //! @return returns a std::string which is the last extension of a filename
    std::string File::GetLastExtension( const std::string &FILENAME)
    {
      // create string "cleaned_filename" and initialize with filename without the path or s_PostFixes
      const std::string cleaned_filename( File::CleanFilename( FILENAME, GetPostFixes()));

      // create size_t "dot_position" and initialize with last occurance of s_ExtensionDelimiter in "stripped_filename"
      const size_t dot_position( cleaned_filename.find_last_of( GetExtensionDelimiter()));

      // true if s_ExtensionDelimiter does not occur in "FILENAME"
      // or if it occurs at the first position (indicating hidden file)
      if( dot_position == std::string::npos || dot_position == 0)
      {
        // there is no extension so return empty string
        return std::string();
      }

      // create size_t "extension_length" and initialize with the difference between "info_position" and "dot_position"
      const size_t extension_length( cleaned_filename.size() - dot_position);

      // create string "extension"; initialize with substring of length "extension_length"
      // use "dot_position" + 1 so that s_ExtensionDelimiter is not included in "extension"
      const std::string extension( cleaned_filename.substr( dot_position + 1, extension_length));

      // return "extension"
      return extension;
    }

    //! @brief CleanFilename removes the path and any s_PostFixes (e.x. '*' or '~') from a filename
    //! @param FILENAME is the filename which is desired to be cleaned
    //! @param POSTFIX the string of possible characters considered as a postfix
    //!        everything after the first occurrence of any of the characters in POSTFIX will be removed
    //! @return returns a string which is FILENAME without the path or any s_PostFixes
    std::string File::CleanFilename( const std::string &FILENAME, const std::string &POSTFIX)
    {
      // create string "pathless_filename" and initialize with filename without the path or s_PostFixes
      const std::string pathless_filename( File::RemovePath( FILENAME));

      // create string "stripped_filename" and initialize with "pathless_filename" after s_PostFixes have been removed
      const std::string cleaned_filename( File::RemoveAfterFirstOfAnyPostfix( pathless_filename, POSTFIX));

      // "cleaned_filename" has no path and everything after first occurance of anything in "POSTFIX" has been removed
      return cleaned_filename;
    }

    //! @brief RemoveLastExtension removes the last extension from a filename
    //! @param FILENAME the filename whose last extension is desired to be removed
    //! @return returns a string which is FILENAME with its last extension removed
    std::string File::RemoveLastExtension( const std::string &FILENAME)
    {
      // create size_t "filename_size" and initialize with size of "FILENAME"
      const size_t filename_size( FILENAME.size());

      //! create string "filename_last_extension" and initialize with the last exension of "FILENAME"
      const std::string filename_last_extension( File::GetLastExtension( FILENAME));

      // true if FILENAME has no extension or has a s_ExtensionDelimiter with nothing
      // (except possibly a postfix) after it. So, need to determine which case is occurring in the if statement
      if( filename_last_extension.size() == 0)
      {
        // create string "pathless_filename" and initialize with filename without the path or s_PostFixes
        const std::string pathless_filename( File::RemovePath( FILENAME));

        // create size_t "dot_position"; initialize with last occurance of GetExtensionDelimiter in "stripped_filename"
        const size_t dot_position( pathless_filename.find_last_of( GetExtensionDelimiter()));

        // true if "dot_position" is zero
        // indicating the only s_ExtensionDelimiter in "pathless_filename" makes it a hidden file
        if( dot_position == 0 || dot_position == std::string::npos)
        {
          // FILENAME really has no extension to remove so return "FILENAME"
          return FILENAME;
        }
        // "dot_position" is not a hidden file indicator
        // so there is a s_ExtensionDelimiter at the end of "FILENAME" with nothing after it
        else
        {
          // create size_t "dot_position_in_original_filename"
          // initialize with the position where the dot occurs in "FILENAME"
          // ( since "pathless_filename" was used to get "dot_position")
          const size_t dot_position_in_original_filename( dot_position + ( filename_size - pathless_filename.size()));

          // create string "last_extensionless_filename" and initialize with substring of "FILENAME" which does not
          // contain anything after the last s_ExtensionDelimiter (i.e. the last extension has been removed)
          const std::string last_extensionless_filename
          (
            FILENAME.substr( 0, filename_size - ( filename_size - ( dot_position_in_original_filename)))
          );

          // return "last_extensionless_filename"
          return last_extensionless_filename;
        }

      }

      // create size_t "extension_position"; initialize with position "filename_last_extension" occurs in "FILENAME"
      const size_t extension_position( FILENAME.rfind( filename_last_extension));

      // create string "stripped_filename" and initialize with substring of "FILENAME" not including the last extension
      const std::string stripped_filename
      (
        FILENAME.substr( 0, filename_size - ( filename_size - extension_position) - 1)
      );

      // return "stripped_filename"
      return stripped_filename;
    }

    //! @brief Remove known compression type extensions from a filename
    //! @param FILENAME is the filename for which the compression extension, if any, should be removed
    //! @return returns a std::string which lacks the compression extension
    std::string File::RemoveCompressionExtension( const std::string &FILENAME)
    {
      std::string extension( GetLastExtension( FILENAME));
      StreamBufferClass compression( GetStreamBufferClasses().GetCompressionFromExtension( extension));
      if( compression != GetStreamBufferClasses().e_Uncompressed)
      {
        return FILENAME.substr( 0, FILENAME.size() - extension.size() - 1);
      }
      return FILENAME;
    }

    //! @brief RemoveFullExtension removes the full extension from a filename
    //! @param FILENAME the filename whose full extension is desired to be removed
    //! @return returns a string which is FILENAME with its full extension removed
    std::string File::RemoveFullExtension( const std::string &FILENAME)
    {
      //! create string "current_filename" and initialize with "FILENAME"
      std::string current_filename( FILENAME);

      //! create string "new_filename" and initialize with "current_filename" minus the last extension
      std::string new_filename( File::RemoveLastExtension( current_filename));

      // true while there are still extensions being removed from "current_filename"
      while( current_filename != new_filename)
      {
        // set "current_filename" to "new_filename"
        current_filename = new_filename;

        // remove the last extension from "current_filename" and set "new_filename" to this string
        new_filename = File::RemoveLastExtension( current_filename);
      }

      // return "new_filename"
      return new_filename;
    }

    //! @brief RemoveCommonComponents removes common directories & filenames from QUERY that also exist in TEMPLATE
    //! @param QUERY the filename to find the non-common components in
    //! @param TEMPLATE another filename with a common prefix/suffix to be removed
    //! @return returns a string which is FILENAME with its full extension removed
    std::string File::RemoveCommonComponents( const std::string &QUERY, const std::string &TEMPLATE)
    {
      if( QUERY == TEMPLATE)
      {
        return std::string();
      }
      storage::Vector< std::string> components_template( util::SplitString( TEMPLATE, "/"));
      storage::Vector< std::string> components_query( util::SplitString( QUERY, "/"));
      storage::Vector< std::string>::const_iterator
        itr_template( components_template.Begin()), itr_template_end( components_template.End()),
        itr_query( components_query.Begin()), itr_query_end( components_query.End());
      size_t n_skipped_start( 0);
      // skip common components
      while( itr_template != itr_template_end && itr_query != itr_query_end && *itr_query == *itr_template)
      {
        ++itr_template;
        ++itr_query;
        ++n_skipped_start;
      }
      // if the query is a prefix of the template
      if( itr_query == itr_query_end)
      {
        return std::string();
      }

      storage::Vector< std::string>::const_reverse_iterator
        itr_rev_template( components_template.ReverseBegin()), itr_rev_template_end( components_template.ReverseEnd()),
        itr_rev_query( components_query.ReverseBegin()), itr_rev_query_end( components_query.ReverseEnd());
      size_t n_skipped_end( 0);
      // skip common components
      while
      (
        itr_rev_template != itr_rev_template_end && itr_rev_query != itr_rev_query_end
        && *itr_rev_query == *itr_rev_template
      )
      {
        ++itr_rev_template;
        ++itr_rev_query;
        ++n_skipped_end;
      }
      return util::Join( "/", storage::Vector< std::string>( itr_query, itr_query_end - n_skipped_end));
    }

    //! @brief static function to convert a Class Identifier into a filename
    //! this converts any string into a filename by replacing the namespace-"::" with "_"
    //! @param CLASS_IDENTIFIER string to convert
    //! @return string for filename
    const std::string File::ConvertClassIdentifierToFilename( const std::string &CLASS_IDENTIFIER)
    {
      // initialize filename with all "::" replaced by "_"
      std::string filename( util::ReplaceString( CLASS_IDENTIFIER, "::", "_"));
      // commas are ok in Windows, spaces are not contained in class identifiers
      // replace all "<" and ">" from templated objects
      std::replace_if( filename.begin(), filename.end(), std::bind2nd( std::equal_to< char>(), '<'), '+');
      std::replace_if( filename.begin(), filename.end(), std::bind2nd( std::equal_to< char>(), '>'), '-');
      // replace all "/" and "\" just because any string can be given to this function
      std::replace_if( filename.begin(), filename.end(), std::bind2nd( std::equal_to< char>(), '/'), '_');
      std::replace_if( filename.begin(), filename.end(), std::bind2nd( std::equal_to< char>(), '\\'), '_');

      return filename;
    }

    //! @brief FilesMatch checks if the text contents (strings) of each line of two files are the same
    //! @param FILE_NAME_A the first file which will be compared
    //! @param FILE_NAME_B the file which will be compared to "FILE_NAME_A"
    //! @return returns a bool if each line of "FILE_NAME_A" is the same as "CORRECT_FILE_NAME"
    bool File::FilesMatch( const std::string &FILE_NAME_A, const std::string &FILE_NAME_B)
    {
      // create IFStream "read"
      IFStream read_a, read_b;

      // open each file
      if( !File::TryOpenIFStream( read_a, FILE_NAME_A))
      {
        BCL_MessageStd( "Could not open " + FILE_NAME_A + " for reading");
        return false;
      }
      if( !File::TryOpenIFStream( read_b, FILE_NAME_B))
      {
        File::CloseClearFStream( read_a);
        BCL_MessageStd( "Could not open " + FILE_NAME_B + " for reading");
        return false;
      }

      // determine whether the streams are equal
      const bool files_match( StreamsMatch( read_a, read_b));

      // close each file
      File::CloseClearFStream( read_a);
      File::CloseClearFStream( read_b);

      return files_match;
    }

    //! @brief FilesMatchWithinAbsoluteTolerance is like FilesMatch, except that strings that are numerical in both
    //! files are considered equal if the numbers they represent are within an absolute tolerance
    //! @param FILE_NAME_A the first file which will be compared
    //! @param FILE_NAME_B the file which will be compared to "FILE_NAME_A"
    //! @param ABSOLUTE_TOLERANCE the absolute tolerance to use, should be ~= precision of #s in the file
    //! @return returns a bool if each line of "FILE_NAME_A" is the same as "CORRECT_FILE_NAME"
    bool File::FilesMatchWithinAbsoluteTolerance
    (
      const std::string &FILE_NAME_A,
      const std::string &FILE_NAME_B,
      const double &ABSOLUTE_TOLERANCE
    )
    {
      // create IFStream "read"
      IFStream read_a, read_b;

      // open each file
      File::MustOpenIFStream( read_a, FILE_NAME_A);
      File::MustOpenIFStream( read_b, FILE_NAME_B);

      // determine whether the streams are equal
      const bool files_match( StreamsMatchWithinAbsoluteTolerance( read_a, read_b, ABSOLUTE_TOLERANCE));

      // close each file
      File::CloseClearFStream( read_a);
      File::CloseClearFStream( read_b);

      return files_match;
    }

    //! @brief FilesMatch checks if the binary contents (bits) of each of two files are the same
    //! @param FILE_NAME_A the first file which will be compared
    //! @param FILE_NAME_B the file which will be compared to "FILE_NAME_A"
    //! @return returns a bool if each bit of "FILE_NAME_A" is the same as "FILE_NAME_B"
    bool File::BinaryFilesMatch( const std::string &FILE_NAME_A, const std::string &FILE_NAME_B)
    {
      // create IFStream "read"
      IFStream read_a, read_b;

      // open each file
      File::MustOpenIFStream( read_a, FILE_NAME_A, std::ios_base::binary);
      File::MustOpenIFStream( read_b, FILE_NAME_B, std::ios_base::binary);

      // determine whether the streams are equal
      const bool files_match( BinaryStreamsMatch( read_a, read_b));

      // close each file
      File::CloseClearFStream( read_a);
      File::CloseClearFStream( read_b);

      return files_match;
    }

    //! @brief StreamsMatch tests whether the text contents (strings) of each whitespace-delimited token of two streams are the same
    //! @param STREAM_A, STREAM_B the streams to compare
    //! @return returns a bool if each token of "STREAM_A" is the same as "STREAM_B"
    bool File::StreamsMatch( std::istream &STREAM_A, std::istream &STREAM_B)
    {
      // make temporary strings to hold the tokens of each stream
      std::string tmp_a, tmp_b;

      // iterate over both streams, stop when either stream goes bad
      while( STREAM_A.good() && STREAM_B.good())
      {
        // read the next token from each stream
        STREAM_A >> tmp_a;
        STREAM_B >> tmp_b;

        // test for equality
        if( tmp_a != tmp_b)
        {
          BCL_MessageStd
          (
            "strings from streams differed; stream a string: " + tmp_a +
            " stream b string: " + tmp_b
          );
          return false;
        }
      }

      // test whether the streams are both good, otherwise they are clearly not equal
      if( STREAM_A.good() != STREAM_B.good())
      {
        BCL_MessageStd( "Only one stream was good for reading");
        return false;
      }

      // streams were identical, return true
      return true;
    }

    //! @brief BinaryStreamsMatch tests whether the contents (bits) of each file are identical
    //! @param STREAM_A, STREAM_B the streams to compare
    //! @return returns a bool if each token of "STREAM_A" is the same as "STREAM_B"
    bool File::BinaryStreamsMatch( std::istream &STREAM_A, std::istream &STREAM_B)
    {
      // test whether the streams are both good, otherwise they are clearly not equal
      if( STREAM_A.good() != STREAM_B.good())
      {
        BCL_MessageStd( "Only one stream was good for reading");
        return false;
      }

      // if the streams are both bad, then they match
      if( STREAM_A.good() == false)
      {
        return true;
      }

      // check stream sizes
      std::streampos stream_size( 0);
      {
        // determine the size of each stream
        const std::istream::pos_type stream_a_begin( STREAM_A.tellg()), stream_b_begin( STREAM_B.tellg());

        // move to the end of the streams
        STREAM_A.seekg( std::istream::off_type( 0), std::ios_base::end);
        STREAM_B.seekg( std::istream::off_type( 0), std::ios_base::end);

        // get the end positions
        const std::istream::pos_type stream_a_end( STREAM_A.tellg()), stream_b_end( STREAM_B.tellg());

        // return the streams to their original positions
        STREAM_A.seekg( std::istream::off_type( stream_a_begin), std::ios_base::beg);
        STREAM_B.seekg( std::istream::off_type( stream_b_begin), std::ios_base::beg);

        // check the sizes
        if( ( stream_a_end - stream_a_begin) != ( stream_b_end - stream_b_begin))
        {
          BCL_MessageStd
          (
            "Stream sizes differed, " + util::Format()( ( stream_a_end - stream_a_begin))
            + " vs " + util::Format()( ( stream_b_end - stream_b_begin))
          );
          return false;
        }

        // stream sizes were identical, set the stream size
        stream_size = std::streampos( stream_a_end - stream_a_begin);
      }

      // create a buffer to read bits from the stream (4k at a time)
      static const std::streamsize buffer_size( 4096);
      static char buffer_a[ buffer_size];
      static char buffer_b[ buffer_size];

      // keep track of the remaining size
      std::streampos bytes_remaining( stream_size);
      for( ; bytes_remaining >= buffer_size; bytes_remaining -= buffer_size)
      {
        // read the next block of data
        STREAM_A.read( buffer_a, buffer_size);
        STREAM_B.read( buffer_b, buffer_size);

        // check whether the blocks differed
        if( memcmp( buffer_a, buffer_b, buffer_size) != 0)
        {
          BCL_MessageStd
          (
            "Streams differed in block "
            + util::Format()( size_t( stream_size - bytes_remaining))
            + " - "
            + util::Format()( size_t( stream_size - bytes_remaining) + size_t( buffer_size) - 1)
          );
          return false;
        }
      }

      // handle data at the end of the array
      if( bytes_remaining > 0)
      {
        // read the next block of data
        STREAM_A.read( buffer_a, bytes_remaining);
        STREAM_B.read( buffer_b, bytes_remaining);

        // check whether the blocks differed
        if( memcmp( buffer_a, buffer_b, bytes_remaining) != 0)
        {
          BCL_MessageStd
          (
            "Streams differed in block "
            + util::Format()( size_t( stream_size - bytes_remaining))
            + " - " + util::Format()( size_t( stream_size) - 1)
          );
          return false;
        }
      }

      return true;
    }

    // anonymous namespace to avoid export of symbols
    namespace
    {
      //! @brief test whether two strings that contain numeric values separated by any other characters, are equal to
      //!        within a specified tolerance
      //! @param A, B strings to compare
      //! @return 0 if the strings are equal; < 0 if there are non-numeric differences between the strings,
      //!         and the max absolute difference between numeric values between the strings otherwise
      double PartiallyNumericStringsMaxDifference( const std::string &A, const std::string &B)
      {
        if( A == B)
        {
          return 0.0;
        }
        if( A.empty() || B.empty())
        {
          return A.empty() == B.empty() ? 0.0 : util::GetUndefined< double>();
        }
        double max_abs_difference( 0.0);
        size_t pos_a( 0), pos_b( 0), sz_a( A.size()), sz_b( B.size());
        while( pos_a < sz_a && pos_b < sz_b)
        {
          // get the length of the next floating point values out of each string
          const size_t len_fp_a( util::LengthOfFloatingPointType( A, pos_a));
          const size_t len_fp_b( util::LengthOfFloatingPointType( B, pos_b));

          // test that both strings either have an floating point string at this point, or both don't
          if( bool( len_fp_a) != bool( len_fp_b))
          {
            BCL_MessageStd
            (
              "strings from streams differed because one stream had a number where the other had a string; "
              "stream a string: " + A.substr( pos_a) +
              " stream b string: " + B.substr( pos_b)
            );
            return util::GetUndefined< double>();
          }
          else if( !len_fp_a)
          {
            // compare characters
            if( A[ pos_a] != B[ pos_b])
            {
              BCL_MessageStd
              (
                "strings from streams differed on non-numeric characters; "
                "stream a string " + A + " @ index " + util::Format()( pos_a) + ": " + A.substr( pos_a, 1) +
                " stream b string " + B + " @ index " + util::Format()( pos_b) + ": " + B.substr( pos_b, 1)
              );
              return util::GetUndefined< double>();
            }
          }
          else
          {
            // strings were numeric, extract the numbers
            const double string_a_fp( util::ConvertStringToNumericalValue< double>( A.substr( pos_a, len_fp_a)));
            const double string_b_fp( util::ConvertStringToNumericalValue< double>( B.substr( pos_b, len_fp_b)));

            // update max abs difference
            max_abs_difference = std::max( max_abs_difference, math::Absolute( string_a_fp - string_b_fp));

            pos_a += len_fp_a;
            pos_b += len_fp_b;
          }
          ++pos_a;
          ++pos_b;
        }
        // skip any ending space
        while( pos_a < sz_a && ( !isprint( A[ pos_a]) || isspace( A[ pos_a])))
        {
          ++pos_a;
        }
        // skip any ending space
        while( pos_b < sz_b && ( !isprint( A[ pos_a]) || isspace( A[ pos_a])))
        {
          ++pos_b;
        }
        if( pos_a < sz_a || pos_b < sz_b)
        {
          BCL_MessageStd
          (
            "String lengths/number of numbers differed; "
            "pos/size a: " + util::Format()( pos_a) + "/" + util::Format()( sz_a) +
            "pos/size b: " + util::Format()( pos_b) + "/" + util::Format()( sz_b) +
            " string a: " + A + "\nB: " + B
          );
          return util::GetUndefined< double>();
        }
        return max_abs_difference;
      }
    }

    //! @brief StreamsMatchWithinAbsoluteTolerance is like StreamsMatch, except that strings that are numerical in
    //! both streams are considered equal if the numbers they represent are within an absolute tolerance
    //! @param FILE_NAME_A, STREAM_B the streams to compare
    //! @param ABSOLUTE_TOLERANCE the absolute tolerance to use, should be ~= precision of #s in the file
    //! @return returns true if each string of "STREAM_A" is the same as "STREAM_B" or within the given numerical tolerance
    bool File::StreamsMatchWithinAbsoluteTolerance
    (
      std::istream &STREAM_A,
      std::istream &STREAM_B,
      const double &ABSOLUTE_TOLERANCE
    )
    {
      // make temporary strings to hold the tokens of each stream
      std::string tmp_a, tmp_b;

      // track the maximum absolute difference in excess of the absolute tolerance
      double max_abs_difference( 0.0);

      // iterate over both streams, stop when either stream goes bad
      while( STREAM_A.good() && STREAM_B.good())
      {
        // read the next token from each stream
        STREAM_A >> tmp_a;
        STREAM_B >> tmp_b;

        // test for equality
        const double abs_difference( PartiallyNumericStringsMaxDifference( tmp_a, tmp_b));
        if( !util::IsDefined( abs_difference))
        {
          return false;
        }
        max_abs_difference = std::max( abs_difference, max_abs_difference);
      }

      // test whether the streams are both good, otherwise they are clearly not equal
      if( STREAM_A.good() != STREAM_B.good())
      {
        BCL_MessageStd( "Only one stream was good for reading");
        return false;
      }

      if( max_abs_difference > ABSOLUTE_TOLERANCE)
      {
        BCL_MessageStd( "Streams differed by up to: " + util::Format()( max_abs_difference));
        return false;
      }

      // streams were identical, return true
      return true;
    }

    //! @brief ObjectsMatchWithinAbsoluteTolerance is like StreamsMatchWithin Absolute Tolerance, except
    //! objects are converted into streams and then StreamMathwithinAbolute tolerance is called
    //! @param OBJECT_A, OBJECT_B the objects to compare
    //! @param ABSOLUTE_TOLERANCE the absolute tolerance to use, should be ~= precision of #s in the file
    //! @return returns true if each string of "STREAM_A" is the same as "STREAM_B" or within the given numerical tolerance
    bool File::ObjectsMatchWithinAbsoluteTolerance
    (
      const util::ObjectInterface &OBJECT_A,
      const util::ObjectInterface &OBJECT_B,
      const double &ABSOLUTE_TOLERANCE
    )
    {
      //Turn both objects into strings
      std::string first_object( util::Format()( OBJECT_A));
      std::string second_object( util::Format()( OBJECT_B));

      // Construct standard string streams
      std::stringstream first_stream( first_object);
      std::stringstream second_stream( second_object);

      // Test stream within tolerance
      return File::StreamsMatchWithinAbsoluteTolerance( first_stream, second_stream, ABSOLUTE_TOLERANCE);
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief ClearStream resets any error flags that have been set on a std::ios to goodbit
    //!        please see <http://www.cplusplus.com/reference/iostream/ios_base/iostate.html> for info on stream states
    //! @param STREAM the std::ios which will have its error flags reset
    void File::ClearStream( std::ios &STREAM)
    {
      // set all error states of "STREAM" to goodbit; goodbit is the default parameter of clear()
      STREAM.clear();
    }

    //! @brief RemovePath removes the path from a filename
    //! @param FILENAME is the filename whose path is desired to be removed
    //! @return string which is FILENAME without its path
    std::string File::RemovePath( const std::string &FILENAME)
    {
      // create size_t "path_end_position" and initialize with the last occurance of \ or /
      const size_t path_end_position( FILENAME.find_last_of( "\\/")); //< need to escape the \ character

      // true if / or \ was not found in "FILENAME"
      if( path_end_position == std::string::npos)
      {
        // "FILENAME" did not include a path so just return "FILENAME"
        return FILENAME;
      }

      // create string "stripped_filename" and initialize with entire string after "path_end_position"
      const std::string stripped_filename( FILENAME.substr( path_end_position + 1, std::string::npos));

      // return "stripped_filename"
      return stripped_filename;
    }

    //! @brief split a given filename into path and filename
    //! @param FILENAME the file name to be separated
    //! @return first is the path, second is the single filename
    storage::VectorND< 2, std::string> File::SplitToPathAndFileName( const std::string &FILENAME)
    {
      // create size_t "path_end_position" and initialize with the last occurance of \ or /
      const size_t path_end_position( FILENAME.find_last_of( "\\/")); //< need to escape the \ character

      // true if / or \ was not found in "FILENAME"
      if( path_end_position == std::string::npos)
      {
        // "FILENAME" did not include a path so just return "FILENAME" with an empty path
        return storage::VectorND< 2, std::string>( "", FILENAME);
      }
      else if( path_end_position == FILENAME.size() - 1)
      {
        // FILENAME did not include a real filename
        return storage::VectorND< 2, std::string>( FILENAME, "");
      }

      // return "stripped_filename"
      return storage::VectorND< 2, std::string>
             (
               FILENAME.substr( 0, path_end_position + 1),
               FILENAME.substr( path_end_position + 1, FILENAME.length() - path_end_position)
             );
    }

    //! @brief RemoveAfterFirstOfAnyPostfix removes everything after any character in a string is found in a filename
    //! @param FILENAME the filename which has a postfix to be removed
    //! @param POSTFIX the string of possible characters considered as a postfix
    //!        everything after the first occurrence of any of the characters in POSTFIX will be removed
    //! @return returns string which is FILENAME without its postfix
    std::string File::RemoveAfterFirstOfAnyPostfix( const std::string &FILENAME, const std::string &POSTFIX)
    {
      // create size_t "post_fix_position" and initialize to the position in FILENAME where any of the characters in
      // "POSTFIX" occurs
      const size_t post_fix_position( FILENAME.find_first_of( POSTFIX));

      // create size_t "filename_size" and initialize to the size of "FILENAME"
      const size_t filename_size( FILENAME.size());

      // create size_t "stripped_filename_length" initialize to the size "FILENAME" will be after postifix is removed
      const size_t stripped_filename_length( filename_size - ( filename_size - post_fix_position));

      // create string "stripped_filename"; initialize with the substring of "FILENAME" which does not contain postfix
      const std::string stripped_filename( FILENAME.substr( 0, stripped_filename_length));

      // return "stripped_filename" which has had the postfix removed
      return stripped_filename;
    }

    //! @brief checks if the string is an absolute path
    //! @param FILENAME the filename which is checked for being an absolute path
    //! @return true if it an absolute path ('/' as first char in linux, {drive_letter}:\ as first in windows)
    bool File::IsAbsolutePath( const std::string &FILENAME)
    {
#if defined(__linux__) || defined(__APPLE__)
      return !FILENAME.empty() && FILENAME[ 0] == '/';
#elif defined(_MSC_VER) || defined(WIN32) || defined(__MINGW32__)
      return FILENAME.length() >= 3 && std::isalpha( FILENAME[ 0]) && FILENAME[ 1] == ':' && ( FILENAME[ 2] == '\\' || FILENAME[ 2] == '/');
#else
      return false;
#endif
    }

    //! @brief makes the string into an absolute path
    //! @param FILENAME the filename which is checked for being an absolute path
    //! @return the file name as an absolute path
    std::string File::MakeAbsolutePath( const std::string &FILENAME)
    {
      if( IsAbsolutePath( FILENAME))
      {
        return FILENAME;
      }
      #if defined(_MSC_VER) || defined(WIN32) || defined(__MINGW32__)
      if( FILENAME.find( '/') != std::string::npos)
      {
        std::string filename_copy( FILENAME);
        std::replace( filename_copy.begin(), filename_copy.end(), '/', '\\');
        return MakeAbsolutePath( filename_copy);
      }
      const char true_path_separator( '\\');
      #else
      const char true_path_separator( '/');
      #endif

      // get the current working directory
      char temp[ FILENAME_MAX];
      std::string cwd( getcwd( temp, FILENAME_MAX) ? std::string( temp) : std::string( ""));

      if( cwd[ cwd.size() - 1] != true_path_separator)
      {
        cwd += true_path_separator;
      }

      // determine whether the current working directory was requested
      if( FILENAME.empty() || FILENAME == ".")
      {
        return cwd;
      }

      // test other variants of the current working directory
      if
      (
        FILENAME.size() == size_t( 2)
        && FILENAME[ 0] == '.'
        && FILENAME[ 0] == true_path_separator
      )
      {
        return cwd;
      }

      // determine how many characters to trim from the original filename
      size_t chars_to_trim( 0);
      if( FILENAME[ 0] == '.' && FILENAME[ 1] == true_path_separator)
      {
        chars_to_trim = 2;
      }
      return cwd + FILENAME.substr( chars_to_trim);
    }

  } // namespace io
} // namespace bcl
