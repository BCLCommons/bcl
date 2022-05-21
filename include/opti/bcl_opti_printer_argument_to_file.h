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

#ifndef BCL_OPTI_PRINTER_ARGUMENT_TO_FILE_H_
#define BCL_OPTI_PRINTER_ARGUMENT_TO_FILE_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_printer_default.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PrinterArgumentToFile
    //! @brief PrinterInterface implementation for printing arguments every nth iteration to a specified file
    //! @details template class PrinterArgumentToFile is a PrinterInterface. It is part of the minimizer framework and
    //!        will write out a given Argument to a file. It is mean for storing arguments / models every n th iteration
    //!        in a specified file.
    //!
    //! @see @link example_opti_printer_argument_to_file.cpp @endlink
    //! @author butkiem1
    //! @date Feb 25, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class PrinterArgumentToFile :
      public PrinterDefault< t_ArgumentType, t_ResultType>
    {

    private:

    //////////
    // data //
    //////////

      //! @brief specific path where model should be stored
      std::string m_FilenamePrefix;

      //! @brief specific path where model should be stored
      std::string m_FilenameExtension;

      //! @brief interval for printing a model every n_th iteration
      size_t m_PrintInterval;

      //! @brief flag to overwrite printed argument
      bool m_OverwritePrintedArgument;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PrinterArgumentToFile() :
        m_FilenamePrefix(),
        m_FilenameExtension(),
        m_PrintInterval(),
        m_OverwritePrintedArgument()
      {
      }

      //! @brief constructor with parameters
      //! @param OUTPUT_FILE_NAME filename for printer output file
      //! @param INTERVAL interval for which the argument should be written out to a file
      //! @param OVERWRITE flag for overwriting argument file
      //! @param WRITE_RESULT_TO_SCREEN whether to write the result to the screen
      PrinterArgumentToFile
      (
        const std::string &OUTPUT_FILE_NAME,
        const size_t &INTERVAL,
        const bool OVERWRITE = false,
        const bool &WRITE_RESULT_TO_SCREEN = false
      ) :
        PrinterDefault< t_ArgumentType, t_ResultType>
        (
          false,
          WRITE_RESULT_TO_SCREEN
        ),
        m_FilenamePrefix
        (
          io::File::SplitToPathAndFileName( OUTPUT_FILE_NAME).First()
          + io::File::RemoveLastExtension( io::File::SplitToPathAndFileName( OUTPUT_FILE_NAME).Second())
        ),
        m_FilenameExtension( io::File::GetExtensionDelimiter() + io::File::GetLastExtension( OUTPUT_FILE_NAME)),
        m_PrintInterval( INTERVAL),
        m_OverwritePrintedArgument( OVERWRITE)
      {
      }

      //! @brief Clone function
      //! @return pointer to new PrinterArgumentToFile< t_ArgumentType, typename t_ResultType>
      PrinterArgumentToFile< t_ArgumentType, t_ResultType> *Clone() const
      {
        return new PrinterArgumentToFile< t_ArgumentType, t_ResultType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief returns the class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetAlias() const
      {
        static const std::string s_alias( "File");
        return s_alias;
      }

      //! @brief return prefix of output file
      const std::string &GetFilenamePrefix() const
      {
        return m_FilenamePrefix;
      }

      //! @brief set prefix for output file
      //! @param PREFIX prefix of output file
      void SetFilenamePrefix( const std::string &PREFIX)
      {
        m_FilenamePrefix = PREFIX;
      }

      //! @brief return extension of output file
      const std::string &GetFilenameExtension() const
      {
        return m_FilenameExtension;
      }

      //! @brief set extension for output file
      //! @param EXTENSION extension of output file
      void SetFilenameExtension( const std::string &EXTENSION)
      {
        m_FilenameExtension = EXTENSION;
      }

      //! @brief return interval of printing output file
      const size_t &GetPrintInterval() const
      {
        return m_PrintInterval;
      }

      //! @brief set interval for printing out argument
      //! @param INTERVAL interval of printing out argument
      void SetPrintInterval( const size_t &INTERVAL)
      {
        m_PrintInterval = INTERVAL;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief virtual print function taking an argument an writes it out for a specific iteration
      //! @param ARGUMENT argument of interest
      //! @param ITERATION iteration number to be printed
      //! @return whether the printing was successful
      bool Print
      (
        const t_ArgumentType &ARGUMENT,
        const size_t ITERATION
      ) const
      {

        // check whether print interval is zero
        if( m_PrintInterval == 0)
        {
          BCL_MessageDbg( "PrintInterval is zero!");
          return false;
        }

        // print only every nth interval
        if( ITERATION % m_PrintInterval != 0)
        {
          return true;
        }

        // open output file
        io::OFStream output_stream;

        std::string print_filename( m_FilenamePrefix);

        // if flag is set overwrite printed argument
        if( m_OverwritePrintedArgument)
        {
          // write to a temporary file, then rename it; this minimizes the chance that the execution dies while
          // writing out the file, thus destroying the previous result
          io::DirectoryEntry print_filename_entry( print_filename + m_FilenameExtension);
          io::DirectoryEntry temp_filename( print_filename + ".tmp");
          io::File::MustOpenOFStream( output_stream, temp_filename.GetFullName());
          io::Serialize::Write( ARGUMENT, output_stream);
          io::File::CloseClearFStream( output_stream);
          BCL_Assert
          (
            temp_filename.Rename( print_filename_entry, true),
            "Could not rename file " + temp_filename.GetFullName() + " to " + print_filename_entry.GetFullName()
          );
        }
        else
        {
          io::File::MustOpenOFStream
          (
            output_stream,
            print_filename + io::File::GetExtensionDelimiter() + util::Format()( ITERATION) + m_FilenameExtension
          );
          // write model to file stream
          io::Serialize::Write( ARGUMENT, output_stream);
          io::File::CloseClearFStream( output_stream);
        }

        // return
        return true;
      }

      //! @brief prints information concerning the approximation process based on the status of the given tracker
      //! @param TRACKER holds the status of the approximation process
      void Print( const Tracker< t_ArgumentType, t_ResultType> &TRACKER) const
      {
        PrinterDefault< t_ArgumentType, t_ResultType>::Print( TRACKER);
        Print( TRACKER.GetCurrent()->First(), TRACKER.GetIteration());
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer serializer;
        serializer.SetClassDescription
        (
          "Prints approximations to a file (iteration-dependent filename) and optionally results to the screen. "
        );
        serializer.AddInitializer
        (
          "prefix",
          "Filename prefix; may include path to a directory. Files will be written out to {prefix}{iteration}.{suffix}",
          io::Serialization::GetAgent( &m_FilenamePrefix),
          "approximation"
        );
        serializer.AddInitializer
        (
          "suffix",
          "Suffix to give each file",
          io::Serialization::GetAgent( &m_FilenameExtension),
          ".model"
        );
        serializer.AddInitializer
        (
          "overwrite",
          "Set to true to allow overwriting of existing files",
          io::Serialization::GetAgent( &m_OverwritePrintedArgument),
          "False"
        );
        serializer.AddInitializer
        (
          "interval",
          "Number of iterations between writing approximations to files",
          io::Serialization::GetAgent( &m_PrintInterval),
          "1"
        );
        serializer.AddInitializer
        (
          "result",
          "True to print the score/result/objective function at each iteration to the screen",
          io::Serialization::GetAgent( &this->GetDisplayResult()),
          "True"
        );

        return serializer;
      }

    }; // template class PrinterArgumentToFile

    // instantiate s_Instance
    template< typename t_ArgumentType, typename t_ResultType>
    const util::SiPtr< const util::ObjectInterface> PrinterArgumentToFile< t_ArgumentType, t_ResultType>::s_Instance
    (
      util::Enumerated< PrintInterface< t_ArgumentType, t_ResultType> >::AddInstance( new PrinterArgumentToFile< t_ArgumentType, t_ResultType>())
    );

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_PRINTER_ARGUMENT_TO_FILE_H_
