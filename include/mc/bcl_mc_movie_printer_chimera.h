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

#ifndef BCL_MC_MOVIE_PRINTER_CHIMERA_H_
#define BCL_MC_MOVIE_PRINTER_CHIMERA_H_

// include the namespace header
#include "bcl_mc.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_mc_movie_printer_interface.h"
#include "storage/bcl_storage_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoviePrinterChimera
    //! @brief creates a movie script for chimera
    //!
    //! @see @link example_mc_movie_printer_chimera.cpp @endlink
    //! @author woetzen
    //! @date Jul 21, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoviePrinterChimera :
      public MoviePrinterInterface
    {

    private:

    //////////
    // data //
    //////////

      std::string m_FirstScript; //!< first portion of the script
      std::string m_MovieScript; //!< the actual movie script with all frames

      std::string m_Path;       //!< path for file
      std::string m_FilePrefix; //!< prefix for file

      storage::Vector< std::string> m_LabelColNames; //!< column names
      storage::Vector< std::string> m_LabelRowNames; //!< row names

      bool  m_RayTracing; //!< enable ray tracing -true by default
      size_t m_Width;     //!< widht of movie frames
      size_t m_Height;    //!< height of movie frames

      static const size_t s_WaitTimes[]; //!< wait times depending on step status

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MoviePrinterChimera();

      //! @brief Clone function
      //! @return pointer to new MoviePrinterInterface
      MoviePrinterChimera *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns defautl extension for movie script file
      //! @return extension for script files
      const std::string &GetScriptFileExtension() const;

      //! @brief initialize with all required information
      //! @param PREFIX prefix as absolute path
      //! @param COL_NAMES the column names to be considered
      //! @param ROW_NAMES the row name to be considered
      //! @param WIDTH width in pixel
      //! @param HEIGHT height in pixel
      //! @param ENABLE_RAY true or false
      void Initialize
      (
        const std::string &PREFIX,
        const storage::Vector< std::string> &COL_NAMES,
        const storage::Vector< std::string> &ROW_NAMES,
        const size_t WIDTH, const size_t HEIGHT,
        const bool ENABLE_RAY
      );

    ////////////////
    // operations //
    ////////////////

      //! @brief reset
      void Reset( const std::string &PREFIX);

      //! @brief add start frame
      //! @param FILE_NAME the filename of the starting frame
      //! @param KEEP_FRAME_UP is start frame permanent (like a density map)
      void SetStartFrame
      (
        const std::string &FILE_NAME,
        const bool KEEP_FRAME_UP
      );

      //! @brief add frame to movie
      //! @param FILE_NAME the filename of the pdb the model was written to
      //! @param STEP_STATUS the mc step status
      //! @param LABEL add label to frame
      void AddFrame
      (
        const std::string &FILE_NAME,
        const opti::StepStatus &STEP_STATUS,
        const storage::Table< double> &LABEL
      );

      //! @brief add final frame to movie
      //! @param FILE_NAME the filename of the pdb the model was written to
      //! @param LABEL add label to frame
      void AddFinalFrame
      (
        const std::string &FILE_NAME,
        const storage::Table< double> &LABEL
      );

      //! @brief set the watermark
      //! @param WATERMARK the watermark
      void SetWaterMark( const std::string &WATERMARK);

      //! @brief create ffmpeg command line
      //! @return commandline to be called with ffmpeg to generate the movie from the generated pngs
      std::string FFMpegCommandLine() const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write the movie script to a file
      //! @param OSTREAM
      std::ostream &WriteScript
      (
        std::ostream &OSTREAM
      ) const;

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief convert table to label
      //! @param TABLE
      //! @return string
      std::string TableToLabel( const storage::Table< double> &TABLE) const;

    }; // class MoviePrinterInterface

  } // namespace mc
} // namespace bcl

#endif // BCL_MC_MOVIE_PRINTER_CHIMERA_H_
