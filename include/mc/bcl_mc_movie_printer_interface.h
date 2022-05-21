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

#ifndef BCL_MC_MOVIE_PRINTER_INTERFACE_H_
#define BCL_MC_MOVIE_PRINTER_INTERFACE_H_

// include the namespace header
#include "bcl_mc.h"

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "opti/bcl_opti_step_status.h"
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoviePrinterInterface
    //! @brief creates a movie script for programs like chimera, pymol, vmd ...
    //!
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date Jul 21, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoviePrinterInterface :
      public util::ObjectInterface
    {

    public:

    //////////
    // data //
    //////////

      //! wait time for final frame
      static size_t GetFinalWaitTime()
      {
        return 50;
      }

      //! @brief get the default bcl watermark
      static const std::string &GetDefaultWatermark();

      //! @brief return command line flag for printing the script necessary to create a movie
      //! @return command line flag for printing the script necessary to create a movie
      static util::ShPtr< command::FlagInterface> &GetFlagMoviePrinter();

      //! @brief return command line parameter for specifying the prefix for where movie files are generated
      //! @return command line parameter for specifying the prefix for where movie files are generated
      static util::ShPtr< command::ParameterInterface> &GetParameterMoviePrefix();

      //! @brief return command line parameter for specifying the type of movie printer that should be used
      //! @return command line parameter for specifying the type of movie printer that should be used
      static util::ShPtr< command::ParameterInterface> &GetParameterMoviePrinterType();

      //! @brief return command line parameter for specifying the width resolution of the movie
      //! @return command line parameter for specifying the width resolution of the movie
      static util::ShPtr< command::ParameterInterface> &GetParameterMovieWidth();

      //! @brief return command line parameter for specifying the height resolution of the movie
      //! @return command line parameter for specifying the height resolution of the movie
      static util::ShPtr< command::ParameterInterface> &GetParameterMovieHeight();

      //! @brief return command line parameter for specifying whether or not all movie frames should be ray traced
      //! @return command line parameter for specifying whether or not all movie frames should be ray traced
      static util::ShPtr< command::ParameterInterface> &GetParameterRayTrace();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new MoviePrinterInterface
      virtual MoviePrinterInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns default extension for movie script file
      //! @return extension for script files
      virtual const std::string &GetScriptFileExtension() const = 0;

      //! @brief initialize with all required information
      //! @param PREFIX prefix as absolute path
      //! @param COL_NAMES the column names to be considered
      //! @param ROW_NAMES the row name to be considered
      //! @param WIDTH width in pixel
      //! @param HEIGHT height in pixel
      //! @param ENABLE_RAY true or false
      virtual void Initialize
      (
        const std::string &PREFIX,
        const storage::Vector< std::string> &COL_NAMES,
        const storage::Vector< std::string> &ROW_NAMES,
        const size_t WIDTH, const size_t HEIGHT,
        const bool ENABLE_RAY
      ) = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief reset
      //! @param PREFIX prefix for the name of the movie file
      virtual void Reset( const std::string &PREFIX) = 0;

      //! @brief add start frame
      //! @param FILE_NAME the filename of the starting frame
      //! @param KEEP_FRAME_UP is start frame permanent (like a density map)
      virtual void SetStartFrame
      (
        const std::string &FILE_NAME,
        const bool KEEP_FRAME_UP
      ) = 0;

      //! @brief add frame to movie
      //! @param FILE_NAME the filename of the pdb the model was written to
      //! @param STEP_STATUS the mc step status
      //! @param LABEL add label to frame
      virtual void AddFrame
      (
        const std::string &FILE_NAME,
        const opti::StepStatus &STEP_STATUS,
        const storage::Table< double> &LABEL
      ) = 0;

      //! @brief add final frame to movie
      //! @param FILE_NAME the filename of the pdb the model was written to
      //! @param LABEL add label to frame
      virtual void AddFinalFrame
      (
        const std::string &FILE_NAME,
        const storage::Table< double> &LABEL
      ) = 0;

      //! @brief set the watermark
      //! @param WATERMARK the watermark
      virtual void SetWaterMark( const std::string &WATERMARK = GetDefaultWatermark()) = 0;

      //! @brief create ffmpeg command line
      //! @return commandline to be called with ffmpeg to generate the movie from the generated pngs
      virtual std::string FFMpegCommandLine() const = 0;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write the movie script to a file
      //! @param OSTREAM
      virtual std::ostream &WriteScript
      (
        std::ostream &OSTREAM
      ) const = 0;

    }; // class MoviePrinterInterface

  } // namespace mc
} // namespace bcl

#endif // BCL_MC_MOVIE_PRINTER_INTERFACE_H_
