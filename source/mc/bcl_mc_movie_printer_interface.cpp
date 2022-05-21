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
#include "command/bcl_command_flag_interface.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "mc/bcl_mc_movie_printers.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

  //////////
  // data //
  //////////

    //! @brief get the default bcl watermark
    const std::string &MoviePrinterInterface::GetDefaultWatermark()
    {
      return GetCopyright();
    }

    //! @brief return command line flag for printing the script necessary to create a movie
    //! @return command line flag for printing the script necessary to create a movie
    util::ShPtr< command::FlagInterface> &MoviePrinterInterface::GetFlagMoviePrinter()
    {
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "movie",
          "\twrite a movie script and change preferences for movie generation"
        )
      );

      // if the flag is initialized for the first time
      if( s_flag->GetParameterList().IsEmpty())
      {
        util::ShPtr< command::FlagStatic> flag( s_flag);
        // insert parameters into flag
        flag->PushBack( GetParameterMoviePrefix());
        flag->PushBack( GetParameterMoviePrinterType());
        flag->PushBack( GetParameterMovieWidth());
        flag->PushBack( GetParameterMovieHeight());
        flag->PushBack( GetParameterRayTrace());
      }

      return s_flag;
    }

    //! @brief return command line parameter for specifying the prefix for where movie files are generated
    //! @return command line parameter for specifying the prefix for where movie files are generated
    util::ShPtr< command::ParameterInterface> &MoviePrinterInterface::GetParameterMoviePrefix()
    {
      static util::ShPtr< command::ParameterInterface> s_parameter
      (
        new command::Parameter
        (
          "movie_prefix",
          "prefix for where movie files are generated",
          ""
        )
      );

      return s_parameter;
    }

    //! @brief return command line parameter for specifying the type of movie printer that should be used
    //! @return command line parameter for specifying the type of movie printer that should be used
    util::ShPtr< command::ParameterInterface> &MoviePrinterInterface::GetParameterMoviePrinterType()
    {
      static util::ShPtr< command::ParameterInterface> s_parameter
      (
        new command::Parameter
        (
          "movieprinter",
          "choice of movie printer for different programs",
          command::ParameterCheckEnumerate< MoviePrinters>(),
          GetMoviePrinters().e_Pymol
        )
      );

      return s_parameter;
    }

    //! @brief return command line parameter for specifying the width resolution of the movie
    //! @return command line parameter for specifying the width resolution of the movie
    util::ShPtr< command::ParameterInterface> &MoviePrinterInterface::GetParameterMovieWidth()
    {
      static util::ShPtr< command::ParameterInterface> s_parameter
      (
        new command::Parameter
        (
          "movie_width",
          "width of movie",
          command::ParameterCheckRanged< size_t>( 480, 1920),
          "720"
        )
      );

      return s_parameter;
    }

    //! @brief return command line parameter for specifying the height resolution of the movie
    //! @return command line parameter for specifying the height resolution of the movie
    util::ShPtr< command::ParameterInterface> &MoviePrinterInterface::GetParameterMovieHeight()
    {
      static util::ShPtr< command::ParameterInterface> s_parameter
      (
        new command::Parameter
        (
          "movie_height",
          "height of movie",
          command::ParameterCheckRanged< size_t>( 320, 1200),
          "480"
        )
      );

      return s_parameter;
    }

    //! @brief return command line parameter for specifying whether or not all movie frames should be ray traced
    //! @return command line parameter for specifying whether or not all movie frames should be ray traced
    util::ShPtr< command::ParameterInterface> &MoviePrinterInterface::GetParameterRayTrace()
    {
      static util::ShPtr< command::ParameterInterface> s_parameter
      (
        new command::Parameter
        (
          "ray",
          "raytrace all frames - 0=no, 1=yes",
          command::ParameterCheckRanged< size_t>( 0, 1),
          "0"
        )
      );

      return s_parameter;
    }

  } // namespace mc
} // namespace bcl
