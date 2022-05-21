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
#include "mc/bcl_mc_movie_printers.h"

// includes from bcl - sorted alphabetically
#include "mc/bcl_mc_movie_printer_chimera.h"
#include "mc/bcl_mc_movie_printer_pymol.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Default constructor
    MoviePrinters::MoviePrinters() :
      e_Chimera( AddEnum( "Chimera", util::ShPtr< MoviePrinterInterface>( new MoviePrinterChimera()))),
      e_Pymol( AddEnum( "Pymol", util::ShPtr< MoviePrinterInterface>( new MoviePrinterPymol())))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MoviePrinters::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get enumerated list of MoviePrinters
    const MoviePrinters &GetMoviePrinters()
    {
      return MoviePrinters::GetEnums();
    }

  } // namespace mc

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< mc::MoviePrinterInterface>, mc::MoviePrinters>;

  } // namespace util
} // namespace bcl
