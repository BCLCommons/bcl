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

#ifndef BCL_MC_MOVIE_PRINTERS_H_
#define BCL_MC_MOVIE_PRINTERS_H_

// include the namespace header
#include "bcl_mc.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_mc_movie_printer_interface.h"
#include "util/bcl_util_enumerate.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoviePrinters
    //! @brief manages the mc printers
    //!
    //! @see @link example_mc_movie_printers.cpp @endlink
    //! @author woetzen
    //! @date November 21, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoviePrinters :
      public util::Enumerate< util::ShPtr< MoviePrinterInterface>, MoviePrinters>
    {
      friend class util::Enumerate< util::ShPtr< MoviePrinterInterface>, MoviePrinters>;

    public:

    //////////
    // data //
    //////////

      MoviePrinter e_Chimera; //!< chimera MoviePrinter
      MoviePrinter e_Pymol;   //!< pymol MoviePrinter

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Default constructor
      MoviePrinters();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    }; // class MoviePrinters

    //! @brief get enumerated list of MoviePrinters
    BCL_API const MoviePrinters &GetMoviePrinters();

  } // namespace mc

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< mc::MoviePrinterInterface>, mc::MoviePrinters>;

  } // namespace util
} // namespace bcl

#endif // BCL_MC_MOVIE_PRINTERS_H_
