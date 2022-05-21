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

#ifndef BCL_QUALITY_SUPERIMPOSE_MEASURES_H
#define BCL_QUALITY_SUPERIMPOSE_MEASURES_H

// include the namespace header
#include "bcl_quality.h"

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_quality_superimpose_interface.h"
#include "util/bcl_util_enumerate.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace quality
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SuperimposeMeasures
    //! @brief SuperimposeMeasures class provides enumeration of different measures for superimposing two given coordinate vectors
    //! @details SuperimposeMeasures enumerates SuperimposeInterface derived class that calculates superimposition of two
    //! given coordinate vectors according to a quality measure and returns the transformation matrix.
    //!
    //! @see @link example_quality_superimpose_measures.cpp @endlink
    //! @author karakam
    //! @date Oct 31, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SuperimposeMeasures :
      public util::Enumerate< util::ShPtr< SuperimposeInterface>, SuperimposeMeasures>
    {
      friend class util::Enumerate< util::ShPtr< SuperimposeInterface>, SuperimposeMeasures>;

    public:

    //////////
    // data //
    //////////

      const SuperimposeMeasure e_RMSD;
      const SuperimposeMeasure e_RMSD_XYSuperimposition;
      const SuperimposeMeasure e_GDT_1A;
      const SuperimposeMeasure e_GDT_2A;
      const SuperimposeMeasure e_GDT_4A;
      const SuperimposeMeasure e_GDT_8A;
      const SuperimposeMeasure e_LCS;
      const SuperimposeMeasure e_MaxSub;
      const SuperimposeMeasure e_NoSuperimpose;

      //! @brief return command line flag for defining the measures to calculate superimposition
      //! @return command line flag for  defining the measure to calculate superimposition
      static util::ShPtr< command::FlagInterface> &GetFlagSuperimposeMeasure();

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct all SuperimposeMeasures
      SuperimposeMeasures();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    }; // class SuperimposeMeasures

    //! @brief function that returns static instance of SuperimposeMeasures
    //! @return static instance of SuperimposeMeasures
    BCL_API SuperimposeMeasures &GetSuperimposeMeasures();

  } // namespace quality

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< quality::SuperimposeInterface>, quality::SuperimposeMeasures>;

  } // namespace util
} // namespace bcl

#endif //BCL_QUALITY_SUPERIMPOSE_MEASURES_H
