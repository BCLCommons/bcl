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

#ifndef BCL_QUALITY_MEASURES_H
#define BCL_QUALITY_MEASURES_H

// include the namespace header
#include "bcl_quality.h"

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_quality_measure_interface.h"
#include "util/bcl_util_enumerate.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace quality
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Measures
    //! @brief Measures class provides enumeration of different quality measures to compare two given coordinate vectors
    //! @details Measures enumerates Interface derived class that compare two given coordinate vectors and returns a quality
    //! measure
    //!
    //! @see @link example_quality_measures.cpp @endlink
    //! @author alexanns, karakam
    //! @date Apr 22, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Measures :
      public util::Enumerate< util::ShPtr< MeasureInterface>, Measures>
    {
      friend class util::Enumerate< util::ShPtr< MeasureInterface>, Measures>;

    public:

    //////////
    // data //
    //////////

      const Measure e_RMSD;
      const Measure e_RMSD_NoSuperimposition;
      const Measure e_RMSD_XYSuperimposition;
      const Measure e_DME;
      const Measure e_DMF_HA;
      const Measure e_DMF_TS;
      const Measure e_LCS;
      const Measure e_GDT_HA;
      const Measure e_GDT_TS;
      const Measure e_GDT_1A;
      const Measure e_GDT_2A;
      const Measure e_GDT_4A;
      const Measure e_GDT_8A;
      const Measure e_MaxSub;
      const Measure e_Zero;

      //! @brief returns the default distance cutoff set for HA (high accuracy)
      //! @return the default distance cutoff set HA
      static const storage::Set< double> &GetDistanceCutoffsHA();

      //! @brief returns the default distance cutoff set for TS (total score)
      //! @return the default distance cutoff set TS
      static const storage::Set< double> &GetDistanceCutoffsTS();

      //! @brief return command line flag for defining the quality measures to be calculated
      //! @return command line flag for defining the quality measures to be calculated
      static util::ShPtr< command::FlagInterface> &GetFlagQualityMeasures();

      //! @brief function to return the list of quality measures in a set defined by the command line flag
      //! @return the list of quality measures in a set defined by the command line flag
      static storage::Set< Measure> GetCommandLineQualityMeasures();

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct all Measures
      Measures();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    }; // class Measures

    //! @brief function that returns the static instance of the Measures class
    //! @return the static instance of the Measures class
    BCL_API
    Measures &GetMeasures();

  } // namespace quality

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< quality::MeasureInterface>, quality::Measures>;

  } // namespace util
} // namespace bcl

#endif //BCL_QUALITY_MEASURES_H
