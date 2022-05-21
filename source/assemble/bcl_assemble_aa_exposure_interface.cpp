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
#include "assemble/bcl_assemble_aa_exposure_interface.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_count.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  /////////////////
  // data access //
  /////////////////

    //! @brief gives the flag allowing user to set the sequence exclusion over the command line
    //! @return the flag allowing user to set the sequence exclusion over the command line
    const util::ShPtr< command::FlagInterface> &AAExposureInterface::GetFlagMinimalSequenceSeparation()
    {
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "minimal_sequence_exclusion",
          "The sequence exclusion to use when calculating exposure from a structure.",
          command::Parameter
          (
            "sequence_exclusion",
            "size_t which is the number of residues to be excluded in sequence when calculating exposures",
            util::Format()( AANeighborCount:: GetDefaultMinimalSequenceSeparation())
          )
        )
      );

      return s_flag;
    }

    //! @brief access to the distance cutoff
    //! @return distance cutoff above which the neighbor does not have influence on the score anymore
    double AAExposureInterface::GetDistanceCutoff() const
    {
      return GetThresholdRange().GetMax();
    }

  } // namespace assemble
} // namespace bcl
