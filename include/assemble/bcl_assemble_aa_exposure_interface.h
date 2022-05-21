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

#ifndef BCL_ASSEMBLE_AA_EXPOSURE_INTERFACE_H_
#define BCL_ASSEMBLE_AA_EXPOSURE_INTERFACE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AAExposureInterface
    //! @brief calculate the exposure for an amino acids surrounded by its neighbors
    //! @details for an AANeighbors list, the surface exposure is calculated
    //!
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date Apr 17, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AAExposureInterface :
      public math::FunctionInterfaceSerializable< AANeighborList, double>
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new AAExposureInterface
      virtual AAExposureInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief gives the flag allowing user to set the sequence exclusion over the command line
      //! @return the flag allowing user to set the sequence exclusion over the command line
      static const util::ShPtr< command::FlagInterface> &GetFlagMinimalSequenceSeparation();

      //! @brief access to the distance cutoff
      //! @return distance cutoff above which the neighbor does not have influence on the exposure anymore
      double GetDistanceCutoff() const;

      //! @brief min and max exposure measure
      //! @return range in which exposure can be
      virtual const math::Range< double> &GetRange() const = 0;

      //! @brief get threshold range
      //! @details threshold are used, to have a continuous function for the exposure measure, instead of a stepwise
      //! @return the default thresholds used
      virtual const math::Range< double> &GetThresholdRange() const = 0;

      //! @brief set the threshold range
      //! @param RANGE the range in which other neighbors that are considered for the exposure do not count full
      virtual void SetThresholdRange( const math::Range< double> &RANGE) = 0;

      //! @brief access to the minimal sequence separation
      //! @return minimal sequence separation in sequence distance
      virtual size_t GetMinimalSequenceSeparation() const = 0;

      //! @brief access to the minimal sequence separation
      //! @param MINIMAL_SEQUENCE_SEPARATION in sequence distance
      virtual void SetMinimalSequenceSeparation( const size_t MINIMAL_SEQUENCE_SEPARATION) = 0;

      //! @brief return the histogram filename where statistics are stored
      //! @return the filename containing the thresholds, sequence separation and histograms for each environment and
      //!         aa type from a databank of proteins
      virtual const std::string &GetHistogramFileName() const = 0;

      //! @brief is direct measure - exposed surface correlates with the measure or actually measures buried surface
      //! @return true, if exposure measure correlates with exposed surface, false if it correlates to buried surface
      virtual bool IsDirect() const = 0;

    }; // class AAExposureInterface

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_AA_EXPOSURE_INTERFACE_H_ 
