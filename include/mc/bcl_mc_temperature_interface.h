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

#ifndef BCL_MC_TEMPERATURE_INTERFACE_H_
#define BCL_MC_TEMPERATURE_INTERFACE_H_

// include the namespace header
#include "bcl_mc.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "opti/bcl_opti_tracker_base.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class TemperatureInterface
    //! @brief Interface class for Temperature classes
    //! @details This class provides the interface for various Monte Carlo related temperature adjustment classes to be derived
    //! from. The only function that has to be overwritten is GetTemperature(), it is the responsibility of that class
    //! to make sure it has access to Tracker information if required.
    //!
    //! @remarks example unnecessary
    //! @author karakam, woetzen
    //! @date 08.14.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API TemperatureInterface :
      public util::SerializableInterface
    {

    public:

    ////////////////////////////////////
    //  construction and destruction  //
    ////////////////////////////////////

      //! @brief virtual copy constructor
      virtual TemperatureInterface *Clone() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief returns last calculated temperature without updating
      //! @return last calculated temperature without updating
      virtual double GetLastCalculatedTemperature() const = 0;

      //! @brief returns temperature
      //! @param TRACKER the current tracker
      //! @return temperature
      virtual double GetTemperature( const opti::TrackerBase &TRACKER) const = 0;

      //! @brief Provide the temperature calculator with the last delta value, if needed
      virtual void TrackDelta( const double &DELTA)
      {
      }

      //! @brief reset this temperature
      virtual void Reset() = 0;

    }; // class TemperatureInterface

  } // namespace mc
} // namespace bcl

#endif //BCL_MC_TEMPERATURE_INTERFACE_H_
