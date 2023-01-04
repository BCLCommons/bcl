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

#ifndef BCL_MC_METROPOLIS_H_
#define BCL_MC_METROPOLIS_H_

// include the namespace header
#include "bcl_mc.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_mc_temperature_default.h"
#include "bcl_mc_temperature_interface.h"
#include "io/bcl_io_serialization.h"
#include "opti/bcl_opti_step_status.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_triplet.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_serializable_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically
#include <math.h>

namespace bcl
{
  namespace mc
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Metropolis
    //! @brief determines whether a specified amount of change in energy should be accepted.
    //! @details Metropolis class utilizes Metropolis criterium, used commonly in Monte Carlo based minimizations, to decide
    //! whether a mutate that leads to a specific change in the energy should be accepted or not. The change is accepted
    //! if the energy has been improved. In other cases, it can still get accepted but with a probability that is
    //! inversely proportional to the amount of increase in the energy.
    //!
    //! @see @link example_mc_metropolis.cpp @endlink
    //! @author karakam
    //! @date Feb 10, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ResultType>
    class Metropolis :
      public util::SerializableInterface
    {

    private:

    //////////
    // data //
    //////////

      //! Temperature adjustment class
      util::Implementation< TemperatureInterface> m_Temperature;

      //! boolean to whether keep track of all temperature history
      bool m_KeepHistory;

      //! Minimal change in the energy that has to happen for improved
      double m_MinimalChange;

      //! history of all temperatures
      mutable storage::List< storage::Pair< opti::StepStatusEnum, double> > m_TemperatureHistory;

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
      Metropolis() :
        m_Temperature(),
        m_KeepHistory(),
        m_MinimalChange(),
        m_TemperatureHistory()
      {
      }

      //! @brief constructor from a TemperatureInterface
      //! @param SP_TEMPERATURE_INTERFACE ShPtr to TemperatureInterface to adjust temperature
      //! @param KEEP_HISTORY boolean to determine whether to keep track of history or not
      //! @param MINIMAL_CHANGE Minimal change in the energy that has to happen for improved
      Metropolis
      (
        const util::ShPtr< TemperatureInterface> &SP_TEMPERATURE_INTERFACE,
        const bool KEEP_HISTORY = false,
        const double MINIMAL_CHANGE = double( 0.0)
      ) :
        m_Temperature( *SP_TEMPERATURE_INTERFACE),
        m_KeepHistory( KEEP_HISTORY),
        m_MinimalChange( MINIMAL_CHANGE),
        m_TemperatureHistory()
      {
      }

      //! @brief Clone function
      //! @return pointer to new Metropolis< t_ArgumentType, t_ResultType>
      Metropolis< t_ResultType> *Clone() const
      {
        return new Metropolis< t_ResultType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get temperature history
      //! @return temperature history
      const storage::List< storage::Pair< opti::StepStatusEnum, double> > &GetTemperatureHistory() const
      {
        return m_TemperatureHistory;
      }

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const
      {
        static const std::string s_name( "Metropolis");
        return s_name;
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer serializer;
        serializer.SetClassDescription
        (
          "Metropolis criterion that decides whether to accept or reject a change based on a simulated"
          "temperature and a score difference."
        );
        serializer.AddInitializer
        (
          "keep history",
          "track the complete temperature history",
          io::Serialization::GetAgent( &m_KeepHistory)
        );
        serializer.AddInitializer
        (
          "minimum change",
          "minimum score change considered as improvement",
          io::Serialization::GetAgent( &m_MinimalChange)
        );
        serializer.AddInitializer
        (
          "temperature control",
          "temperature control to adjust the temperature during the simulation",
          io::Serialization::GetAgent( &m_Temperature)
        );

        return serializer;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief reset this metropolis
      void Reset()
      {
        // reset members
        m_Temperature->Reset();
        m_TemperatureHistory.Reset();
      }

    ///////////////
    // operators //
    ///////////////

      // TODO: add detailed description of Metropolis stuff
      //! @brief function to evaluate the delta and decide on whether it should be accepted or rejected
      //! @param PREVIOUS_RESULT
      //! @param CURRENT_RESULT
      //! @param TRACKER current tracker base
      //! @return StepStatus; e_Accepted or e_Rejected
      opti::StepStatus Evaluate
      (
        const t_ResultType &PREVIOUS_RESULT,
        const t_ResultType &CURRENT_RESULT,
        const opti::TrackerBase &TRACKER
      )
      {
        // calculate the delta between two arguments
        double delta( CURRENT_RESULT - PREVIOUS_RESULT);

        // consider if opti is toward maximum
        if( TRACKER.GetImprovementType() == opti::e_LargerIsBetter)
        {
          delta = -delta;
        } // without this, assumes opti::e_SmallerIsBetter

        // Rejected by default
        opti::StepStatus step_status( opti::e_Rejected);

        // Accept the change if difference is larger than the minimal change( set to 0 by default)
        if( delta < -m_MinimalChange)
        {
          // improved
          step_status = opti::e_Improved;
        }
        else
        {
          // if the min_change was not 0 and this value was between 0 and min change
          if( m_MinimalChange != double( 0.0) && delta < double( 0.0))
          {
            // set delta to 0
            delta = 0.0;
          }

          // track the delta
          m_Temperature->TrackDelta( delta);

          // otherwise accept with a certain probability between 0 and 1 that is dependent on the amount of
          // change vs current temperature
          const double temperature
          (
            std::max( m_Temperature->GetTemperature( TRACKER), std::numeric_limits< double>::epsilon())
          );
          if( random::GetGlobalRandom().Double() < exp( -delta / temperature))
          {
            // improved
            step_status = opti::e_Accepted;
          }
        }

        // if history if being kept
        if( m_KeepHistory)
        {
          // insert it into the history
          m_TemperatureHistory.PushBack
          (
            storage::Pair< opti::StepStatusEnum, double>
            (
              step_status, m_Temperature->GetLastCalculatedTemperature()
            )
          );
        }

        // else reject it
        return step_status;
      }

    //////////////////////
    // input and output //
    //////////////////////

    }; // template class Metropolis

    // instantiate s_Instance
    template< typename t_ResultType>
    const util::SiPtr< const util::ObjectInterface> Metropolis< t_ResultType>::s_Instance
    (
      GetObjectInstances().AddInstance( new Metropolis< t_ResultType>())
    );

  } // namespace mc
} // namespace bcl

#endif // BCL_MC_METROPOLIS_H_
