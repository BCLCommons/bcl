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

#ifndef BCL_MC_TEMPERATURE_ACCEPTED_H_
#define BCL_MC_TEMPERATURE_ACCEPTED_H_

// include the namespace header
#include "bcl_mc.h"

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_mc_temperature_interface.h"
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_si_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class TemperatureAccepted
    //! @brief temperature adjustment according to the ratio of accepted steps
    //! @details This class adjusts the temperature according to the ratio of accepted steps. It tries the match the
    //! current ratio of accepted steps to the one calculated between given start ratio and end ratio and using the
    //! number of steps as a linear predictor. It only tries to adjust the temperature every Nth step where N is
    //! specified by the user. If the ratio is lower than the expected one, it increases the temperature, otherwise it
    //! decreases the temperature. For cases where the actual ratio is close to expected ratio ( within 0.1) it
    //! adjusts the temperature by multiplying/diving by a small coefficient, while for cases where the difference
    //! is larger, the temperature is adjusted by multiplying/diving by a larger coefficient.
    //!
    //! @see @link example_mc_temperature_accepted.cpp @endlink
    //! @author karakam, woetzen, fischea
    //! @date 08.14.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API TemperatureAccepted :
      public TemperatureInterface
    {

    public:

    //////////
    // data //
    //////////

      //! @brief return command line flag for setting temperature adjustment based on accepted steps ratio
      //! @return command line flag for setting temperature adjustment based on accepted steps ratio
      static util::ShPtr< command::FlagInterface> &GetFlagTemperature();

      //! @brief return command line parameter for setting the start fraction
      //! @return command line parameter for setting the start fraction
      static util::ShPtr< command::ParameterInterface> &GetParameterStartFraction();

      //! @brief return command line parameter for setting the end fraction
      //! @return command line parameter for setting the end fraction
      static util::ShPtr< command::ParameterInterface> &GetParameterEndFraction();

      //! @brief return command line parameter for setting the end fraction
      //! @return command line parameter for setting the end fraction
      static util::ShPtr< command::ParameterInterface> &GetParameterStartTemperature();

      //! @brief return command line parameter for setting the end fraction
      //! @return command line parameter for setting the end fraction
      static util::ShPtr< command::ParameterInterface> &GetParameterUpdateInterval();

    private:

    //////////
    // data //
    //////////

      //! pair of last iteration number and corresponding temperature
      mutable storage::Pair< size_t, double> m_Temperature;

      //! starting fraction
      double m_StartFraction;

      //! starting fraction
      double m_EndFraction;

      //! start temperature
      double m_StartTemperature;

      //! magnitude of temperature adjustment as scaling
      double m_Delta;

      //! length of interval between adjustments
      size_t m_UpdateInterval;

      //! list of previous delta values
      storage::Vector< double> m_PreviousDeltas;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    ////////////////////////////////////
    //  construction and destruction  //
    ////////////////////////////////////

      //! @brief default constructor
      TemperatureAccepted();

      //! @brief constructor from number iterations, the other paramaters are taken from the commandline arguments
      //! @param NUMBER_OF_ITERATIONS number of iterations
      TemperatureAccepted( const size_t NUMBER_OF_ITERATIONS);

      //! @brief constructor from starting and ending temperature
      //! @param START_FRACTION starting fraction
      //! @param END_FRACTION ending fraction
      //! @param NUMBER_OF_ITERATIONS number of iterations
      //! @param START_TEMPERATURE starting temperature
      //! @param UPDATE_INTERVAL interval length between each update
      TemperatureAccepted
      (
        const double START_FRACTION,
        const double END_FRACTION,
        const size_t NUMBER_OF_ITERATIONS,
        const double START_TEMPERATURE = double( 1.0),
        const size_t UPDATE_INTERVAL = 10
      );

      //! @brief virtual copy constructor
      TemperatureAccepted *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief returns last calculated temperature without updating
      //! @return last calculated temperature without updating
      double GetLastCalculatedTemperature() const
      {
        return m_Temperature.Second();
      }

      //! @brief return current temperature
      //! @param TRACKER the current tracker
      //! @return current temperature
      double GetTemperature( const opti::TrackerBase &TRACKER) const;

      //! @brief Provide the temperature calculator with the last delta value, if needed
      void TrackDelta( const double &DELTA)
      {
        m_PreviousDeltas.PushBack( DELTA);
      }

      //! @brief reset temperature
      void Reset();

    private:

    ////////////////
    // operations //
    ////////////////

      //! @brief calculate new temperature
      //! @param TRACKER the current tracker
      void UpdateTemperature( const opti::TrackerBase &TRACKER) const;

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream
      //! @param INDENT indentation
      //! @return ostream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class TemperatureAccepted

  } // namespace mc
} // namespace bcl

#endif //BCL_MC_TEMPERATURE_ACCEPTED_H_
