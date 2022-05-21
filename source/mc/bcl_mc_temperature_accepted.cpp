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
#include "mc/bcl_mc_temperature_accepted.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> TemperatureAccepted::s_Instance
    (
      util::Enumerated< TemperatureInterface>::AddInstance( new TemperatureAccepted())
    );

    //! @brief return command line flag for setting temperature adjustment based on accepted steps ratio
    //! @return command line flag for setting temperature adjustment based on accepted steps ratio
    util::ShPtr< command::FlagInterface> &TemperatureAccepted::GetFlagTemperature()
    {
      // initialize static instance of the flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "mc_temperature_fraction",
          "\tmodify the start and end fractions of accepted steps for temperature adjustments for monte carlo minimization"
        )
      );

      // if the flag is initialized for the first time
      if( s_flag->GetParameterList().IsEmpty())
      {
        util::ShPtr< command::FlagStatic> flag( s_flag);
        // insert parameters into flag
        flag->PushBack( GetParameterStartFraction());
        flag->PushBack( GetParameterEndFraction());
        flag->PushBack( GetParameterStartTemperature());
        flag->PushBack( GetParameterUpdateInterval());
      }

      // end
      return s_flag;
    }

    //! @brief return command line parameter for setting the start fraction
    //! @return command line parameter for setting the start fraction
    util::ShPtr< command::ParameterInterface> &TemperatureAccepted::GetParameterStartFraction()
    {
      // initialize static instance of this parameter
      static util::ShPtr< command::ParameterInterface> s_parameter
      (
        new command::Parameter
        (
          "start_ratio",
          "\tthis a fraction of accepted steps at the start of the algorithm",
          command::ParameterCheckRanged< double>( 0.0, 1.0),
          "0.5"
        )
      );

      // end
      return s_parameter;
    }

    //! @brief return command line parameter for setting the end fraction
    //! @return command line parameter for setting the end fraction
    util::ShPtr< command::ParameterInterface> &TemperatureAccepted::GetParameterEndFraction()
    {
      // initialize static instance of this parameter
      static util::ShPtr< command::ParameterInterface> s_parameter
      (
        new command::Parameter
        (
          "end_ratio",
          "\tthis a fraction of accepted steps at the end of the algorithm",
          command::ParameterCheckRanged< double>( 0.0, 1.0),
          "0.2"
        )
      );

      // end
      return s_parameter;
    }

    //! @brief return command line parameter for setting the end fraction
    //! @return command line parameter for setting the end fraction
    util::ShPtr< command::ParameterInterface> &TemperatureAccepted::GetParameterStartTemperature()
    {
      // initialize static instance of this parameter
      static util::ShPtr< command::ParameterInterface> s_parameter
      (
        new command::Parameter
        (
          "start_temperature", "\tthis the starting temperature ", "1"
        )
      );

      // end
      return s_parameter;
    }

    //! @brief return command line parameter for setting the end fraction
    //! @return command line parameter for setting the end fraction
    util::ShPtr< command::ParameterInterface> &TemperatureAccepted::GetParameterUpdateInterval()
    {
      // initialize static instance of this parameter
      static util::ShPtr< command::ParameterInterface> s_parameter
      (
        new command::Parameter
        (
          "update_interval", "\tnumber of steps between each update", "10"
        )
      );

      // end
      return s_parameter;
    }

  ////////////////////////////////////
  //  construction and destruction  //
  ////////////////////////////////////

    //! @brief default constructor
    TemperatureAccepted::TemperatureAccepted() :
      m_Temperature(),
      m_StartFraction(),
      m_EndFraction(),
      m_StartTemperature(),
      m_Delta(),
      m_UpdateInterval()
    {
    }

    //! @brief constructor from number iterations, the other paramaters are taken from the commandline arguments
    //! @param NUMBER_OF_ITERATIONS number of iterations
    TemperatureAccepted::TemperatureAccepted
    (
      const size_t NUMBER_OF_ITERATIONS
    ) :
      m_Temperature( 0, GetParameterStartTemperature()->GetNumericalValue< double>()),
      m_StartFraction( GetParameterStartFraction()->GetNumericalValue< double>()),
      m_EndFraction( GetParameterEndFraction()->GetNumericalValue< double>()),
      m_StartTemperature( GetParameterStartTemperature()->GetNumericalValue< double>()),
      m_Delta( 0.0),
      m_UpdateInterval( GetParameterUpdateInterval()->GetNumericalValue< size_t>())
    {
    }

    //! @brief constructor from starting and ending temperature
    //! @param START_FRACTION starting fraction
    //! @param END_FRACTION ending fraction
    //! @param NUMBER_OF_ITERATIONS number of iterations
    //! @param START_TEMPERATURE starting temperature
    //! @param UPDATE_INTERVAL interval length between each update
    TemperatureAccepted::TemperatureAccepted
    (
      const double START_FRACTION,
      const double END_FRACTION,
      const size_t NUMBER_OF_ITERATIONS,
      const double START_TEMPERATURE,
      const size_t UPDATE_INTERVAL
    ) :
      m_Temperature( 0, START_TEMPERATURE),
      m_StartFraction( START_FRACTION),
      m_EndFraction( END_FRACTION),
      m_StartTemperature( START_TEMPERATURE),
      m_Delta( 0.0),
      m_UpdateInterval( UPDATE_INTERVAL)
    {
      //check that start fraction is in interval (0, 1], end fraction in interval [0, 1) and that start > end
      if( m_StartFraction > 0.0 && m_StartFraction <= 1.0 && m_EndFraction > 0.0 && m_EndFraction <= 1.0)
      {
        // calculate the delta as in the linear case
        // the /2 is required since not the correct fraction but the overall fraction of accepted steps is monitored
        m_Delta = ( m_EndFraction - m_StartFraction) / NUMBER_OF_ITERATIONS;
      }
    }

    //! @brief virtual copy constructor
    TemperatureAccepted *TemperatureAccepted::Clone() const
    {
      return new TemperatureAccepted( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &TemperatureAccepted::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &TemperatureAccepted::GetAlias() const
    {
      static const std::string s_name( "AcceptedTemperatureControl");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer TemperatureAccepted::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        "Temperature control that adjusts the temperature to achieve a given fraction of accepted versus"
        "rejected mutations."
      );
      serializer.AddInitializer
      (
        "start fraction",
        "start fraction of accepted versus rejected mutations",
        io::Serialization::GetAgent( &m_StartFraction),
        "0.5"
      );
      serializer.AddInitializer
      (
        "end fraction",
        "end fraction of accepted versus rejected mutations",
        io::Serialization::GetAgent( &m_EndFraction),
        "0.2"
      );
      serializer.AddInitializer
      (
        "start temperature",
        "temperature at the beginning of the simulation",
        io::Serialization::GetAgent( &m_StartTemperature),
        "1"
      );
      serializer.AddInitializer
      (
        "delta",
        "magnitude of temperature change",
        io::Serialization::GetAgent( &m_Delta),
        "0.00015"
      );
      serializer.AddInitializer
      (
        "update interval",
        "length of the interval between updates",
        io::Serialization::GetAgent( &m_UpdateInterval),
        "10"
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief reset temperature
    void TemperatureAccepted::Reset()
    {
      // reset members
      m_Temperature = storage::Pair< size_t, double>( 0, m_StartTemperature);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief return current temperature
    //! @param TRACKER the current tracker
    //! @return current temperature
    double TemperatureAccepted::GetTemperature( const opti::TrackerBase &TRACKER) const
    {
      // if the iteration number does not match the counter
      if( TRACKER.GetIteration() - m_Temperature.First() >= m_UpdateInterval)
      {
        // update temperature
        UpdateTemperature( TRACKER);
      }

      // return temperature
      return m_Temperature.Second();
    }

    //! @brief calculate new temperature
    //! @param TRACKER the base tracker
    void TemperatureAccepted::UpdateTemperature( const opti::TrackerBase &TRACKER) const
    {
      // update only every m_ConstantIntervalLength' step
      if( TRACKER.GetIteration() - m_Temperature.First() < m_UpdateInterval)
      {
        return;
      }

      // update the iteration number
      m_Temperature.First() = TRACKER.GetIteration();

      // calculate target fraction
      const double target_fraction( m_StartFraction + m_Temperature.First() * m_Delta);

      const size_t points_to_consider
      (
        std::min( size_t( std::max( 1.0 / std::max( target_fraction, 1.0e-6), 10.0)), m_UpdateInterval)
      );
      const size_t nth_point_desired
      (
        std::max( size_t( target_fraction * points_to_consider), size_t( 1)) - 1
      );
      // copy the vector of deltas under consideration
      storage::Vector< double> deltas_copy( m_PreviousDeltas.End() - points_to_consider, m_PreviousDeltas.End());
      std::nth_element( deltas_copy.Begin(), deltas_copy.Begin() + nth_point_desired, deltas_copy.End());
      static const double ln2( std::log( 2.0));
      m_Temperature.Second() = 0.7 * m_Temperature.Second() + 0.3 * deltas_copy( nth_point_desired) / ln2;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &TemperatureAccepted::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Temperature, ISTREAM);
      io::Serialize::Read( m_StartFraction, ISTREAM);
      io::Serialize::Read( m_EndFraction, ISTREAM);
      io::Serialize::Read( m_StartTemperature, ISTREAM);
      io::Serialize::Read( m_Delta, ISTREAM);
      io::Serialize::Read( m_UpdateInterval, ISTREAM);

      // end
      return ISTREAM;
    };

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &TemperatureAccepted::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Temperature, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_StartFraction, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_EndFraction, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_StartTemperature, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Delta, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_UpdateInterval, OSTREAM, INDENT);

      // end
      return OSTREAM;
    };

  } // namespace mc
} // namespace bcl
