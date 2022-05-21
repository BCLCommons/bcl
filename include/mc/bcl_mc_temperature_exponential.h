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

#ifndef BCL_MC_TEMPERATURE_EXPONENTIAL_H_
#define BCL_MC_TEMPERATURE_EXPONENTIAL_H_

// include the namespace header
#include "bcl_mc.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_mc_temperature_interface.h"
#include "util/bcl_util_si_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class TemperatureExponential
    //! @brief Allows exponential adjustment of the temperature
    //! @details TemperatureInterface derived class that allows changing the temperature exponentially, from a given
    //! start temperature to a given end temperature and number of total expected number of iterations
    //!
    //! @see @link example_mc_temperature_exponential.cpp @endlink
    //! @author karakam, woetzen, fischea
    //! @date 08.14.2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API TemperatureExponential :
      public TemperatureInterface
    {

    private:

    //////////
    // data //
    //////////

      //! pair of last iteration number and corresponding temperature
      mutable size_t m_LastIterationNumber;

      //! magnitude of temperature adjustment as scaling
      double m_Scale;

      //! starting temperature
      double m_StartTemperature;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    ////////////////////////////////////
    //  construction and destruction  //
    ////////////////////////////////////

      //! @brief default constructor
      TemperatureExponential();

      //! @brief constructor from a starting and ending temperature and total number of iterations
      //! @param START_TEMPERATURE starting temperature
      //! @param END_TEMPERATURE ending temperature
      //! @param NUMBER_OF_ITERATIONS number of iterations
      TemperatureExponential
      (
        const double START_TEMPERATURE,
        const double END_TEMPERATURE,
        const size_t NUMBER_OF_ITERATIONS
      );

      //! @brief virtual copy constructor
      TemperatureExponential *Clone() const;

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
      double GetLastCalculatedTemperature() const;

      //! @brief return current temperature
      //! @param TRACKER the current tracker
      //! @return current temperature
      double GetTemperature( const opti::TrackerBase &TRACKER) const
      {
        // update the iteration number and the corresponding temperature
        m_LastIterationNumber = TRACKER.GetIteration();

        // return temperature
        return GetLastCalculatedTemperature();
      }

      //! @brief reset this temperature
      void Reset();

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

    }; // class TemperatureExponential

  } // namespace mc
} // namespace bcl

#endif //BCL_MC_TEMPERATURE_EXPONENTIAL_H_
