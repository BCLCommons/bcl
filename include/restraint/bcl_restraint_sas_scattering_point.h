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

#ifndef BCL_RESTRAINT_SAS_SCATTERING_POINT_H_
#define BCL_RESTRAINT_SAS_SCATTERING_POINT_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SasScatteringPoint
    //! @brief Wrapper Class of SaxsScatteringData and SansScatteringData
    //! @details Class to explicitly list intensity scattering angle and error
    //!
    //! @see @link example_restraint_saxs_scattering_point.cpp @endlink
    //! @author putnamdk
    //! @date Jul 28, 2011
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SasScatteringPoint :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! @brief Qvalue represents the scattering angle from Saxs Data
      double m_Qvalue;

      //! @brief Intensity represents the intensity for a given scattering angle for Saxs Data
      double m_Intensity;

      //! @brief Error is the measurement error
      double m_Error;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SasScatteringPoint();

      //! @brief constructor from given input data
      //! @param QVALUE - scattering angle
      //! @param INTENSITY - Saxs intensity
      //! @param ERROR - measurement error
      //! @param COMPUTED_INTENSITY - computed intensity from Debye Formula
      SasScatteringPoint
      (
        const double QVALUE,
        const double INTENSITY,
        const double MEASUREMENT_ERROR
      );

      //! @brief Clone function
      //! @return pointer to new SaxsScatteringPoint
      SasScatteringPoint *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief Accessor Function to Private data variable
      //! @return the Intensity value
      const double &GetIntensity() const
      {
        return m_Intensity;
      }

      //! @brief Mutator function to change Intensity
      //! @param INTENSITY pass by const reference to allow scaling factor adjustment
      void SetIntensity( const double &INTENSITY)
      {
        m_Intensity = INTENSITY;
      }

      //! @brief Accessor Function to Private data variable
      //! @return The scattering angle (q-value)
      const double &GetQvalue() const
      {
        return m_Qvalue;
      }

      //! @brief Accessor Function to Private data variable
      //! @return The error value
      const double &GetError() const
      {
        return m_Error;
      }

      //! @brief Mutator function to change error
      //! @param ERROR_VALUE For simulated data, error must be computed.
      void SetError( const double &ERROR_VALUE)
      {
        m_Error = ERROR_VALUE;
      }

      //! @brief Boolean function to return state of error value
      //! @return return true if the error is defined
      bool IsErrorDefined() const;

    ////////////////
    // operators  //
    ////////////////

      //! @brief compare two SaxsScatteringPoint
      //! @param POINT the point to compare to this point
      //! @return if *this and POINT have identical data
      bool operator ==( const SasScatteringPoint &POINT) const
      {
        return
            m_Qvalue == POINT.m_Qvalue &&
            m_Intensity == POINT.m_Intensity &&
            m_Error == POINT.m_Error;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class SasScatteringPoint
  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_SAS_SCATTERING_POINT_H_
