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

#ifndef BCL_RESTRAINT_SAS_DISTANCE_DENSITY_POINT_H_
#define BCL_RESTRAINT_SAS_DISTANCE_DENSITY_POINT_H_

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
    //! @class SasDistanceDensityPoint
    //! @brief Storage Class to store a transformed data point from I(q) to P(r)
    //! @details Class to explicitly list Rvalue, density and error
    //!
    //! @see @link example_restraint_saxs_distance_density_point.cpp @endlink
    //! @author putnamdk
    //! @date May 23, 2014
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SasDistanceDensityPoint :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! @brief Qvalue represents the scattering angle from Saxs Data
      double m_Rvalue;

      //! @brief Intensity represents the intensity for a given scattering angle for Saxs Data
      double m_Density;

      //! @brief Error is the measurement error
      double m_Error;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SasDistanceDensityPoint();

      //! @brief constructor from given input data
      //! @param QVALUE - scattering angle
      //! @param INTENSITY - Saxs intensity
      //! @param ERROR - measurement error
      //! @param COMPUTED_INTENSITY - computed intensity from Debye Formula
      SasDistanceDensityPoint
      (
        const double RVALUE,
        const double DENSITY,
        const double MEASUREMENT_ERROR
      );

      //! @brief Clone function
      //! @return pointer to new SasDistanceDensityPoint
      SasDistanceDensityPoint *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief Accessor Function to Private data variable
      //! @return the Density value
      const double &GetDensity() const
      {
        return m_Density;
      }

      //! @brief Mutator function to change Density
      //! @param DENSITY value to set for given point
      void SetDensity( const double &DENSITY)
      {
        m_Density = DENSITY;
      }

      //! @brief Accessor Function to Private data variable
      //! @return The distance in angstrom (r-value)
      const double &GetRvalue() const
      {
        return m_Rvalue;
      }

      //! @brief Mutator function to set Rvalue
      //! @param RVALUE value to set for given point
      void SetRvalue( const double &RVALUE)
      {
        m_Rvalue = RVALUE;
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

    ////////////////
    // operators  //
    ////////////////

      //! @brief compare two SasDistanceDensityPoint
      //! @param POINT the point to compare to this point
      //! @return if *this and POINT have identical data
      bool operator ==( const SasDistanceDensityPoint &POINT) const
      {
        return
          m_Rvalue == POINT.m_Rvalue &&
          m_Density == POINT.m_Density &&
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

    }; // class SasDistanceDensityPoint
  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_SAS_DISTANCE_DENSITY_POINT_H_
