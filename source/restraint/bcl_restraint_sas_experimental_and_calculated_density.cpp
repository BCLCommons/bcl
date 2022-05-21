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
#include "restraint/bcl_restraint_sas_experimental_and_calculated_density.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "math/bcl_math.h"
#include "math/bcl_math_cubic_spline.h"
#include "math/bcl_math_cubic_spline_damped.h"
#include "restraint/bcl_restraint_sas_analysis.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SasExperimentalAndCalculatedDensity::s_Instance
    (
      GetObjectInstances().AddInstance( new SasExperimentalAndCalculatedDensity())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SasExperimentalAndCalculatedDensity::SasExperimentalAndCalculatedDensity() :
      m_ExperimentalDensity(),
      m_CalculatedDensity()
    {
    }

    //!@brief copy constructor
    SasExperimentalAndCalculatedDensity::SasExperimentalAndCalculatedDensity( const SasExperimentalAndCalculatedDensity &RHS)
    {
      m_ExperimentalDensity = RHS.GetExperimentalDensity();
      m_ExperimentalDmax  = RHS.GetExperimentalDensity().GetDmax();
      m_CalculatedDensity = RHS.GetCalculatedDensity();
      m_CalculatedDmax = RHS.GetCalculatedDensity().GetDmax();
    }

    //! @brief constructor from given input data
    SasExperimentalAndCalculatedDensity::SasExperimentalAndCalculatedDensity
    (
      const SasDensityData &EXPERIMENTAL,
      const SasDensityData &CALCULATED
    ) :
      m_ExperimentalDensity( EXPERIMENTAL),
      m_CalculatedDensity( CALCULATED)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SasExperimentalAndCalculatedDensity
    SasExperimentalAndCalculatedDensity *SasExperimentalAndCalculatedDensity::Clone() const
    {
      return new SasExperimentalAndCalculatedDensity( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SasExperimentalAndCalculatedDensity::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  /////////////////
  // operations  //
  /////////////////

    //! @brief scale the calculated and experimental data
    //! @param SCALING_FACTOR - the value to scale all of the Intensity values by
    void SasExperimentalAndCalculatedDensity::ScaleExperimentalDensity( const double &EXP_SCALING_FACTOR)
    {
      // Scale Experimental and Calculated Data
      m_ExperimentalDensity = SasAnalysis::ScaleDensityData( m_ExperimentalDensity, EXP_SCALING_FACTOR);
    }

    //! @brief scale the calculated and experimental data
    //! @param SCALING_FACTOR - the value to scale all of the Intensity values by
    void SasExperimentalAndCalculatedDensity::ScaleCalculatedDensity( const double &CAL_SCALING_FACTOR)
    {
      // Scale Experimental and Calculated Data
      m_CalculatedDensity = SasAnalysis::ScaleDensityData( m_CalculatedDensity, CAL_SCALING_FACTOR);
    }

    //! @brief align the max peak of the experimental and calculated curves
    void SasExperimentalAndCalculatedDensity::ShiftDensity()
    {
      double cal_hx_max( m_CalculatedDensity.GetHxmax());
      double exp_hx_max( m_ExperimentalDensity.GetHxmax());

      double shift( exp_hx_max - cal_hx_max);

      m_CalculatedDensity = SasAnalysis::ShiftDensity( m_CalculatedDensity, shift);

      // Compute spline for experimental Density
      const size_t data_size( m_ExperimentalDensity.GetDensitySize());

      // initialize math vectors with size of data set
      linal::Vector< double> data_values( data_size);

      const double delta( m_ExperimentalDensity.GetBinSize());

      // populate math vectors with values from SAXS_DATA
      size_t pos( 0);

      for
      (
        storage::Vector< SasDistanceDensityPoint>::const_iterator
         data_itr( m_ExperimentalDensity.Begin()),
         data_itr_end( m_ExperimentalDensity.End());
        data_itr != data_itr_end;
        ++data_itr, ++pos
      )
      {
        data_values( pos) = data_itr->GetDensity();
      }

      // Create CublicSpline objects for the dataset
      math::CubicSplineDamped data_function;

      // Train the data set
      data_function.Train
      (
        m_ExperimentalDensity.Begin()->GetRvalue(),
        delta,
        data_values
      );

      // iterate over calculated data using the splines from the experimental data
      for
      (
        storage::Vector< SasDistanceDensityPoint>::iterator
         cal_data_itr( m_CalculatedDensity.Begin()),
         exp_data_itr( m_ExperimentalDensity.Begin()),
         cal_data_itr_end( m_CalculatedDensity.End());
        cal_data_itr != cal_data_itr_end;
        ++cal_data_itr, ++exp_data_itr
      )
      {
        if( cal_data_itr->GetRvalue() < m_ExperimentalDensity.GetDmax())
        {
          exp_data_itr->SetDensity( data_function( cal_data_itr->GetRvalue()));
          //exp_data_itr->SetRvalue ( cal_data_itr->GetRvalue());
        }
        else
        {
          exp_data_itr->SetDensity( 0.0);
          //exp_data_itr->SetRvalue ( cal_data_itr->GetRvalue());
        }
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write to filename in three columns : R-value, density, experimental error
    //! @param filename to write to
    std::ostream &SasExperimentalAndCalculatedDensity::WriteToOstream( std::ostream &OSTREAM) const
    {
      OSTREAM << "R_Value Experimental_Density Experimental_Error Computed_Density" << '\n';

      // iterate over experimental and calculated data to get values
      for
      (
        storage::Vector< SasDistanceDensityPoint>::const_iterator
          exp_data_itr( GetExperimentalDensity().Begin()),
          cal_data_itr( GetCalculatedDensity().Begin()),
          exp_data_itr_end( GetExperimentalDensity().End());
        exp_data_itr != exp_data_itr_end;
        ++exp_data_itr, ++cal_data_itr
      )
      {
        // write the data to the ofstream
        OSTREAM << exp_data_itr->GetRvalue()    << ' '
                << exp_data_itr->GetDensity()   << ' '
                << exp_data_itr->GetError()     << ' '
                << cal_data_itr->GetDensity()   << '\n';
      }
      return OSTREAM;
    }

    //! @brief write to filename in three columns : R-value, experimental density, and calculated density
    //! @param filename to write to
    void SasExperimentalAndCalculatedDensity::WriteToGnuplotFileName( const std::string &FILENAME) const
    {
      io::OFStream write;
      io::File::MustOpenOFStream( write, FILENAME);

      // write the data
      WriteToOstream( write);

      // close the ofstream
      io::File::CloseClearFStream( write);
    }

    //! @brief write to filename in four columns : Q-value, experimental intensity, error and calculated intensity
    //! @param filename to write to
    void SasExperimentalAndCalculatedDensity::WriteToFileName( const std::string &FILENAME) const
    {
      io::OFStream write;
      io::File::MustOpenOFStream( write, FILENAME);

      // write the data
      WriteToOstream( write);

      // close the ofstream
      io::File::CloseClearFStream( write);
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SasExperimentalAndCalculatedDensity::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ExperimentalDensity, ISTREAM);
      io::Serialize::Read( m_CalculatedDensity, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SasExperimentalAndCalculatedDensity::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ExperimentalDensity, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_CalculatedDensity, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
