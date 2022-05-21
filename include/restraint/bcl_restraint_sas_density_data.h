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

#ifndef BCL_RESTRAINT_SAS_DENSITY_DATA_H_
#define BCL_RESTRAINT_SAS_DENSITY_DATA_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_restraint_handler_base.h"
#include "bcl_restraint_sas_distance_density_point.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "command/bcl_command_flag_interface.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SasDensityData
    //! @brief Stores SAXS raw data transformed into a distance distribution function
    //! @details Container and processing class for saxs density (r, density, error) to be used for analysis
    //!
    //! @see @link example_restraint_sas_density_data.cpp @endlink
    //! @author putnamdk
    //! @date July 3, 2014
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SasDensityData :
      public HandlerBase< SasDensityData>
    {

    private:

    //////////
    // data //
    //////////

      //! storage for density distribution data, < r-value, density, error>
      storage::Vector< SasDistanceDensityPoint> m_DensityDistributionData;

      //! width of the bin
      double m_BinSize;

      //! How many bins the data has
      size_t m_BinNumber;

      //! maximum dimension of the set
      double m_Dmax;

      //! maximum height of the set
      double m_Hmax;

      //! x-coordinate of the maximum height
      double m_Hxmax;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      typedef storage::Vector< SasDistanceDensityPoint>::iterator       density_iterator;
      typedef storage::Vector< SasDistanceDensityPoint>::const_iterator const_density_iterator;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor from extension
      SasDensityData( const std::string &EXTENSION = ".pofr");

      //! @brief constructor from given input data
      explicit SasDensityData
      (
        const storage::Vector< SasDistanceDensityPoint> &DENSITY_DATA
      );

      //! @brief constructor from a histogram
      explicit SasDensityData
      (
        const math::Histogram &DENSITY_HISTOGRAM,
        const double &DMAX
      );

      //! @brief Clone function
      //! @return pointer to new SasDensityData
      SasDensityData *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the density data
      //! @return the density data
      const storage::Vector< SasDistanceDensityPoint> &GetDensityData() const
      {
        return m_DensityDistributionData;
      }

      const double &GetBinSize() const
      {
        return m_BinSize;
      }

      const size_t &GetBinNumber() const
      {
        return m_BinNumber;
      }

      const double &GetDmax() const
      {
        return m_Dmax;
      }

      const double &GetHmax() const
      {
        return m_Hmax;
      }

      const double &GetHxmax() const
      {
        return m_Hxmax;
      }

      const double ComputeHmax() const;

      const double ComputeHxmax() const;

      //! @brief return iterator on begin for density data
      //! @return iterator pointing to the beginning of the container, i.e. the first element
      density_iterator Begin()
      {
        return m_DensityDistributionData.Begin();
      }

      //! @brief return const_iterator on begin
      //! @return const_iterator pointing to the beginning of the container, i.e. the first element
      const_density_iterator Begin() const
      {
        return m_DensityDistributionData.Begin();
      }

      //! @brief return iterator on end
      //! @return iterator pointing to the end of the container, i.e. behind the last element
      density_iterator End()
      {
        return m_DensityDistributionData.End();
      }

      //! @brief return const_iterator on end
      //! @return const_iterator pointing to the end of the container, i.e. behind the last element
      const_density_iterator End() const
      {
        return m_DensityDistributionData.End();
      }

      //! @brief Get SasDistanceDensityPoint at a specific position
      //! @return SasDistanceDensityPoint at a specific position
      const SasDistanceDensityPoint &GetDensityLocation( int LOCATION) const
      {
        return m_DensityDistributionData( LOCATION);
      }

      //! @brief Get Size of SaxsDensityDistributionData
      //! @return size of SaxsDensityDistributionData
      const size_t GetDensitySize() const
      {
        return m_DensityDistributionData.GetSize();
      }

      //! @brief pushback function to add object to m_Data vector
      //! @param DATAPOINT_OBJECT DataPoint values Q, I, and Error
      void PushBackDensity( const SasDistanceDensityPoint &DATAPOINT_OBJECT);

      //! @param VALUE to set binsize to
      void SetBinSize( const double &BIN_SIZE);

      //! @param VALUE to set binsize to
      void SetBinNumber( const size_t &BIN_NUMBER);

      //! @param VALUE to set binsize to
      void SetDmax( const double &DMAX);

      //! @param VALUE to set hmax to
      void SetHmax( const double &HMAX);

      //! @param VALUE to set hxmax to
      void SetHxmax( const double &HXMAX);

    /////////////////
    // operations  //
    /////////////////

      //! @brief preallocate density data memory
      //! @param SIZE size to preallocate
      void AllocateDensityMemory( const size_t &SIZE);

    private:

      //! @brief read experimental data from Gnom
      //! @param ISTREAM input data stream
      void ReadGnomData( std::istream &ISTREAM);

    //////////////////////
    // input and output //
    //////////////////////

    public:

      //! @brief reads pofr restraints from an input stream
      //! @param input stream to read the restraints from
      //! @return the read in restraints
      SasDensityData ReadRestraints( std::istream &ISTREAM) const;

      //! @brief read experimental data from BCL needs to be public function
      //! @param ISTREAM input data stream
      void ReadBCLModel( std::istream &ISTREAM);

      //! @brief reads in the member data from a formatted file containing 3 columns:
      //! @brief scattering angle q ( 4*Pi*sin(theta)/lambda) where lambda is measured in Angstroms, I is the intensity
      //! @brief at a given a value, E is the experimental error.
      //! @brief Crysol generated files must used the paramater /dro 0.0.  The algorithm does not support adding the
      //! @brief hydration layer around the molecule.
      //! @param ISTREAM input stream
      //! @param FORMAT the file format to use for reading
      //! @return istream which was read from
      std::istream &ReadFromDataFile
      (
        std::istream &ISTREAM
      );

      //! @brief writes out the member data from a formatted file containing 3 columns
      //! @param OSTREAM output stream
      //! @return ostream which was written to
      std::ostream &WriteToDataFile( std::ostream &OSTREAM) const;

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

    }; // class SasDensityData

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_SAS_DENSITY_DATA_H_
