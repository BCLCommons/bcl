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

#ifndef BCL_RESTRAINT_SAS_SCATTERING_DATA_H_
#define BCL_RESTRAINT_SAS_SCATTERING_DATA_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_restraint_handler_base.h"
#include "bcl_restraint_sas_scattering_point.h"
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
    //! @class SasScatteringData
    //! @brief Stores SAXS raw data
    //! @details Container and processing class for saxs data (q, intensity, error) to be used during folding runs
    //!
    //! @see @link example_restraint_saxs_scattering_data.cpp @endlink
    //! @author putnamdk
    //! @date April 26, 2013
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SasScatteringData :
      public HandlerBase< SasScatteringData>
    {

    private:

    //////////
    // data //
    //////////

      //! storage for saxs data, < q-value, intensity, error>
      storage::Vector< SasScatteringPoint> m_Data;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      typedef storage::Vector< SasScatteringPoint>::iterator       scattering_iterator;
      typedef storage::Vector< SasScatteringPoint>::const_iterator const_scattering_iterator;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor from extension
      SasScatteringData( const std::string &EXTENSION = ".saxs");

      //! @brief constructor from given input data
      explicit SasScatteringData
      (
        const storage::Vector< SasScatteringPoint> &INIT_DATA
      );

      //! @brief Clone function
      //! @return pointer to new SasScatteringData
      SasScatteringData *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns short name for the class used over the command line
      //! @return short name for the class used over the command line
      const std::string &GetAlias() const;

      //! @brief returns the scattering data
      //! @return the scattering data
      const storage::Vector< SasScatteringPoint> &GetScatteringData() const
      {
        return m_Data;
      }

      //! @brief return iterator on begin
      //! @return iterator pointing to the beginning of the container, i.e. the first element
      scattering_iterator Begin()
      {
        return m_Data.Begin();
      }

      //! @brief return const_iterator on begin
      //! @return const_iterator pointing to the beginning of the container, i.e. the first element
      const_scattering_iterator Begin() const
      {
        return m_Data.Begin();
      }

      //! @brief return iterator on end
      //! @return iterator pointing to the end of the container, i.e. behind the last element
      scattering_iterator End()
      {
        return m_Data.End();
      }

      //! @brief return const_iterator on end
      //! @return const_iterator pointing to the end of the container, i.e. behind the last element
      const_scattering_iterator End() const
      {
        return m_Data.End();
      }

      //! @brief Get SaxsScatteringPoint at a specific position
      //! @return SaxsScatteringPoint at a specific position
      const SasScatteringPoint &GetScatteringLocation( int LOCATION) const
      {
        return m_Data( LOCATION);
      }

      //! @brief Get Size of SasScatteringData
      //! @return size of SasScatteringData
      const size_t GetScatteringSize() const
      {
        return m_Data.GetSize();
      }

      //! @brief bool test for error values
      //! @return true if error is defined for all values of the dataset
      const bool IsErrorDefined() const;

      //! @brief pushback function to add object to m_Data vector
      //! @param DATAPOINT_OBJECT DataPoint values Q, I, and Error
      void PushBackScattering( const SasScatteringPoint &DATAPOINT_OBJECT);

    /////////////////
    // operations  //
    /////////////////

      //! @brief preallocate scattering data memory
      //! @param SIZE size to preallocate
      void AllocateScatteringMemory( const size_t &SIZE);

    private:

      //! @brief read experimental data from BCL
      //! @param ISTREAM input data stream
      void ReadBCLProcessedData( std::istream &ISTREAM);

      //! @brief read experimental data from Crysol
      //! @param ISTREAM input data stream
      void ReadCrysolData( std::istream &ISTREAM);

      //! @brief read experimental data
      //! @param ISTREAM input data stream
      void ReadExperimentalData( std::istream &ISTREAM);

      //! @brief read experimental data from Gnom
      //! @param ISTREAM input data stream
      //! @param INTENSITY_COLUMN column in gnom file to read intensity from
      void ReadGnomData( std::istream &ISTREAM, size_t &INTENSITY_COLUMN, bool USE_EXTRAPOLATION = false);

    //////////////////////
    // input and output //
    //////////////////////

    public:

      //! @brief read experimental data from BCL needs to be public function
      //! @param ISTREAM input data stream
      // void ReadBCLModel( std::istream &ISTREAM);

      //! @brief reads in the member data from a formatted file containing 3 columns:
      //! @brief scattering angle q ( 4*Pi*sin(theta)/lambda) where lambda is measured in Angstroms, I is the intensity
      //! @brief at a given a value, E is the experimental error.
      //! @brief Crysol generated files must used the paramater /dro 0.0.  The algorithm does not support adding the
      //! @brief hydration layer around the molecule.
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &ReadFromDataFile( std::istream &ISTREAM);

      //! @brief reads saxs restraints from an input stream
      //! @param input stream to read the restraints from
      //! @return the read in restraints
      SasScatteringData ReadRestraints( std::istream &ISTREAM) const;

      //! @brief reads the computed SAXS profile from the fit density curve from gnom
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &ReadFitFromGnom( std::istream &ISTREAM);

      //! @brief writes out the member data from a formatted file containing 3 columns
      //! @param OSTREAM output stream
      //! @return ostream which was written to
      std::ostream &WriteToDataFile( std::ostream &OSTREAM, bool HEADER = true) const;

      //! @brief writes out the member data to a specified file name
      //! @param FILENAME the name of the output file
      void WriteToDataFileName( const std::string &FILENAME, const bool &HEADER) const;

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

    }; // class SasScatteringData

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_SAS_SCATTERING_DATA_H_
