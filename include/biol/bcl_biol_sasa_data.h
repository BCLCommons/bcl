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

#ifndef BCL_BIOL_SASA_DATA_H_
#define BCL_BIOL_SASA_DATA_H_

// include the namespace header
#include "bcl_biol.h"
#include "bcl_biol_sasa_point.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_interface.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SasaData
    //! @brief Stores Solvent Accessible Surface Area data from MSMS
    //! @details Container class for sasa data
    //! @details (Atom Number, Solvent Excluded Surface (SES), Solvent Accessible Surface (SAS)
    //!
    //! @see @link example_biol_sasa_data.cpp @endlink
    //! @author putnamdk
    //! @date April 9, 2013
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SasaData :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! storage for sasa data, < Atom Number, Solvent Excluded Surface, Solvent Accessible Surfacer>
      storage::Vector< SasaPoint> m_SasaData;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      typedef storage::Vector< SasaPoint>::iterator       iterator;
      typedef storage::Vector< SasaPoint>::const_iterator const_iterator;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SasaData();

      //! @brief constructor from given input data
      explicit SasaData( const storage::Vector< SasaPoint> &INIT_DATA);

      //! @brief Clone function
      //! @return pointer to new SasaData
      SasaData *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief flag for specifying input data format
      //! @return command line flag for specifying input data format
      static util::ShPtr< command::FlagInterface> &GetFlagInputFormat();

      //! @brief returns the storage data
      //! @return the storage data
      const storage::Vector< SasaPoint> &GetData() const
      {
        return m_SasaData;
      }

      //! @brief return iterator on begin
      //! @return iterator pointing to the beginning of the container, i.e. the first element
      iterator Begin()
      {
        return m_SasaData.Begin();
      }

      //! @brief return const_iterator on begin
      //! @return const_iterator pointing to the beginning of the container, i.e. the first element
      const_iterator Begin() const
      {
        return m_SasaData.Begin();
      }

      //! @brief return iterator on end
      //! @return iterator pointing to the end of the container, i.e. behind the last element
      iterator End()
      {
        return m_SasaData.End();
      }

      //! @brief return const_iterator on end
      //! @return const_iterator pointing to the end of the container, i.e. behind the last element
      const_iterator End() const
      {
        return m_SasaData.End();
      }

      //! @brief pushback function to add object to m_Data vector
      //! @param DATAPOINT_OBJECT DataPoint values Q, I, and Error
      void PushBack( const SasaPoint &DATAPOINT_OBJECT);

    /////////////////
    // operations  //
    /////////////////

      //! @brief preallocate memory
      //! @param SIZE size to preallocate
      void AllocateMemory( const size_t &SIZE);

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief reads in the member data from a formatted file containing 3 columns:
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &ReadFromDataFile( std::istream &ISTREAM);

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

    }; // class SasaData

  } // namespace biol
} // namespace bcl

#endif // BCL_BIOL_SASA_DATA_H_
