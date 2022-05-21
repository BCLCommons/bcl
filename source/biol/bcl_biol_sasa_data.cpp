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
#include "biol/bcl_biol_sasa_data.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "math/bcl_math_cubic_spline.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SasaData::s_Instance
    (
      GetObjectInstances().AddInstance( new SasaData())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SasaData::SasaData() :
      m_SasaData()
    {
    }

    //! @brief constructor from given input data
    SasaData::SasaData( const storage::Vector< SasaPoint> &INIT_DATA) :
      m_SasaData( INIT_DATA)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SasaData
    SasaData *SasaData::Clone() const
    {
      return new SasaData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SasaData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief pushback function to add object to Dataset vector
    //! @param DATAPOINT_OBJECT DataPoint values Q, I, and Error
    void SasaData::PushBack( const SasaPoint &DATAPOINT_OBJECT)
    {
      m_SasaData.PushBack( DATAPOINT_OBJECT);
    }

  /////////////////
  // operations  //
  /////////////////

    //! @brief preallocate memory
    //! @param SIZE size to preallocate
    void SasaData::AllocateMemory( const size_t &SIZE)
    {
      m_SasaData.AllocateMemory( SIZE);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads in the member data from a formatted file containing 3 columns:
    //! @param ISTREAM input stream
    //! @param FORMAT the file format to use for reading
    //! @return istream which was read from
    std::istream &SasaData::ReadFromDataFile( std::istream &ISTREAM)
    {
      // build a string to hold the line information
      std::string read_line;

      // skip the first line
      std::getline( ISTREAM, read_line);

      // while the end of the file is not reached
      while( std::getline( ISTREAM, read_line).good())
      {
        // split the string
        const storage::Vector< std::string> split_line( util::SplitString( read_line, " "));

        // push back the data
        m_SasaData.PushBack
        (
          SasaPoint
            (
              util::ConvertStringToNumericalValue< double>( split_line( 0)),    // Atom Number
              util::ConvertStringToNumericalValue< double>( split_line( 1)),    // Solvent Excluded Surface
              util::ConvertStringToNumericalValue< double>( split_line( 2))     // Solvent Accessible Surface
            )
         );
        }

      size_t size( m_SasaData.GetSize());

      BCL_Assert( size != 0, "The number of Sasa values are: " + util::Format()( size));

      // return the stream
      return ISTREAM;
    }

    //! @brief writes out the member data from a formatted file containing 3 columns
    //! @param OSTREAM output stream
    //! @return ostream which was written to
    std::ostream &SasaData::WriteToDataFile( std::ostream &OSTREAM) const
    {
      // initialize format
      const util::Format format;

      // iterate over m_SasaData
      for
      (
        storage::Vector< SasaPoint>::const_iterator itr( m_SasaData.Begin()), itr_end( m_SasaData.End());
         itr != itr_end; ++itr
      )
      {
        // write the data
        OSTREAM << format( itr->GetAtomNumber()) << '\t'
                << format( itr->GetSolventExcludedSurface()) << '\t'
                << format( itr->GetSolventAccessibleSurface()) << '\n';
      }

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SasaData::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SasaData, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SasaData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SasaData, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
