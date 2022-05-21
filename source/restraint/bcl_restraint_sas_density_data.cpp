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
#include "linal/bcl_linal_vector.h"
#include "math/bcl_math_histogram.h"
#include "restraint/bcl_restraint_sas_analysis.h"
#include "restraint/bcl_restraint_sas_density_data.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SasDensityData::s_Instance
    (
      util::Enumerated< HandlerBase< SasDensityData> >::AddInstance
      (
        new SasDensityData()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SasDensityData::SasDensityData( const std::string &EXTENSION) :
      HandlerBase< SasDensityData>( EXTENSION),
      m_DensityDistributionData(),
      m_BinSize( 0),
      m_BinNumber( 0),
      m_Dmax( 0),
      m_Hmax( 0),
      m_Hxmax( 0)
    {
    }

    //! @brief constructor from given input data
    SasDensityData::SasDensityData
    (
      const storage::Vector< SasDistanceDensityPoint> &DENSITY_DATA
    ) :
      HandlerBase< SasDensityData>( ".pofr"),
      m_DensityDistributionData( DENSITY_DATA),
      m_Dmax( 0),
      m_Hmax( 0),
      m_Hxmax( 0)
    {
      storage::Vector< SasDistanceDensityPoint>::const_iterator data_itr( DENSITY_DATA.Begin());
      data_itr++;

      m_BinSize = data_itr->GetRvalue();
      m_BinNumber = DENSITY_DATA.GetSize();
      m_Dmax = SasAnalysis::ComputeDmax( *this);
      m_Hmax = SasAnalysis::FindDensitymax( *this);
      m_Hxmax = SasAnalysis::FindxDensitymax( *this);
    }

    //! @brief constructor from a histogram
    SasDensityData::SasDensityData
    (
      const math::Histogram &DENSITY_HISTOGRAM,
      const double &DMAX
    ) :
      m_Dmax( DMAX)
    {
      linal::Vector< double> binning( DENSITY_HISTOGRAM.GetBinning());
      linal::Vector< double> counts( DENSITY_HISTOGRAM.GetHistogram());

      double h_max( 0.0);
      double hx_max( 0.0);

      // Convert Histogram to SasDensityData
      for
      (
        const double *x( binning.Begin()), *x_end( binning.End()), *y( counts.Begin()), *y_end( counts.End());
        x != x_end && y != y_end;
        ++x, ++y
      )
      {
        // 0.5 was subtracted from x to align the bin on the left boundary to match experimental data
        m_DensityDistributionData.PushBack( SasDistanceDensityPoint( *x - 0.5, *y, util::GetUndefined< double>()));
        if( *y > h_max)
        {
          h_max = *y;
          hx_max = *x - 0.5;
        }
      }

      m_BinSize = DENSITY_HISTOGRAM.GetBinSize();
      m_BinNumber = DENSITY_HISTOGRAM.GetNumberOfBins();
      m_Hmax = h_max;
      m_Hxmax = hx_max;
    }

    //! @brief Clone function
    //! @return pointer to new SasDensityData
    SasDensityData *SasDensityData::Clone() const
    {
      return new SasDensityData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SasDensityData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief pushback function to add object to Dataset vector
    //! @param DATAPOINT_OBJECT  //! @brief pushback function to add object to m_Data vectorDataPoint values R, P(R), and Error
    void SasDensityData::PushBackDensity( const SasDistanceDensityPoint &DENSITY_POINT_OBJECT)
    {
      m_DensityDistributionData.PushBack( DENSITY_POINT_OBJECT);
    }

    //! @param VALUE to set binsize to
    void SasDensityData::SetBinSize( const double &BIN_SIZE)
    {
      m_BinSize = BIN_SIZE;
    }

    //! @param VALUE to set binsize to
    void SasDensityData::SetBinNumber( const size_t &BIN_NUMBER)
    {
      m_BinNumber = BIN_NUMBER;
    }

    //! @param VALUE to set binsize to
    void SasDensityData::SetDmax( const double &DMAX)
    {
      m_Dmax = DMAX;
    }

    //! @param VALUE to set Hmax to
    void SasDensityData::SetHmax( const double &HMAX)
    {
      m_Hmax = HMAX;
    }

    //! @param VALUE to set Hmax to
    void SasDensityData::SetHxmax( const double &HXMAX)
    {
      m_Hxmax = HXMAX;
    }

    const double SasDensityData::ComputeHmax() const
    {
      return SasAnalysis::FindDensitymax( *this);
    }

    const double SasDensityData::ComputeHxmax() const
    {
      return SasAnalysis::FindxDensitymax( *this);
    }

  /////////////////
  // operations  //
  /////////////////

    //! @brief preallocate Density Data memory
    //! @param SIZE size to preallocate
    void SasDensityData::AllocateDensityMemory( const size_t &SIZE)
    {
      m_DensityDistributionData.AllocateMemory( SIZE);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read experimental data from BCLModel
    //! @param ISTREAM input data stream
    void SasDensityData::ReadBCLModel( std::istream &ISTREAM)
    {
      // String to hold line data
      std::string read_line;

      // skip the header line
      std::getline( ISTREAM, read_line);

      // iterate over each line in the file
      while( std::getline( ISTREAM, read_line).good())
      {
        // split the string
        const storage::Vector< std::string> split_line( util::SplitString( read_line, " "));

        double p( util::ConvertStringToNumericalValue< double>( split_line( 0)));
        double density( util::ConvertStringToNumericalValue< double>( split_line( 1)));
        double error( util::ConvertStringToNumericalValue< double>( split_line( 2)));

        // push back the data p, p(r), Error(p)
        m_DensityDistributionData.PushBack( SasDistanceDensityPoint( p, density, error));
      } // close while loop
    } // close function

    //! @brief read experimental data from Gnom
    //! @param ISTREAM input data stream
    void SasDensityData::ReadGnomData( std::istream &ISTREAM)
    {

      std::string read_line;

      // read first line of input file
      std::getline( ISTREAM, read_line);

      // split the string based on white space as a delimeter
      storage::Vector< std::string> Line( util::SplitString( read_line, " "));

      // get the size of the first line
      size_t line_size( Line.GetSize());

      // paramaters to read through gnom filetype
      bool readDDF( false);

      // while the end of the file is not reached
      while( std::getline( ISTREAM, read_line).good())
      {
        // split the string
        const storage::Vector< std::string> split_line( util::SplitString( read_line, " "));

        line_size = split_line.GetSize();
        if( line_size == 0)
        {
          continue;
        }
        if( split_line( 0) == "R" && split_line( 1) == "P(R)")
        {
          readDDF = true;
        }
        else if( split_line( 0) == "Reciprocal")
        {
          readDDF = false;
        }
        else if( readDDF && line_size == 3)
        {
          // push back the data
          m_DensityDistributionData.PushBack
          (
            SasDistanceDensityPoint
            (
              util::ConvertStringToNumericalValue< double>( split_line( 0)),  // R-value
              util::ConvertStringToNumericalValue< double>( split_line( 1)),  // P(R)-value
              util::ConvertStringToNumericalValue< double>( split_line( 2))   // error
            )
          );
        }
      } // close while loop

      storage::Vector< SasDistanceDensityPoint>::const_iterator data_itr( m_DensityDistributionData.Begin());
      data_itr++;

      m_BinSize = data_itr->GetRvalue();
      m_BinNumber = m_DensityDistributionData.GetSize();
      m_Dmax = SasAnalysis::ComputeDmax( *this);
      m_Hmax = SasAnalysis::FindDensitymax( *this);
      m_Hxmax = SasAnalysis::FindxDensitymax( *this);
    } // close the function

    //! @brief reads saxs restraints from an input stream
    //! @brief input stream to read the restraints from
    //! @return the read in restraints
    SasDensityData SasDensityData::ReadRestraints( std::istream &ISTREAM) const
    {
      SasDensityData data;
      data.ReadFromDataFile( ISTREAM);
      return data;
    }

    //! @brief reads in the member data from a formatted file containing 3 columns:
    //! @brief scattering angle q ( 4*Pi*sin(theta)/lambda) where lambda is measured in Angstroms, I is the intensity
    //! @brief at a given a value, E is the experimental error.
    //! @brief Crysol generated files must used the paramater /dro 0.0.  The algorithm does not support adding the
    //! @brief hydration layer around the molecule.
    //! @param ISTREAM input stream
    //! @param FORMAT the file format to use for reading
    //! @return istream which was read from
    std::istream &SasDensityData::ReadFromDataFile( std::istream &ISTREAM)
    {
      // build a string to hold the line information
      std::string read_line;

      enum Types { gnom, unknown};
      Types filetype( unknown);

      // read first line of input file
      std::getline( ISTREAM, read_line);

      // split the string based on white space as a delimeter
      storage::Vector< std::string> Line( util::SplitString( read_line, " "));

      // get the size of the first line
      size_t line_size( Line.GetSize());

      size_t number_blank_lines( 0);

      // ignore empty lines at the first of the file
      while( line_size <= 1)
      {
        number_blank_lines++;
        std::getline( ISTREAM, read_line);
        Line = util::SplitString( read_line, " ");
        line_size = Line.GetSize();

        if( number_blank_lines == 1000)
        {
          BCL_Exit( "Can not read file passed to SasDensityData::ReadFromDataFile, greater than 1000 blank lines", -1);
        }

      }

      // Verify the file is a GNOM file
      if( line_size >= 5 && Line( 1) == "G" && Line( 2) == "N" && Line( 3) == "O" && Line( 4) == "M")
      {
        BCL_MessageStd( "Filetype is gnom:");
        filetype = gnom;
      }
      else
      {
        BCL_Exit( "Could not read file; last read line was: " + read_line, -1);
      }

      switch( filetype)
      {
        case gnom:
          ReadGnomData( ISTREAM);
          break;
        case unknown:
        default:
          break;
      }

      size_t size( 0);

      // make sure there are r values for analysis
      size = m_DensityDistributionData.GetSize();

      BCL_Assert( size != 0, "The number of Q values are: " + util::Format()( size));

      // return the stream
      return ISTREAM;
    }

    //! @brief writes out the member data from a formatted file containing 3 columns
    //! @param OSTREAM output stream
    //! @return ostream which was written to
    std::ostream &SasDensityData::WriteToDataFile( std::ostream &OSTREAM) const
    {
      // initialize format
      const util::Format format;

      OSTREAM << format( " BCL DENSITY INPUT PARAMETERS: RValue Density Error") << '\n';

      // iterate over m_DensityDistributionData
      for
      (
        storage::Vector< SasDistanceDensityPoint>::const_iterator density_itr
        (
          m_DensityDistributionData.Begin()),
         density_itr_end( m_DensityDistributionData.End()
        );
        density_itr != density_itr_end; ++density_itr
      )
      {
        // write the data
        OSTREAM << format( density_itr->GetRvalue()) << ' '
               << format( density_itr->GetDensity()) << ' '
               << format( density_itr->GetError()) << '\n';
      }

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SasDensityData::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_DensityDistributionData, ISTREAM);
      io::Serialize::Read( m_BinSize, ISTREAM);
      io::Serialize::Read( m_BinNumber, ISTREAM);
      io::Serialize::Read( m_Dmax, ISTREAM);
      io::Serialize::Read( m_Hmax, ISTREAM);
      io::Serialize::Read( m_Hxmax, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SasDensityData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_DensityDistributionData, OSTREAM, INDENT);
      io::Serialize::Write( m_BinSize, OSTREAM, INDENT);
      io::Serialize::Write( m_BinNumber, OSTREAM, INDENT);
      io::Serialize::Write( m_Dmax, OSTREAM, INDENT);
      io::Serialize::Write( m_Hmax, OSTREAM, INDENT);
      io::Serialize::Write( m_Hxmax, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl
