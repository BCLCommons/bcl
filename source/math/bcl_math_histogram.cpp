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
#include "math/bcl_math_histogram.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_function_interface.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    //class Histogram

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Histogram::s_Instance
    (
      GetObjectInstances().AddInstance( new Histogram())
    );

    //! @brief string to indicate the lower boundary in output
    //! @return string to indicate the lower boundary in output
    const std::string &Histogram::GetLowerBoundaryString()
    {
      static std::string s_lower_boundary_string( "...<");
      return s_lower_boundary_string;
    }

    //! @brief string to indicate a bin
    //! @return string to indicate a bin
    const std::string &Histogram::GetBinString()
    {
      static std::string s_bin_string( "<..>");
      return s_bin_string;
    }

    //! @brief string to indicate the upper boundary in output
    //! @return string to indicate the upper boundary in output
    const std::string &Histogram::GetUpperBoundaryString()
    {
      static std::string s_upper_boundary_string( ">...");
      return s_upper_boundary_string;
    }

    //! @brief flag for setting the minimum value of the histogram
    //! @return flag for setting the minimum value of the histogram
    const util::ShPtr< command::FlagInterface> Histogram::GetFlagMin()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "histogram_min",
          "The minimum value to be binned by the histogram.",
          command::Parameter( "minimum", "double which is the minimum value to be binned by the histogram", "0")
        )
      );

      return s_flag;
    }

    //! @brief flag for setting the bin size of the histogram
    //! @return flag for setting the bin size of the histogram
    const util::ShPtr< command::FlagInterface> Histogram::GetFlagBinSize()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "histogram_binsize",
          "The size of the bins in the histogram.",
          command::Parameter( "binsize", "double which is the size of the bins in the histogram", "1")
        )
      );

      return s_flag;
    }

    //! @brief flag for setting the number of bins of the histogram
    //! @return flag for setting the number of bins of the histogram
    const util::ShPtr< command::FlagInterface> Histogram::GetFlagNumberOfBins()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "histogram_numbins",
          "The number of the bins in the histogram.",
          command::Parameter( "number of bins", "size_t which is the number of bins in the histogram", "10")
        )
      );

      return s_flag;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor - initializes all to 0
    Histogram::Histogram() :
      m_LowerUpperBoundaries( 0.0, 0.0),
      m_BinSize( 0.0),
      m_Histogram( 0),
      m_LowerUpperBoundariesCounts( 0.0, 0.0)
    {
    }

    //! @brief construct histogram from starting value MIN, BIN_SIZE and NUMBER_OF_BINS
    //! @param MIN the minimal value representing the left boundary
    //! @param BIN_SIZE the size (width) of one bin
    //! @param NUMBER_OF_BINS the number of bin in the histogram
    Histogram::Histogram
    (
      const double MIN,
      const double BIN_SIZE,
      const size_t NUMBER_OF_BINS
    ) :
      m_LowerUpperBoundaries( MIN, MIN + BIN_SIZE * NUMBER_OF_BINS),
      m_BinSize( BIN_SIZE),
      m_Histogram( NUMBER_OF_BINS),
      m_LowerUpperBoundariesCounts( 0, 0)
    {
    }

    //! @brief construct a histogram given a map of data with bins being the key and counts being the value
    //! @param DATA histogram data with key being the bin and the value being the counts for that bin
    Histogram::Histogram( const storage::Map< double, double> &DATA) :
      m_LowerUpperBoundaries( 0.0, 0.0),
      m_BinSize( 0.0),
      m_Histogram( DATA.GetSize()),
      m_LowerUpperBoundariesCounts( 0.0, 0.0)
    {
      // make sure the data is not empty
      if( DATA.IsEmpty())
      {
        return;
      }

      // get the information about the histogram contained in DATA
      const double smallest_bin( DATA.Begin()->first);
      const double largest_bin( DATA.ReverseBegin()->first);
      const size_t num_bins( DATA.GetSize());
      const double bin_size( num_bins > 1 ? ( largest_bin - smallest_bin) / double( num_bins - 1) : 1.0);

      // set the lower and upper values which are going to be one half the bin size above and below the smallest and
      // largest bins in DATA
      m_LowerUpperBoundaries.First() = smallest_bin - bin_size / 2.0;
      m_LowerUpperBoundaries.Second() = largest_bin + bin_size / 2.0;

      BCL_MessageDbg( "lower upper boundaries " + util::Format()( m_LowerUpperBoundaries));

      m_BinSize = bin_size;

      double *ptr( m_Histogram.Begin());

      // fill m_Histogram
      for
      (
        storage::Map< double, double>::const_iterator itr( DATA.Begin()), itr_end( DATA.End()); itr != itr_end;
        ++itr, ++ptr
      )
      {
        ( *ptr) = itr->second;
      }

      BCL_MessageDbg( "histogram " + util::Format()( *this));
    }

    //! copy constructor
    Histogram *Histogram::Clone() const
    {
      return new Histogram( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Histogram::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! get number of bins
    size_t Histogram::GetNumberOfBins() const
    {
      return m_Histogram.GetSize();
    }

    //! get binsize
    double const &Histogram::GetBinSize() const
    {
      return m_BinSize;
    }

    //! GetBoundaries
    const storage::VectorND< 2, double> &Histogram::GetBoundaries() const
    {
      return m_LowerUpperBoundaries;
    }

    //! GetBoundaries counts
    const storage::VectorND< 2, double> &Histogram::GetBoundariesCounts() const
    {
      return m_LowerUpperBoundariesCounts;
    }

    //! GetHistogram
    const linal::Vector< double> &Histogram::GetHistogram() const
    {
      return m_Histogram;
    }

    //! @brief a vector of the bin coordinates of the histogram
    //! @return the center of each bin in a vector
    linal::Vector< double> Histogram::GetBinning() const
    {
      return linal::FillVector< double>( m_Histogram.GetSize(), m_LowerUpperBoundaries.First() + 0.5 * m_BinSize, m_BinSize);
    }

    //! get sum of all counts
    double Histogram::GetSumOfAllCounts() const
    {
      return m_Histogram.Sum() + m_LowerUpperBoundariesCounts.First() + m_LowerUpperBoundariesCounts.Second();
    }

    //! This functions returns the sum of counts between [MIN, MAX)
    double Histogram::GetCountsInBetween( const double MIN, const double MAX) const
    {
      double sum( 0);

      // add counts below lower limit
      if( MIN < m_LowerUpperBoundaries.First())
      {
        sum += m_LowerUpperBoundariesCounts.First();
      }

      // add counts above higher limit
      if( MAX >= m_LowerUpperBoundaries.Second())
      {
        sum += m_LowerUpperBoundariesCounts.Second();
      }

      // add counts within the given interval
      for
      (
        const double
          *ptr( m_Histogram.Begin() + size_t( ( MIN - m_LowerUpperBoundaries.First()) / m_BinSize)),
          *ptr_end( m_Histogram.Begin() + size_t( ( MAX - m_LowerUpperBoundaries.First()) / m_BinSize));
        ptr != ptr_end;
        ++ptr
      )
      {
        sum += *ptr;
      }

      return sum;
    }

    //! @brief set the count for a bin
    //! @param INDEX index of the bin
    //! @param COUNT the new count
    void Histogram::SetCount( const size_t INDEX, const double COUNT)
    {
      m_Histogram( INDEX) = COUNT;
    }

  ////////////////
  // operations //
  ////////////////

    //! calculate the densitymap from the data and parameters
    Histogram &Histogram::CalculateHistogram( const storage::Vector< double> &VALUES_VECTOR)
    {
      //reset all counts to zero
      Reset();

      //count
      for
      (
        std::vector< double>::const_iterator itr( VALUES_VECTOR.Begin()), itr_end( VALUES_VECTOR.End());
        itr != itr_end; ++itr
      )
      {
        PushBack( *itr);
      }
      return *this;
    }

    //! pushback a value to its position in the histogram
    //! @param VALUE the value to be considered to derive the bin
    //! @param WEIGHT the increment for the bin that belongs to VALUE - 1 by default
    void Histogram::PushBack( const double VALUE, const double WEIGHT)
    {
      //checks that value and weight is defined
      if( !util::IsDefined( VALUE) || !util::IsDefined( WEIGHT))
      {
        return;
      }

      // lower as lower boundary
      if( VALUE < m_LowerUpperBoundaries.First())
      {
        m_LowerUpperBoundariesCounts.First() += WEIGHT;
        return;
      }
      else if( VALUE > m_LowerUpperBoundaries.Second())
      {
        m_LowerUpperBoundariesCounts.Second() += WEIGHT;
        return;
      }

      // determine the nominal bin for this value
      const size_t nominal_bin( ( VALUE - m_LowerUpperBoundaries.First()) / m_BinSize);

      // higher then upper boundary; add to outer counts bin
      if( nominal_bin >= m_Histogram.GetSize())
      {
        if( nominal_bin)
        {
          m_Histogram( nominal_bin - 1) += WEIGHT;
        }
        else
        {
          m_LowerUpperBoundariesCounts.First() += WEIGHT / 2.0;
          m_LowerUpperBoundariesCounts.Second() += WEIGHT / 2.0;
        }
      }
      else
      {
        // within one of the bins
        m_Histogram( nominal_bin) += WEIGHT;
      }
    }

    //! @brief combine this with a given histogram by adding up all counts
    //! all parameters have to be identical for this operation to work
    //! @param HISTOGRAM histogram to be added to this one
    //! @return true if it was successful - that is, when boundaries and bin sizes do match
    bool Histogram::Combine( const Histogram &HISTOGRAM)
    {
      // allow combining, if this is empty
      if( IsEmpty())
      {
        *this = HISTOGRAM;
        return true;
      }

      // check that the parameters agree
      if
      (
           m_LowerUpperBoundaries.First() != HISTOGRAM.m_LowerUpperBoundaries.First()
        || m_LowerUpperBoundaries.Second() != HISTOGRAM.m_LowerUpperBoundaries.Second()
        || m_BinSize != HISTOGRAM.m_BinSize
      )
      {
        BCL_MessageCrt
        (
          "combining histograms with different parameters:\n"
          "left boundary:  " + util::Format()( m_LowerUpperBoundaries.First()) + " != " + util::Format()( HISTOGRAM.m_LowerUpperBoundaries.First()) + "\n" +
          "right boundary: " + util::Format()( m_LowerUpperBoundaries.Second()) + " != " + util::Format()( HISTOGRAM.m_LowerUpperBoundaries.Second()) + "\n" +
          "bin size:       " + util::Format()( m_BinSize) + " != " + util::Format()( HISTOGRAM.m_BinSize)
        );
        return false;
      }

      // add the boundary counts and the histogram
      m_LowerUpperBoundariesCounts.First()  += HISTOGRAM.m_LowerUpperBoundariesCounts.First();
      m_LowerUpperBoundariesCounts.Second() += HISTOGRAM.m_LowerUpperBoundariesCounts.Second();
      m_Histogram                           += HISTOGRAM.m_Histogram;

      // end
      return true;
    }

    //! @brief add a pseudo-count to histogram and boundaries
    //! @param VALUE pseudo-count to be added
    void Histogram::AddPseudoCount( const double VALUE)
    {
      // add pseudo count to boundaries
      m_LowerUpperBoundariesCounts.First() += VALUE;
      m_LowerUpperBoundariesCounts.Second() += VALUE;
      // add pseudo count to boundaries
      m_Histogram += VALUE;
    }

    //! @brief reset all counts
    void Histogram::Reset()
    {
      m_Histogram = 0.0;
      m_LowerUpperBoundariesCounts = storage::VectorND< 2, double>( 0, 0);
    }

    //! @brief checks if there is any count
    //! @return
    bool Histogram::IsEmpty() const
    {
      return GetSumOfAllCounts() == double( 0.0);
    }

    //! @brief determine index of last bin with count == 0 from the back
    //! @param COUNT_THRESHOLD if count > COUNT_THRESHOLD, it is considered to contain information - default = 0.0
    //! @return index of last bin with counts in it - number of bins, if right boundary has counts, 0 if nothing has count
    size_t Histogram::GetIndexOfLastInformationContainingBin( const double COUNT_THRESHOLD) const
    {
      // if the boundary count is filled, the last index with information is the last bin
      if( m_LowerUpperBoundariesCounts.Second() > COUNT_THRESHOLD)
      {
        return GetNumberOfBins();
      }

      size_t last_index_information( 0);
      for( size_t i( 0), number_bins( GetNumberOfBins()); i < number_bins; ++i)
      {
        if( m_Histogram( i) > COUNT_THRESHOLD)
        {
          last_index_information = i;
        }
      }

      // end
      return last_index_information;
    }

    //! @brief determine index of first bin with count != 0 from the front
    //! @param COUNT_THRESHOLD if count > COUNT_THRESHOLD, it is considered to contain information - default = 0.0
    //! @return index of first bin with counts in it - number of bins, 0 if all have counts, number of bins, if all have no counts
    size_t Histogram::GetIndexOfFirstInformationContainingBin( const double COUNT_THRESHOLD) const
    {
      // if the boundary count is filled, the last index with information is the last bin
      if( m_LowerUpperBoundariesCounts.First() > COUNT_THRESHOLD)
      {
        return 0;
      }

      for( size_t i( 0), number_bins( GetNumberOfBins()); i < number_bins; ++i)
      {
        if( m_Histogram( i) > COUNT_THRESHOLD)
        {
          return i;
        }
      }

      // end
      return GetNumberOfBins();
    }

    //! @brief remove the bins before the given index
    //! counts are added to left boundary
    //! @param FIRST_INDEX index of bin before which histogram is cut
    void Histogram::RemoveBinsBeforeIndex( const size_t FIRST_INDEX)
    {
      // assert that the index is smaller than the number of bins
      BCL_Assert
      (
        FIRST_INDEX <= GetNumberOfBins(),
        "Cannot remove bins before index, " + util::Format()( FIRST_INDEX) +
        " since it is larger than the number of bins"
      );

      // add the counts of the to be removed bins to the boundary
      for( size_t i( 0); i != FIRST_INDEX; ++i)
      {
        m_LowerUpperBoundariesCounts.First() += m_Histogram( i);
      }

      // remove the bins
      m_Histogram = linal::Vector< double>( m_Histogram.Begin() + FIRST_INDEX, m_Histogram.End());

      // set new lower boundary
      m_LowerUpperBoundaries.First() = m_LowerUpperBoundaries.First() + double( FIRST_INDEX) * m_BinSize;
    }

    //! @brief reset the bins before the given index
    //! counts in bins and at boundary are set to 0.0
    //! @param FIRST_INDEX index of bin before which histogram is reset
    void Histogram::ResetBinsBeforeIndex( const size_t FIRST_INDEX)
    {
      // assert that the index is smaller than the number of bins
      BCL_Assert
      (
        FIRST_INDEX <= GetNumberOfBins(),
        "Cannot reset bins before index, " + util::Format()( FIRST_INDEX) +
        " since it is larger than the number of bins"
      );

      // reset counts of the bins up to first index
      for( size_t i( 0); i != FIRST_INDEX; ++i)
      {
        m_Histogram( i) = 0.0;
      }

      // set new lower boundary count
      m_LowerUpperBoundariesCounts.First() = 0.0;
    }

    //! @brief remove the bins after the given index
    //! will add the counts from the removed bins to the boundary counts
    void Histogram::RemoveBinsAfterIndex( const size_t LAST_INDEX)
    {
      // do nothing if last index is beyond the number of bins
      if( LAST_INDEX >= GetNumberOfBins())
      {
        return;
      }

      // add the counts of the to be removed bins to the boundary
      for( size_t i( LAST_INDEX + 1), number_bins( GetNumberOfBins()); i != number_bins; ++i)
      {
        m_LowerUpperBoundariesCounts.Second() += m_Histogram( i);
      }

      // remove the bins
      m_Histogram = linal::Vector< double>( LAST_INDEX + 1, m_Histogram.Begin());

      // set new upper boundary
      m_LowerUpperBoundaries.Second() = m_LowerUpperBoundaries.First() + m_Histogram.GetSize() * m_BinSize;
    }

    //! @brief reset the bins after the given index
    //! counts in bins and at boundary are set to 0.0
    //! @param LAST_INDEX index of bin after which histogram is reset
    void Histogram::ResetBinsAfterIndex( const size_t LAST_INDEX)
    {
      // do nothing if last index is beyond the number of bins
      if( LAST_INDEX >= GetNumberOfBins())
      {
        return;
      }

      // reset the counts of the bins to the boundary
      for( size_t i( LAST_INDEX + 1), number_bins( GetNumberOfBins()); i != number_bins; ++i)
      {
        m_Histogram( i) = 0.0;
      }

      // set new upper boundary count
      m_LowerUpperBoundariesCounts.Second() = 0.0;
    }

    //! @brief normalize the counts in the histogram
    //! each bin is divided by the total number of counts (even the boundary counts)
    void Histogram::Normalize()
    {
      // total sum
      double total( m_Histogram.Sum());
      total += m_LowerUpperBoundariesCounts.First();
      total += m_LowerUpperBoundariesCounts.Second();

      // normalize only if total is larger than 0
      if( total == double( 0))
      {
        return;
      }

      // normalize
      m_LowerUpperBoundariesCounts.First() /= total;
      m_LowerUpperBoundariesCounts.Second() /= total;
      m_Histogram /= total;
    }

    //! @brief calculate the mean
    //! @return weighted mean of x
    double Histogram::CalculateMean() const
    {
      double weighted_mean( 0);
      double total( 0);

      double current_bin_center( m_LowerUpperBoundaries.First() + 0.5 * m_BinSize);
      // iterate over counts
      for
      (
        const double *count( m_Histogram.Begin()), *count_end( m_Histogram.End());
        count != count_end;
        ++count, current_bin_center += m_BinSize
      )
      {
        weighted_mean += current_bin_center * ( *count);
        total += *count;
      }

      // end
      return weighted_mean / total;
    }

    //! @brief calculate the standard deviation
    //! @return standard deviation
    double Histogram::CalculateSD() const
    {
      double weighted_square_norm( 0);
      double weighted_mean( 0);
      double total( 0);

      double current_bin_center( m_LowerUpperBoundaries.First() + 0.5 * m_BinSize);
      // iterate over counts
      for
      (
        const double *count( m_Histogram.Begin()), *count_end( m_Histogram.End());
        count != count_end;
        ++count, current_bin_center += m_BinSize
      )
      {
        weighted_mean += current_bin_center * ( *count);
        total += *count;
        weighted_square_norm += Sqr( current_bin_center) * ( *count);
      }

      // calculate mean
      weighted_mean /= total;
      // calculate square norm
      weighted_square_norm /= total;

      // calculate standard deviation and return
      return Sqrt( Absolute( weighted_square_norm - Sqr( weighted_mean)));
    }

    //! @brief extends the histogram in the lower or upper direction
    //! @param NUM_BINS_LOWER the number of bins in the lower direction that the histogram should be extended
    //! @param LOWER_BIN_VALUES value that will be assigned to all of the newly extended bins in the lower direction
    //! @param FLOOR the desired lowest value for the lower boundary. bins won't be extended lower than this
    //! @param NUM_BINS_UPPER the number of bins in the upper direction that the histogram should be extended
    //! @param UPPER_BIN_VALUES value that will be assigned to all of the newly extended bins in the upper direction
    //! @param CEILING the desired highest value for the upper boundary. bins won't be extended higher than this
    //! @return bool indicating of extension was successful. Might fail if there are boundary counts since extension
    //!         might put the boundary counts into a bin but don't know from boundary counts alone.
    bool Histogram::ExtendBoundaries
    (
      const size_t NUM_BINS_LOWER,
      const double LOWER_BIN_VALUES,
      const double FLOOR,
      const size_t NUM_BINS_UPPER,
      const double UPPER_BIN_VALUES,
      const double CEILING
    )
    {
      // can't work if there are counts in the upper or lower boundaries
      if( m_LowerUpperBoundariesCounts.First() != 0 || m_LowerUpperBoundariesCounts.Second() != 0)
      {
        BCL_MessageStd( "can't extend boundaries of histogram with boundary counts");
        return false;
      }

      // make sure the floor and ceiling are outside of the boundaries
      if( FLOOR > m_LowerUpperBoundaries.First() || CEILING < m_LowerUpperBoundaries.Second())
      {
        BCL_MessageStd( "floor or ceiling within boundaries of histogram");
        return false;
      }

      // check how many bins in lower and upper direction the histogram can be extended
      const size_t lower_bin_room( Absolute( m_LowerUpperBoundaries.First() - FLOOR)    / m_BinSize);
      const size_t upper_bin_room( Absolute( CEILING - m_LowerUpperBoundaries.Second()) / m_BinSize);

      // set the number of bins to extend in lower direction
      const size_t lower_additional_bins( std::min( lower_bin_room, NUM_BINS_LOWER));

      // set the number of bins to extend in upper direction
      const size_t upper_additional_bins( std::min( upper_bin_room, NUM_BINS_UPPER));

      // set the upper and lower boundaries
      m_LowerUpperBoundaries.First() = m_LowerUpperBoundaries.First() - lower_additional_bins * m_BinSize;
      m_LowerUpperBoundaries.Second() = m_LowerUpperBoundaries.Second() + upper_additional_bins * m_BinSize;

      // add the new lower bins
      linal::Vector< double> new_values( lower_additional_bins + m_Histogram.GetSize() + upper_additional_bins);
      new_values.ReplaceElements( 0, linal::Vector< double>( lower_additional_bins, LOWER_BIN_VALUES));
      new_values.ReplaceElements( lower_additional_bins, m_Histogram);
      new_values.ReplaceElements
      (
        lower_additional_bins + m_Histogram.GetSize(),
        linal::Vector< double>( upper_additional_bins, UPPER_BIN_VALUES)
      );

      // set m_Histogram
      m_Histogram = new_values;

      return true;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write into std::ostream x and y values vertically
    //! @param OSTREAM output stream
    //! @param FORMAT_BINNING format for the bin middle coordinate
    //! @param FORMAT_VALUES format for the written values (counts)
    //! @param INDENT the indentation
    //! @return the ostream written to
    std::ostream &Histogram::WriteHorizontally
    (
      std::ostream &OSTREAM,
      const util::Format &FORMAT_BINNING,
      const util::Format &FORMAT_VALUES,
      const size_t INDENT
    ) const
    {
      //write header line
      io::Serialize::InsertIndent( OSTREAM, INDENT);
      OSTREAM << "\t\t" << GetLowerBoundaryString() << '\t';
      for( size_t i( 0); i < m_Histogram.GetSize(); ++i)
      {
        OSTREAM << GetBinString() << '\t';
      }
      OSTREAM << GetUpperBoundaryString() << '\n';
      io::Serialize::InsertIndent( OSTREAM, INDENT);
      OSTREAM << "center\t\t";

      //write x values
      //write lower boundary
      OSTREAM << FORMAT_BINNING( m_LowerUpperBoundaries.First()) << '\t';

      //write middle point of each bin
      for( double i( 0.5); i < m_Histogram.GetSize(); ++i)
      {
        OSTREAM << FORMAT_BINNING( m_LowerUpperBoundaries.First() + i * m_BinSize) << '\t';
      }

      //write upper boundary
      OSTREAM << FORMAT_BINNING( m_LowerUpperBoundaries.Second()) << '\n';
      io::Serialize::InsertIndent( OSTREAM, INDENT);
      OSTREAM << "counts\t\t";

      //write all counts
      OSTREAM << FORMAT_VALUES( m_LowerUpperBoundariesCounts.First()) << '\t';
      for( const double *ptr( m_Histogram.Begin()), *ptr_end( m_Histogram.End()); ptr != ptr_end; ++ptr)
      {
        OSTREAM << FORMAT_VALUES( *ptr) << '\t';
      }
      OSTREAM << FORMAT_VALUES( m_LowerUpperBoundariesCounts.Second()) << '\n';

      return OSTREAM;
    }

    //! write into std::ostream x and y values vertically
    std::ostream &Histogram::WriteVertically( std::ostream &OSTREAM, const util::Format &FORMAT_BINNING, const util::Format &FORMAT_VALUES) const
    {
      //write header line
      OSTREAM << "center\tcounts" << '\n';

      //write lower boundary
      OSTREAM << GetLowerBoundaryString() << "\t\t" << FORMAT_BINNING( m_LowerUpperBoundaries.First()) << '\t' << FORMAT_VALUES( m_LowerUpperBoundariesCounts.First()) << '\n';

      //write middle point of each bin and the according count
      double i( 0.5);
      for( const double *ptr( m_Histogram.Begin()), *ptr_end( m_Histogram.End()); ptr != ptr_end; ++ptr, ++i)
      {
        OSTREAM << GetBinString() << "\t\t"
               << FORMAT_BINNING( m_LowerUpperBoundaries.First() + i * m_BinSize) << '\t'
               << FORMAT_VALUES( *ptr) << '\n';
      }

      //write upper boundary
      OSTREAM << GetUpperBoundaryString() << "\t\t" << FORMAT_BINNING( m_LowerUpperBoundaries.Second()) << '\t' << FORMAT_VALUES( m_LowerUpperBoundariesCounts.Second()) << '\n';

      return OSTREAM;
    }

    //! read horizontally written Histogram from std::istream
    std::istream &Histogram::Read( std::istream &ISTREAM)
    {
      //return
      return ReadHorizontally( ISTREAM);
    }

    //! @brief write the data of the histogram in gnuplot format
    //! @param OSTREAM the stream to write gnuplot script to
    //! @param BINNING the data of the histgram that is going to be printed in gnuplot format
    //! @param COUNTS the counts that correspond to the bins
    //! @return the stream it was written to
    std::ostream &Histogram::WriteGnuPlotHeatMapFormatted
    (
      std::ostream &OSTREAM, const linal::Vector< double> &BINNING, const linal::Vector< double> &COUNTS
    )
    {
      for
      (
        const double *x( BINNING.Begin()), *x_end( BINNING.End()), *y( COUNTS.Begin()), *y_end( COUNTS.End());
        x != x_end && y != y_end;
        ++x, ++y
      )
      {
        OSTREAM << '\n';
        for( size_t i( 0); i < 2; ++i)
        {
          OSTREAM << i << '\t' << *x << '\t' << *y << '\n';
        }
      }

      // indicate end of data for heatmap
      OSTREAM << "e\n";

      // end
      return OSTREAM;
    }

    //! @brief write the data of the histogram in gnuplot format for making a line plot
    //! @param OSTREAM the stream to write gnuplot script to
    //! @param BINNING the data of the histgram that is going to be printed in gnuplot format
    //! @param COUNTS the counts that correspond to the bins
    //! @return the stream it was written to
    std::ostream &Histogram::WriteGnuPlotLinePlotFormatted
    (
      std::ostream &OSTREAM, const linal::Vector< double> &BINNING, const linal::Vector< double> &COUNTS
    )
    {
      // write out the data to file
      for
      (
        const double *x( BINNING.Begin()), *x_end( BINNING.End()), *y( COUNTS.Begin()), *y_end( COUNTS.End());
        x != x_end && y != y_end;
        ++x, ++y
      )
      {
        OSTREAM << *x << '\t' << *y << '\n';
      }

      // indicate end of data for gnuplot
      OSTREAM << "e\n";

      // end
      return OSTREAM;
    }

    //! @brief generate a linear gnuplot
    //! @param OSTREAM the stream to write gnuplot script to
    //! @param TITLE title for the linear map
    //! @return the stream it was written to
    std::ostream &Histogram::WriteLinearGnuplot( std::ostream &OSTREAM, const std::string &TITLE) const
    {
      // write comments
      OSTREAM << "# BCL generated linear plot from histogram\n";

      OSTREAM << "set terminal png transparent enhanced # size 2160,640 \n";
      OSTREAM << "set output \"" << TITLE << ".png\"\n";
      OSTREAM << "set encoding iso\n";
      OSTREAM << "set view map\n";
      OSTREAM << "set title \"" << TITLE << "\"\n";
      OSTREAM << "unset key\n\n";

      OSTREAM << "set xlabel \"x\"\n";
      OSTREAM << "set xrange [" << m_LowerUpperBoundaries.First() << ":" << m_LowerUpperBoundaries.Second() << "]\n";
      OSTREAM << "set autoscale y\n";
      OSTREAM << "plot '-' using 1:2 with lines\n";
      OSTREAM << "# number x values " << m_Histogram.GetSize() << " binning " << m_BinSize << '\n';
      OSTREAM << '\n';

      // write histogram data and return OSTREAM
      return WriteGnuPlotLinePlotFormatted( OSTREAM, GetBinning(), m_Histogram);
    }

    //! @brief read horizontally written Histogram from std::istream
    //! @param ISTREAM inout stream to read from
    //! @return the input stream read from
    std::istream &Histogram::ReadHorizontally( std::istream &ISTREAM)
    {
      // check whether file is valid Histogram file
      std::string identify;
      ISTREAM >> identify;
      BCL_Assert( identify == GetLowerBoundaryString(), "BCL_HISTOGRAM is not written horizontally");

      size_t number_of_bins( 0);

      // count the number of bins
      do
      {
        ISTREAM >> identify;
        ++number_of_bins;
      } while( !ISTREAM.eof() && identify != "center");

      // decrement for the boundary bins
      number_of_bins -= 2;

      //read lower boundary
      ISTREAM >> m_LowerUpperBoundaries.First();
      do
      {
        ISTREAM >> identify;
        if( util::IsNumerical( identify))
        {
          m_LowerUpperBoundaries.Second() = util::ConvertStringToNumericalValue< double>( identify);
        }
      } while( !ISTREAM.eof() && identify != "counts");

      //set histogram to new size
      m_Histogram = linal::Vector< int>( number_of_bins);

      //read count of lower boundary
      ISTREAM >> m_LowerUpperBoundariesCounts.First();

      for( double *ptr( m_Histogram.Begin()), *ptr_end( m_Histogram.End()); ptr != ptr_end; ++ptr)
      {
        ISTREAM >> *ptr;
      }

      // read count of upper boundary
      ISTREAM >> m_LowerUpperBoundariesCounts.Second();

      // calculate bin size
      m_BinSize = ( m_LowerUpperBoundaries.Second() - m_LowerUpperBoundaries.First()) / number_of_bins;

      //return
      return ISTREAM;
    }

  } // namespace math
} // namespace bcl
