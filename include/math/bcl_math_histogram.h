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

#ifndef BCL_MATH_HISTOGRAM_H_
#define BCL_MATH_HISTOGRAM_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_math_function_interface_serializable.h"
#include "linal/bcl_linal_vector.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Histogram
    //! @brief is a class to collect statistics.
    //! @details it has a certain number of equidistant bins and boundaries, and collects counts in those bins,
    //! and counts for values to the left and right boundary of the interval.
    //!
    //! @see @link example_math_histogram.cpp @endlink
    //! @author woetzen
    //! @date 11.05.2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Histogram :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      storage::VectorND< 2, double> m_LowerUpperBoundaries;       //!< min and max value in the histogram
      double                        m_BinSize;                    //!< GetSize of one bin
      linal::Vector< double>        m_Histogram;                  //!< calculate histogram
      storage::VectorND< 2, double> m_LowerUpperBoundariesCounts; //!< counts that are at the border of the histogram

      //! @brief string to indicate the lower boundary in output
      //! @return string to indicate the lower boundary in output
      static const std::string &GetLowerBoundaryString();

      //! @brief string to indicate a bin
      //! @return string to indicate a bin
      static const std::string &GetBinString();

      //! @brief string to indicate the upper boundary in output
      //! @return string to indicate the upper boundary in output
      static const std::string &GetUpperBoundaryString();

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief flag for setting the minimum value of the histogram
      //! @return flag for setting the minimum value of the histogram
      static const util::ShPtr< command::FlagInterface> GetFlagMin();

      //! @brief flag for setting the bin size of the histogram
      //! @return flag for setting the bin size of the histogram
      static const util::ShPtr< command::FlagInterface> GetFlagBinSize();

      //! @brief flag for setting the number of bins of the histogram
      //! @return flag for setting the number of bins of the histogram
      static const util::ShPtr< command::FlagInterface> GetFlagNumberOfBins();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor - initializes all to 0
      Histogram();

      //! @brief construct histogram from starting value MIN, BIN_SIZE and NUMBER_OF_BINS
      //! @param MIN the minimal value representing the left boundary
      //! @param BIN_SIZE the size (width) of one bin
      //! @param NUMBER_OF_BINS the number of bin in the histogram
      Histogram
      (
        const double MIN,
        const double BIN_SIZE,
        const size_t NUMBER_OF_BINS
      );

      //! @brief construct a histogram given a map of data with bins being the key and counts being the value
      //! @param DATA histogram data with key being the bin and the value being the counts for that bin
      Histogram( const storage::Map< double, double> &DATA);

      //! copy constructor
      Histogram *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief access number of bins
      //! @return the number of bins (without the bins below or above the boundaries)
      size_t GetNumberOfBins() const;

      //! @brief access binsize
      //! @return the width of a bin
      double const &GetBinSize() const;

      //! @brief GetBoundaries
      //! @return first is the left (lower) boundary threshold, second the right (upper) boundary threshold
      const storage::VectorND< 2, double> &GetBoundaries() const;

      //! @brief GetBoundaries counts
      //! @return first is the left hand (lower) boundary count, second is the right hand (upper) boundary count
      const storage::VectorND< 2, double> &GetBoundariesCounts() const;

      //! @brief GetHistogram
      //! @return linal::Vector< double> that contains the counts for the bins within the boundaries
      const linal::Vector< double> &GetHistogram() const;

      //! @brief binning
      //! @return Vector with coordinates representing the middle of all bins
      linal::Vector< double> GetBinning() const;

      //! @brief get sum of all counts
      //! @return the sum of the boundary counts and plus the sum of the counts in all bins
      double GetSumOfAllCounts() const;

      //! This functions returns the sum of counts between [MIN, MAX)
      //! @param MIN left bound of interval
      //! @param MAX rigth bound of interval
      //! @return the sum of counts in the given interval
      double GetCountsInBetween( const double MIN, const double MAX) const;

      //! @brief set the count for a bin
      //! @param INDEX index of the bin
      //! @param COUNT the new count
      void SetCount( const size_t INDEX, const double COUNT);

    ////////////////
    // operations //
    ////////////////

      //! @brief calculate the historgam from a vector of values
      //! @param VALUES_VECTOR vector of values
      //! @return a reference to this histograms
      Histogram &CalculateHistogram( const storage::Vector< double> &VALUES_VECTOR);

      //! @brief pushback a value to its position in the histogram
      //! @param VALUE the value to be considered to derive the bin
      //! @param WEIGHT the increment for the bin that belongs to VALUE - 1 by default
      void PushBack( const double VALUE, const double WEIGHT = double( 1));

      //! @brief combine this with a given histogram by adding up all counts
      //! all parameters have to be identical for this operation to work
      //! @param HISTOGRAM histogram to be added to this one
      //! @return true if it was successful
      bool Combine( const Histogram &HISTOGRAM);

      //! @brief add a pseudo count
      //! add a pseudo count to all bins and to the boundaries of the histogram
      //! @param VALUE the value to be added to each bin and to the boundaries
      void AddPseudoCount( const double VALUE);

      //! @brief reset all counts
      void Reset();

      //! @brief checks if there is any count
      //! @return true for empty histogram (sum of all counts is 0)
      bool IsEmpty() const;

      //! @brief is range defined
      //! @return returns true, if the statistic is defined, false if it is not (e.g. ConsiderValue was never called)
      bool IsDefined() const
      {
        return !IsEmpty();
      }

      //! @brief determine index of last bin with count == 0 from the back
      //! @param COUNT_THRESHOLD if count > COUNT_THRESHOLD, it is considered to contain information - default = 0.0
      //! @return index of last bin with counts in it - number of bins, if right boundary has counts, 0 if nothing has count
      size_t GetIndexOfLastInformationContainingBin( const double COUNT_THRESHOLD = 0.0) const;

      //! @brief determine index of first bin with count != 0 from the front
      //! @param COUNT_THRESHOLD if count > COUNT_THRESHOLD, it is considered to contain information - default = 0.0
      //! @return index of first bin with counts in it - number of bins, 0 if all have counts, number of bins, if all have no counts
      size_t GetIndexOfFirstInformationContainingBin( const double COUNT_THRESHOLD = 0.0) const;

      //! @brief remove the bins before the given index
      //! counts are added to left boundary
      //! @param FIRST_INDEX index of bin before which histogram is cut
      void RemoveBinsBeforeIndex( const size_t FIRST_INDEX);

      //! @brief reset the bins before the given index
      //! counts in bins and at boundary are set to 0.0
      //! @param FIRST_INDEX index of bin before which histogram is reset
      void ResetBinsBeforeIndex( const size_t FIRST_INDEX);

      //! @brief remove the bins after the given index
      //! counts are added to right boundary
      //! @param LAST_INDEX index of bin after which histogram is cut
      void RemoveBinsAfterIndex( const size_t LAST_INDEX);

      //! @brief reset the bins after the given index
      //! counts in bins and at boundary are set to 0.0
      //! @param LAST_INDEX index of bin after which histogram is reset
      void ResetBinsAfterIndex( const size_t LAST_INDEX);

      //! @brief normalize the counts in the histogram
      //! each bin is divided by the total number of counts (even the boundary counts)
      void Normalize();

      //! @brief calculate the mean
      //! @return weighted mean of x
      double CalculateMean() const;

      //! @brief calculate the standard deviation
      //! @return standard deviation
      double CalculateSD() const;

      //! @brief extends the histogram in the lower or upper direction
      //! @param NUM_BINS_LOWER the number of bins in the lower direction that the histogram should be extended
      //! @param LOWER_BIN_VALUES value that will be assigned to all of the newly extended bins in the lower direction
      //! @param FLOOR the desired lowest value for the lower boundary. bins won't be extended lower than this
      //! @param NUM_BINS_UPPER the number of bins in the upper direction that the histogram should be extended
      //! @param UPPER_BIN_VALUES value that will be assigned to all of the newly extended bins in the upper direction
      //! @param CEILING the desired highest value for the upper boundary. bins won't be extended higher than this
      //! @return bool indicating of extension was successful. Might fail if there are boundary counts since extension
      //!         might put the boundary counts into a bin but don't know from boundary counts alone.
      bool ExtendBoundaries
      (
        const size_t NUM_BINS_LOWER,
        const double LOWER_BIN_VALUES,
        const double FLOOR,
        const size_t NUM_BINS_UPPER,
        const double UPPER_BIN_VALUES,
        const double CEILING
      );

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write into std::ostream x and y values horizontally
      //! @param OSTREAM output stream
      //! @param FORMAT_BINNING format for the bin middle coordinate
      //! @param FORMAT_VALUES format for the written values (counts)
      //! @param INDENT the indentation
      //! @return the ostream written to
      std::ostream &WriteHorizontally
      (
        std::ostream &OSTREAM,
        const util::Format &FORMAT_BINNING = util::Format().W( 8).FFP( 3),
        const util::Format &FORMAT_VALUES = util::Format().W( 8).FFP( 3),
        const size_t INDENT = 0
      ) const;

      //! @brief write into std::ostream x and y values vertically
      //! @param OSTREAM output stream
      //! @param FORMAT_BINNING format for the bin middle coordinate
      //! @param FORMAT_VALUES format for the written values (counts)
      //! @return the ostream written to
      std::ostream &WriteVertically
      (
        std::ostream &OSTREAM,
        const util::Format &FORMAT_BINNING = util::Format().W( 8).FFP( 3),
        const util::Format &FORMAT_VALUES = util::Format().W( 8).FFP( 3)
      ) const;

      //! @brief write the data of the histogram in gnuplot format
      //! @param OSTREAM the stream to write gnuplot script to
      //! @param BINNING the data of the histgram that is going to be printed in gnuplot format
      //! @param COUNTS the counts that correspond to the bins
      //! @return the stream it was written to
      static std::ostream &WriteGnuPlotHeatMapFormatted
      (
        std::ostream &OSTREAM, const linal::Vector< double> &BINNING, const linal::Vector< double> &COUNTS
      );

      //! @brief write the data of the histogram in gnuplot format for making a line plot
      //! @param OSTREAM the stream to write gnuplot script to
      //! @param BINNING the data of the histgram that is going to be printed in gnuplot format
      //! @param COUNTS the counts that correspond to the bins
      //! @return the stream it was written to
      static std::ostream &WriteGnuPlotLinePlotFormatted
      (
        std::ostream &OSTREAM, const linal::Vector< double> &BINNING, const linal::Vector< double> &COUNTS
      );

      //! @brief generate a linear gnuplot
      //! @param OSTREAM the stream to write gnuplot script to
      //! @param TITLE title for the linear map
      //! @return the stream it was written to
      std::ostream &WriteLinearGnuplot( std::ostream &OSTREAM, const std::string &TITLE) const;

      //! @brief read horizontally written Histogram from std::istream
      //! @param ISTREAM inout stream to read from
      //! @return the input stream read from
      std::istream &ReadHorizontally( std::istream &ISTREAM);

    protected:

      //! write into std::ostream ( horizontally)
      //! @param OSTREAM output stream
      //! @param INDENT the indentation
      //! @return the ostream written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return WriteHorizontally( OSTREAM, util::Format().W( 8).FFP( 3), util::Format().W( 8).FFP( 3), INDENT);
      }

      //! @brief read horizontally written Histogram from std::istream
      //! @param ISTREAM inout stream to read from
      //! @return the input stream read from
      std::istream &Read( std::istream &ISTREAM);

    }; // class Histogram

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_HISTOGRAM_H_
