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
#include "math/bcl_math_histogram_2d.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_operations.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

  //////////
  // data //
  //////////

    //! @brief string to indicate the center
    //! @return string to indicate the center
    const std::string &Histogram2D::GetCenterString()
    {
      static std::string s_center_string( "center");
      return s_center_string;
    }

    //! @brief string to indicate counts
    //! @return string to indicate counts
    const std::string &Histogram2D::GetCountsString()
    {
      static std::string s_counts( "counts");
      return s_counts;
    }

    //! @brief string to indicate the left boundary in output
    //! @return string to indicate the left boundary in output
    const std::string &Histogram2D::GetLeftBoundaryString()
    {
      static std::string s_left_boundary_string( "...<");
      return s_left_boundary_string;
    }

    //! @brief string to indicate a bin
    //! @return string to indicate a bin
    const std::string &Histogram2D::GetBinString()
    {
      static std::string s_bin_string( "<..>");
      return s_bin_string;
    }

    //! @brief string to indicate the right boundary in output
    //! @return string to indicate the right boundary in output
    const std::string &Histogram2D::GetRightBoundaryString()
    {
      static std::string s_right_boundary_string( ">...");
      return s_right_boundary_string;
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Histogram2D::s_Instance
    (
      GetObjectInstances().AddInstance( new Histogram2D())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor - initializes everything to 0
    Histogram2D::Histogram2D() :
      m_LowerUpperBoundariesX( 0.0, 0.0),
      m_LowerUpperBoundariesY( 0.0, 0.0),
      m_BinSizeXY( 0.0, 0.0),
      m_Histogram( 0, 0)
    {
    }

    //! construct Histogram from X and Y starting values MIN_X_Y, from Pair BINSIZE_X_Y and from NUMBER_OF_BINS_X_Y
    Histogram2D::Histogram2D
    (
      const storage::VectorND< 2, double> &MIN_X_Y,
      const storage::VectorND< 2, double> &BINSIZE_X_Y,
      const storage::VectorND< 2, size_t> &NUMBER_OF_BINS_X_Y
    ) :
      m_LowerUpperBoundariesX( MIN_X_Y.First(), MIN_X_Y.First() + BINSIZE_X_Y.First() * NUMBER_OF_BINS_X_Y.First()),
      m_LowerUpperBoundariesY( MIN_X_Y.Second(), MIN_X_Y.Second() + BINSIZE_X_Y.Second() * NUMBER_OF_BINS_X_Y.Second()),
      m_BinSizeXY( BINSIZE_X_Y),
      m_Histogram( NUMBER_OF_BINS_X_Y.Second() + 2, NUMBER_OF_BINS_X_Y.First() + 2)
    {
    }

    //! copy constructor
    Histogram2D *Histogram2D::Clone() const
    {
      return new Histogram2D( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Histogram2D::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! get number of bins X
    size_t Histogram2D::GetNumberOfBinsX() const
    {
      return m_Histogram.GetNumberCols() - 2;
    }

    //! get number of bins Y
    size_t Histogram2D::GetNumberOfBinsY() const
    {
      return m_Histogram.GetNumberRows() - 2;
    }

    //! get binsize XY
    const storage::VectorND< 2, double> &Histogram2D::GetBinSizeXY() const
    {
      return m_BinSizeXY;
    }

    //! GetBoundaries X
    const storage::VectorND< 2, double> &Histogram2D::GetBoundariesX() const
    {
      return m_LowerUpperBoundariesX;
    }

    //! GetBoundaries Y
    const storage::VectorND< 2, double> &Histogram2D::GetBoundariesY() const
    {
      return m_LowerUpperBoundariesY;
    }

    //! GetBoundaries counts X
    storage::VectorND< 2, linal::Vector< double> > Histogram2D::GetBoundariesCountsX() const
    {
      return storage::VectorND< 2, linal::Vector< double> >( m_Histogram.GetCol( 0), m_Histogram.GetCol( m_Histogram.GetNumberCols() - 1));
    }

    //! GetBoundaries counts Y
    storage::VectorND< 2, linal::Vector< double> > Histogram2D::GetBoundariesCountsY() const
    {
      return storage::VectorND< 2, linal::Vector< double> >( m_Histogram.GetRow( 0), m_Histogram.GetRow( m_Histogram.GetNumberRows() - 1));
    }

    //! GetHistogram
    linal::Matrix< double> Histogram2D::GetHistogram() const
    {
      return m_Histogram.CreateSubMatrix( m_Histogram.GetNumberRows() - 2, m_Histogram.GetNumberCols() - 2, 1, 1);
    }

    //! get sum of all counts
    double Histogram2D::GetSumOfAllCounts() const
    {
      return m_Histogram.Sum();
    }

    //! GetBinning
    storage::VectorND< 2, linal::Vector< double> > Histogram2D::GetBinningXY() const
    {
      storage::VectorND< 2, linal::Vector< double> > binning( linal::Vector< double>( m_Histogram.GetNumberCols() - 2, double( 0.0)), linal::Vector< double>( m_Histogram.GetNumberRows() - 2, double( 0.0)));

      //binning in X direction
      double i( 0.5);
      for( double *ptr( binning.First().Begin()), *ptr_end( binning.First().End()); ptr != ptr_end; ++ptr, ++i)
      {
        *ptr = m_LowerUpperBoundariesX.First() + i * m_BinSizeXY.First();
      }

      //binning in Y direction
      i = 0.5;
      for( double *ptr( binning.Second().Begin()), *ptr_end( binning.Second().End()); ptr != ptr_end; ++ptr, ++i)
      {
        *ptr = m_LowerUpperBoundariesY.First() + i * m_BinSizeXY.Second();
      }

      //return
      return binning;
    }

  ////////////////
  // operations //
  ////////////////

    //! calculate the Histogram2D from a VALUES_VECTOR
    Histogram2D &Histogram2D::CalculateHistogram
    (
      const storage::Vector< storage::VectorND< 2, double> > &VALUES_VECTOR
    )
    {
      //reset all counts to zero
      Reset();

      //count
      for
      (
        std::vector< storage::VectorND< 2, double> >::const_iterator
          itr( VALUES_VECTOR.Begin()), itr_end( VALUES_VECTOR.End());
        itr != itr_end; ++itr
      )
      {
        PushBack( *itr);
      }

      return *this;
    }

    //! reset all counts
    void Histogram2D::Reset()
    {
      m_Histogram.SetZero();
    }

    //! pushback a pair of values to the right position in the histogram
    void Histogram2D::PushBack
    (
      const storage::VectorND< 2, double> &PAIR_OF_VALUES,
      const double &WEIGHT
    )
    {
      if( !util::IsDefined( PAIR_OF_VALUES.First()) || !util::IsDefined( PAIR_OF_VALUES.Second())) return;

      size_t bin_x( 0), bin_y( 0);

      //determine field in x-direction
      if( PAIR_OF_VALUES.First() < m_LowerUpperBoundariesX.First())
      {
        bin_x = 0;
      }
      else if( PAIR_OF_VALUES.First() >= m_LowerUpperBoundariesX.Second())
      {
        bin_x = m_Histogram.GetNumberCols() - 1;
      }
      else
      {
        bin_x = size_t( ( PAIR_OF_VALUES.First() - m_LowerUpperBoundariesX.First()) / m_BinSizeXY.First() + 1);
      }

      //determine field in y-direction
      if( PAIR_OF_VALUES.Second() < m_LowerUpperBoundariesY.First())
      {
        bin_y = 0;
      }
      else if( PAIR_OF_VALUES.Second() >= m_LowerUpperBoundariesY.Second())
      {
        bin_y = m_Histogram.GetNumberRows() - 1;
      }
      else
      {
        bin_y = size_t( ( PAIR_OF_VALUES.Second() - m_LowerUpperBoundariesY.First()) / m_BinSizeXY.Second() + 1);
      }

      m_Histogram( bin_y, bin_x) += WEIGHT;
    }

    //! @brief combine this with a given histogram2d by adding up all counts
    //! all parameters have to be identical for this operation to work
    //! @param HISTOGRAM_2D histogram to be added to this one
    //! @return true if it was successful
    bool Histogram2D::Combine( const Histogram2D &HISTOGRAM_2D)
    {
      // allow combining, if this is empty
      if( IsEmpty())
      {
        *this = HISTOGRAM_2D;
        return true;
      }

      // check that the parameters agree
      if
      (
           m_LowerUpperBoundariesX.First()  != HISTOGRAM_2D.m_LowerUpperBoundariesX.First()
        || m_LowerUpperBoundariesX.Second() != HISTOGRAM_2D.m_LowerUpperBoundariesX.Second()
        || m_LowerUpperBoundariesY.First()  != HISTOGRAM_2D.m_LowerUpperBoundariesY.First()
        || m_LowerUpperBoundariesY.Second() != HISTOGRAM_2D.m_LowerUpperBoundariesY.Second()
        || m_BinSizeXY.First() != HISTOGRAM_2D.m_BinSizeXY.First()
        || m_BinSizeXY.Second() != HISTOGRAM_2D.m_BinSizeXY.Second()
      )
      {
        return false;
      }

      // add the histogram
      m_Histogram += HISTOGRAM_2D.m_Histogram;

      // end
      return true;
    }

    //! checks if there is any count
    bool Histogram2D::IsEmpty() const
    {
      return GetSumOfAllCounts() == 0;
    }

    //! @brief normalize the counts in the histogram
    //! each bin is divided by the total number of counts (even the boundary counts)
    void Histogram2D::Normalize()
    {
      // create conts double "total_counts" and initialize with result of "GetSumOfAllCounts"
      const double total_counts( GetSumOfAllCounts());

      // if total_counts is 0 then return since no normalization needed
      if( !total_counts)
      {
        return;
      }

      // normalize
      m_Histogram /= total_counts;
    }

    //! @brief normalize the row counts in the histogram
    //! @detail each cell is divided by the maximum number of counts in the row
    void Histogram2D::NormalizeRows()
    {
      for( size_t row_index( 0), row_count( m_Histogram.GetNumberRows()); row_index < row_count; ++row_index)
      {
        // copy row from histogram and calculate the count
        linal::Vector< double> row( m_Histogram.GetRow( row_index));
        const double max_count( row.Max());

        // iterate over each element dividing by the total count of this column
        for( linal::Vector< double>::iterator itr( row.Begin()), itr_end( row.End()); itr != itr_end; ++itr)
        {
          *itr /= max_count;
        }

        // replace the row in the histogram with the column-normalized one
        m_Histogram.ReplaceRow( row_index, row);
      }
    }

    //! @brief normalize the column counts in the histogram
    //! each cell is divided by the total number of counts in the column
    void Histogram2D::NormalizeY()
    {
      for( size_t col( 0), col_count( m_Histogram.GetNumberCols()); col < col_count; ++col)
      {
        // copy column from histogram and calculate the count
        linal::Vector< double> col_vec( m_Histogram.GetCol( col));
        const double total_count_in_current_col( col_vec.Sum());

        // if the count if larger than 1, normalize this column
        if( total_count_in_current_col > 1)
        {
          // iterate over each element dividing by the total count of this column
          for( linal::Vector< double>::iterator itr( col_vec.Begin()), itr_end( col_vec.End()); itr != itr_end; ++itr)
          {
            *itr /= total_count_in_current_col;
          }

          // replace the row in the histogram with the column-normalized one
          m_Histogram.ReplaceCol( col, col_vec);
        }
      }
    }

    //! @brief Normalize by another histogram 2d, which represents background counts
    //! @param BACKGROUND background counts histogram
    void Histogram2D::NormalizeByBackground( const Histogram2D &HIST)
    {
      BCL_Assert
      (
        m_LowerUpperBoundariesX == HIST.m_LowerUpperBoundariesX
        && m_LowerUpperBoundariesY == HIST.m_LowerUpperBoundariesY
        && m_BinSizeXY == HIST.m_BinSizeXY,
        "Normalization could not occur because histograms were of different size or binning"
      );
      for( size_t i( 0), nr( HIST.m_Histogram.GetNumberRows()); i < nr; ++i)
      {
        auto this_row( m_Histogram.GetRow( i));
        auto bg_row( HIST.m_Histogram.GetRow( i));
        auto row_threshold( std::max( this_row.Sum() / 100.0 / double( bg_row.GetSize()), 10.5));
        auto itr_bg( bg_row.Begin());
        for
        (
          auto itr_this( this_row.Begin()), itr_this_end( this_row.End());
          itr_this != itr_this_end;
          ++itr_this, ++itr_bg
        )
        {
          // normalize each non-zero bin
          if( *itr_this > row_threshold || *itr_this < -row_threshold)
          {
            *itr_this /= *itr_bg;
          }
          else
          {
            *itr_this = 0.0;
          }
        }
      }
    }

    //! @brief calculate the logarithm of the column counts + 1 in the histogram
    void Histogram2D::Log()
    {
      for( size_t col( 0), col_count( m_Histogram.GetNumberCols()); col < col_count; ++col)
      {
        linal::Vector< double> col_vec( m_Histogram.GetCol( col));
        for( linal::Vector< double>::iterator itr( col_vec.Begin()), itr_end( col_vec.End()); itr != itr_end; ++itr)
        {
          *itr = std::log10( *itr + 1.0);
        }
        m_Histogram.ReplaceCol( col, col_vec);
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! write into std::ostream x and y values horizontally
    std::ostream &Histogram2D::Write
    (
      std::ostream &OSTREAM,
      const util::Format &FORMAT_BINNING,
      const util::Format &FORMAT_VALUES
    ) const
    {
      //write header line
      OSTREAM << GetCenterString() << '\t' << GetLeftBoundaryString() << '\t';
      for( size_t i( 0); i < m_Histogram.GetNumberCols() - 2; ++i)
      {
        OSTREAM << GetBinString() << '\t';
      }
      OSTREAM << GetRightBoundaryString() << '\n' << '\t' << GetCenterString() << '\t' << GetCountsString() << '\t';

      //write x values
      //write lower boundary X
      OSTREAM << FORMAT_BINNING( m_LowerUpperBoundariesX.First()) << '\t';

      //write middle point of each bin in X direction
      double i( 0.5);
      for( i = 0.5; i < m_Histogram.GetNumberCols() - 2; ++i)
      {
        OSTREAM << FORMAT_BINNING( m_LowerUpperBoundariesX.First() + i * m_BinSizeXY.First()) << '\t';
      }

      //write upper boundary X
      OSTREAM << FORMAT_BINNING( m_LowerUpperBoundariesX.Second()) << '\n';

      //write lower boundary Y
      OSTREAM << '\t' << GetLeftBoundaryString() << '\t' << FORMAT_BINNING( m_LowerUpperBoundariesY.First()) << '\t';

      //write lower boundary counts
      for( const double *ptr( m_Histogram[ 0]); ptr != m_Histogram[ 0] + m_Histogram.GetNumberCols(); ++ptr)
      {
        OSTREAM << FORMAT_VALUES( *ptr) << '\t';
      }

      OSTREAM << '\n';

      //write all counts
      i = 0.5;
      for( size_t row( 1); row < m_Histogram.GetNumberRows() - 1; ++row, ++i)
      {
        OSTREAM << '\t' << GetBinString() << '\t' << FORMAT_BINNING( m_LowerUpperBoundariesY.First() + i * m_BinSizeXY.Second()) << '\t';
        for( const double *ptr( m_Histogram[ row]); ptr != m_Histogram[ row] + m_Histogram.GetNumberCols(); ++ptr)
        {
          OSTREAM << FORMAT_VALUES( *ptr) << '\t';
        }
        OSTREAM << '\n';
      }

      //write upper boundary Y
      OSTREAM << '\t' << GetRightBoundaryString() << '\t' << FORMAT_BINNING( m_LowerUpperBoundariesY.Second()) << '\t';

      //write upper boundary counts
      for
      (
        const double *ptr( m_Histogram[ m_Histogram.GetNumberRows() - 1]);
        ptr != m_Histogram[ m_Histogram.GetNumberRows() - 1] + m_Histogram.GetNumberCols(); ++ptr
      )
      {
        OSTREAM << FORMAT_VALUES( *ptr) << '\t';
      }

      OSTREAM << '\n';

      return OSTREAM;
    }

    //! read Histogram2D from std::istream
    std::istream &Histogram2D::Read( std::istream &ISTREAM)
    {
      // check whether file is valid Histogram2D file
      std::string identify;
      ISTREAM >> identify;
      BCL_Assert( identify == GetCenterString(), "unexpected string. Should be \"" + GetCenterString() + "\" instead of " + identify);
      ISTREAM >> identify;
      BCL_Assert( identify == GetLeftBoundaryString(), "unexpected string. Should be \"" + GetLeftBoundaryString() + "\" instead of " + identify);

      size_t number_cols( 1), number_rows( 1);
      //count number of cols
      while( !ISTREAM.eof() && identify != GetRightBoundaryString())
      {
        ISTREAM >> identify;
        number_cols++;
        if( identify != GetRightBoundaryString())
        {
          BCL_Assert( identify == GetBinString(), "unexpected string. Should be \"" + GetBinString() + "\"");
        }
      }

      ISTREAM >> identify;
      BCL_Assert( identify == GetCenterString(), "unexpected string. Should be \"" + GetCenterString() + "\"");
      ISTREAM >> identify;
      BCL_Assert( identify == GetCountsString(), "unexpected string. Should be \"" + GetCountsString() + "\"");

      //read lower boundary X
      ISTREAM >> m_LowerUpperBoundariesX.First();

      //read upper boundary X
      while( !ISTREAM.eof())
      {
        ISTREAM >> identify;
        if( identify != GetLeftBoundaryString())
        {
          m_LowerUpperBoundariesX.Second() = util::ConvertStringToNumericalValue< double>( identify);
        }
        else
        {
           break;
        }
      }

      //read lower boundary Y
      ISTREAM >> m_LowerUpperBoundariesY.First();

      //read values
      storage::Vector< double> values;

      std::string tmp;
      while( !ISTREAM.eof())
      {
        do
        {
          ISTREAM >> tmp;
          if( util::IsNumerical( tmp))
          {
            values.PushBack( util::ConvertStringToNumericalValue< double>( tmp));
          }

        } while( !ISTREAM.eof() && tmp != GetBinString() && tmp != GetRightBoundaryString());
        number_rows++;
        if( tmp != GetRightBoundaryString())
        {
          ISTREAM >> tmp;
          continue;
        }
        ISTREAM >> m_LowerUpperBoundariesY.Second();
        break;
      }

      for( size_t i = 0; i < number_cols && !ISTREAM.eof(); ++i)
      {
        ISTREAM >> tmp;
        values.PushBack( util::ConvertStringToNumericalValue< double>( tmp));
      }

      m_Histogram = linal::Matrix< double>( number_rows, number_cols, values);

      m_BinSizeXY.First() = ( m_LowerUpperBoundariesX.Second() - m_LowerUpperBoundariesX.First()) / double( number_cols - 2);
      m_BinSizeXY.Second() = ( m_LowerUpperBoundariesY.Second() - m_LowerUpperBoundariesY.First()) / double( number_rows - 2);

      //return
      return ISTREAM;
    }

  } // namespace math
} // namespace bcl
