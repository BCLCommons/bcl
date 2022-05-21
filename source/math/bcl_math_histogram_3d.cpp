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
#include "math/bcl_math_histogram_3d.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_operations.h"
#include "linal/bcl_linal_vector_operations.h"
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
    const std::string &Histogram3D::GetCenterString()
    {
      static std::string s_center_string( "center");
      return s_center_string;
    }

    //! @brief string to indicate counts
    //! @return string to indicate counts
    const std::string &Histogram3D::GetCountsString()
    {
      static std::string s_counts( "counts");
      return s_counts;
    }

    //! @brief string to indicate the left boundary in output
    //! @return string to indicate the left boundary in output
    const std::string &Histogram3D::GetLeftBoundaryString()
    {
      static std::string s_left_boundary_string( "...<");
      return s_left_boundary_string;
    }

    //! @brief string to indicate a bin
    //! @return string to indicate a bin
    const std::string &Histogram3D::GetBinString()
    {
      static std::string s_bin_string( "<..>");
      return s_bin_string;
    }

    //! @brief string to indicate the right boundary in output
    //! @return string to indicate the right boundary in output
    const std::string &Histogram3D::GetRightBoundaryString()
    {
      static std::string s_right_boundary_string( ">...");
      return s_right_boundary_string;
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Histogram3D::s_Instance
    (
      GetObjectInstances().AddInstance( new Histogram3D())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor - initializes everything to 0
    Histogram3D::Histogram3D() :
      m_LowerUpperBoundariesX( 0.0, 0.0),
      m_LowerUpperBoundariesY( 0.0, 0.0),
      m_LowerUpperBoundariesZ( 0.0, 0.0),
      m_BinSizeXYZ( 0.0, 0.0, 0.0),
      m_Histogram( 0, 0, 0, 0.0)
    {
    }

    //! construct Histogram from X and Y starting values MIN_X_Y, from Pair BINSIZE_X_Y and from NUMBER_OF_BINS_X_Y
    Histogram3D::Histogram3D
    (
      const storage::VectorND< 3, double> &MIN_X_Y_Z,
      const storage::VectorND< 3, double> &BINSIZE_X_Y_Z,
      const storage::VectorND< 3, size_t> &NUMBER_OF_BINS_X_Y_Z,
      const double &INITIAL_VAL
    ) :
      m_LowerUpperBoundariesX( MIN_X_Y_Z.First(), MIN_X_Y_Z.First() + BINSIZE_X_Y_Z.First() * NUMBER_OF_BINS_X_Y_Z.First()),
      m_LowerUpperBoundariesY( MIN_X_Y_Z.Second(), MIN_X_Y_Z.Second() + BINSIZE_X_Y_Z.Second() * NUMBER_OF_BINS_X_Y_Z.Second()),
      m_LowerUpperBoundariesZ( MIN_X_Y_Z.Third(), MIN_X_Y_Z.Third() + BINSIZE_X_Y_Z.Third() * NUMBER_OF_BINS_X_Y_Z.Third()),
      m_BinSizeXYZ( BINSIZE_X_Y_Z),
      m_Histogram( NUMBER_OF_BINS_X_Y_Z.First(), NUMBER_OF_BINS_X_Y_Z.Second(), NUMBER_OF_BINS_X_Y_Z.Third(), INITIAL_VAL)
    {
      SetupBinning();
    }

    //! copy constructor
    Histogram3D *Histogram3D::Clone() const
    {
      return new Histogram3D( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Histogram3D::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! get number of bins X
    size_t Histogram3D::GetNumberOfBinsX() const
    {
      return m_Histogram.NumberLayers();
    }

    //! get number of bins Y
    size_t Histogram3D::GetNumberOfBinsY() const
    {
      return m_Histogram.GetNumberRows();
    }

    //! get number of bins Y
    size_t Histogram3D::GetNumberOfBinsZ() const
    {
      return m_Histogram.GetNumberCols();
    }

    //! get binsize XY
    const storage::VectorND< 3, double> &Histogram3D::GetBinSizeXYZ() const
    {
      return m_BinSizeXYZ;
    }

    //! GetBoundaries X
    const storage::VectorND< 2, double> &Histogram3D::GetBoundariesX() const
    {
      return m_LowerUpperBoundariesX;
    }

    //! GetBoundaries Y
    const storage::VectorND< 2, double> &Histogram3D::GetBoundariesY() const
    {
      return m_LowerUpperBoundariesY;
    }

    //! GetBoundaries Y
    const storage::VectorND< 2, double> &Histogram3D::GetBoundariesZ() const
    {
      return m_LowerUpperBoundariesZ;
    }

    //! GetHistogram
    const Tensor< double> &Histogram3D::GetHistogram() const
    {
      return m_Histogram;
    }

    //! GetHistogram
    Tensor< double> &Histogram3D::GetChangeableHistogram()
    {
      return m_Histogram;
    }

    //! get sum of all counts
    double Histogram3D::GetSumOfAllCounts() const
    {
      return m_Histogram.GetValues().Sum();
    }

    //! GetBinning
    storage::VectorND< 3, linal::Vector< double> > Histogram3D::GetBinningXYZ() const
    {
      return storage::VectorND< 3, linal::Vector< double> >( m_BinningX, m_BinningY, m_BinningZ);
    }

  ////////////////
  // operations //
  ////////////////

    //! reset all counts
    void Histogram3D::Reset()
    {
      m_Histogram.GetValues() = 0.0;
    }

    //! pushback a pair of values to the right position in the histogram
    void Histogram3D::PushBack
    (
      const double &X,
      const double &Y,
      const double &Z,
      const double &WEIGHT
    )
    {
      if( !WEIGHT || !util::IsDefined( X) || !util::IsDefined( Y) || !util::IsDefined( Z))
      {
        return;
      }

      size_t bin_x( 0), bin_y( 0), bin_z( 0);

      //determine field in x-direction
      if( X < m_LowerUpperBoundariesX.First())
      {
        bin_x = 0;
      }
      else if( X >= m_LowerUpperBoundariesX.Second())
      {
        bin_x = m_Histogram.NumberLayers() - 1;
      }
      else
      {
        bin_x = size_t( ( X - m_LowerUpperBoundariesX.First()) / m_BinSizeXYZ.First());
      }

      //determine field in y-direction
      if( Y < m_LowerUpperBoundariesY.First())
      {
        bin_y = 0;
      }
      else if( Y >= m_LowerUpperBoundariesY.Second())
      {
        bin_y = m_Histogram.GetNumberRows() - 1;
      }
      else
      {
        bin_y = size_t( ( Y - m_LowerUpperBoundariesY.First()) / m_BinSizeXYZ.Second());
      }

      //determine field in z-direction
      if( Z < m_LowerUpperBoundariesZ.First())
      {
        bin_z = 0;
      }
      else if( Z >= m_LowerUpperBoundariesZ.Second())
      {
        bin_z = m_Histogram.GetNumberCols() - 1;
      }
      else
      {
        bin_z = size_t( ( Z - m_LowerUpperBoundariesZ.First()) / m_BinSizeXYZ.Third());
      }

      m_Histogram( bin_x, bin_y, bin_z) += WEIGHT;
    }

    //! Use tri-linear interpolation to obtain a value at an arbitrary point
    //! @param X, Y, Z the coordinates of interest
    double Histogram3D::Interpolate( const double &X, const double &Y, const double &Z) const
    {
      if( !util::IsDefined( X) || !util::IsDefined( Y) || !util::IsDefined( Z))
      {
        return util::GetUndefined< double>();
      }

      size_t bin_x_lo( 0), bin_x_hi( 0), bin_y_lo( 0), bin_y_hi( 0), bin_z_lo( 0), bin_z_hi( 0);

      //determine field in x-direction
      if( X < m_BinningX( 0))
      {
        bin_x_lo = 0;
        bin_x_hi = 1;
      }
      else if( X >= m_BinningX.Last())
      {
        bin_x_lo = m_Histogram.NumberLayers() - 2;
        bin_x_hi = m_Histogram.NumberLayers() - 1;
      }
      else
      {
        bin_x_lo = size_t( ( X - m_BinningX( 0)) / m_BinSizeXYZ.First());
        bin_x_hi = std::min( bin_x_lo + 1, m_Histogram.NumberLayers() - 1);
      }

      //determine field in y-direction
      if( Y < m_BinningY( 0))
      {
        bin_y_lo = 0;
        bin_y_hi = 1;
      }
      else if( Y >= m_BinningY.Last())
      {
        bin_y_lo = m_Histogram.GetNumberRows() - 2;
        bin_y_hi = m_Histogram.GetNumberRows() - 1;
      }
      else
      {
        bin_y_lo = size_t( ( Y - m_BinningY( 0)) / m_BinSizeXYZ.Second());
        bin_y_hi = std::min( bin_y_lo + 1, m_Histogram.GetNumberRows() - 1);
      }

      //determine field in z-direction
      if( Z < m_BinningZ( 0))
      {
        bin_z_lo = 0;
        bin_z_hi = 1;
      }
      else if( Z >= m_BinningZ.Last())
      {
        bin_z_lo = m_Histogram.GetNumberCols() - 2;
        bin_z_hi = m_Histogram.GetNumberCols() - 1;
      }
      else
      {
        bin_z_lo = size_t( ( Z - m_BinningZ( 0)) / m_BinSizeXYZ.Third());
        bin_z_hi = std::min( bin_z_lo + 1, m_Histogram.GetNumberCols() - 1);
      }

      const double x_to_lo( ( X - m_BinningX( bin_x_lo)));
      const double x_to_hi( ( m_BinningX( bin_x_hi) - X));
      const double y_to_lo( ( Y - m_BinningY( bin_y_lo)));
      const double y_to_hi( ( m_BinningY( bin_y_hi) - Y));
      const double z_to_lo( ( Z - m_BinningZ( bin_z_lo)));
      const double z_to_hi( ( m_BinningZ( bin_z_hi) - Z));

      const double c00
      (
        m_Histogram( bin_x_lo, bin_y_lo, bin_z_lo) * x_to_hi
        + m_Histogram( bin_x_hi, bin_y_lo, bin_z_lo) * x_to_lo
      );
      const double c01
      (
        m_Histogram( bin_x_lo, bin_y_lo, bin_z_hi) * x_to_hi
        + m_Histogram( bin_x_hi, bin_y_lo, bin_z_hi) * x_to_lo
      );
      const double c10
      (
        m_Histogram( bin_x_lo, bin_y_hi, bin_z_lo) * x_to_hi
        + m_Histogram( bin_x_hi, bin_y_hi, bin_z_lo) * x_to_lo
      );
      const double c11
      (
        m_Histogram( bin_x_lo, bin_y_hi, bin_z_hi) * x_to_hi
        + m_Histogram( bin_x_hi, bin_y_hi, bin_z_hi) * x_to_lo
      );

      const double c0( c00 * y_to_hi + c10 * y_to_lo);
      const double c1( c01 * y_to_hi + c11 * y_to_lo);

      return ( c0 * z_to_hi + c1 * z_to_lo) / ( m_BinSizeXYZ.First() * m_BinSizeXYZ.Second() * m_BinSizeXYZ.Third());
    }

    //! Get the value for the nearest bin
    //! @param X, Y, Z the coordinates of interest
    double Histogram3D::Value( const double &X, const double &Y, const double &Z) const
    {
      if( !util::IsDefined( X) || !util::IsDefined( Y) || !util::IsDefined( Z))
      {
        return util::GetUndefined< double>();
      }

      size_t bin_x( 0), bin_y( 0), bin_z( 0);

      //determine field in x-direction
      if( X < m_BinningX( 0) || !m_BinSizeXYZ.First())
      {
        bin_x = 0;
      }
      else if( X >= m_LowerUpperBoundariesX.Second())
      {
        bin_x = m_Histogram.NumberLayers() - 1;
      }
      else
      {
        bin_x = size_t( ( X - m_LowerUpperBoundariesX.First()) / m_BinSizeXYZ.First());
      }

      //determine field in y-direction
      if( Y < m_BinningY( 0) || !m_BinSizeXYZ.Second())
      {
        bin_y = 0;
      }
      else if( Y >= m_LowerUpperBoundariesY.Second())
      {
        bin_y = m_Histogram.GetNumberRows() - 1;
      }
      else
      {
        bin_y = size_t( ( Y - m_LowerUpperBoundariesY.First()) / m_BinSizeXYZ.Second());
      }

      //determine field in z-direction
      if( Z < m_BinningZ( 0) || !m_BinSizeXYZ.Third())
      {
        bin_z = 0;
      }
      else if( Z >= m_LowerUpperBoundariesZ.Second())
      {
        bin_z = m_Histogram.GetNumberCols() - 1;
      }
      else
      {
        bin_z = size_t( ( Z - m_LowerUpperBoundariesZ.First()) / m_BinSizeXYZ.Third());
      }
      return m_Histogram( bin_x, bin_y, bin_z);
    }

    //! @brief combine this with a given histogram2d by adding up all counts
    //! all parameters have to be identical for this operation to work
    //! @param HISTOGRAM_3D histogram to be added to this one
    //! @return true if it was successful
    bool Histogram3D::Combine( const Histogram3D &HISTOGRAM_3D)
    {
      // allow combining, if this is empty
      if( IsEmpty())
      {
        *this = HISTOGRAM_3D;
        return true;
      }

      // check that the parameters agree
      if
      (
           m_LowerUpperBoundariesX  != HISTOGRAM_3D.m_LowerUpperBoundariesX
        || m_LowerUpperBoundariesY  != HISTOGRAM_3D.m_LowerUpperBoundariesY
        || m_LowerUpperBoundariesZ  != HISTOGRAM_3D.m_LowerUpperBoundariesZ
        || m_BinSizeXYZ             != HISTOGRAM_3D.m_BinSizeXYZ
      )
      {
        return false;
      }

      // add the histogram
      linal::VectorReference< double> ref( m_Histogram.GetValues());
      ref += HISTOGRAM_3D.m_Histogram.GetValues();

      // end
      return true;
    }

    //! checks if there is any count
    bool Histogram3D::IsEmpty() const
    {
      return GetSumOfAllCounts() == 0;
    }

    //! @brief normalize the counts in the histogram
    //! each bin is divided by the total number of counts (even the boundary counts)
    void Histogram3D::Normalize()
    {
      m_Histogram.GetValues().SetToSum( 1.0);
    }

    //! @brief Normalize by another histogram 3d, which represents background counts
    //! @param BACKGROUND background counts histogram
    void Histogram3D::NormalizeByBackground( const Histogram3D &HIST)
    {
      BCL_Assert
      (
        m_LowerUpperBoundariesX == HIST.m_LowerUpperBoundariesX
        && m_LowerUpperBoundariesY == HIST.m_LowerUpperBoundariesY
        && m_LowerUpperBoundariesZ == HIST.m_LowerUpperBoundariesZ
        && m_BinSizeXYZ == HIST.m_BinSizeXYZ,
        "Normalization could not occur because histograms were of different size or binning"
      );
      auto itr_bg( HIST.m_Histogram.GetValues().Begin());
      for
      (
        auto itr_this( m_Histogram.Begin()), itr_this_end( m_Histogram.End());
        itr_this != itr_this_end;
        ++itr_this, ++itr_bg
      )
      {
        // normalize each non-zero bin
        if( *itr_this && *itr_bg)
        {
          *itr_this /= *itr_bg;
        }
        else
        {
          *itr_this = 0.0;
        }
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! write into std::ostream x and y values horizontally
    std::ostream &Histogram3D::Write
    (
      std::ostream &OSTREAM,
      const util::Format &FORMAT_BINNING,
      const util::Format &FORMAT_VALUES
    ) const
    {
      // Write as a simple X Y Z Value file (readable in VisIt)
      OSTREAM << "X Y Z value\n";
      for( size_t x_i( 0), x_mx( m_Histogram.NumberLayers()); x_i < x_mx; ++x_i)
      {
        const double x( m_LowerUpperBoundariesX.First() + m_BinSizeXYZ.First() * ( 0.5 + double( x_i)));
        for( size_t y_i( 0), y_mx( m_Histogram.GetNumberRows()); y_i < y_mx; ++y_i)
        {
          const double y( m_LowerUpperBoundariesY.First() + m_BinSizeXYZ.Second() * ( 0.5 + double( y_i)));
          for( size_t z_i( 0), z_mx( m_Histogram.GetNumberCols()); z_i < z_mx; ++z_i)
          {
            const double z( m_LowerUpperBoundariesZ.First() + m_BinSizeXYZ.Third() * ( 0.5 + double( z_i)));
            OSTREAM << FORMAT_BINNING( x) << ' ' << FORMAT_BINNING( y) << ' ' << FORMAT_BINNING( z) << ' '
                    << FORMAT_VALUES( m_Histogram( x_i, y_i, z_i)) << '\n';
          }
        }
      }

      return OSTREAM;
    }

    //! read Histogram3D from std::istream
    std::istream &Histogram3D::Read( std::istream &ISTREAM)
    {
      storage::Vector< double> values;
      storage::Vector< double> x_values;
      storage::Vector< double> y_values;
      storage::Vector< double> z_values;
      std::string tmp;
      std::getline( ISTREAM, tmp);
      if( tmp.empty())
      {
        std::getline( ISTREAM, tmp);
      }
      BCL_Assert( util::StartsWith( tmp, "X Y Z "), "Histogram 3D corrupted. Should have started with X Y Z");
      storage::Vector< double> tmp_values;
      while( std::getline( ISTREAM, tmp))
      {
        tmp_values = util::SplitStringToNumerical< double>( tmp, " ");
        if( tmp_values.IsEmpty())
        {
          break;
        }
        BCL_Assert
        (
          tmp_values.GetSize() == 4,
          "Histogram3D Corrupted. Should have had 4 values on each line, instead had: " + tmp
        );
        values.PushBack( tmp_values( 3));
        if( x_values.IsEmpty())
        {
          x_values.PushBack( tmp_values( 0));
          y_values.PushBack( tmp_values( 1));
          z_values.PushBack( tmp_values( 2));
        }
        else if( tmp_values( 0) > x_values.LastElement())
        {
          x_values.PushBack( tmp_values( 0));
        }
        else if( tmp_values( 1) > y_values.LastElement())
        {
          y_values.PushBack( tmp_values( 1));
        }
        else if( tmp_values( 2) > z_values.LastElement())
        {
          z_values.PushBack( tmp_values( 2));
        }
      }
      m_BinSizeXYZ( 0) = x_values.GetSize() > size_t( 1)
                         ? ( x_values.LastElement() - x_values.FirstElement()) / double( x_values.GetSize() - 1)
                         : 0.0;
      m_BinSizeXYZ( 1) = y_values.GetSize() > size_t( 1)
                         ? ( y_values.LastElement() - y_values.FirstElement()) / double( y_values.GetSize() - 1)
                         : 0.0;
      m_BinSizeXYZ( 2) = z_values.GetSize() > size_t( 1)
                         ? ( z_values.LastElement() - z_values.FirstElement()) / double( z_values.GetSize() - 1)
                         : 0.0;
      m_LowerUpperBoundariesX.First() = x_values( 0) - m_BinSizeXYZ( 0) / 2.0;
      m_LowerUpperBoundariesY.First() = y_values( 0) - m_BinSizeXYZ( 1) / 2.0;
      m_LowerUpperBoundariesZ.First() = z_values( 0) - m_BinSizeXYZ( 2) / 2.0;
      m_LowerUpperBoundariesX.Second() = x_values.LastElement() + m_BinSizeXYZ( 0) / 2.0;
      m_LowerUpperBoundariesY.Second() = y_values.LastElement() + m_BinSizeXYZ( 1) / 2.0;
      m_LowerUpperBoundariesZ.Second() = z_values.LastElement() + m_BinSizeXYZ( 2) / 2.0;
      if( values.GetSize())
      {
        m_Histogram = Tensor< double>( x_values.GetSize(), y_values.GetSize(), z_values.GetSize(), &values( 0));
      }
      else
      {
        m_Histogram = Tensor< double>();
      }
      SetupBinning();

      //return
      return ISTREAM;
    }

    //! @brief setup bins according to the boundaries
    void Histogram3D::SetupBinning()
    {
      m_BinningX = linal::Vector< double>( m_Histogram.NumberLayers(), double( 0.0));
      m_BinningY = linal::Vector< double>( m_Histogram.GetNumberRows(), double( 0.0));
      m_BinningZ = linal::Vector< double>( m_Histogram.GetNumberCols(), double( 0.0));

      //binning in X direction
      double i( 0.5);
      for( double *ptr( m_BinningX.Begin()), *ptr_end( m_BinningX.End()); ptr != ptr_end; ++ptr, ++i)
      {
        *ptr = m_LowerUpperBoundariesX.First() + i * m_BinSizeXYZ.First();
      }

      //binning in Y direction
      i = 0.5;
      for( double *ptr( m_BinningY.Begin()), *ptr_end( m_BinningY.End()); ptr != ptr_end; ++ptr, ++i)
      {
        *ptr = m_LowerUpperBoundariesY.First() + i * m_BinSizeXYZ.Second();
      }

      //binning in Z direction
      i = 0.5;
      for( double *ptr( m_BinningZ.Begin()), *ptr_end( m_BinningZ.End()); ptr != ptr_end; ++ptr, ++i)
      {
        *ptr = m_LowerUpperBoundariesZ.First() + i * m_BinSizeXYZ.Third();
      }
    }

  } // namespace math
} // namespace bcl
