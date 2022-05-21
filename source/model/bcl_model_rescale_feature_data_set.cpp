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
#include "model/bcl_model_rescale_feature_data_set.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "math/bcl_math_running_average_sd.h"
#include "model/bcl_model_feature_data_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> RescaleFeatureDataSet::s_Instance
    (
      GetObjectInstances().AddInstance( new RescaleFeatureDataSet())
    );

    //! @brief Type as string
    //! @param TYPE the type
    //! @return the string for the kernel
    const std::string &RescaleFeatureDataSet::GetTypeName( const Type &TYPE)
    {
      static const std::string s_names[] =
      {
        "None",
        "MinMax",
        "AveStd",
        GetStaticClassName< Type>()
      };
      return s_names[ TYPE];
    }

    //! global static function that allows altering the default null value; note that undefined is interpreted as the
    //! middle of the output range
    float &RescaleFeatureDataSet::GetDefaultValueForEmptyRangedColumns()
    {
      static float s_default_value( util::GetUndefined< float>());
      return s_default_value;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RescaleFeatureDataSet::RescaleFeatureDataSet() :
      m_RescaleRanges(),
      m_Range()
    {
    }

    //! @brief construct from feature result data set
    //! @param DATA the feature data set from which to obtain the ranges
    //! @param RANGE features will be rescaled to this range
    //! @param TYPE type of rescaling to perform
    RescaleFeatureDataSet::RescaleFeatureDataSet
    (
      const FeatureDataSetInterface< float> &DATA,
      const math::Range< float> &RANGE,
      const Type &TYPE
    ) :
      m_RescaleRanges(),
      m_Range( RANGE)
    {
      // get the ranges from the underlying matrix
      if( TYPE == e_AveStd)
      {
        GetInputRangesFromMatrixAveStd( DATA.GetMatrix());
      }
      else if( TYPE == e_MinMax)
      {
        GetInputRangesFromMatrix( DATA.GetMatrix());
      }
      else
      {
        // non-rescaler
        m_RescaleFunction.Resize( DATA.GetFeatureSize(), math::LinearFunction( 1.0, 0.0));
        m_DescaleFunction.Resize( DATA.GetFeatureSize(), math::LinearFunction( 1.0, 0.0));
      }
    }

    //! @brief construct from a matrix const interface
    //! @param DATA the matrix from which to obtain the ranges
    //! @param RANGE features will be rescaled to this range
    //! @param TYPE type of rescaling to perform
    RescaleFeatureDataSet::RescaleFeatureDataSet
    (
      const linal::MatrixConstInterface< float> &DATA,
      const math::Range< float> &RANGE,
      const Type &TYPE
    ) :
      m_RescaleRanges(),
      m_Range( RANGE)
    {
      if( TYPE == e_AveStd)
      {
        GetInputRangesFromMatrixAveStd( DATA);
      }
      else if( TYPE == e_MinMax)
      {
        GetInputRangesFromMatrix( DATA);
      }
      else
      {
        m_RescaleFunction.Resize( DATA.GetNumberCols(), math::LinearFunction( 1.0, 0.0));
        m_DescaleFunction.Resize( DATA.GetNumberCols(), math::LinearFunction( 1.0, 0.0));
      }
    }

    //! @brief construct from from ranges to to range
    //! @param FROM_RANGES the original ranges
    //! @param TO_RANGE the rescaled range
    RescaleFeatureDataSet::RescaleFeatureDataSet
    (
      const RescaleFeatureDataSet &RESCALING,
      const util::ShPtrVector< math::FunctionInterfaceSerializable< double, double> > &MODIFIERS
    ) :
      m_RescaleRanges( RESCALING.m_RescaleRanges),
      m_Range( RESCALING.m_Range)
    {
      BCL_Assert( MODIFIERS.GetSize() == m_RescaleRanges.GetSize(), "Require as many modifiers as rescaling ranges");
      for( size_t col( 0), n_cols( MODIFIERS.GetSize()); col < n_cols; ++col)
      {
        m_RescaleRanges( col)
          = math::Range< float>
            (
              ( *MODIFIERS( col))( m_RescaleRanges( col).GetMin()),
              ( *MODIFIERS( col))( m_RescaleRanges( col).GetMax())
            );
      }
      ConstructScalingFunctions();
    }

    //! @brief Clone function
    //! @return pointer to new RescaleFeatureDataSet
    RescaleFeatureDataSet *RescaleFeatureDataSet::Clone() const
    {
      return new RescaleFeatureDataSet( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RescaleFeatureDataSet::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief rescale a value for a particular column using the appropriate rescale function
    //! @param COLUMN column that the value logically belongs to
    //! @param VALUE value to rescale
    //! @return the rescaled value
    float RescaleFeatureDataSet::RescaleValue( const size_t &COLUMN, const float &VALUE) const
    {
      return m_RescaleFunction( COLUMN)( VALUE);
    }

    //! @brief descale a value for a particular column using the appropriate descale function
    //! @param COLUMN column that the value logically belongs to
    //! @param VALUE value to descale
    //! @return the descaled value
    float RescaleFeatureDataSet::DescaleValue( const size_t &COLUMN, const float &VALUE) const
    {
      return m_DescaleFunction( COLUMN)( VALUE);
    }

    //! @brief scales back to pre-scaled values
    //! @param FEATURE the feature data set interface to scale back to original values
    //! @return FeatureDataSet of de-scaled values
    FeatureDataSet< float> RescaleFeatureDataSet::DeScale( const FeatureDataSetInterface< float> &FEATURE) const
    {
      // create a matrix copy
      linal::Matrix< float> mat( FEATURE.GetMatrix());

      // descale it
      DeScaleMatrix( mat);

      // create a new data set out of it
      return FeatureDataSet< float>( mat);
    }

    //! @brief scales back to pre-rescaled values
    //! @param MATRIX the matrix to descale
    void RescaleFeatureDataSet::DeScaleMatrix( linal::MatrixInterface< float> &MATRIX) const
    {
      // check for identity rescaling
      if( m_DescaleFunction.IsEmpty() || m_RescaleRanges.IsEmpty())
      {
        return;
      }

      // make sure the matrix has the right number of columns, if it is not empty
      BCL_Assert
      (
        MATRIX.GetNumberRows() == size_t( 0) || MATRIX.GetNumberCols() == m_RescaleRanges.GetSize(),
        "DeScaleMatrix given matrix with incorrect # columns!"
      );

      // get matrix and dimensionality info
      const size_t number_features( MATRIX.GetNumberCols());
      const size_t number_data_pts( MATRIX.GetNumberRows());

      // scale back to original values
      for( size_t row_id( 0); row_id < number_data_pts; ++row_id)
      {
        storage::Vector< math::LinearFunction>::const_iterator itr_func( m_DescaleFunction.Begin());
        for( float *itr( MATRIX[ row_id]), *itr_end( itr + number_features); itr != itr_end; ++itr, ++itr_func)
        {
          *itr = ( *itr_func)( *itr);
        }
      }
    }

    //! @brief scale a MatrixInterface
    //! @param MATRIX the matrix to rescale
    void RescaleFeatureDataSet::RescaleMatrix( linal::MatrixInterface< float> &MATRIX) const
    {
      // check for identity rescaling
      if( m_RescaleFunction.IsEmpty())
      {
        return;
      }
      BCL_Assert
      (
        m_RescaleFunction.GetSize() == MATRIX.GetNumberCols(),
        "Wrong number of rescaling columns! " + util::Format()( m_RescaleFunction.GetSize())
        + " vs " + util::Format()( MATRIX.GetNumberCols())
      );

      // get matrix and dimensionality info
      const size_t number_features( MATRIX.GetNumberCols());
      const size_t number_data_pts( MATRIX.GetNumberRows());

      // scale back to original values
      for( size_t row_id( 0); row_id < number_data_pts; ++row_id)
      {
        storage::Vector< math::LinearFunction>::const_iterator itr_func( m_RescaleFunction.Begin());
        for( float *itr( MATRIX[ row_id]), *itr_end( itr + number_features); itr != itr_end; ++itr, ++itr_func)
        {
          *itr = ( *itr_func)( *itr);
        }
      }
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief scale a FeatureDataSet
    //! @param FEATURE the FeatureDataSetInterface to rescale
    //! @return rescaled FeatureDataSet to one range
    FeatureDataSet< float> RescaleFeatureDataSet::operator()( const FeatureDataSetInterface< float> &FEATURE) const
    {
      // get matrix and dimensionality info
      linal::Matrix< float> mat( FEATURE.GetMatrix());
      RescaleMatrix( mat);

      // return scaled feature data set
      return FeatureDataSet< float>( mat);
    }

    //! @brief test for inequality
    //! @param OTHER other rescaling function
    //! @return true if the rescaling functions are inequal
    bool RescaleFeatureDataSet::operator !=( const RescaleFeatureDataSet &OTHER) const
    {
      return !( m_Range == OTHER.m_Range) || !( m_RescaleRanges == OTHER.m_RescaleRanges);
    }

    //! @brief test for equality
    //! @param OTHER other rescaling function
    //! @return true if the rescaling functions are equal
    bool RescaleFeatureDataSet::operator ==( const RescaleFeatureDataSet &OTHER) const
    {
      return m_Range == OTHER.m_Range && m_RescaleRanges == OTHER.m_RescaleRanges;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RescaleFeatureDataSet::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_RescaleRanges, ISTREAM);
      io::Serialize::Read( m_Range, ISTREAM);
      ConstructScalingFunctions();

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &RescaleFeatureDataSet::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write
      io::Serialize::Write( m_RescaleRanges, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Range, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief determine m_RescaleRanges from the given input matrix
    //! @param DATA the matrix from which to obtain the ranges
    void RescaleFeatureDataSet::GetInputRangesFromMatrix( const linal::MatrixConstInterface< float> &DATA)
    {
      m_RescaleRanges.Reset();

      // pre-allocating vector memory
      const size_t number_data_pts( DATA.GetNumberRows());
      m_RescaleRanges.AllocateMemory( DATA.GetNumberCols());

      // find ranges of each column
      float tmp( 0);
      for( size_t col( 0), number_features( DATA.GetNumberCols()); col < number_features; ++col)
      {
        // set initial min and max
        float min( std::numeric_limits< float>::infinity());
        float max( -std::numeric_limits< float>::infinity());

        // iterate through all values in that feature col
        // no need to worry about nans here, b/c they will always fail both the < and > comparison
        for( size_t row( 0); row < number_data_pts; ++row)
        {
          // set as min or max if qualifies
          tmp = DATA( row, col);
          if( util::IsDefined( tmp))
          {
            if( tmp < min)
            {
              min = tmp;
            }
            else if( tmp > max)
            {
              max = tmp;
            }
          }
        }
        // add to range vector
        if( min < max)
        {
          m_RescaleRanges.PushBack( math::Range< float>( min, max));
        }
        else
        {
          m_RescaleRanges.PushBack( math::Range< float>( 0.0, 0.0));
        }
      }
      ConstructScalingFunctions();
    }

    //! @brief determine m_RescaleRanges from the given input matrix
    //! @param DATA the matrix from which to obtain the ranges
    void RescaleFeatureDataSet::GetInputRangesFromMatrixAveStd( const linal::MatrixConstInterface< float> &DATA)
    {
      m_RescaleRanges.Reset();

      // pre-allocating vector memory
      m_RescaleRanges.Resize( DATA.GetNumberCols());

      math::RunningAverageSD< linal::Vector< float> > ave_std;

      // compute averages and std for whole dataset
      for( size_t row( 0), n_rows( DATA.GetNumberRows()); row < n_rows; ++row)
      {
        ave_std += DATA.GetRow( row);
      }

      // find ranges of each column
      const linal::Vector< float> &ave( ave_std.GetAverage());
      const linal::Vector< float> &std( ave_std.GetStandardDeviation());

      for( size_t col( 0), number_features( DATA.GetNumberCols()); col < number_features; ++col)
      {
        // THIS LINE OF CODE IS CRITICAL; otherwise floating point numbers smaller than float min may be written out,
        // but they cannot be read back in; at least on gcc 4.7.  This is likely a compiler/library bug
        const float average( math::Absolute( ave( col)) < std::numeric_limits< float>::min() ? 0.0 : ave( col));

        // determine whether the column has any variance
        if( !util::IsDefined( std( col)))
        {
          m_RescaleRanges( col) = math::Range< float>( average, average);
        }
        else
        {
          m_RescaleRanges( col)
            = math::Range< float>
              (
                average - float( 2.0) * std( col),
                average + float( 2.0) * std( col)
              );
        }
      }
      ConstructScalingFunctions();
    }

    //! @brief determine m_RescaleRanges from the given input matrix
    //! @param DATA the matrix from which to obtain the ranges
    //! @param FRACTION (0,1) the minimum fraction of each column's values that should be covered by the chosen range
    void RescaleFeatureDataSet::GetInputRangesFromMatrixRange
    (
      const linal::MatrixConstInterface< float> &DATA,
      const float &FRACTION
    )
    {
      m_RescaleRanges.Reset();

      // pre-allocating vector memory
      m_RescaleRanges.Resize( DATA.GetNumberCols());

      // create a temporary vector for holding each column
      linal::Vector< float> column_storage( DATA.GetNumberRows());

      // determine the sorted feature range desired to be between -1 and 1 for each column
      const size_t sorted_feature_range_min( DATA.GetNumberRows() * ( 1.0 - FRACTION) / 2.0);
      const size_t sorted_feature_range_max( ( DATA.GetNumberRows() - 1) * FRACTION / 2.0);

      // for each column
      for( size_t col( 0), number_features( DATA.GetNumberCols()); col < number_features; ++col)
      {
        // load all data from the column into the vector
        for( size_t row( 0), n_rows( DATA.GetNumberRows()); row < n_rows; ++row)
        {
          column_storage( row) = DATA( row, col);
        }

        // use the nth-element algorithm to discover the actual value range that would cover 80% of the data
        linal::Vector< float>::iterator min_range_itr( column_storage.Begin() + sorted_feature_range_min);
        std::nth_element( column_storage.Begin(), min_range_itr, column_storage.End());
        float min_range_value( *min_range_itr);
        linal::Vector< float>::iterator max_range_itr( column_storage.Begin() + sorted_feature_range_max);
        std::nth_element( min_range_itr, max_range_itr, column_storage.End());
        float max_range_value( *max_range_itr);
        if( min_range_value == max_range_value)
        {
          // probably a binary or nearly binary variable; take the end values instead
          linal::Vector< float>::iterator itr_min( column_storage.Begin()), itr_max( column_storage.End() - 1);
          std::nth_element( column_storage.Begin(), itr_min, min_range_itr);
          std::nth_element( max_range_itr, itr_max, column_storage.End());
          min_range_value = *itr_min;
          max_range_value = *itr_max;
          BCL_MessageCrt( "Reverting to default scaling on descriptor: " + util::Format()( col));
        }

        m_RescaleRanges( col) = math::Range< float>( min_range_value, max_range_value);
      }

      ConstructScalingFunctions();
    }

    //! @brief construct rescaling and descaling functions
    void RescaleFeatureDataSet::ConstructScalingFunctions()
    {
      const float midpoint
      (
        util::IsDefined( GetDefaultValueForEmptyRangedColumns())
        ? GetDefaultValueForEmptyRangedColumns()
        : m_Range.GetMiddle()
      );
      const size_t number_columns( m_RescaleRanges.GetSize());
      m_RescaleFunction.Resize( number_columns);
      m_DescaleFunction.Resize( number_columns);
      for( size_t column( 0); column < number_columns; ++column)
      {
        // check range, if ~0, set all values to midpoint
        if( m_RescaleRanges( column).GetWidth() < std::numeric_limits< float>::epsilon())
        {
          m_RescaleFunction( column) = math::LinearFunction( 0.0, midpoint);
          m_DescaleFunction( column) = math::LinearFunction( 0.0, m_RescaleRanges( column).GetMiddle());
        }
        else
        {
          m_RescaleFunction( column) = math::LinearFunction( m_RescaleRanges( column), m_Range);
          m_DescaleFunction( column) = math::LinearFunction( m_Range, m_RescaleRanges( column));
        }
      }
    }

  } // namespace model
} // namespace bcl
