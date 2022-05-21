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

// include header of this class
#include "bcl_model_feature_data_set.h"
// includes from bcl - sorted alphabetically
#include "bcl_model_feature_reference.h"
#include "bcl_model_rescale_feature_data_set.h"
#include "linal/bcl_linal_matrix.h"
#include "linal/bcl_linal_matrix_const_reference.h"
#include "linal/bcl_linal_matrix_reference.h"
#include "linal/bcl_linal_vector.h"
#include "math/bcl_math_range_set.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    template< typename t_DataType>
    FeatureDataSet< t_DataType>::FeatureDataSet() :
      m_Data(),
      m_Rescale(),
      m_FeatureLabels()
    {
    }

    //! @brief constructor from interface
    template< typename t_DataType>
    FeatureDataSet< t_DataType>::FeatureDataSet
    (
      const FeatureDataSetInterface< t_DataType> &DATA_SET
    ) :
      m_Data( DATA_SET.GetMatrix()),
      m_Rescale( DATA_SET.GetScaling().HardCopy()),
      m_FeatureLabels( DATA_SET.GetFeatureLabelSet())
    {
    }

    //! @brief constructor from number of examples and feature size
    //! @param EXAMPLES number of distinct entities in the dataset (aka number rows or number features)
    //! @param FEATURE_SIZE size of each feature in the data set (aka number cols or feature size)
    //! @param FILL_VALUE value to fill the matrix with
    template< typename t_DataType>
    FeatureDataSet< t_DataType>::FeatureDataSet
    (
      const size_t &EXAMPLES,
      const size_t &FEATURE_SIZE,
      const t_DataType &FILL_VALUE
    ) :
      m_Data( EXAMPLES, FEATURE_SIZE, FILL_VALUE),
      m_Rescale()
    {
    }

    //! @brief constructor from matrix
    //! @param MATRIX the matrix to construct from
    template< typename t_DataType>
    FeatureDataSet< t_DataType>::FeatureDataSet( const linal::MatrixConstInterface< t_DataType> &MATRIX) :
      m_Data( MATRIX),
      m_Rescale()
    {
    }

    //! @brief constructor from matrix and existing rescaling
    //! @param MATRIX the matrix to construct from
    //! @param RESCALING existing rescaling
    template< typename t_DataType>
    FeatureDataSet< t_DataType>::FeatureDataSet
    (
      const linal::MatrixConstInterface< t_DataType> &MATRIX,
      const RescaleFeatureDataSet &RESCALING
    ) :
      m_Data( MATRIX),
      m_Rescale( RESCALING.Clone())
    {
    }

    //! @brief constructor from matrix and set of columns
    //! @param MATRIX the matrix to construct from
    //! @param COLS the colums to copy
    template< typename t_DataType>
    FeatureDataSet< t_DataType>::FeatureDataSet
    (
      const linal::MatrixConstInterface< t_DataType> &MATRIX,
      const math::RangeSet< size_t> &COLS
    )
    {
      // only use the range set if ranges were actually given
      if( MATRIX.GetNumberOfElements() == 0 || COLS.IsEmpty())
      {
        m_Data = linal::Matrix< t_DataType>( MATRIX);
      }

      // create a vector of column indices that will be kept
      storage::Vector< size_t> cols_to_keep;
      for
      (
        storage::Set< math::Range< size_t> >::const_iterator
          itr_range( COLS.GetRanges().Begin()), itr_range_end( COLS.GetRanges().End());
        itr_range != itr_range_end;
        ++itr_range
      )
      {
        cols_to_keep.AllocateMemory( cols_to_keep.GetSize() + itr_range->GetWidth());
        for
        (
          size_t col( itr_range->GetMin()), end_col( itr_range->GetMax());
          col <= end_col;
          ++col
        )
        {
          cols_to_keep.PushBack( col);
        }
      }

      BCL_Assert( !cols_to_keep.IsEmpty(), "Invalid or empty range given");

      BCL_Assert
      (
        cols_to_keep.LastElement() < MATRIX.GetNumberCols(),
        "Out of bounds column " + util::Format()( cols_to_keep.LastElement())
      );

      // create the matrix to be of appropriate dimensions
      m_Data = linal::Matrix< t_DataType>( MATRIX.GetNumberRows(), cols_to_keep.GetSize());

      // copy the data elements into m_Data
      for( size_t row( 0), n_rows( MATRIX.GetNumberRows()); row < n_rows; ++row)
      {
        const t_DataType *src_row_ptr( MATRIX[ row]);
        t_DataType *data_row_ptr( m_Data[ row]);
        for( size_t col( 0), n_cols( cols_to_keep.GetSize()); col < n_cols; ++col)
        {
          data_row_ptr[ col] = src_row_ptr[ cols_to_keep( col)];
        }
      }
    }

    //! @brief Clone function
    //! @return pointer to new FeatureDataSet< t_DataType>
    template< typename t_DataType>
    FeatureDataSet< t_DataType> *FeatureDataSet< t_DataType>::Clone() const
    {
      return new FeatureDataSet< t_DataType>( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_DataType>
    const std::string &FeatureDataSet< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get size of
    //! @return the number of data items in the feature
    template< typename t_DataType>
    size_t FeatureDataSet< t_DataType>::GetFeatureSize() const
    {
      return m_Data.GetNumberCols();
    }

    //! @brief number of features
    //! @return the number of features
    template< typename t_DataType>
    size_t FeatureDataSet< t_DataType>::GetNumberFeatures() const
    {
      return m_Data.GetNumberRows();
    }

    //! @brief get a ShPtr to the feature label set, which indicates what each column represents
    //! @return a shptr to the feature label set
    template< typename t_DataType>
    const util::ShPtr< FeatureLabelSet> &FeatureDataSet< t_DataType>::GetFeatureLabelSet() const
    {
      return m_FeatureLabels;
    }

    //! @brief set a ShPtr to the feature label set, which indicates what each column represents
    //! @param NEW_SET shptr to the new feature label set
    template< typename t_DataType>
    void FeatureDataSet< t_DataType>::SetFeatureLabelSet( const util::ShPtr< FeatureLabelSet> &NEW_SET)
    {
      m_FeatureLabels = NEW_SET;
    }

    //! @brief set a ShPtr to the feature label set, which indicates what each column represents
    template< typename t_DataType>
    void FeatureDataSet< t_DataType>::SetFeatureLabelSet( const FeatureLabelSet &NEW_SET)
    {
      m_FeatureLabels = util::ShPtr< FeatureLabelSet>( NEW_SET.Clone());
    }

    //! @brief get the owning matrix
    //! @return matrixconstinterface
    template< typename t_DataType>
    linal::MatrixConstReference< t_DataType> FeatureDataSet< t_DataType>::GetMatrix() const
    {
      return linal::MatrixConstReference< t_DataType>( m_Data);
    }

    //! @brief get owning matrix at POS of size LENGTH
    //! @param POS position from which to start the matrix to be returned
    //! @param LENGTH the length in rows of the matrix to be returned
    //! @return the matrixconstinterface from POS of size LENGTH rows
    template< typename t_DataType>
    linal::MatrixConstReference< t_DataType> FeatureDataSet< t_DataType>::GetMatrix( const size_t POS, const size_t LENGTH) const
    {
      return linal::MatrixConstReference< t_DataType>( LENGTH, m_Data.GetNumberCols(), operator[]( POS));
    }

    //! @brief get a changeable reference to the matrix
    //! @return matrix reference
    template< typename t_DataType>
    linal::MatrixReference< t_DataType> FeatureDataSet< t_DataType>::GetMatrix()
    {
      return linal::MatrixReference< t_DataType>( m_Data);
    }

    //! @brief test whether this feature data set is rescaled
    //! @return true if the dataset is rescaled
    template< typename t_DataType>
    bool FeatureDataSet< t_DataType>::IsRescaled() const
    {
      return m_Rescale.IsDefined();
    }

    //! @brief get a pointer to the rescaling object
    //! @return the rescaling object
    template< typename t_DataType>
    const util::ShPtr< RescaleFeatureDataSet> &FeatureDataSet< t_DataType>::GetScaling() const
    {
      return m_Rescale;
    }

    //! @brief test whether this feature data set has the same scaling as another
    //! @param OTHER the other feature dataset
    //! @return true if this feature data set has the same scaling as OTHER
    template< typename t_DataType>
    bool FeatureDataSet< t_DataType>::HasSameScaling( const FeatureDataSetInterface< t_DataType> &OTHER) const
    {
      if( m_Rescale == OTHER.GetScaling())
      {
        return true;
      }
      return m_Rescale.IsDefined() == OTHER.GetScaling().IsDefined()
             && ( !m_Rescale.IsDefined() || *m_Rescale == *OTHER.GetScaling());
    }

  ////////////////
  // operations //
  ////////////////

    namespace
    {
      // anonymous namespace to prevent export of functions

      //! @brief scales back to pre-rescaled values
      //! @param RESCALER the rescaling object
      //! @param MATRIX the matrix to descale
      //! @return FeatureDataSet of scaled values
      void Descale
      (
        util::ShPtr< RescaleFeatureDataSet> &RESCALER,
        linal::Matrix< float> &MATRIX
      )
      {
        if( RESCALER.IsDefined())
        {
          RESCALER->DeScaleMatrix( MATRIX);
          RESCALER = util::ShPtr< RescaleFeatureDataSet>();
        }
      }

      //! @brief trivial implementation for char
      //! @param RESCALER the rescaling object
      //! @param MATRIX the matrix to descale
      //! @return FeatureDataSet of scaled values
      void Descale
      (
        util::ShPtr< RescaleFeatureDataSet> &RESCALER,
        linal::Matrix< char> &MATRIX
      )
      {
        BCL_Assert( !RESCALER.IsDefined(), "Character matrix should never have a rescaling!");
      }

      //! @brief rescale a matrix given a rescaler
      //! @param CURRENT_RESCALER the current rescaling object
      //! @param MATRIX the current matrix
      void RescaleImpl
      (
        util::ShPtr< RescaleFeatureDataSet> &CURRENT_RESCALER,
        linal::Matrix< float> &MATRIX
      )
      {
        CURRENT_RESCALER->RescaleMatrix( MATRIX);
      }

      //! @brief rescale a matrix given a rescaler
      //! @param CURRENT_RESCALER the current rescaling object
      //! @param MATRIX the current matrix
      void RescaleImpl
      (
        util::ShPtr< RescaleFeatureDataSet> &CURRENT_RESCALER,
        linal::Matrix< char> &MATRIX
      )
      {
        BCL_Exit( "Character matrices cannot be rescaled!", -1);
      }

      //! @brief rescale to the target range, works even if presently scaled to a different range
      //! @param MATRIX the current matrix
      //! @param TO_RANGE the to range
      //! @param SCALING the scaling
      util::ShPtr< RescaleFeatureDataSet> CreateRescale
      (
        linal::Matrix< float> &MATRIX,
        const math::Range< float> &TO_RANGE,
        const RescaleFeatureDataSet::Type &SCALING
      )
      {
        return util::ShPtr< RescaleFeatureDataSet>
               (
                 new RescaleFeatureDataSet( MATRIX, TO_RANGE, SCALING)
               );
      }

      //! @brief Trivial implementation for character ranges
      //! @param MATRIX the current matrix
      //! @param TO_RANGE the to range
      //! @param SCALING the scaling
      util::ShPtr< RescaleFeatureDataSet> CreateRescale
      (
        linal::Matrix< char> &MATRIX,
        const math::Range< float> &TO_RANGE,
        const RescaleFeatureDataSet::Type &SCALING
      )
      {
        BCL_Assert( SCALING == RescaleFeatureDataSet::e_None, "Character matrix should never have a rescaling!");
        return util::ShPtr< RescaleFeatureDataSet>();
      }

    }

    //! @brief scales back to pre-rescaled values
    //! @return FeatureDataSet of scaled values
    template< typename t_DataType>
    FeatureDataSet< t_DataType> &FeatureDataSet< t_DataType>::DeScale()
    {
      // trivial case; non-rescaled data
      Descale( m_Rescale, m_Data);
      return *this;
    }

    //! @brief rescale to the target range, works even if presently scaled to a different range
    //! @param SCALING type of scaling to use
    template< typename t_DataType>
    FeatureDataSet< t_DataType> &FeatureDataSet< t_DataType>::Rescale
    (
      const math::Range< float> &TO_RANGE,
      const RescaleFeatureDataSet::Type &SCALING
    )
    {
      // handle existing scaling
      if( IsRescaled())
      {
        if
        (
          m_Rescale->GetRange().GetMin() == TO_RANGE.GetMin()
          && m_Rescale->GetRange().GetMax() == TO_RANGE.GetMax()
        )
        {
          // already have the correct rescaling, return
          return *this;
        }
        else
        {
          Descale( m_Rescale, m_Data);
        }
      }

      m_Rescale = CreateRescale( m_Data, TO_RANGE, SCALING);
      RescaleImpl( m_Rescale, m_Data);
      return *this;
    }

    //! @brief rescale to the target range, works even if presently scaled to a different range
    template< typename t_DataType>
    FeatureDataSet< t_DataType> &FeatureDataSet< t_DataType>::Rescale( const RescaleFeatureDataSet &TO_RANGE)
    {
      // handle existing scaling
      if( IsRescaled())
      {
        if( *m_Rescale == TO_RANGE)
        {
          // already have the correct rescaling, return
          return *this;
        }
        else
        {
          BCL_MessageStd( "Calling rescale on an already rescaled function!");
          DeScale();
        }
      }

      m_Rescale = util::ShPtr< RescaleFeatureDataSet>( TO_RANGE.Clone());
      RescaleImpl( m_Rescale, m_Data);
      return *this;
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    template< typename t_DataType>
    std::istream &FeatureDataSet< t_DataType>::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Data, ISTREAM);
      io::Serialize::Read( m_Rescale, ISTREAM);
      io::Serialize::Read( m_FeatureLabels, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    template< typename t_DataType>
    std::ostream &FeatureDataSet< t_DataType>::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Data, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Rescale, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_FeatureLabels, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace model
} // namespace bcl
