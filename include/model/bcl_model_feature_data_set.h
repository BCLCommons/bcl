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

#ifndef BCL_MODEL_FEATURE_DATA_SET_H_
#define BCL_MODEL_FEATURE_DATA_SET_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_model_feature_data_set_interface.h"
#include "bcl_model_feature_reference.h"
#include "bcl_model_rescale_feature_data_set.h"
#include "linal/bcl_linal_matrix.h"
#include "linal/bcl_linal_matrix_const_reference.h"
#include "linal/bcl_linal_matrix_reference.h"
#include "type/bcl_type_compare.h"
#include "type/bcl_type_enable_if.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FeatureDataSet
    //! @brief the data set containing features
    //!
    //! @tparam t_DataType can be float, float, int, complex, etc...
    //!
    //! @see @link example_model_feature_data_set.cpp @endlink
    //! @author loweew, woetzen
    //! @date Sep 28, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class FeatureDataSet :
      public FeatureDataSetInterface< t_DataType>
    {

    private:

    //////////
    // data //
    //////////

      //! the linal::matrix data contained in the feature data set
      linal::Matrix< t_DataType> m_Data;

      //! ShPtr to indicate current rescaling; undefined unless currently rescaled
      util::ShPtr< RescaleFeatureDataSet> m_Rescale;

      //! ShPtr to the feature label set, which indicates what each column represents
      util::ShPtr< FeatureLabelSet> m_FeatureLabels;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FeatureDataSet();

      //! @brief constructor from interface
      FeatureDataSet
      (
        const FeatureDataSetInterface< t_DataType> &DATA_SET
      );

      //! @brief constructor from number of examples and feature size
      //! @param EXAMPLES number of distinct entities in the dataset (aka number rows or number features)
      //! @param FEATURE_SIZE size of each feature in the data set (aka number cols or feature size)
      //! @param FILL_VALUE value to fill the matrix with
      FeatureDataSet
      (
        const size_t &EXAMPLES,
        const size_t &FEATURE_SIZE,
        const t_DataType &FILL_VALUE
      );

      //! @brief constructor from matrix
      //! @param MATRIX the matrix to construct from
      FeatureDataSet( const linal::MatrixConstInterface< t_DataType> &MATRIX);

      //! @brief constructor from matrix and existing rescaling
      //! @param MATRIX the matrix to construct from
      //! @param RESCALING existing rescaling
      FeatureDataSet
      (
        const linal::MatrixConstInterface< t_DataType> &MATRIX,
        const RescaleFeatureDataSet &RESCALING
      );

      //! @brief constructor from matrix and set of columns
      //! @param MATRIX the matrix to construct from
      //! @param COLS the colums to copy
      FeatureDataSet
      (
        const linal::MatrixConstInterface< t_DataType> &MATRIX,
        const math::RangeSet< size_t> &COLS
      );

      //! @brief Clone function
      //! @return pointer to new FeatureDataSet< t_DataType>
      FeatureDataSet< t_DataType> *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get size of
      //! @return the number of data items in the feature
      size_t GetFeatureSize() const;

      //! @brief number of features
      //! @return the number of features
      size_t GetNumberFeatures() const;

      //! @brief get a ShPtr to the feature label set, which indicates what each column represents
      //! @return a shptr to the feature label set
      const util::ShPtr< FeatureLabelSet> &GetFeatureLabelSet() const;

      //! @brief set a ShPtr to the feature label set, which indicates what each column represents
      //! @param NEW_SET shptr to the new feature label set
      void SetFeatureLabelSet( const util::ShPtr< FeatureLabelSet> &NEW_SET);

      //! @brief set a ShPtr to the feature label set, which indicates what each column represents
      void SetFeatureLabelSet( const FeatureLabelSet &NEW_SET);

      //! @brief get the owning matrix
      //! @return matrixconstinterface
      linal::MatrixConstReference< t_DataType> GetMatrix() const;

      //! @brief get owning matrix at POS of size LENGTH
      //! @param POS position from which to start the matrix to be returned
      //! @param LENGTH the length in rows of the matrix to be returned
      //! @return the matrixconstinterface from POS of size LENGTH rows
      linal::MatrixConstReference< t_DataType> GetMatrix( const size_t POS, const size_t LENGTH) const;

      //! @brief get a changeable reference to the matrix
      //! @return matrix reference
      linal::MatrixReference< t_DataType> GetMatrix();

      //! @brief get the owning matrix
      //! @return the internally held matrix
      linal::Matrix< t_DataType> &GetRawMatrix()
      {
        return m_Data;
      }

      //! @brief get the owning matrix
      //! @return the interally held matrix
      const linal::Matrix< t_DataType> &GetRawMatrix() const
      {
        return m_Data;
      }

      //! @brief bool indicating whether or not feature has ownership
      //! @return bool - if true, then it has ownership; if it's false, the GetMatrix will return an ownership matrix that's actually a copy
      bool HasOwnership() const
      {
        return true;
      }

      //! @brief test whether this feature data set is rescaled
      //! @return true if the dataset is rescaled
      bool IsRescaled() const;

      //! @brief get a pointer to the rescaling object
      //! @return the rescaling object
      const util::ShPtr< RescaleFeatureDataSet> &GetScaling() const;

      //! @brief test whether this feature data set has the same scaling as another
      //! @param OTHER the other feature dataset
      //! @return true if this feature data set has the same scaling as OTHER
      bool HasSameScaling( const FeatureDataSetInterface< t_DataType> &OTHER) const;

    ////////////////
    // operations //
    ////////////////

      //! @brief scales back to pre-rescaled values
      //! @return FeatureDataSet of scaled values
      FeatureDataSet< t_DataType> &DeScale();

      //! @brief rescale to the target range, works even if presently scaled to a different range
      //! @param SCALING type of scaling to use
      FeatureDataSet< t_DataType> &Rescale
      (
        const math::Range< float> &TO_RANGE,
        const RescaleFeatureDataSet::Type &SCALING = RescaleFeatureDataSet::e_MinMax
      );

      //! @brief scales back to pre-rescaled values
      //! @return FeatureDataSet of scaled values
      FeatureDataSet< t_DataType> &Rescale( const RescaleFeatureDataSet &TO_RANGE);

    ///////////////
    // operators //
    ///////////////

      //! @brief pointer to POS
      //! @param POS position in array of pointer
      //! @return const pointer to POS in range containing all elements of Matrix
      const t_DataType *operator[]( const size_t POS) const
      {
        BCL_Assert( POS < m_Data.GetNumberRows(), "Trying to access out of range of FeatureDataSet");
        return m_Data.Begin() + POS * m_Data.GetNumberCols();
      }

      //! @brief non-owning pointer to position
      //! @param POS position in array of pointer
      //! @return returns a feature reference / non-owning pointer
      FeatureReference< t_DataType> operator()( const size_t POS) const
      {
        return FeatureReference< t_DataType>( m_Data.GetNumberCols(), operator[]( POS));
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // template class FeatureDataSet

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API FeatureDataSet< float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API FeatureDataSet< char>;

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_FEATURE_DATA_SET_H_
