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

#ifndef BCL_MODEL_FEATURE_DATA_REFERENCE_H_
#define BCL_MODEL_FEATURE_DATA_REFERENCE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_feature_data_set_interface.h"
#include "bcl_model_feature_reference.h"
#include "linal/bcl_linal_matrix.h"
#include "linal/bcl_linal_matrix_const_interface.h"
#include "linal/bcl_linal_matrix_const_reference.h"
#include "linal/bcl_linal_matrix_operations.h"
#include "linal/bcl_linal_vector.h"
#include "math/bcl_math_range_set.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FeatureDataReference
    //! @brief the data set containing features
    //!
    //! @tparam t_DataType can be float, float, int, complex, etc...
    //!
    //! @see @link example_model_feature_data_reference.cpp @endlink
    //! @author mendenjl
    //! @date Apr 06, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class FeatureDataReference :
      public FeatureDataSetInterface< t_DataType>
    {

    private:

    //////////
    // data //
    //////////

      //! the linal::matrix data contained in the feature data set
      linal::MatrixConstReference< t_DataType> m_Data;

      //! ShPtr to the feature label set, which indicates what each column represents
      util::ShPtr< FeatureLabelSet> m_FeatureLabels;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FeatureDataReference() :
        m_Data()
      {
      }

      //! @brief constructor from matrix
      //! @param MATRIX the matrix to construct from
      FeatureDataReference
      (
        const linal::MatrixConstInterface< t_DataType> &MATRIX,
        const util::ShPtr< FeatureLabelSet> &FEATURE_LABEL_SET =
          util::ShPtr< FeatureLabelSet>()
      ) :
        m_Data( MATRIX),
        m_FeatureLabels( FEATURE_LABEL_SET)
      {
      }

      //! @brief Clone function
      //! @return pointer to new FeatureDataReference< t_DataType>
      FeatureDataReference< t_DataType> *Clone() const
      {
        return new FeatureDataReference< t_DataType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get size of
      //! @return the number of data items in the feature
      size_t GetFeatureSize() const
      {
        return m_Data.GetNumberCols();
      }

      //! @brief number of features
      //! @return the number of features
      size_t GetNumberFeatures() const
      {
        return m_Data.GetNumberRows();
      }

      //! @brief get the owning matrix
      //! @return matrixconstinterface
      linal::MatrixConstReference< t_DataType> GetMatrix() const
      {
        return m_Data;
      }

      //! @brief get owning matrix at POS of size LENGTH
      //! @param POS position from which to start the matrix to be returned
      //! @param LENGTH the length in rows of the matrix to be returned
      //! @return the matrixconstinterface from POS of size LENGTH rows
      linal::MatrixConstReference< t_DataType> GetMatrix( const size_t POS, const size_t LENGTH) const
      {
        return linal::MatrixConstReference< t_DataType>( LENGTH, m_Data.GetNumberCols(), operator[]( POS));
      }

      //! @brief bool indicating whether or not feature has ownership
      //! @return bool - if true, then it has ownership; if it's false, the GetMatrix will return an ownership matrix that's actually a copy
      bool HasOwnership() const
      {
        return false;
      }

      //! @brief get a pointer to the rescaling object
      //! @return the rescaling object
      const util::ShPtr< RescaleFeatureDataSet> &GetScaling() const
      {
        static util::ShPtr< RescaleFeatureDataSet> s_scaling;
        return s_scaling;
      }

      //! @brief get a ShPtr to the feature label set, which indicates what each column represents
      //! @return a shptr to the feature label set
      const util::ShPtr< FeatureLabelSet> &GetFeatureLabelSet() const
      {
        return m_FeatureLabels;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief pointer to POS
      //! @param POS position in array of pointer
      //! @return const pointer to POS in range containing all elements of Matrix
      const t_DataType *operator[]( const size_t POS) const
      {
        BCL_Assert( POS < m_Data.GetNumberRows(), "Trying to access out of range of FeatureDataReference");
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
      std::istream &Read( std::istream &ISTREAM)
      {
        // read members
        io::Serialize::Read( m_Data, ISTREAM);

        // return the stream
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_Data, OSTREAM, INDENT);

        // return the stream
        return OSTREAM;
      }

    }; // template class FeatureDataReference

    // instantiate s_Instance
    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> FeatureDataReference< t_DataType>::s_Instance
    (
      GetObjectInstances().AddInstance( new FeatureDataReference< t_DataType>())
    );

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_FEATURE_DATA_REFERENCE_H_
