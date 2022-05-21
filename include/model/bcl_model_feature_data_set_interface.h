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

#ifndef BCL_MODEL_FEATURE_DATA_SET_INTERFACE_H_
#define BCL_MODEL_FEATURE_DATA_SET_INTERFACE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_model_feature_label_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FeatureDataSetInterface
    //! @brief this interface is for feature data sets
    //!
    //! @remarks example unnecessary
    //! @author loweew, woetzen, butkiem1, mendenjl
    //! @date 09/15/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class FeatureDataSetInterface :
      public util::ObjectInterface
    {

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief get size of
      //! @return the number of data items in the feature
      virtual size_t GetFeatureSize() const = 0;

      //! @brief number of features
      //! @return the number of features
      virtual size_t GetNumberFeatures() const = 0;

      //! @brief get the non-owning matrix
      //! @return matrixconstinterface
      virtual linal::MatrixConstReference< t_DataType> GetMatrix() const = 0;

      //! @brief get non-owning matrix at POS of size LENGTH
      //! @param POS position from which to start the matrix to be returned
      //! @param LENGTH the length in rows of the matrix to be returned
      //! @return the matrixconstinterface from POS of size LENGTH rows
      virtual linal::MatrixConstReference< t_DataType> GetMatrix( const size_t POS, const size_t LENGTH) const = 0;

      //! @brief get a ShPtr to the feature label set, which indicates what each column represents
      //! @return a shptr to the feature label set
      virtual const util::ShPtr< FeatureLabelSet> &GetFeatureLabelSet() const = 0;

      //! @brief bool indicating whether or not feature has ownership
      //! @return bool - if true, then it has ownership; if it's false, the GetMatrix will return an ownership matrix that's actually a copy
      virtual bool HasOwnership() const = 0;

      //! @brief test whether this feature data set is rescaled
      //! @return true if the dataset is rescaled
      virtual bool IsRescaled() const
      {
        return false;
      }

      //! @brief get a pointer to the rescaling object
      //! @return the rescaling object
      virtual const util::ShPtr< RescaleFeatureDataSet> &GetScaling() const = 0;

    ///////////////
    // operators //
    ///////////////

      //! @brief pointer to POS
      //! @param POS position in array of pointer
      //! @return const pointer to POS in range containing all elements of Matrix
      virtual const t_DataType *operator[]( const size_t POS) const = 0;

      //! @brief non-owning pointer to position
      //! @param POS position in array of pointer
      //! @return returns a feature reference / non-owning pointer
      virtual FeatureReference< t_DataType> operator()( const size_t POS) const = 0;

    }; // template FeatureDataSetInterface

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_FEATURE_DATA_SET_INTERFACE_H_
