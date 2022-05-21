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

#ifndef BCL_MODEL_KAPPA_NEAREST_NEIGHBOR_H_
#define BCL_MODEL_KAPPA_NEAREST_NEIGHBOR_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_model_feature_data_set.h"
#include "bcl_model_interface.h"
#include "descriptor/bcl_descriptor_dataset.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class KappaNearestNeighbor
    //! @brief performs the kappa nearest neighbor unsupervised training algorithm
    //!        (lazy training algorithm).
    //!
    //! @details reference: Zheng W and Tropsha A (2000) Novel variable selection quantitative structure–property
    //!                     relationship approach based on the k-nearest-neighbor principle.
    //!                     J Chem Inf Comput Sci 40: 185–194
    //!
    //! @see @link example_model_kappa_nearest_neighbor.cpp @endlink
    //! @author loweew
    //! @date 03/02/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API KappaNearestNeighbor :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      //! data to train the kNN on
      util::ShPtr< descriptor::Dataset> m_TrainingData;

      //! kappa value indicating number of nearest neighbors to compare
      size_t m_Kappa;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////
    // data //
    //////////

      //! the default input range
      static const math::Range< float> s_DefaultInputRange;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      KappaNearestNeighbor();

      //! @brief constructor from training data, kappa, and rescale
      //! @param TRAINING_DATA data to train the k-NN on
      //! @param KAPPA number of nearest neighbors to consider
      KappaNearestNeighbor
      (
        const util::ShPtr< descriptor::Dataset> &TRAINING_DATA,
        const size_t KAPPA
      );

      //! @brief copy constructor
      KappaNearestNeighbor *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief get the output feature size for this model
      //! @return the output feature size for this model
      size_t GetNumberOutputs() const
      {
        return m_TrainingData.IsDefined() ? m_TrainingData->GetResultSize() : size_t( 0);
      }

      //! @brief sets training data
      //! @param TRAINING_DATA data to be classified
      void SetTrainingData
      (
        util::ShPtr< descriptor::Dataset> &TRAINING_DATA
      )
      {
        m_TrainingData = TRAINING_DATA;
      }

      //! @brief gets training data
      //! @return training data as descriptor::Dataset
      const util::ShPtr< descriptor::Dataset> &GetTrainingData()
      {
        return m_TrainingData;
      }

      //! @brief sets kappa
      //! @param KAPPA number of neighbors to consider
      void SetKappa( size_t KAPPA)
      {
        m_Kappa = KAPPA;
      }

      //! @brief gets kappa value
      //! @return size_t of kappa which is the number of nearest neighbors to consider
      size_t GetKappa() const
      {
        return m_Kappa;
      }

      //! @brief Set the scaling of a feature set according to the model
      //! @param FEATURES feature set of interest
      //! @note this allows for external classes that own a dataset to ensure that a new dataset is never created
      //!       when operator() is called
      void Rescale( FeatureDataSet< float> &FEATURE) const;

      //! @brief performs the kappa nearest neighbor algorithm
      //! @param QUERY_VECTOR math::Vector< float> of features
      //!        and known output of query point
      //! @param RESULT_STORAGE pointer to vector or matrix where result should be stored; must be preallocated
      static storage::Vector< storage::Pair< float, size_t> > FindWithoutRescaling
      (
        const FeatureReference< float> &QUERY_VECTOR,
        const size_t KAPPA,
        const FeatureDataSet< float> &TRAINING_FEATURES
      );

      //! @brief performs the kappa nearest neighbor algorithm
      //! @param QUERY_VECTORS math::Vector< float> of features
      //!        and known output of query point
      //! @param RESULT_STORAGE pointer to vector or matrix where result should be stored; must be preallocated
      static storage::Vector< storage::Triplet< float, size_t, size_t> > FindWithoutRescaling
      (
        const FeatureDataSet< float> &QUERY_VECTORS,
        const size_t KAPPA,
        const FeatureDataSet< float> &TRAINING_FEATURES
      );

    //////////////
    // operator //
    //////////////

      //! @brief predict result with model using a NOT rescaled feature vector
      //! @param FEATURE not rescaled feature vector
      //! @return predicted result vector using a model
      FeatureDataSet< float> PredictWithoutRescaling( const FeatureDataSetInterface< float> &FEATURE) const;

      //! @brief predict result with model using a rescaled feature vector
      //! @param FEATURE normalized or rescaled feature vector
      //! @return predicted result vector using a model
      FeatureDataSet< float> operator()( const FeatureDataSetInterface< float> &FEATURE) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! read KappaNearestNeighbor from std::istream
      std::istream &Read( std::istream &ISTREAM);

      //! write KappaNearestNeighbor into std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    private:

      //! @brief performs the kappa nearest neighbor algorithm
      //! @param QUERY_VECTOR FeatureReference< float> of features
      //! @param RESULT_STORAGE pointer to vector or matrix where result should be stored; must be preallocated
      void EuclideanDistance
      (
        const FeatureReference< float> &QUERY_VECTOR,
        float *RESULT_STORAGE
      ) const;

    }; // class KappaNearestNeighbor

  } // namespace model
} // namespace bcl

#endif //BCL_MODEL_KAPPA_NEAREST_NEIGHBOR_H_
