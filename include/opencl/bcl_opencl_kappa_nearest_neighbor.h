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

#ifndef BCL_OPENCL_KAPPA_NEAREST_NEIGHBOR_H_
#define BCL_OPENCL_KAPPA_NEAREST_NEIGHBOR_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_opencl_euclidean_distance.h"
#include "bcl_opencl_insertion_sort.h"
#include "bcl_opencl_matrix.h"
#include "bcl_opencl_model_interface.h"
#include "descriptor/bcl_descriptor_dataset.h"
#include "model/bcl_model_feature_data_set.h"
#include "model/bcl_model_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class KappaNearestNeighbor
    //! @brief performs the kappa nearest neighbor unsupervised training algorithm
    //!        (lazy training algorithm) - optimized for the gpu
    //!
    //! @details reference: Zheng W and Tropsha A (2000) Novel variable selection quantitative structure–property
    //!                     relationship approach based on the k-nearest-neighbor principle.
    //!                     J Chem Inf Comput Sci 40: 185–194
    //!
    //! @see @link example_opencl_kappa_nearest_neighbor.cpp @endlink
    //! @author loweew
    //! @date Apr 12, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API KappaNearestNeighbor :
      public ModelInterface
    {

    private:

    //////////
    // data //
    //////////

      //! data to train the kNN on
      util::ShPtr< descriptor::Dataset> m_ReferenceData;

      //! ref features on device
      Matrix< float> m_ReferenceFeatures;

      //! ref results on device
      Matrix< float> m_ReferenceResults;

      //! query features on device
      Matrix< float> m_QueryFeatures;

      //! kappa value indicating number of nearest neighbors to compare
      size_t m_Kappa;

      //! queue
      CommandQueue m_Queue;

      //! program
      cl::Program m_Program;

      //! euclidean distance
      EuclideanDistance< float> m_GpuDist;

      //! insertion sort
      InsertionSort    m_GpuSort;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////
    // data //
    //////////

      //! the default input range
      static const math::Range< float> s_DefaultInputRange;

      //! the default output range
      static const math::Range< float> s_DefaultOutputRange;
      static const char *s_CLCompilerOptions;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      KappaNearestNeighbor();

      //! @brief default constructor
      KappaNearestNeighbor( const CommandQueue &QUEUE);

      //! @brief constructor from training data, kappa, and rescale
      //! @param REFERENCE_DATA data to train the k-NN on
      //! @param KAPPA number of nearest neighbors to consider
      KappaNearestNeighbor
      (
        util::ShPtr< descriptor::Dataset> &REFERENCE_DATA,
        const size_t KAPPA,
        const CommandQueue &QUEUE
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

      //! @brief sets training data
      //! @param REFERENCE_DATA data to be classified
      void SetReferenceData
      (
        util::ShPtr< descriptor::Dataset> &REFERENCE_DATA
      )
      {
        m_ReferenceData = REFERENCE_DATA;
      }

      //! @brief gets training data
      //! @return training data as descriptor::Dataset
      util::ShPtr< descriptor::Dataset> &GetReferenceData()
      {
        return m_ReferenceData;
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

      //! @brief get the output feature size for this model
      //! @return the output feature size for this model
      size_t GetNumberOutputs() const
      {
        return m_ReferenceResults.GetNumberCols() - m_ReferenceResults.GetColPadding();
      }

    //////////////
    // operator //
    //////////////

      //! @brief Set the scaling of a feature set according to the model
      //! @param FEATURES feature set of interest
      //! @note this allows for external classes that own a dataset to ensure that a new dataset is never created
      //!       when operator() is called
      void Rescale( model::FeatureDataSet< float> &FEATURE) const;

      //! @brief predict result with model using a NOT rescaled feature vector
      //! @param FEATURE not rescaled feature vector
      //! @return predicted result vector using a model
      model::FeatureDataSet< float> PredictWithoutRescaling( const model::FeatureDataSetInterface< float> &FEATURE) const;

      //! @brief predict result with model using a rescaled feature vector
      //! @param FEATURE normalized or rescaled feature vector
      //! @return predicted result vector using a model
      model::FeatureDataSet< float> operator()( const model::FeatureDataSetInterface< float> &FEATURE) const;

      //! @brief predict result with model using a rescaled feature vector
      //! @param FEATURE normalized or rescaled Matrix
      //! @return predicted result vector using a model
      Matrix< float> operator()( const Matrix< float> &FEATURE) const;

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

      //! @brief sets the reference data on the device
      void SetReferenceDataOnDevice();

      //! @brief calculates the weighted results based on distance
      //! @param SORTED_DISTANCES the kappa closest distances
      //! @param SORTED_INDECES the index corresponding to that distance for the result
      //! @return the output matrix of results
      Matrix< float> CalculateWeightedResults
      (
        const Matrix< float> &SORTED_DISTANCES,
        const Matrix< int> &SORTED_INDECES
      ) const;

      //! @brief compile programs for given precision
      //! @param PRECISION float or double
      //! @return ERROR error that occured, CL_SUCCESS if no error
      cl_int CompilePrograms( const util::CPPDataTypes::Types &PRECISION);

    }; // class KappaNearestNeighbor

  } // namespace opencl
} // namespace bcl

#endif //BCL_OPENCL_KAPPA_NEAREST_NEIGHBOR_H_
