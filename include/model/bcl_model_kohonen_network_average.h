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

#ifndef BCL_MODEL_KOHONEN_NETWORK_AVERAGE_H_
#define BCL_MODEL_KOHONEN_NETWORK_AVERAGE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_feature_data_set.h"
#include "bcl_model_interface.h"
#include "bcl_model_kohonen_node.h"
#include "descriptor/bcl_descriptor_dataset.h"
#include "linal/bcl_linal_vector.h"
#include "sched/bcl_sched_mutex.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class KohonenNetworkAverage
    //! @brief takes the average output of the training data won by each node.
    //!
    //! @see @link example_model_kohonen_network_average.cpp @endlink
    //! @author lemmonwa, mendenjl
    //! @date 7/19/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API KohonenNetworkAverage :
      public Interface
    {

    //////////
    // data //
    //////////

    public:

      //! add self framework enums
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    protected:

      //! the actual map, consisting of reference vectors
      storage::Vector< KohonenNode> m_Codebook;

      //! the dimensions of the map
      linal::Vector< double> m_MapDimensions;

      //! rescale input to transfer function input range
      util::ShPtr< RescaleFeatureDataSet> m_RescaleInput;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief default constructor
      KohonenNetworkAverage()
      {
      }

      //! @brief constructor from parameters
      //! @param MAP_DIMENSIONS the number and dimensions of the map. Ex. < x=40, y=20 >
      //! @param RESCALE method to rescale features and results
      //! @param NORMALIZED_DATA data that can be used to set the initial positions of the nodes
      //!        defaults to empty, so all nodes have identical feature vectors of 0
      KohonenNetworkAverage
      (
        const linal::Vector< double> &MAP_DIMENSIONS,
        const util::ShPtr< RescaleFeatureDataSet> &RESCALE = util::ShPtr< RescaleFeatureDataSet>()
      ) :
        m_MapDimensions( MAP_DIMENSIONS),
        m_RescaleInput( RESCALE)
      {
      }

      //! @brief Clone function
      //! @return pointer to new KohonenNetworkAverage
      KohonenNetworkAverage *Clone() const
      {
        return new KohonenNetworkAverage( *this);
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

      //! @brief get the output feature size for this model
      //! @return the output feature size for this model
      size_t GetNumberOutputs() const
      {
        return !m_Codebook.IsEmpty() ? m_Codebook.FirstElement().GetResultVector().GetSize() : size_t( 0);
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief returns map dimensions
      //! @return the dimensionality of the kohonen network, i.e. 40 X 20
      const linal::Vector< double> &GetMapDimensions() const
      {
        return m_MapDimensions;
      }

      //! @brief set the nodes up according to the map dimensions
      //! @param NORMALIZED_DATA data that can be used to set the initial positions of the nodes
      void InitializeNodes( const descriptor::Dataset &NORMALIZED_DATA);

      //! @brief reset sets the weight of all nodes to 0
      void Reset();

      //! @brief gets the code book
      //! @return return constant reference to map distance function
      const storage::Vector< KohonenNode> &GetCodeBook() const
      {
        return m_Codebook;
      }

      //! @brief gets the code book
      //! @return return constant reference to map distance function
      storage::Vector< KohonenNode> &GetCodeBook()
      {
        return m_Codebook;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief determines the winning node indices of a normalized data set
      //! @param NORMALIZED_DATA all the training data, already normalized
      //! @param RANGE the range of training data to determine the winning indices for
      //! @param ORDER order in which to visit
      //! @param PREVIOUS_WINNERS previous rounds node winners, will be updated
      //! @note this function is threaded for performance
      void GetWinningNodeIndices
      (
        const FeatureDataSetInterface< float> &NORMALIZED_DATA,
        const math::Range< size_t> &RANGE,
        const storage::Vector< size_t> &ORDER,
        storage::Vector< size_t> &PREVIOUS_WINNERS,
        storage::Vector< sched::Mutex> &MUTECES
      ) const;

      //! @brief gets the winning node
      //! @param VECTOR the sample of data from the training data, already normalized
      //! @param HINT a hint for the index of the winning node; doubles speed on average if set to previous winner
      //! @return the winning node, paired with the error
      size_t GetIndexOfWinningNode
      (
        const linal::VectorConstInterface< float> &VECTOR,
        const size_t &HINT = size_t( 0)
      ) const;

      //! @brief Set the scaling of a feature set according to the model
      //! @param FEATURES feature set of interest
      //! @note this allows for external classes that own a dataset to ensure that a new dataset is never created
      //!       when operator() is called
      void Rescale( FeatureDataSet< float> &FEATURE) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief predict result with model using a NOT rescaled feature vector
      //! @param FEATURE not rescaled feature vector
      //! @return predicted result vector using a model
      FeatureDataSet< float> PredictWithoutRescaling
      (
        const FeatureDataSetInterface< float> &FEATURE
      ) const;

      //! @brief predict result with model using a rescaled feature vector
      //! @param FEATURE normalized or rescaled feature vector
      //! @return predicted result vector using a model
      FeatureDataSet< float> operator()( const FeatureDataSetInterface< float> &FEATURE) const;

      //! @brief adds each node from the other network to this network. The other network is probably a partial network.
      //! @param OTHER another Kohonen network to add from
      //! @return the resulting network
      KohonenNetworkAverage &operator +=( const KohonenNetworkAverage &OTHER);

      //! @brief adds each node from the other network to this network with a given weight / alpha
      //! @param OTHER another Kohonen network to add from
      //! @param WEIGHT the weight to give each node from OTHER
      //! @return the resulting network
      KohonenNetworkAverage &AddNetworkWithWeight( const KohonenNetworkAverage &OTHER, const float &WEIGHT);

      //! @brief give centroid values to nodes with no associated nearest features
      KohonenNetworkAverage &FixEmptyNodes();

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief thread worker function for finding the winning indices for a data set with the network
      //! @param DATA the data to operate on
      //! @param RANGE the designated range of the training data to work on
      //! @param ORDER the order in which to visit the nodes
      //! @param WINNING_NODES_AND_MUTECES reference to vector where winning indices are stored, along with mutexes (one per winning node)
      void GetWinningNodeIndicesHelper
      (
        const FeatureDataSetInterface< float> &DATA,
        const storage::Pair< math::Range< size_t>, util::SiPtr< const storage::Vector< size_t> > > &RANGE_AND_ORDER,
        storage::Pair
        <
          util::SiPtr< storage::Vector< size_t> >,
          util::SiPtr< storage::Vector< sched::Mutex> >
        > &WINNING_NODES_AND_MUTECES
      ) const;

    }; // class KohonenNetworkAverage

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_KOHONEN_NETWORK_AVERAGE_H_
