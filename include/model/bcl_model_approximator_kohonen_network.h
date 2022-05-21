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

#ifndef BCL_MODEL_APPROXIMATOR_KOHONEN_NETWORK_H_
#define BCL_MODEL_APPROXIMATOR_KOHONEN_NETWORK_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_approximator_base.h"
#include "bcl_model_kohonen_network_average.h"
#include "bcl_model_objective_function_wrapper.h"
#include "bcl_model_training_schedule.h"
#include "math/bcl_math_range.h"
#include "sched/bcl_sched_mutex.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ApproximatorKohonenNetwork
    //! @brief trains a kohonen network after the manner of the Helsinky code, SOM_PAK,
    //!
    //! @see @link example_model_approximator_kohonen_network.cpp @endlink
    //! @author lemmonwa, mendenjl
    //! @date Jul 19, 2013; Sep 03, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ApproximatorKohonenNetwork :
      public ApproximatorBase
    {

    //////////
    // data //
    //////////

    public:

      //! static instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    ///////////
    // enums //
    ///////////

      //! all possible neighbor adapt functions
      enum NeighborKernel
      {
        e_Gaussian,
        e_Bubble,
        s_NumberNeighborhoodKernels
      };

      //! @brief Kernel as string
      //! @param KERNEL the kernel
      //! @return the string for the kernel
      static const std::string &GetKernelName( const NeighborKernel &KERNEL);

      //! @brief NeighborKernelEnum enum I/O helper
      typedef util::WrapperEnum< NeighborKernel, &GetKernelName, s_NumberNeighborhoodKernels> NeighborKernelEnum;

      //! all possible neighbor adapt functions
      enum Initialization
      {
        e_FirstVectors,           //!< Initial vectors taken from the start of the training data
        e_RandomlyChosenVectors,  //!< Initial vectors randomly chosen from the training data
        e_RandomlyChosenElements, //!< Initial vectors composed of randomly chosen values for each element of the vector
        e_KMeans,                 //!< Initial vectors composed of samples from KMeans-based clustering
        s_NumberInitializers
      };

      //! @brief Initializer as string
      //! @param INITIALIZER the initialization
      //! @return the string for the initializer
      static const std::string &GetInitializationName( const Initialization &INITIALIZER);

      //! @brief Initialization enum I/O helper
      typedef util::WrapperEnum< Initialization, &GetInitializationName, s_NumberInitializers> InitializationEnum;

    protected:

      //! the actual map, consisting of reference vectors
      KohonenNetworkAverage m_Network;

      //! length of training
      size_t m_Length;

      //! index of current iteration
      size_t m_CurrentIteration;

      //! radius (for SOM), max distance to neighbors
      float m_Radius;

      //! current radius
      float m_CurrentRadius;

      //! # of steps between node updates
      //! 1 leads to sequential iteration (slow, can't thread, but most accurate)
      //! 0 -> dataset size gives batch iteration (fast, can thread, least accurate)
      //! intermediate values ( moderate values, e.g. 256, threadable, typically best tradoff)
      size_t m_UpdateEveryNthFeature;

      //! the neighborhood kernel to use
      NeighborKernelEnum m_NeighborKernel;

      //! initialization method
      InitializationEnum m_Initializer;

      //! storage of last mapped node for each input data point
      //! This is used to speed up the search for the closest node
      storage::Vector< size_t> m_LastClosestNodes;

      //! bool - whether to always call Track in the tracker. Derived classes often do not want this feature
      bool m_NoTrack;

      //! scaling type
      RescaleFeatureDataSet::TypeEnum m_RescaleType;

      //! Training schedule; used to allow for balancing and shuffling
      TrainingSchedule m_Schedule;

      //! Cutoff; used for balancing
      float m_Cutoff;

      //! Thread-related variables

      //! number of threads for distributing nodes
      size_t m_NumberThreadsNodes;

      //! muteces; one for each data point
      storage::Vector< sched::Mutex> m_Muteces;

      //! Thread-local variables

      //! ranges of the data for each thread to operate over
      storage::Vector< math::Range< size_t> > m_DataSetRanges;

      //! ranges of nodes in the network for each thread to operate over when accumulating results and performing the
      //! neighbor adapt
      storage::Vector< math::Range< size_t> > m_NodeRanges;

      //! new network to replace the old one during batch iterate
      storage::Vector< KohonenNetworkAverage> m_NewNetworks;
      storage::Vector< KohonenNetworkAverage> m_NewNetworksSwap;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ApproximatorKohonenNetwork();

      //! @brief constructor from all necessary parameters
      //! @param MAP_DIMENSIONS dimensions of the map
      //! @param INITIAL_LENGTH how many iterations to train for (if applicable).
      //! @param INITAL_RADIUS the initial neighborhood radius
      //! @param OBJECTIVE_FUNCTION the objective function from the approximator framework
      //! @param UPDATE_EVERY_NTH_FEATURE update the nodes after seeing this many features
      //! @param NEIGHBOR_KERNEL the neighbor kernel type
      ApproximatorKohonenNetwork
      (
        const linal::Vector< double> &MAP_DIMENSIONS,
        const size_t &INITIAL_LENGTH,
        const float &INITAL_RADIUS,
        const util::ShPtr< ObjectiveFunctionWrapper> &OBJECTIVE_FUNCTION,
        const size_t &UPDATE_EVERY_NTH_FEATURE,
        const NeighborKernel &NEIGHBOR_KERNEL
      );

      //! @brief clone function
      //! @return pointer to new ApproximatorKohonenNetwork
      ApproximatorKohonenNetwork *Clone() const
      {
        return new ApproximatorKohonenNetwork( *this);
      }

    /////////////////
    // data access //
    /////////////////

    public:

      //! @brief get training data set for a specific iterate in approximator framework
      //! @return dataset interface with training data
      const util::ShPtr< descriptor::Dataset> &
      GetTrainingData() const
      {
        return m_TrainingData;
      }

      //! @brief get neighbor adapt function being used
      //! @return a reference to the neighbor  adapt function
      const NeighborKernel &GetNeighborKernel() const
      {
        return m_NeighborKernel;
      }

      //! @brief set training data set for a specific iterate in approximator framework
      //! @param DATA training data set
      void SetTrainingData
      (
        util::ShPtr< descriptor::Dataset> &DATA
      );

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief construct a model from the current iterate
      //! @return shptr to the new model interface
      util::ShPtr< Interface> GetCurrentModel() const;

      //! @brief returns the current approximation
      //! @return current argument result pair
      const util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> > GetCurrentApproximation() const;

      //! @brief conducts the next approximation step and stores the approximation
      void Next();

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

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief evaluates whether the approximation can continue
      //! @return true, if the approximation can continue - otherwise false
      bool CanContinue() const
      {
        return true;
      }

      //! @brief thread helper function for iterate, works on a subset of the training data, computes partial
      //!        adjustments to the network, then adds them to the final result.
      //! @param RANGE_ID the partition of the training data to work on
      //! @param START_DATA_NUMBER the index in the training data for the current epoch
      void TrainThread
      (
        const size_t &RANGE_ID,
        const size_t &START_DATA_NUMBER,
        const size_t &END_DATA_NUMBER
      );

      //! @brief setup the data ranges and node ranges for training
      //! This should only be called while setting the training data or after reading
      void SetupThreadRanges();

      //! @brief spreads the adaptation from a particular node to other nodes within RADIUS
      //! @param NETWORK the network to adapt, usually a copy to avoid conflicting changes
      //! @param NODE_TO_SPREAD is the node whose feature and result vectors should be propagated
      //! @param RADIUS the current neighborhood radius
      void AdaptNeighbors
      (
        KohonenNetworkAverage &NETWORK,
        const KohonenNode &NODE_TO_SPREAD,
        const float &RADIUS
      ) const;

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class ApproximatorKohonenNetwork

  } // namespace model

} // namespace bcl

#endif // BCL_MODEL_APPROXIMATOR_KOHONEN_NETWORK_H_

