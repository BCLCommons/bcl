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

#ifndef BCL_MODEL_TRAINING_SCHEDULE_H_
#define BCL_MODEL_TRAINING_SCHEDULE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically
#include "descriptor/bcl_descriptor.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_model_objective_function_interface.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class TrainingSchedule
    //! @brief Storage class that stores the dataset labels used for monitoring and independently assessing a
    //!        model and the final objective function result and parity, and generically all meta-data associated with
    //!        training of the model, except the feature and result descriptors.
    //!
    //!
    //! @see @link example_model_training_schedule.cpp @endlink
    //! @author mendenjl
    //! @date Sep 23, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API TrainingSchedule :
      public util::SerializableInterface
    {
    private:

    //////////
    // data //
    //////////

      //! Bool: whether to balance the data for classification tasks
      bool m_Balance;

      //! true to shuffle feature / result ordering after every run through the data (iteration)
      bool m_Shuffle;

      //! For balancing: max number of repeated features
      size_t m_BalanceMaxRepeatedFeatures;

      //! For balancing: this is the ratio the rarer class (positives/negatives) will be oversampled to
      float m_BalanceMaxOversampling;

      //! result class for each feature.  Used only with peer-based varieties of input dropout
      storage::Vector< size_t> m_ResultClass;

      //! order in which to visit the features (modified from consecutive if m_Shuffle is on)
      storage::Vector< size_t> m_Order;

      //! features that have the same result class
      storage::Vector< storage::Vector< size_t> > m_PeerFeatures;

      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from members
      TrainingSchedule
      (
        const bool &SHUFFLE = false,
        const bool &BALANCE = false,
        const size_t &MAX_REPEATS = std::numeric_limits< size_t>::max(),
        const float &TARGET_OVERSAMPLING = 1.0
      );

      //! @brief Clone function
      //! @return pointer to new TrainingSchedule
      TrainingSchedule *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief get whether this class will balance the order vector
      bool IsBalanced() const
      {
        return m_Balance;
      }

      //! @brief get whether this class will balance the order vector
      bool IsShuffled() const
      {
        return m_Shuffle;
      }

      //! @brief get maximum number of times a feature will be repeated
      size_t GetMaxRepeats() const
      {
        return m_BalanceMaxRepeatedFeatures;
      }

      //! @brief get maximum/target rate (0-1) of balancing relative to the dominant class
      float GetMaxBalancingRate() const
      {
        return m_BalanceMaxOversampling;
      }

      //! @brief get the number of unique result classes
      size_t GetNumberOfClasses() const
      {
        return m_PeerFeatures.GetSize();
      }

      //! @brief get the result class for every feature
      const storage::Vector< size_t> &GetClasses() const
      {
        return m_ResultClass;
      }

      //! @brief get member ids for each class
      const storage::Vector< storage::Vector< size_t> > &GetClassMembers() const
      {
        return m_PeerFeatures;
      }

      //! @brief get the size of the overbalanced array
      size_t GetSize() const
      {
        return m_Order.GetSize();
      }

      //! @brief get the order vector
      const storage::Vector< size_t> &GetOrder() const
      {
        return m_Order;
      }

      //! @brief retrieve the index of the feature to visit ORDER_INDEX-th
      size_t operator()( const size_t &ORDER_INDEX) const
      {
        return m_Order( ORDER_INDEX);
      }

      //! @brief Perform any end-of-iteration operations necessary (e.g. shuffling)
      void Next()
      {
        if( m_Shuffle)
        {
          m_Order.Shuffle();
        }
      }

      //! @brief Setup the balancer
      //! @param RESULTS results feature data set
      //! @param CUTOFF cutoff from the objective function that separates classes
      void Setup
      (
        const FeatureDataSetInterface< float> &RESULTS,
        const float &CUTOFF
      );

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class TrainingSchedule

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_TRAINING_SCHEDULE_H_
