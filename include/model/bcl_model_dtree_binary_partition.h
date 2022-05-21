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

#ifndef BCL_MODEL_DTREE_BINARY_PARTITION_H_
#define BCL_MODEL_DTREE_BINARY_PARTITION_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_feature_result_and_state.h"
#include "linal/bcl_linal_vector.h"
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DtreeBinaryPartition
    //! @brief data class containing information about a partition of the data
    //!
    //! @see @link example_model_dtree_binary_partition.cpp @endlink
    //! @author mendenjl
    //! @date May 09, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DtreeBinaryPartition :
      public util::ObjectInterface
    {

    public:

    //////////
    // enum //
    //////////

      // data retrievable from this object via GetData
      enum Data
      {
        e_SplitRating,
        e_InitialNumIncorrect,
        e_RatingTimesInitialNumIncorrect,
        e_InitialIncorrectPlusFinalCorrect,
        s_NumberData
      };

      //! @brief Data as string
      //! @param DATA the data of interest
      //! @return the string for the data
      static const std::string &GetDataName( const Data &DATA);

      //! @brief NeighborKernelEnum enum I/O helper
      typedef util::WrapperEnum< Data, &GetDataName, s_NumberData> DataEnum;

    private:

    //////////
    // data //
    //////////

      size_t m_FeatureIndex; //!< feature index on which to split
      float  m_SplitValue;   //!< value at which to split for the feature at that index
      float  m_SplitRating;  //!< rating of the split; higher values imply better splits
      size_t m_Size;         //!< total # features affected by this partition
      size_t m_InitialIncorrect; //!< number of nodes incorrectly classified before using this partition
      size_t m_FinalIncorrect;   //!< total # of features incorrect after classification with this partition

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      DtreeBinaryPartition();

      //! @brief constructor from members
      DtreeBinaryPartition
      (
        const size_t &FEATURE_INDEX,
        const float  &SPLIT_VALUE,
        const float  &SPLIT_RATING,
        const size_t &SIZE = 0,
        const size_t &NUM_INCORRECT = 0
      );

      //! @brief Clone function
      //! @return pointer to new DtreeBinaryPartition
      DtreeBinaryPartition *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the features
      //! @return the features
      const size_t &GetFeatureIndex() const
      {
        return m_FeatureIndex;
      }

      //! @brief get the split value
      //! @return the split value
      const float &GetSplitValue() const
      {
        return m_SplitValue;
      }

      //! @brief get the result state
      //! @return the result state
      const float &GetSplitRating() const
      {
        return m_SplitRating;
      }

      //! @brief get the total # of features affected by the partition
      //! @return the total # of features affected by the partition
      const size_t &GetSize() const
      {
        return m_Size;
      }

      //! @brief get the # of incorrect features prior to the partition
      //! @return the # of incorrect features prior to the partition
      const size_t &GetInitialNumIncorrect() const
      {
        return m_InitialIncorrect;
      }

      //! @brief get the # of incorrect features after the partition
      //! @return the # of incorrect features after the partition
      const size_t &GetFinalNumIncorrect() const
      {
        return m_FinalIncorrect;
      }

      //! @brief get the value specified by the enum
      //! @param DATA the data to retrieve from the partition
      //! @return the value specified by the enum
      float GetData( const Data &DATA) const;

      //! @brief add incorrect estimate info and total dataset size to a binary partition
      //! @param DATA_SIZE size of the dataset
      //! @param STATE_COUNTS total state counts (from determine total state counts)
      void DetermineInitialAccuracy( const size_t &DATA_SIZE, const linal::Vector< size_t> &STATE_COUNTS);

      //! @brief determine the final accuracy of the partition
      //! @param DATA the initial dataset
      void DetermineAccuracy( const storage::Vector< FeatureResultAndState> &DATA);

    ///////////////
    // operators //
    ///////////////

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

      //! @brief estimate # incorrect
      //! @param SIZE the total # of elements represented by the state counts vector
      //! @param STATE_COUNTS the count of features in each state
      static size_t EstimateNumberIncorrect( const size_t &SIZE, const linal::Vector< size_t> &STATE_COUNTS);

    }; // class DtreeBinaryPartition

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_DTREE_BINARY_PARTITION_H_
