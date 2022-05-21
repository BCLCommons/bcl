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

#ifndef BCL_DESCRIPTOR_SEQUENCE_SEGMENT_STATISTICS_H_
#define BCL_DESCRIPTOR_SEQUENCE_SEGMENT_STATISTICS_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_sequence.h"
#include "bcl_descriptor_segment_finder.h"
#include "bcl_descriptor_segment_info.h"
#include "bcl_descriptor_window.h"
#include "bcl_descriptor_window_weighting_interface.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SequenceSegmentStatistics
    //! @brief Computes information about the segments in a sequence
    //!        A segment is a region of consecutive elements of a sequence which have the same values returned by a user
    //!        defined conditional descriptor.  The user can also select a descriptor to take mean/sd for over this
    //!        segment if they desire. Information about adjacent segments, as well as overall sequence statistics may
    //!        also be obtained. This class is similar to WindowConditional average but differs in that it looks at
    //!        segments, which are not relative to any particular amino acid, rather than length-restricted windows.
    //!
    //! @tparam t_DataType type of the argument for which the description is going to be generated
    //!
    //! @see @link example_descriptor_sequence_segment_statistics.cpp @endlink
    //! @author mendenjl
    //! @date Mar 06, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class SequenceSegmentStatistics :
      public BaseSequence< t_DataType, float>
    {

    private:

    //////////
    // data //
    //////////

      //! descriptor upon which the averaging is conditional; The average continues only so long as this descriptor
      //! returns the same value(s)
      util::Implementation< Base< t_DataType, float> > m_ConditionalDescriptor;

      //! descriptor for which to compute statistics
      util::Implementation< Base< t_DataType, float> > m_Descriptor;

      //! statistics that will be calculated on the sequence segments
      storage::Vector< SegmentFinder::StatisticEnum> m_Statistics;
      storage::Vector< size_t>                       m_StatisticsSizes; //!< # of values returned by each statistic
      size_t                                         m_StatisticsSizePerCondition; //!< Number of values collected per condition

      //! set of values returned by the condition to collect segment statistics for
      //! Note that the overall statistics will always be returned, whether or not m_Conditions is empty
      storage::Vector< linal::Vector< float> > m_Conditions;

    public:

    //////////
    // data //
    //////////

      //! instances of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SequenceSegmentStatistics();

      //! @brief Clone function
      //! @return pointer to new SequenceSegmentStatistics
      SequenceSegmentStatistics< t_DataType> *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief function to return derived-class-held implementations to this interface
      //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
      //! implementations
      iterate::Generic< Base< t_DataType, float> > GetInternalDescriptors();

    private:

      //! @brief calculate the descriptors
      //! @param STORAGE storage for the descriptor
      void Calculate
      (
        linal::VectorReference< float> &STORAGE
      );

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const
      {
        return e_PreferCache;
      }

    }; // class SequenceSegmentStatistics

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API SequenceSegmentStatistics< biol::AABase>;
    BCL_EXPIMP_TEMPLATE template class BCL_API SequenceSegmentStatistics< char>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_SEQUENCE_SEGMENT_STATISTICS_H_
