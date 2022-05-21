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

#ifndef BCL_DESCRIPTOR_WINDOW_SEGMENT_STATISTICS_H_
#define BCL_DESCRIPTOR_WINDOW_SEGMENT_STATISTICS_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"
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
    //! @class WindowSegmentStatistics
    //! @brief Computes information about the current and nearby segments
    //!        A segment is a region of consecutive elements of a sequence which have the same values returned by a user
    //!        defined conditional descriptor.  The user can also select a descriptor to take mean/sd for over this
    //!        segment if they desire. Information about adjacent segments, as well as overall sequence statistics may
    //!        also be obtained. This class is similar to WindowConditional average but differs in that it looks at
    //!        segments, which are not relative to any particular amino acid, rather than length-restricted windows.
    //!
    //! @tparam t_DataType type of the argument for which the description is going to be generated
    //!
    //! @see @link example_descriptor_window_segment_statistics.cpp @endlink
    //! @author mendenjl
    //! @date Mar 05, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class WindowSegmentStatistics :
      public BaseElement< t_DataType, float>
    {

    private:

    //////////
    // data //
    //////////

      //! descriptor upon which the averaging is conditional; The average continues only so long as this descriptor
      //! returns the same value(s)
      util::Implementation< Base< t_DataType, float> > m_ConditionalDescriptor;

      //! descriptor for which to compute statistics
      util::Implementation< Base< t_DataType, float> > m_ValueDescriptor;

      //! maximum number of window segments in each direction to include in the description; 0 means only the current
      //! segment
      //! This allows obtaining the averages and lengths not only for the current segments, but also other nearby segments
      size_t m_SegmentsRadius;

      //! if set, ignore regions with undefined condition values (default is to consider them separate regions)
      bool m_IgnoreUndefinedConditionSegments;

      //! Position of each element in the sequence in the segment finder
      //! This vector is / used only when ignoring undefined condition segments, otherwise the position is that of the
      //! element in the sequence
      storage::Vector< size_t> m_SegmentFinderIndex;

      //! statistics that will be calculated on all elements of the segment-window
      storage::Vector< SegmentInfo::StatisticEnum> m_Statistics;
      storage::Vector< size_t> m_StatisticsSizes;

      //! segment finder
      SegmentFinder m_SegmentFinder;

      //! number of unconditional statistics values returned per iterator
      size_t m_SegmentWindowValueCount;

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
      WindowSegmentStatistics();

      //! @brief Clone function
      //! @return pointer to new WindowSegmentStatistics
      WindowSegmentStatistics< t_DataType> *Clone() const;

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

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      void SetObjectHook();

      //! @brief run the segment finder over the sequence
      void RunSegmentFinder();

      //! @brief calculate the descriptors
      //! @param ELEMENT the element of the sequence of interest
      //! @param STORAGE storage for the descriptor
      void Calculate
      (
        const iterate::Generic< const t_DataType> &ELEMENT,
        linal::VectorReference< float> &STORAGE
      );

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

    }; // class WindowSegmentStatistics

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API WindowSegmentStatistics< biol::AABase>;
    BCL_EXPIMP_TEMPLATE template class BCL_API WindowSegmentStatistics< char>;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_WINDOW_SEGMENT_STATISTICS_H_
