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

#ifndef BCL_DESCRIPTOR_WINDOW_ALIGNMENT_TYPE_H_
#define BCL_DESCRIPTOR_WINDOW_ALIGNMENT_TYPE_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util_wrapper_enum.h"

// includes from bcl - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @file bcl_descriptor_window_alignment_type.h
    //! @brief enum wrapper indicating how windows are aligned relative to the element of the window under consideration
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Mar 12, 2013
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! Alignment of the window used for the periodogram
    enum WindowAlignmentType
    {
      e_Left,   //!< Window considers this element and m_Size prior elements
      e_Right,  //!< Window considers this element and the next m_Size elements
      e_Center, //!< Window considers this element and the m_Size prior and next elements
      e_JufoCenter, //!< Window alignment for old jufo; used for testing the old jufo only
      s_NumberOfWindowAlignments
    };

    //! @brief Alignment as string
    //! @param ALIGNMENT the alignment
    //! @return the Alignment as string
    const std::string &GetWindowAlignmentString( const WindowAlignmentType &ALIGNMENT);

    //! SSEInfoTypeEnum simplifies the usage of the SSEInfoType enum of this class
    typedef util::WrapperEnum< WindowAlignmentType, &GetWindowAlignmentString, s_NumberOfWindowAlignments>
      WindowAlignmentEnum;

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_WINDOW_ALIGNMENT_TYPE_H_
