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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "chemistry/bcl_chemistry_fragment_feed_interface.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_feed.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  /////////////////
  // data access //
  /////////////////

    //! @brief add this feed to the enumerated list
    //! @return this
    FragmentFeedInterface *FragmentFeedInterface::AddFeed()
    {
      // so that static initialization order is irrelevant, put this fragment feed at the start or beginning of the list
      // depending on static class name
      if
      (
        FragmentFeed::GetFeeds().IsEmpty()
        || FragmentFeed::GetFeeds().FirstElement()->GetClassIdentifier() > GetClassIdentifier()
      )
      {
        // add this feeder to the enumerated feeds; should only be done for s_Instance variables
        FragmentFeed::GetFeeds().PushFront( util::OwnPtr< FragmentFeedInterface>( this->Clone()));

        // add the initialization flag
        FragmentFeed::GetParameters().InsertElements( 0, GetParameter(), 1);
      }
      else
      {
        // add this feeder to the enumerated feeds; should only be done for s_Instance variables
        FragmentFeed::GetFeeds().PushBack( util::OwnPtr< FragmentFeedInterface>( this->Clone()));

        // add the initialization flag
        FragmentFeed::GetParameters().InsertElements( 1, GetParameter(), 1);
      }

      GetObjectInstances().AddInstance( this);

      return this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FragmentFeedInterface::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &FragmentFeedInterface::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
