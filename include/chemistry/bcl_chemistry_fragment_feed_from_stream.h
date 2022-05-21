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

#ifndef BCL_CHEMISTRY_FRAGMENT_FEED_FROM_STREAM_H_
#define BCL_CHEMISTRY_FRAGMENT_FEED_FROM_STREAM_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_fragment_feed_interface.h"
#include "io/bcl_io_ifstream.h"
#include "sdf/bcl_sdf_mdl_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FragmentFeedFromStream
    //! @brief Iteratively loads molecules from an arbitrary source
    //!
    //! @see @link example_chemistry_fragment_feed_from_stream.cpp @endlink
    //! @author mendenjl
    //! @date May 20, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentFeedFromStream :
      public FragmentFeedInterface
    {

    //////////
    // data //
    //////////

    private:

      std::istream &m_Input; //!< Input stream
      std::istream::pos_type m_InitialPosition; //!< Initial position in the stream
      size_t        m_Size;  //!< Number of molecules in the current file; only calculated if requested

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from stream
      FragmentFeedFromStream( std::istream &INPUT_STREAM);

      //! @brief Copy constructor
      FragmentFeedFromStream( const FragmentFeedFromStream &FEED);

      //! @brief destructor
      ~FragmentFeedFromStream();

      //! @brief Clone function
      //! @return pointer to new MolecularConformationShared
      FragmentFeedFromStream *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get the parameter for use with initialization
      //! @return the static list of initialization flags
      command::Parameter GetParameter() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief open the feed from the given location (file or db origin)
      //! @param LOCATION location to open
      //! return true if the given location existed and contained at least one molecule
      bool Open( const std::string &LOCATION);

      //! @brief skip the current molecule
      //! @return false if the end of the source was reached
      bool Skip();

      //! @brief get the size of the currently-open source
      //! @param MAX_TO_CONSIDER maximum number to return; do not read more than this
      //! @return size of source that is currently open, or MAX_TO_CONSIDER, whichever is smaller
      size_t GetSize( const size_t &MAX_TO_CONSIDER = std::numeric_limits< size_t>::max());

      //! @brief operator ++ (prefix, e.g. ++a)
      //! @param HANDLER to read the molecule into
      //! @return true if the next molecule could be retrieved
      bool RetrieveNextMolecule( sdf::MdlHandler &HANDLER);

    }; // class FragmentFeedFromStream

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_FRAGMENT_FEED_FROM_STREAM_H_
