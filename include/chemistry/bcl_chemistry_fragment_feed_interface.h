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

#ifndef BCL_CHEMISTRY_FRAGMENT_FEED_INTERFACE_H_
#define BCL_CHEMISTRY_FRAGMENT_FEED_INTERFACE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "command/bcl_command_parameter.h"
#include "sdf/bcl_sdf_mdl_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FragmentFeedInterface
    //! @brief Iteratively loads molecules from an arbitrary source
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date May 16, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentFeedInterface :
      public util::ObjectInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new MolecularConformationShared
      virtual FragmentFeedInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief get the parameter for use with initialization
      //! @return the static list of initialization flags
      virtual command::Parameter GetParameter() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief open the feed from the given location (file or db origin)
      //! @param LOCATION location to open
      //! @return true if the given location existed
      virtual bool Open( const std::string &LOCATION) = 0;

      //! @brief skip the current molecule
      //! @return false if the end of the source was reached
      virtual bool Skip() = 0;

      //! @brief get the size of the currently-open source
      //! @param MAX_TO_CONSIDER maximum number to return; do not read more than this
      //! @return size of source that is currently open, or MAX_TO_CONSIDER, whichever is smaller
      virtual size_t GetSize( const size_t &MAX_TO_CONSIDER = std::numeric_limits< size_t>::max()) = 0;

      //! @brief operator ++ (prefix, e.g. ++a)
      //! @param HANDLER to read the molecule into
      //! @return true if the next molecule could be retrieved
      virtual bool RetrieveNextMolecule( sdf::MdlHandler &HANDLER) = 0;

    ////////////////////
    // helper methods //
    ////////////////////

    protected:

      //! @brief add this feed to the enumerated list, so that it will be used by FragmentFeed
      //! @return this
      FragmentFeedInterface *AddFeed();

    //////////////////////
    // input and output //
    //////////////////////

    private:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class FragmentFeedInterface

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_FRAGMENT_FEED_INTERFACE_H_
