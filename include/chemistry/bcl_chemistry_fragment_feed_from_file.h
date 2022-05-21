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

#ifndef BCL_CHEMISTRY_FRAGMENT_FEED_FROM_FILE_H_
#define BCL_CHEMISTRY_FRAGMENT_FEED_FROM_FILE_H_

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
    //! @class FragmentFeedFromFile
    //! @brief Iteratively loads molecules from an arbitrary source
    //!
    //! @see @link example_chemistry_fragment_feed_from_file.cpp @endlink
    //! @author mendenjl
    //! @date May 16, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentFeedFromFile :
      public FragmentFeedInterface
    {

    //////////
    // data //
    //////////

    private:

      io::IFStream   m_Input;    //!< Input stream
      std::string    m_Filename; //!< Filename
      size_t         m_Size;     //!< Number of molecules in the current file; only calculated if requested

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FragmentFeedFromFile() :
        m_Size( 0)
      {
      }

      //! @brief Copy constructor
      FragmentFeedFromFile( const FragmentFeedFromFile &FEED);

      //! @brief destructor
      ~FragmentFeedFromFile();

      //! @brief Clone function
      //! @return pointer to new MolecularConformationShared
      FragmentFeedFromFile *Clone() const;

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

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief count the number of molecules in a file
      //! @param FILENAME name of the file
      //! @param MAX maximum number to look for
      static size_t CountMolecules
      (
        const std::string &FILENAME,
        const size_t &MAX = std::numeric_limits< size_t>::max()
      );

    }; // class FragmentFeedFromFile

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_FRAGMENT_FEED_FROM_FILE_H_ 
