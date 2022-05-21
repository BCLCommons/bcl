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
#include "chemistry/bcl_chemistry_fragment_feed_from_stream.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FragmentFeedFromStream::FragmentFeedFromStream( std::istream &INPUT_STREAM) :
      m_Input( INPUT_STREAM),
      m_InitialPosition( INPUT_STREAM.tellg()),
      m_Size( 0)
    {
    }

    //! @brief Copy constructor
    FragmentFeedFromStream::FragmentFeedFromStream( const FragmentFeedFromStream &FEED) :
      m_Input( FEED.m_Input),
      m_InitialPosition( FEED.m_InitialPosition),
      m_Size( FEED.m_Size)
    {
      m_Input.clear();
      m_Input.seekg( m_InitialPosition, std::ios::beg);
      m_Input.rdbuf()->pubseekoff( m_InitialPosition, std::ios_base::beg, std::ios_base::in);
    }

    //! @brief destructor
    FragmentFeedFromStream::~FragmentFeedFromStream()
    {
    }

    //! @brief Clone function
    //! @return pointer to new MolecularConformationShared
    FragmentFeedFromStream *FragmentFeedFromStream::Clone() const
    {
      return new FragmentFeedFromStream( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &FragmentFeedFromStream::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the parameter for use with initialization
    //! @return the static list of initialization flags
    command::Parameter FragmentFeedFromStream::GetParameter() const
    {
      return command::Parameter();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief open the feed from the given location (file or db origin)
    //! @param LOCATION location to open
    //! return true if the given location existed and contained at least one molecule
    bool FragmentFeedFromStream::Open( const std::string &LOCATION)
    {
      m_Input.clear();
      m_Input.seekg( m_InitialPosition, std::ios_base::beg);
      m_Input.rdbuf()->pubseekoff( m_InitialPosition, std::ios_base::beg, std::ios_base::in);
      m_Size = util::GetUndefined< size_t>();
      return true;
    }

    //! @brief skip the current molecule
    //! @return false if the end of the source was reached
    bool FragmentFeedFromStream::Skip()
    {
      std::string temp;
      size_t n_read( 0);
      while( m_Input.good() && std::getline( m_Input, temp))
      {
        ++n_read;
        if( sdf::MdlHandler::IsTerminalLine( temp))
        {
          return true;
        }
      }
      // minimum sdf length if 5 lines for an empty mol (3-description, one header, M  END)
      // mol-files do not end in $$$$ but only have one molecule, so this
      // check allows us to still recognize potentially valid mol-files
      return n_read >= 5;
    }

    //! @brief get the size of the currently-open source
    //! @param MAX_TO_CONSIDER maximum number to return; do not read more than this
    //! @return size of source that is currently open, or MAX_TO_CONSIDER, whichever is smaller
    size_t FragmentFeedFromStream::GetSize( const size_t &MAX_TO_CONSIDER)
    {
      if( !util::IsDefined( m_Size))
      {
        m_Size = 0;

        // determine the size by skipping all the molecules
        while( m_Size < MAX_TO_CONSIDER && Skip())
        {
          ++m_Size;
        }
        m_Input.clear();
        m_Input.seekg( m_InitialPosition, std::ios_base::beg);
        m_Input.rdbuf()->pubseekoff( m_InitialPosition, std::ios_base::beg, std::ios_base::in);
      }
      return m_Size;
    }

    //! @brief operator ++ (prefix, e.g. ++a)
    //! @param HANDLER to read the molecule into
    //! @return true if the next molecule could be retrieved
    bool FragmentFeedFromStream::RetrieveNextMolecule( sdf::MdlHandler &HANDLER)
    {
      if( !m_Input.good())
      {
        return false;
      }
      HANDLER.ReadFromSDF( m_Input);
      return HANDLER.IsValid();
    }

  } // namespace chemistry
} // namespace bcl
