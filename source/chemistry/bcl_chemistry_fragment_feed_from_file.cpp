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
#include "chemistry/bcl_chemistry_fragment_feed_from_file.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_parameter_check_file_existence.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> FragmentFeedFromFile::s_Instance
    (
      ( new FragmentFeedFromFile())->AddFeed()
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Copy constructor
    FragmentFeedFromFile::FragmentFeedFromFile( const FragmentFeedFromFile &FEED) :
      m_Input(),
      m_Filename( FEED.m_Filename),
      m_Size( FEED.m_Size)
    {
      if( m_Size)
      {
        Open( m_Filename);
      }
    }

    //! @brief destructor
    FragmentFeedFromFile::~FragmentFeedFromFile()
    {
      io::File::CloseClearFStream( m_Input);
    }

    //! @brief Clone function
    //! @return pointer to new MolecularConformationShared
    FragmentFeedFromFile *FragmentFeedFromFile::Clone() const
    {
      return new FragmentFeedFromFile( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &FragmentFeedFromFile::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the parameter for use with initialization
    //! @return the static list of initialization flags
    command::Parameter FragmentFeedFromFile::GetParameter() const
    {
      return
        command::Parameter
        (
          "filenames",
          "files containing molecules in sdf format",
          command::ParameterCheckFileExistence()
        );
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief open the feed from the given location (file or db origin)
    //! @param LOCATION location to open
    //! return true if the given location existed and contained at least one molecule
    bool FragmentFeedFromFile::Open( const std::string &LOCATION)
    {
      io::File::CloseClearFStream( m_Input);
      if( !io::File::TryOpenIFStream( m_Input, LOCATION))
      {
        m_Size = 0;
        return false;
      }
      else if( LOCATION != m_Filename)
      {
        m_Filename = LOCATION;
        m_Size = util::GetUndefined< size_t>();
      }
      return true;
    }

    //! @brief skip the current molecule
    //! @return false if the end of the source was reached
    bool FragmentFeedFromFile::Skip()
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
    size_t FragmentFeedFromFile::GetSize( const size_t &MAX_TO_CONSIDER)
    {
      if( !util::IsDefined( m_Size))
      {
        m_Size = CountMolecules( m_Filename, MAX_TO_CONSIDER);
      }
      return m_Size;
    }

    //! @brief operator ++ (prefix, e.g. ++a)
    //! @param HANDLER to read the molecule into
    //! @return true if the next molecule could be retrieved
    bool FragmentFeedFromFile::RetrieveNextMolecule( sdf::MdlHandler &HANDLER)
    {
      if( !m_Input.good())
      {
        BCL_MessageStd( "FragmentFeedFromFile: bad input stream");
        io::File::CloseClearFStream( m_Input);
        return false;
      }
      HANDLER.ReadFromSDF( m_Input);
      return HANDLER.IsValid() || m_Input.good();
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief count the number of molecules in a file
    //! @param FILENAME name of the file
    //! @param MAX maximum number to look for
    size_t FragmentFeedFromFile::CountMolecules( const std::string &FILENAME, const size_t &MAX)
    {
      // make another fragment feed
      FragmentFeedFromFile size_tester;
      // open the same file
      size_tester.Open( FILENAME);

      size_t size( 0);

      // determine the size by skipping all the molecules
      while( size < MAX && size_tester.Skip())
      {
        ++size;
      }
      return size;
    }

  } // namespace chemistry
} // namespace bcl
