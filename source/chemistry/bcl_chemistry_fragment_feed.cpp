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
#include "chemistry/bcl_chemistry_fragment_feed.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_feed_from_file.h"
#include "chemistry/bcl_chemistry_fragment_feed_from_stream.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_thunk_job.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    //! @brief access the static list of feeds
    //! @return the static list of feeds
    storage::List< util::OwnPtr< FragmentFeedInterface> > &FragmentFeed::GetFeeds()
    {
      static storage::List< util::OwnPtr< FragmentFeedInterface> > s_feeds;
      return s_feeds;
    }

    //! @brief access the static list of parameters
    //! @return the static list of parameters
    storage::Vector< command::Parameter> &FragmentFeed::GetParameters()
    {
      static storage::Vector< command::Parameter> s_parameters;
      return s_parameters;
    }

    //! @brief access the static list of flags for a given feed name
    //! @param NAME the name of the feed; different feeds have the same parameter but different prefixs
    //! @return the static list of flags for a given feed name
    const util::ShPtrVector< command::FlagInterface> &FragmentFeed::GetFlags( const std::string &NAME)
    {
      static storage::Map< std::string, util::ShPtrVector< command::FlagInterface> > s_flag_map;
      storage::Map< std::string, util::ShPtrVector< command::FlagInterface> >::iterator
        itr( s_flag_map.Find( NAME));
      if( itr != s_flag_map.End())
      {
        return itr->second;
      }

      // construct new flags with this name as the prefix
      util::ShPtrVector< command::FlagInterface> new_flags;
      for
      (
        storage::Vector< command::Parameter>::const_iterator
          itr_param( GetParameters().Begin()), itr_param_end( GetParameters().End());
        itr_param != itr_param_end;
        ++itr_param
      )
      {
        new_flags.PushBack
        (
          util::ShPtr< command::FlagInterface>
          (
            new command::FlagDynamic
            (
              NAME + "_" + itr_param->GetName(),
              NAME + " " + itr_param->GetDescription(),
              *itr_param
            )
          )
        );
      }
      new_flags.PushBack
      (
        util::ShPtr< command::FlagInterface>
        (
          new command::FlagStatic
          (
            NAME + "_start",
            "index (0-offset) of first molecule to load for " + NAME,
            command::Parameter
            (
              "start",
              "index (0-offset) of first molecule to load",
              command::ParameterCheckRanged< size_t>(),
              "0"
            )
          )
        )
      );
      new_flags.PushBack
      (
        util::ShPtr< command::FlagInterface>
        (
          new command::FlagStatic
          (
            NAME + "_max",
            "Specify the maximum number of molecules to be loaded for " + NAME,
            command::Parameter
            (
              "number",
              "number of molecules",
              command::ParameterCheckRanged< size_t>( 1, std::numeric_limits< size_t>::max()),
              util::Format()( std::numeric_limits< size_t>::max())
            )
          )
        )
      );
      return s_flag_map[ NAME] = new_flags;
    }

    sched::Mutex FragmentFeed::s_LoadingMutex;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FragmentFeed::FragmentFeed( const std::string &NAME) :
      m_MoleculeIndex( 0),
      m_ReadMoleculeIndex( 0),
      m_MaxToRead( std::numeric_limits< size_t>::max()),
      m_HPref( sdf::GetCommandLineHydrogensPref()),
      m_Name( NAME),
      m_Neutralize( sdf::GetCommandLineNeutralizationPref()),
      m_NThreads( 0),
      m_MolInThreadBlockId( 0),
      m_NumberMolsPerGroupRead( 1),
      m_ShowStatus( true)
    {
      Restart();
    }

    //! @brief constructor from pre-fix
    FragmentFeed::FragmentFeed( const std::string &NAME, const sdf::HydrogenHandlingPref &PREF) :
      m_MoleculeIndex( 0),
      m_ReadMoleculeIndex( 0),
      m_MaxToRead( std::numeric_limits< size_t>::max()),
      m_HPref( PREF),
      m_Name( NAME),
      m_Neutralize( sdf::GetCommandLineNeutralizationPref()),
      m_NThreads( 0),
      m_MolInThreadBlockId( 0),
      m_NumberMolsPerGroupRead( 1),
      m_ShowStatus( true)
    {
      Restart();
    }

    //! @brief constructor from h-preference and filenames
    //! @param PREF hydrogen handling preference
    //! @param FILENAMES filenames for this fragment feed
    //! @param MAX_TO_READ maximum # of molecules to read
    FragmentFeed::FragmentFeed
    (
      const storage::Vector< std::string> &FILENAMES,
      const sdf::HydrogenHandlingPref &PREF,
      const size_t &MAX_TO_READ,
      const sdf::NeutralizationPref &NEUTRALIZE
    ) :
      m_MoleculeIndex( 0),
      m_ReadMoleculeIndex( 0),
      m_MaxToRead( MAX_TO_READ),
      m_HPref( PREF),
      m_Name(),
      m_Neutralize( NEUTRALIZE == sdf::e_CmdLine ? sdf::GetCommandLineNeutralizationPref() : sdf::NeutralizationPrefEnum( NEUTRALIZE)),
      m_NThreads( 0),
      m_MolInThreadBlockId( 0),
      m_NumberMolsPerGroupRead( 1),
      m_ShowStatus( false)
    {
      // load the initializers
      m_Feeds.Reset();

      // get the initializers for this feed
      storage::Vector< std::string> initializers( FILENAMES);

      // make a reference on the feed
      util::OwnPtr< FragmentFeedInterface> feed( new FragmentFeedFromFile);

      // walk through each initializer, check if it can be opened
      for
      (
        storage::Vector< std::string>::const_iterator
          itr_init( initializers.Begin()), itr_init_end( initializers.End());
        itr_init != itr_init_end;
        ++itr_init
      )
      {
        if( !feed->Open( *itr_init))
        {
          // skip initializers that could not be opened
          BCL_MessageStd( "FragmentFeed: couldn't open " + *itr_init + " for reading");
          continue;
        }

        m_Feeds.PushBack( storage::Pair< util::OwnPtr< FragmentFeedInterface>, std::string>( feed, *itr_init));
      }

      m_ItrFeed = m_Feeds.Begin();

      // check for empty feed
      m_MoleculeIndex = 0;
      if( m_Feeds.IsEmpty())
      {
        return;
      }

      // open the first feed
      m_ItrFeed->First()->Open( m_ItrFeed->Second());

      --m_MoleculeIndex; // decrement index, to account for 0 offset
      SetupThreads();

      // load in the next molecule, which is really the first molecule
      operator++();
    }

    //! @brief constructor from a stream, h-preference, and maximum # to read
    //! @param STREAM input stream for this fragment feed
    //! @param PREF hydrogen handling preference
    //! @param MAX_TO_READ maximum # of molecules to read
    //! @param INPUT_START the first molecule # to read
    FragmentFeed::FragmentFeed
    (
      std::istream &STREAM,
      const sdf::HydrogenHandlingPref &PREF,
      const size_t &MAX_TO_READ,
      const size_t &INPUT_START,
      const sdf::NeutralizationPref &NEUTRALIZE
    ) :
      m_MoleculeIndex( 0),
      m_ReadMoleculeIndex( 0),
      m_MaxToRead( MAX_TO_READ),
      m_HPref( PREF),
      m_Name(),
      m_Neutralize( NEUTRALIZE == sdf::e_CmdLine ? sdf::GetCommandLineNeutralizationPref() : sdf::NeutralizationPrefEnum( NEUTRALIZE)),
      m_NThreads( 0),
      m_MolInThreadBlockId( 0),
      m_NumberMolsPerGroupRead( 1),
      m_ShowStatus( false)
    {
      // load the initializers
      m_Feeds.Reset();

      // make a reference on the feed
      util::OwnPtr< FragmentFeedInterface> feed( new FragmentFeedFromStream( STREAM));

      if( !feed->Open( ""))
      {
        // skip initializers that could not be opened
        return;
      }

      m_Feeds.PushBack( storage::Pair< util::OwnPtr< FragmentFeedInterface>, std::string>( feed, ""));
      m_ItrFeed = m_Feeds.Begin();

      // check for empty feed
      m_MoleculeIndex = 0;
      if( m_Feeds.IsEmpty())
      {
        return;
      }

      // skip any undesired molecules at the beginning
      while( m_MoleculeIndex < INPUT_START && m_ItrFeed->First()->Skip())
      {
        ++m_MoleculeIndex;
      }
      m_MoleculeIndex = 0;

      --m_MoleculeIndex; // decrement index, to account for 0 offset
      SetupThreads();

      // load in the next molecule, which is really the first molecule
      operator++();
    }

    //! @brief virtual desctructor
    FragmentFeed::~FragmentFeed()
    {
      // join all existing threads
      CleanUp();
    }

    //! @brief Clone function
    //! @return pointer to new FragmentFeeder
    FragmentFeed *FragmentFeed::Clone() const
    {
      return new FragmentFeed( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &FragmentFeed::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief virtual function that will be called if cleanable is registered
    void FragmentFeed::CleanUp()
    {
      // join all existing threads
      for( size_t thread( 0), tsize( m_StandardizerThreads.GetSize()); thread < tsize; ++thread)
      {
        sched::GetScheduler().Join( m_StandardizerThreads( thread));
      }
    }

    //! @brief add command line flags for reading molecules from all available sources
    //! @param NAME name of feed
    void FragmentFeed::AddFlags( command::Command &CMD, const std::string &NAME)
    {
      const util::ShPtrVector< command::FlagInterface> &flags( GetFlags( NAME));
      for
      (
        util::ShPtrVector< command::FlagInterface>::const_iterator itr( flags.Begin()), itr_end( flags.End());
        itr != itr_end;
        ++itr
      )
      {
        CMD.AddFlag( *itr);
      }
    }

    //! @brief restart the feed by re-reading the command line flags
    void FragmentFeed::Restart()
    {
      const util::ShPtrVector< command::FlagInterface> &flags( GetFlags( m_Name));

      // load data from the flags
      size_t start_index( flags( flags.GetSize() - 2)->GetFirstParameter()->GetNumericalValue< size_t>());
      m_MaxToRead = flags.LastElement()->GetFirstParameter()->GetNumericalValue< size_t>();

      // reduce max to read so that start_index + m_MaxToRead <= std::numeric_limits< size_t>()
      m_MaxToRead = std::min( m_MaxToRead, std::numeric_limits< size_t>::max() - start_index);

      if( GetFeeds().IsEmpty())
      {
        return;
      }

      util::GetLogger().LogStatus( "Counting molecules for " + m_Name);

      // load the initializers
      m_Feeds.Reset();
      size_t total_seen( 0); // track the total number of molecules seen in the sources considered
      storage::List< util::OwnPtr< FragmentFeedInterface> >::iterator itr_feed( GetFeeds().Begin());
      for
      (
        size_t feed_number( 0), number_feeds( GetFeeds().GetSize());
        feed_number < number_feeds && total_seen < m_MaxToRead;
        ++feed_number, ++itr_feed
      )
      {
        // get the initializers for this feed
        storage::Vector< std::string> initializers( flags( feed_number)->GetStringList());

        // make a reference on the feed
        util::OwnPtr< FragmentFeedInterface> &feed( *itr_feed);

        // walk through each initializer, check if it can be opened
        for
        (
          storage::Vector< std::string>::const_iterator
            itr_init( initializers.Begin()), itr_init_end( initializers.End());
          itr_init != itr_init_end;
          ++itr_init
        )
        {
          if( !feed->Open( *itr_init))
          {
            BCL_MessageCrt
            (
              "Could not open feed of type " + feed->GetClassIdentifier() + " with initializer " + *itr_init
            );
            // skip initializers that could not be opened
            continue;
          }

          // check whether the first desired index is high enough that this source will be skipped anyway
          if( start_index >= feed->GetSize( m_MaxToRead + start_index))
          {
            // skip the feed but subtract its size from the start index
            start_index -= feed->GetSize( m_MaxToRead + start_index);
            continue;
          }

          if( m_Feeds.IsEmpty())
          {
            // first source after starting loading; add molecules after start index to total seen
            total_seen += feed->GetSize() - start_index;
          }
          else
          {
            total_seen += feed->GetSize();
          }
          m_Feeds.PushBack( storage::Pair< util::OwnPtr< FragmentFeedInterface>, std::string>( feed, *itr_init));
          if( total_seen >= m_MaxToRead)
          {
            total_seen = m_MaxToRead;
            break;
          }
        }
      }

      m_MaxToRead = std::min( total_seen, m_MaxToRead);
      util::GetLogger().LogStatus( "Found " + util::Format()( m_MaxToRead) + " molecules");

      m_ItrFeed = m_Feeds.Begin();

      // check for empty feed
      m_MoleculeIndex = 0;
      if( m_Feeds.IsEmpty())
      {
        return;
      }

      // open the first feed
      m_ItrFeed->First()->Open( m_ItrFeed->Second());

      // get to the first desired molecule
      for( size_t molecules_skipped( 0); molecules_skipped < start_index; ++molecules_skipped)
      {
        m_ItrFeed->First()->Skip();
      }

      --m_MoleculeIndex; // decrement index, to account for 0 offset

      SetupThreads();

      // load in the next molecule, which is really the first molecule
      operator++();
    }

    //! @return the data member
    FragmentFeed::reference FragmentFeed::operator *()
    {
      return m_ThreadMolecules( m_MolInThreadBlockId);
    }

    //! @return a pointer to the data member
    FragmentFeed::pointer FragmentFeed::operator ->()
    {
      return &m_ThreadMolecules( m_MolInThreadBlockId);
    }

    //! @brief operator ++ (prefix, e.g. ++a)
    //! @return a reference to the iterator after incrementing
    FragmentFeed &FragmentFeed::operator ++()
    {
      if( ++m_MoleculeIndex >= m_MaxToRead)
      {
        WriteStatus();
        return *this;
      }
      if( ++m_MolInThreadBlockId >= m_NumberMolsPerGroupRead)
      {
        LoadNextBlock();
      }
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FragmentFeed::Read( std::istream &ISTREAM)
    {
      // there is no need to read or write this class; it is read in from the users input on the command line
      // additionally, it cannot be read because it contains siptrs and iterators
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &FragmentFeed::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // there is no need to read or write this class; it is read in from the users input on the command line
      // additionally, it cannot be read because it contains siptrs and iterators
      return OSTREAM;
    }

    //! @brief write the status message
    void FragmentFeed::WriteStatus() const
    {
      if( m_ShowStatus)
      {
        const size_t status_bar_length( 20);

        // determine progress percent
        const size_t percent( float( m_MoleculeIndex) * 100.0 / float( m_MaxToRead));

        // determine number of stars in the status bar
        const size_t number_stars( percent * status_bar_length / 100);

        const std::string status
        (
          "["
          + std::string( number_stars, '*')
          + std::string( status_bar_length - number_stars, ' ')
          + "] "
          + util::Format()( percent) + "% "
          + util::Format()( m_MoleculeIndex) + "/" + util::Format()( m_MaxToRead)
          + " molecules read"
        );
        util::GetLogger().LogStatus( status);
      }
    }

    //! @brief read the next molecule from any stream into HANDLER
    //! @param HANDLER Mdl handler to read into
    //! @return read molecule index (undefined on failure)
    size_t FragmentFeed::ThreadReadNextMolecule( sdf::MdlHandler &HANDLER)
    {
      m_ReadingMutex.Lock();
      if( m_ReadMoleculeIndex == m_MaxToRead || m_ReadMoleculeIndex - m_MoleculeIndex >= m_NumberMolsPerGroupRead)
      {
        m_ReadingMutex.Unlock();
        return util::GetUndefined< size_t>();
      }
      size_t read_index( m_ReadMoleculeIndex);
      // try to retrieve the next molecule from this source
      if( m_ItrFeed->First()->RetrieveNextMolecule( HANDLER))
      {
        // success, return true
        //HANDLER.SetNumberMoleculesRead( m_ReadMoleculeIndex);
        ++m_ReadMoleculeIndex;
        m_ReadingMutex.Unlock();
        return read_index;
      }

      // reached the last molecule in this source, open next source.
      // Because empty feeds are discarded, there is no need for a loop here
      if( m_ShowStatus)
      {
        // open using the next initializer, check that there are still feeds and that they can be opened
        // this assert will only be thrown if the file or db was truncated after Restart() was last called
        BCL_Assert
        (
          ++m_ItrFeed != m_Feeds.End()
          && m_ItrFeed->First()->Open( m_ItrFeed->Second())
          && m_ItrFeed->First()->RetrieveNextMolecule( HANDLER),
          "Could not retrieve molecule from source: " + m_ItrFeed->Second()
        );
        //HANDLER.SetNumberMoleculesRead( m_ReadMoleculeIndex);
        ++m_ReadMoleculeIndex;
        m_ReadingMutex.Unlock();
        return read_index;
      }

      // status not shown, instead we must check whether the end was reached since sizes are not calculated ahead
      // of time unless status is shown
      if( m_ItrFeed == m_Feeds.End() || ++m_ItrFeed == m_Feeds.End())
      {
        m_MaxToRead = m_ReadMoleculeIndex;
        m_ReadingMutex.Unlock();
        return util::GetUndefined< size_t>();
      }
      else if( !m_ItrFeed->First()->Open( m_ItrFeed->Second()))
      {
        while( ++m_ItrFeed != m_Feeds.End() && !m_ItrFeed->First()->Open( m_ItrFeed->Second()))
        {
        }
        if( m_ItrFeed == m_Feeds.End())
        {
          // last feed
          m_MaxToRead = m_ReadMoleculeIndex;
          m_ReadingMutex.Unlock();
          return util::GetUndefined< size_t>();
        }
      }
      if( !m_ItrFeed->First()->RetrieveNextMolecule( HANDLER))
      {
        m_MaxToRead = m_ReadMoleculeIndex;
        m_ReadingMutex.Unlock();
        return util::GetUndefined< size_t>();
      }
      //HANDLER.SetNumberMoleculesRead( m_ReadMoleculeIndex);
      ++m_ReadMoleculeIndex;
      m_ReadingMutex.Unlock();
      return read_index;
    }

    //! @brief standardize a particular molecule # in the thread block
    void FragmentFeed::Standardize()
    {
      sdf::MdlHandler handler;
      size_t next_mol_to_read( ThreadReadNextMolecule( handler));
      while( util::IsDefined( next_mol_to_read))
      {
        // try to get the next molecule id to load
        m_ThreadMolecules( next_mol_to_read - m_MoleculeIndex) = sdf::FragmentFactory::MakeFragment( handler, m_HPref, m_Neutralize);
        next_mol_to_read = ThreadReadNextMolecule( handler);
      }
    }

    //! @brief Load next block of molecules and launch threads to standardize them
    void FragmentFeed::LoadNextBlock()
    {
      m_MolInThreadBlockId = 0;
      m_ReadMoleculeIndex = m_MoleculeIndex;
      for( size_t i( 0); i < m_NThreads; ++i)
      {
        // launch a thread to handle this job
        sched::GetScheduler().RunJob( m_StandardizerThreads( i));
      }
      // wait until all molecules have been standardized, then finish
      CleanUp();
      WriteStatus();
    }

    //! @brief setup the data members associated with threads
    void FragmentFeed::SetupThreads()
    {
      // join all existing threads
      m_NThreads = std::max( size_t( 1), std::min( sched::GetNumberCPUs(), m_MaxToRead));
      m_NumberMolsPerGroupRead = std::min( 100 * m_NThreads, m_MaxToRead);
      m_StandardizerThreads.Resize( m_NThreads);
      m_ThreadMolecules.Resize( m_NumberMolsPerGroupRead);
      m_MolInThreadBlockId = m_NumberMolsPerGroupRead;
      for( size_t i( 0); i < m_NThreads; ++i)
      {
        // create the thread job
        m_StandardizerThreads( i) =
          util::ShPtr< sched::JobInterface>
          (
            new sched::ThunkJob< FragmentFeed, void>
            (
              1,
              *this,
              &FragmentFeed::Standardize,
              sched::JobInterface::e_READY,
              NULL
            )
          );
      }
    }

  } // namespace chemistry
} // namespace bcl
