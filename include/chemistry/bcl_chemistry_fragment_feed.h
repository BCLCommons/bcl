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

#ifndef BCL_CHEMISTRY_FRAGMENT_FEED_H_
#define BCL_CHEMISTRY_FRAGMENT_FEED_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_fragment_complete.h"
#include "bcl_chemistry_fragment_feed_interface.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_interface.h"
#include "sched/bcl_sched_job_interface.h"
#include "sched/bcl_sched_mutex.h"
#include "sdf/bcl_sdf_mdl_handler.h"
#include "sdf/bcl_sdf_molecule_reading_pref.h"
#include "util/bcl_util_cleanable_interface.h"
#include "util/bcl_util_sh_ptr_vector.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FragmentFeed
    //! @brief Iterates through a series of locations (filenames or db origins) with molecules.  In this way, it is not
    //! necessary to store an entire fragment ensemble in memory.  Pointers to the fragment desired go out of scope
    //! immediately after calling operator++
    //!
    //! @see @link example_chemistry_fragment_feed.cpp @endlink
    //! @author mendenjl
    //! @date May 16, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentFeed :
      public util::ObjectInterface,
      public util::CleanableInterface
    {

    /////////////
    // friends //
    /////////////

      friend class FragmentFeedInterface; //!< accesses GetFlags and GetFeeds

    //////////
    // data //
    //////////

    private:

      size_t                    m_MoleculeIndex; //!< molecule index currently on (=0 at start index)
      size_t                    m_ReadMoleculeIndex; //!< Highest molecule index that has been read
      size_t                    m_MaxToRead;     //!< Maximum # of molecules to read
      sdf::HydrogenHandlingPref m_HPref;         //!< hydrogen handling pref
      std::string               m_Name;          //!< Name of this feed; prefix for all flags
      sdf::NeutralizationPref   m_Neutralize;    //!< Neutralization Preference

      //! list of feed sources and initializers
      storage::List< storage::Pair< util::OwnPtr< FragmentFeedInterface>, std::string> > m_Feeds;

      //! Corresponding iterator
      storage::List< storage::Pair< util::OwnPtr< FragmentFeedInterface>, std::string> >::iterator m_ItrFeed;

      //! Threading information
      size_t m_NThreads; //!< Number of threads
      sched::Mutex m_ReadingMutex; //!< Mutex to protect access to m_ItrFeed and m_Feeds
      size_t m_MolInThreadBlockId; //!< Index of current molecule in the standardizer threads
      size_t m_NumberMolsPerGroupRead; //!< Number of molecules read per iteration
      util::ShPtrVector< sched::JobInterface> m_StandardizerThreads; //!< Jobs for threads used to standardize molecules
      storage::Vector< FragmentComplete>      m_ThreadMolecules; //!< Fragments, one for each thread

      // bool: whether or not to show the status messages as molecules are loaded
      bool m_ShowStatus;

      //! @brief access the static list of feeds
      //! @return the static list of feeds
      static storage::List< util::OwnPtr< FragmentFeedInterface> > &GetFeeds();

      //! @brief access the static list of parameters
      //! @return the static list of parameters
      static storage::Vector< command::Parameter> &GetParameters();

      //! @brief access the static list of flags for a given feed name
      //! @param NAME the name of the feed; different feeds have the same parameter but different prefixs
      //! @return the static list of flags for a given feed name
      static const util::ShPtrVector< command::FlagInterface> &GetFlags( const std::string &NAME);

      //! @brief restrict the number of loading feeds that can be loaded from at once
      static sched::Mutex s_LoadingMutex;

    public:

    //////////////
    // typedefs //
    //////////////

      // These typedefs are common for all iterators
      typedef const FragmentComplete& reference;         //!< type returned by operator*
      typedef const FragmentComplete* pointer;           //!< type returned by operator->
      typedef FragmentComplete  value_type;              //!< type needed to hold the result of an operator*
      typedef std::input_iterator_tag iterator_category; //!< Type of iterator

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FragmentFeed( const std::string &NAME = "input");

      //! @brief constructor from h-preference
      FragmentFeed( const std::string &NAME, const sdf::HydrogenHandlingPref &PREF);

      //! @brief constructor from h-preference and filenames
      //! @param PREF hydrogen handling preference
      //! @param FILENAMES filenames for this fragment feed
      //! @param MAX_TO_READ maximum # of molecules to read
      FragmentFeed
      (
        const storage::Vector< std::string> &FILENAMES,
        const sdf::HydrogenHandlingPref &PREF,
        const size_t &MAX_TO_READ = std::numeric_limits< size_t>::max(),
        const sdf::NeutralizationPref &NEUTRALIZE = sdf::e_CmdLine
      );

      //! @brief constructor from a stream, h-preference, and maximum # to read
      //! @param STREAM input stream for this fragment feed
      //! @param PREF hydrogen handling preference
      //! @param MAX_TO_READ maximum # of molecules to read
      //! @param INPUT_START the first molecule # to read
      FragmentFeed
      (
        std::istream &STREAM,
        const sdf::HydrogenHandlingPref &PREF,
        const size_t &MAX_TO_READ = std::numeric_limits< size_t>::max(),
        const size_t &INPUT_START = size_t( 0),
        const sdf::NeutralizationPref &NEUTRALIZE = sdf::e_CmdLine
      );

      //! @brief virtual desctructor
      virtual ~FragmentFeed();

      //! @brief Clone function
      //! @return pointer to new FragmentFeeder
      FragmentFeed *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get the name of this feed
      //! @return the name of this feed
      const std::string &GetName() const
      {
        return m_Name;
      }

      //! @brief virtual function that will be called if cleanable is registered
      void CleanUp();

      //! @brief determine how many elements are between the iterator and the first element of the range
      //! @return a long indicating the iterators distance from Begin()
      //! @note return type has to be a long because iterator could be before Begin()
      //! @note O(1) for vectors, O( size of range) for everything else
      size_t GetPosition() const
      {
        return m_MoleculeIndex;
      }

      //! @brief test whether this iterator is at the end of its range
      //! @return true iff the iterator is not at the end
      bool NotAtEnd() const
      {
        return m_MoleculeIndex < m_MaxToRead;
      }

      //! @brief add command line flags for reading molecules from all available sources
      //! @param CMD command to which flags should be appended
      //! @param NAME name of feed
      static void AddFlags( command::Command &CMD, const std::string &NAME = "input");

      //! @brief restart the feed by re-reading the command line flags
      void Restart();

    ///////////////
    // operators //
    ///////////////

      //! @return the data member
      reference operator *();

      //! @return a pointer to the data member
      pointer operator ->();

      //! @brief operator ++ (prefix, e.g. ++a)
      //! @return a reference to the iterator after incrementing
      FragmentFeed &operator ++();

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

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief write the status message
      void WriteStatus() const;

      //! @brief standardize a particular molecule # in the thread block
      void Standardize();

      //! @brief Load next block of molecules and launch threads to standardize them
      void LoadNextBlock();

      //! @brief setup the data members associated with threads
      void SetupThreads();

      //! @brief read the next molecule from any stream into HANDLER
      //! @param HANDLER Mdl handler to read into
      //! @return read molecule index (undefined on failure)
      size_t ThreadReadNextMolecule( sdf::MdlHandler &HANDLER);
    }; // class FragmentFeed

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_FRAGMENT_FEED_H_
