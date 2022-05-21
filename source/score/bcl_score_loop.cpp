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
#include "score/bcl_score_loop.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
#include "command/bcl_command_command_state.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_histogram.h"
#include "score/bcl_score_energy_distribution.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> Loop::s_Instance
    (
      util::Enumerated< math::BinaryFunctionInterfaceSerializable< assemble::SSE, assemble::SSE, double> >::AddInstance
      (
        new Loop()
      )
    );

  /////////////////
  // data access //
  /////////////////

    //! @brief returns default file where loop distance histograms and in consequence energy potentials are read from
    //! @return default file where the statistics and in consequence the energy potentials are read from
    const std::string &Loop::GetDefaultHistogramFilename()
    {
      static const std::string s_default_histogram_filename( "loop.histograms");
      return s_default_histogram_filename;
    }

    //! @brief returns default file where the loop distance table is stored
    //! @return default file where the loop distance table data is stored
    const std::string &Loop::GetDefaultTableFilename()
    {
      static const std::string s_default_table_filename( "loop.tbl");
      return s_default_table_filename;
    }
    //! @brief returns default scheme
    //! @return default scheme
    const std::string &Loop::GetDefaultScheme()
    {
      static const std::string s_default_scheme( "loop");
      return s_default_scheme;
    }

    //! @brief returns default maximum loop length in residues to be considered
    //! @return default maximum loop length in residues to be considered
    size_t Loop::GetDefaultMaxLoopLength()
    {
      static const size_t s_default_max_loop_length( 25);
      return s_default_max_loop_length;
    }

    //! @brief get the name of the object
    //! @return the name of the object
    const std::string &Loop::GetAlias() const
    {
      static const std::string s_name( "Loop");
      return s_name;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer Loop::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "scores the loop length between two SSEs");
      serializer.AddInitializer
      (
        "histogram filename",
        "path to where the statistics and energy potentials are read from",
        io::Serialization::GetAgent( &m_HistogramFileName)
       );
      serializer.AddInitializer
      (
        "max loop length",
        "max number of residues in loop",
        io::Serialization::GetAgent( &m_MaxLoopLength)
       );

      return serializer;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a specified histogram file
    //! @param HISTOGRAM_FILENAME filename of the histogram to be used
    //! @param SCHEME scheme to be used
    //! @param MAX_LOOP_LENGTH maximum number of residues in a loop to be considered
    Loop::Loop
    (
      const std::string &HISTOGRAM_FILENAME,
      const std::string &SCHEME,
      const size_t MAX_LOOP_LENGTH
    ) :
      m_HistogramFileName( HISTOGRAM_FILENAME),
      m_Scheme( SCHEME),
      m_MaxLoopLength( MAX_LOOP_LENGTH),
      m_EnergyFunctions()
    {
      // read the histogram file and store the energy functions, except for the static instance
      if( !command::CommandState::GetGlobalCommandState().IsInStaticInitialization())
      {
        ReadEnergyVector();
      }
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new Loop copied from this one
    Loop *Loop::Clone() const
    {
      return new Loop( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Loop::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &Loop::GetScheme() const
    {
      return m_Scheme;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief score the loop length and distance
    //! @param SEQ_EUC_DISTANCE pair of sequence and euclidean distance
    //! @return score for that combination of sequence and euclidean distance
    double Loop::Score( const storage::Pair< size_t, double> &SEQ_EUC_DISTANCE) const
    {
      // if loop is longer than max loop length, return energy for longer loops which is derived from the distance over
      // log of the sequence distance
      if( SEQ_EUC_DISTANCE.First() >= m_MaxLoopLength)
      {
        return m_EnergyFunctionLongLoops( NormalizeDistance( SEQ_EUC_DISTANCE));
      }

      // return the function value for the pair of residue distance and euclidean distance
      return m_EnergyFunctions( SEQ_EUC_DISTANCE.First())( SEQ_EUC_DISTANCE.Second());
    }

    //! @brief calculate Sequence distance between aas and euclidean distance between ends of bodies of two SSEs
    //! @brief SSE_A first SSE of interest
    //! @brief SSE_B first SSE of interest
    //! @return a pair which first member is the length of the loop end aminoacids and the second is the euclidean distance
    storage::Pair< size_t, double>
    Loop::SequenceAndEuclideanDistance
    (
      const assemble::SSE &SSE_A,
      const assemble::SSE &SSE_B
    )
    {
      // sequences of different chains
      if( SSE_A.GetChainID() != SSE_B.GetChainID())
      {
        BCL_MessageVrb( "passing two different chains");
        return storage::Pair< size_t, double>( util::GetUndefined< size_t>(), util::GetUndefined< double>());
      }

      // order pair, so that first's Sequence first AA has the smaller Sequence ID
      storage::VectorND< 2, util::SiPtr< const assemble::SSE> > ordered_pair( SSE_A, SSE_B);

      // if the first SSE comes after the second SSE
      if( ordered_pair.First()->GetFirstAA()->GetSeqID() > ordered_pair.Second()->GetFirstAA()->GetSeqID())
      {
        std::swap( ordered_pair.First(), ordered_pair.Second());
      }

      const size_t seq_distance
      (
        biol::CalculateSequenceDistance( *ordered_pair.Second(), *ordered_pair.First())
      );

      if( !util::IsDefined( seq_distance))
      {
        BCL_MessageCrt( "undefined distance between two sequences");
        return storage::Pair< size_t, double>( util::GetUndefined< size_t>(), util::GetUndefined< double>());
      }

      storage::VectorND< 2, util::SiPtr< const assemble::SSEGeometryInterface> > ordered_geometries
      (
        util::SiPtr< const assemble::SSEGeometryInterface>( ordered_pair.First()),
        util::SiPtr< const assemble::SSEGeometryInterface>( ordered_pair.Second())
      );

      // get the corresponding fragments and store them
      util::SiPtr< const assemble::SSEGeometryInterface> frag_first_sse
      (
        ordered_pair.First()->GetSSEGeometries().IsEmpty() ?
          util::SiPtr< const assemble::SSEGeometryInterface>( ordered_pair.First()) :
          util::SiPtr< const assemble::SSEGeometryInterface>( ordered_pair.First()->GetSSEGeometries().LastElement())
      );
      util::SiPtr< const assemble::SSEGeometryInterface> frag_second_sse
      (
        ordered_pair.Second()->GetSSEGeometries().IsEmpty() ?
          util::SiPtr< const assemble::SSEGeometryInterface>( ordered_pair.Second()) :
          util::SiPtr< const assemble::SSEGeometryInterface>( ordered_pair.Second()->GetSSEGeometries().FirstElement())
      );

      // update the fragments if they are defined
      // otherwise, such as SSEs that are too short to generate a fragment, stick with the SSE
      if( frag_first_sse.IsDefined())
      {
        ordered_geometries.First() = frag_first_sse;
      }
      if( frag_second_sse.IsDefined())
      {
        ordered_geometries.Second() = frag_second_sse;
      }

      // return undefined if either of the bodies is not defined
      if( !ordered_geometries.First()->IsDefined() || !ordered_geometries.Second()->IsDefined())
      {
        BCL_MessageDbg
        (
          "passing at least one undefined SSE for loop end distance calculation " +
          ordered_geometries.First()->GetIdentification() + " vs " + ordered_geometries.Second()->GetIdentification()
        );
        return storage::Pair< size_t, double>( util::GetUndefined< size_t>(), util::GetUndefined< double>());
      }

      return storage::Pair< size_t, double>
      (
        seq_distance,
        ( ordered_geometries.Second()->BeginOfZ() - ordered_geometries.First()->EndOfZ()).Norm()
      );

    }

    //! @brief normalizes the distance by the sequence length
    //! @param SEQ_EUC_DISTANCE pair of sequence and euclidean distance
    //! @return normalized distance
    double Loop::NormalizeDistance
    (
      const storage::Pair< size_t, double> &SEQ_EUC_DISTANCE
    )
    {
      return SEQ_EUC_DISTANCE.Second() / std::log( double( SEQ_EUC_DISTANCE.First() + 2));
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief the loop length and distance between the two given sses
    //! @brief SSE_A first SSE of interest
    //! @brief SSE_B first SSE of interest
    //! @return potential for loop
    double
    Loop::operator()
    (
      const assemble::SSE &SSE_A,
      const assemble::SSE &SSE_B
    ) const
    {
      // check that sses are from same chain - and neither of them is a coil
      if( SSE_A.GetChainID() != SSE_B.GetChainID())
      {
        BCL_MessageVrb( "cannot score loops between different chains!");
        return double( 0.0);
      }

      // check that both sses are helix or strand
      if( SSE_A.GetType() > biol::GetSSTypes().STRAND || SSE_B.GetType() > biol::GetSSTypes().STRAND)
      {
        BCL_MessageVrb( "can only score loops between SSE types HELIX or STRAND!");
        return double( 0.0);
      }

      //calculate the sequence and euclidean distance
      const storage::Pair< size_t, double> seq_euc_distance( SequenceAndEuclideanDistance( SSE_A, SSE_B));

      // if the distance is not defined, return a zero score
      if( !util::IsDefined( seq_euc_distance.First()) || !util::IsDefined( seq_euc_distance.Second()))
      {
        return 0;
      }

      // score this combination of seq and euclidean distance
      return Score( seq_euc_distance);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write the Scheme and the Function value for the ARGUMENT to the STREAM
    //! @param SSE_A first SSE of interest
    //! @param SSE_B first SSE of interest
    //! @param OSTREAM std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &
    Loop::WriteDetailedSchemeAndValues
    (
      const assemble::SSE &SSE_A,
      const assemble::SSE &SSE_B,
      std::ostream &OSTREAM
    ) const
    {
      //calculate sequence and euclidean distance
      const storage::Pair< size_t, double> seq_euc_distance( SequenceAndEuclideanDistance( SSE_A, SSE_B));

      //write sstype, Angel to membrane plane first and last amino acid and value to the STREAM
      OSTREAM << SSE_A.GetIdentification() << '\t'
              << SSE_B.GetIdentification() << '\t'
              << seq_euc_distance.First() << '\t'
              << seq_euc_distance.Second() << '\t'
              << Score( seq_euc_distance) << '\n';

      //end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculate the maximum loop euclidian distance from the histogram files for each loop sequence length within a certain percentile
    //! @param HISTOGRAM_FILE filename of histogram to be used
    //! @param MAX_NR_LOOP_RESIDUES maximum number of loop residues
    //! @param FRACTION fraction of counts that needs to be below this threshold
    storage::Vector< double> Loop::CalculateMaximumObservedDistances
    (
      const std::string &HISTOGRAM_FILENAME,
      const size_t MAX_NR_LOOP_RESIDUES,
      const double FRACTION
    )
    {
      // make sure gievn fraction is valid
      BCL_Assert( FRACTION >= 0.0 && FRACTION <= 1.0, "The given fraction should be in range [0,1]");

      // initialize vector to be returned
      storage::Vector< double> max_distances_vector( MAX_NR_LOOP_RESIDUES + 1, util::GetUndefinedDouble());

      // initialize histograms vector
      storage::Vector< math::Histogram> histograms;

      // open stream and read
      io::IFStream read;
      io::File::MustOpenIFStream( read, HISTOGRAM_FILENAME);
      read >> histograms;

      // clear the stream
      io::File::CloseClearFStream( read);

      // make sure there are enough histograms
      BCL_Assert
      (
        histograms.GetSize() >= MAX_NR_LOOP_RESIDUES + 1,
        "There are only " + util::Format()( histograms.GetSize()) +
        " histograms in the given file, but requested number is " + util::Format()( MAX_NR_LOOP_RESIDUES)
      );

      // iterate over the histograms
      for( size_t hist_no( 0); hist_no <= MAX_NR_LOOP_RESIDUES; ++hist_no)
      {
        BCL_MessageVrb( "looking at loops of length: " + util::Format()( hist_no));

        // create ref on this histogram
        math::Histogram &histogram( histograms( hist_no));

        // get the total counts and calculated the threshold count sum
        const size_t total_counts( histogram.GetSumOfAllCounts());
        const size_t threshold_count( total_counts * FRACTION);
        size_t last_pos( histogram.GetIndexOfLastInformationContainingBin());

        BCL_MessageVrb
        (
          "total_counts: " + util::Format()( total_counts) +
          "\tthreshold_count: " + util::Format()( threshold_count) +
          "\tlast_pos: " + util::Format()( last_pos)
        );

        // initialize counts_sum
        size_t counts_sum( total_counts);

        // while the counts_sum is still larger than the threshold count
        while( counts_sum >= threshold_count && last_pos > 0)
        {
          // decrease the sum of this bin and move one position back
          if( last_pos == histogram.GetNumberOfBins())
          {
            counts_sum -= histogram.GetBoundariesCounts().Second();
          }
          else
          {
            counts_sum -= histogram.GetHistogram()( last_pos);
          }
          --last_pos;
        }

        // calculate the distance for the right boundary for this bin
        const double max_distance
        (
          histogram.GetBoundaries().First() + ( histogram.GetBinSize() * ( last_pos + 2))
        );

        // print the position
        BCL_MessageVrb( "last_pos: " + util::Format()( last_pos) + "distance: " + util::Format()( max_distance));

        max_distances_vector( hist_no) = max_distance;
      }

      // end
      return max_distances_vector;
    }

    //! @brief set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read from this class
    //! @param ERROR_STREAM stream with which to write errors
    bool Loop::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERROR_STREAM
     )
    {
      ReadEnergyVector();
      return true;
    }

    //! @brief read energy distribution for scoring pairs of AASequences
    void Loop::ReadEnergyVector()
    {
      // read file with according histogram2D
      io::IFStream read;
      io::File::MustOpenIFStream( read, Score::AddHistogramPath( m_HistogramFileName));

      // read the histogram
      storage::Vector< math::Histogram> looplength_distance_histogram;
      read >> looplength_distance_histogram;

      // read the histograms for distance over log seq distance
      storage::Vector< math::Histogram> looplength_log_distance_histogram;
      read >> looplength_log_distance_histogram;

      // clear the stream
      io::File::CloseClearFStream( read);

      // create the energy histogram and store it in m_EnergyFunction
      BCL_Assert( m_MaxLoopLength + 1 <= looplength_distance_histogram.GetSize(), "not enough histograms given");

      // only calculate the energy distributions for individual loop lengths up to the max loop length
      looplength_distance_histogram.Resize( m_MaxLoopLength + 1);
      m_EnergyFunctions = EnergyDistribution::LoopLengthDistancePotential( looplength_distance_histogram);

      BCL_Assert( m_MaxLoopLength + 1 < looplength_log_distance_histogram.GetSize(), "the given max loop length is too big");

      math::Histogram combined_hist;
      // iterate over all energy functions for distance over log loop length
      for
      (
        storage::Vector< math::Histogram>::const_iterator
          itr( looplength_log_distance_histogram.Begin() + m_MaxLoopLength + 1), itr_end( looplength_log_distance_histogram.End());
        itr < itr_end;
        ++itr
      )
      {
        math::Histogram current( *itr);
        current.Normalize();
        combined_hist.Combine( current);
      }

      // derive distribution from combined histograms to score loops that are larger than the max loop length
      m_EnergyFunctionLongLoops = EnergyDistribution::EnergyfunctionFromHistogram( combined_hist, 0.0001);
    }

  } // namespace score
} // namespace bcl
