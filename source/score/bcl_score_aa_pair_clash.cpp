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
#include "score/bcl_score_aa_pair_clash.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_voxel_grid_aa.h"
#include "biol/bcl_biol_aa_base.h"
#include "command/bcl_command_command_state.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_histogram.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! @brief returns default file where the statistics and in consequence the energy potentials are read from
    //! @return default file where the statistics and in consequence the energy potentials are read from
    const std::string &AAPairClash::GetDefaultHistogramFilename()
    {
      // static string
      static const std::string s_default_histogram_filename( "aa_distances_0.05.histograms");

      // end
      return s_default_histogram_filename;
    }

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &AAPairClash::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "aaclash");

      // end
      return s_default_scheme;
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AAPairClash::s_Instance
    (
      util::Enumerated< AAPairDistanceInterface>::AddInstance( new AAPairClash)
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a specified histogram file
    //! @param SIGMOID_WIDTH width of the sigmoid repulsive term, before it reaches 1.0
    //! @param HISTOGRAM_FILENAME filename of the histogram to be used
    //! @param SCHEME scheme to be used
    AAPairClash::AAPairClash
    (
      const double SIGMOID_WIDTH, // = 1.0,
      const std::string &HISTOGRAM_FILENAME, // = GetDefaultHistogramFilename(),
      const std::string &SCHEME // = GetDefaultScheme()
    ) :
      m_ShortestObservedDistance(),
      m_ShortestObservedDistanceMatrix( biol::AATypes::s_NumberStandardAATypes, biol::AATypes::s_NumberStandardAATypes, double( 0.0)),
      m_SigmoidWidth( SIGMOID_WIDTH),
      m_HistogramFileName( HISTOGRAM_FILENAME),
      m_Scheme( SCHEME),
      m_DistanceCutoff( 0.0)
    {
      // read the histogram file and store the distances
      ReadDistanceMap();
    }

    //! @brief constructor from a specified histogram file
    //! @param SIGMOID_WIDTH width of the sigmoid repulsive term, before it reaches 1.0
    //! @param MINIMAL_SEQUENCE_SEPARATION minimal sequence separation
    //! @param HISTOGRAM_FILENAME filename of the histogram to be used
    //! @param SCHEME scheme to be used
    AAPairClash::AAPairClash
    (
      const double SIGMOID_WIDTH,
      const size_t MINIMAL_SEQUENCE_SEPARATION,
      const std::string &HISTOGRAM_FILENAME,
      const std::string &SCHEME
    ) :
      m_ShortestObservedDistance(),
      m_ShortestObservedDistanceMatrix( biol::AATypes::s_NumberStandardAATypes, biol::AATypes::s_NumberStandardAATypes, double( 0.0)),
      m_SigmoidWidth( SIGMOID_WIDTH),
      m_HistogramFileName( HISTOGRAM_FILENAME),
      m_Scheme( SCHEME),
      m_DistanceCutoff( 0.0)
    {
      // read the histogram file and store the distances
      ReadDistanceMap();
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new AAPairClash object that is copied from this one
    AAPairClash *AAPairClash::Clone() const
    {
      return new AAPairClash( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AAPairClash::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief create a matrix of all amino acid pair min distances
    //! @return matrix with minimal distances observed for every amino acid pair
    const linal::Matrix< double> &AAPairClash::GetShortestObservedDistanceMatrix() const
    {
      // end
      return m_ShortestObservedDistanceMatrix;
    }

  ///////////////
  // operators //
  ///////////////

    //! Get the closest observed distance between two AAs in the PDB
    double AAPairClash::GetClosestDistance( const biol::AAType &A, const biol::AAType &B) const
    {
      return A->IsNaturalAminoAcid() && B->IsNaturalAminoAcid() ? m_ShortestObservedDistanceMatrix( A, B) : 0.0;
    }

    //! @brief calculate amino acid pairing potential for given amino acid pair
    //! @param AMINO_ACID_A first amino acid of interest
    //! @param AMINO_ACID_B second amino acid of interest
    //! @return amino acid pairing potential for given amino acid pair
    double AAPairClash::operator()
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B
    ) const
    {
      // calculate the CB Distance
      const double distance
      (
        biol::FirstSidechainAtomDistance( AMINO_ACID_A, AMINO_ACID_B)
      );

      // return the clash score
      return operator()( AMINO_ACID_A, AMINO_ACID_B, distance);
    }

    //! @brief calculate amino acid pairing potential for given amino acid pair
    //! @param AMINO_ACID_A first amino acid of interest
    //! @param AMINO_ACID_B second amino acid of interest
    //! @param DISTANCE distance between the amino acid pair
    //! @return amino acid pairing potential for given amino acid pair
    double AAPairClash::operator()
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B,
      const double DISTANCE
    ) const
    {
      // if distance is not defined return 0
      if( !util::IsDefined( DISTANCE))
      {
        return double( 0.0);
      }

      if( AMINO_ACID_A.GetType()->IsNaturalAminoAcid() && AMINO_ACID_B.GetType()->IsNaturalAminoAcid())
      {
        if( AMINO_ACID_A.GetSeqID() != AMINO_ACID_B.GetSeqID() || AMINO_ACID_A.GetChainID() != AMINO_ACID_B.GetChainID())
        {
          return CalculateRepulsiveTerm
                 (
                   DISTANCE,
                   m_ShortestObservedDistanceMatrix( AMINO_ACID_A.GetType(), AMINO_ACID_B.GetType())
                 );
        }
      }

      // shortest observed distance for this aa pair type
      const storage::Map< storage::Pair< biol::AAType, biol::AAType>, double>::const_iterator
      itr
      (
        m_ShortestObservedDistance.Find
        (
          storage::Pair< biol::AAType, biol::AAType>( AMINO_ACID_A.GetType(), AMINO_ACID_B.GetType())
        )
      );

      if( itr == m_ShortestObservedDistance.End())
      {
        return double( 0.0);
      }

      // calculate the actual repulsive term
      return CalculateRepulsiveTerm( DISTANCE, itr->second);
    }

    //! @brief calculate clash score for a protein model
    //! @param MODEL the protein model of interest
    //! @return amino acid pairing potential for given protein
    double AAPairClash::operator()( const assemble::ProteinModel &MODEL) const
    {
      assemble::VoxelGridAA s_voxel_grid( m_DistanceCutoff);

      auto aas( MODEL.GetAminoAcids());
      s_voxel_grid.SetObjects( aas);
      double score( 0.0);

      storage::Vector< storage::Triplet< util::SiPtr< const biol::AABase>, util::SiPtr< const biol::AABase>, double> >
        res( s_voxel_grid.GetNeighbors( m_DistanceCutoff));
      for( auto itr( res.Begin()), itr_end( res.End()); itr != itr_end; ++itr)
      {
        const auto &triplet( *itr);
        const double closest_distance( GetClosestDistance( triplet.First()->GetType(), triplet.Second()->GetType()));
        score += CalculateRepulsiveTerm( triplet.Third(), closest_distance);
      }
      return score;
    }

    //! @brief calculate clash score for SSEs
    //! @param SSE_A, SSE_B the SSEs to check for clashes
    //! @return clash score
    double AAPairClash::operator()( const assemble::SSE &SSE_A, const assemble::SSE &SSE_B) const
    {
      assemble::VoxelGridAA voxel_grid_a( m_DistanceCutoff);
      voxel_grid_a.SetObjects( util::SiPtrVector< const biol::AABase>( SSE_A.Begin(), SSE_A.End()));
      double score( 0.0);
      for( auto itr_b( SSE_B.Begin()), itr_b_end( SSE_B.End()); itr_b != itr_b_end; ++itr_b)
      {
        storage::Vector< storage::Pair< util::SiPtr< const biol::AABase>, double> >
          res( voxel_grid_a.GetNeighbors( **itr_b, m_DistanceCutoff));
        for( auto itr( res.Begin()), itr_end( res.End()); itr != itr_end; ++itr)
        {
          const auto &pr( *itr);
          const double closest_distance( GetClosestDistance( ( *itr_b)->GetType(), pr.First()->GetType()));
          score += CalculateRepulsiveTerm( pr.Second(), closest_distance);
        }
      }
      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
    //! @param AMINO_ACID_A first amino acid of interest
    //! @param AMINO_ACID_B second amino acid of interest
    //! @param OSTREAM the std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &
    AAPairClash::WriteDetailedSchemeAndValues
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B,
      std::ostream &OSTREAM
    ) const
    {
      // calculate the aa distance
      const double aa_distance( biol::FirstSidechainAtomDistance( AMINO_ACID_A, AMINO_ACID_B));

      // if distance is not defined return 0
      if( !util::IsDefined( aa_distance))
      {
        return OSTREAM;
      }

      // write Scheme
      OSTREAM << AMINO_ACID_A.GetSeqID() << '\t'
              << AMINO_ACID_A.GetType()->GetThreeLetterCode() << '\t'
              << AMINO_ACID_B.GetSeqID() << '\t'
              << AMINO_ACID_B.GetType()->GetThreeLetterCode() << '\t'
              << aa_distance << '\t'
              << operator()( AMINO_ACID_A, AMINO_ACID_B) << '\n';

      // end
      return OSTREAM;
    }

    //! @brief return parameters for data members that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AAPairClash::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Scores whether two amino acids clash");
      serializer.AddInitializer
      (
        "sigmoid width",
        "width of the sigmoidal function to be used",
        io::Serialization::GetAgent( &m_SigmoidWidth),
        "1.0"
      );
      serializer.AddInitializer
      (
        "histogram file name",
        "path to file where the statistics and in consequence the energy potentials are read from",
        io::Serialization::GetAgent( &m_HistogramFileName),
        GetDefaultHistogramFilename()
      );

      return serializer;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief repulsive term from distance and shortest observed distance
    //! @param DISTANCE the actual distance between amino acids of interest
    //! @param SHORTEST_OBSERVED_DISTANCE shortest distance observed
    double AAPairClash::CalculateRepulsiveTerm
    (
      const double DISTANCE,
      const double SHORTEST_OBSERVED_DISTANCE
    ) const
    {
      const double difference( SHORTEST_OBSERVED_DISTANCE - DISTANCE);

      // no repulsion
      if( difference <= 0.0)
      {
        return 0.0;
      }
      // full repulsion
      else if( difference >= m_SigmoidWidth)
      {
        return 1.0;
      }

      return math::WeightBetweenZeroAndPi( ( ( m_SigmoidWidth - difference) / m_SigmoidWidth) * math::g_Pi);
    }

    //! @brief read map of amino acid pair observed distanced from histogram file
    void AAPairClash::ReadDistanceMap()
    {
      m_DistanceCutoff = 0.0;
      // during static initialization, this histogram path is not available. Moreover, the static instance does not need
      // the histograms available
      if( command::CommandState::GetGlobalCommandState().IsInStaticInitialization())
      {
        return;
      }
      // read file with all histograms for each pair of aa types
      io::IFStream read;
      io::File::MustOpenIFStream( read, Score::AddHistogramPath( m_HistogramFileName));

      // initialize temporary strings to read
      std::string tmp_a, tmp_b;

      // number of natural pairs seen
      size_t n_natural_seen( 0);

      // while reading the aatype pair and until reaching the end of the file
      while( read >> tmp_a >> tmp_b && !read.eof())
      {
        // read the string and convert it to first aa type
        storage::Pair< biol::AAType, biol::AAType> aa_type_pair
        (
          biol::GetAATypes().AATypeFromOneLetterCode( tmp_a[0]),
          biol::GetAATypes().AATypeFromOneLetterCode( tmp_b[0])
        );

        // abort if any of the aatypes if unknown
        if
        (
          aa_type_pair.First() == biol::GetAATypes().e_Undefined ||
          aa_type_pair.Second() == biol::GetAATypes().e_Undefined
        )
        {
          // alert user and break
          BCL_MessageCrt
          (
            "undefined AAType found in the histogram " + util::Format()( aa_type_pair)
          );
          break;
        }

        // read the histogram
        math::Histogram current_histogram;
        read >> current_histogram;

        // identify the first non 0 bin
        const size_t shortest_observed_index( current_histogram.GetIndexOfFirstInformationContainingBin( 1.5));
        const double shortest_observed_distance
        (
          shortest_observed_index * current_histogram.GetBinSize() + current_histogram.GetBoundaries().First()
        );

        if( aa_type_pair.First()->IsNaturalAminoAcid() && aa_type_pair.Second()->IsNaturalAminoAcid())
        {
          m_ShortestObservedDistanceMatrix( aa_type_pair.First(), aa_type_pair.Second()) = shortest_observed_distance;
          m_ShortestObservedDistanceMatrix( aa_type_pair.Second(), aa_type_pair.First()) = shortest_observed_distance;
          ++n_natural_seen;
        }
        else
        {
          // insert into map
          m_ShortestObservedDistance[ aa_type_pair] = shortest_observed_distance;
          std::swap( aa_type_pair.First(), aa_type_pair.Second());
          m_ShortestObservedDistance[ aa_type_pair] = shortest_observed_distance;
        }

        // update distance cutoff
        m_DistanceCutoff = std::max( m_DistanceCutoff, shortest_observed_distance);
      }
      BCL_Assert( n_natural_seen == size_t( 210), "Not all natural AAs were in the distance map!");

      // close the stream
      io::File::CloseClearFStream( read);
    }

    //! @brief set the members of this object from the given label
    //! @param LABEL the label containing members that should be read of this class
    //! @return ERROR_STREAM stream with which to write errors
    bool AAPairClash::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      ReadDistanceMap();
      return true;
    }

    //! @brief get the static instance of this class
    const AAPairClash &AAPairClash::GetInstance()
    {
      static AAPairClash s_aa_clash;
      return s_aa_clash;
    }
  } // namespace score
} // namespace bcl
