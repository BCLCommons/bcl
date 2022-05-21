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
#include "score/bcl_score_aa_pair_distance.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_voxel_grid_aa.h"
#include "biol/bcl_biol_aa_base.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_gnuplot_heatmap.h"
#include "math/bcl_math_histogram.h"
#include "score/bcl_score_energy_distribution.h"
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
    const std::string &AAPairDistance::GetDefaultHistogramFilename()
    {
      // static string
      static const std::string s_default_histogram_filename( "aa_distances.histograms");

      // end
      return s_default_histogram_filename;
    }

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &AAPairDistance::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "aadist");

      // end
      return s_default_scheme;
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &AAPairDistance::GetAlias() const
    {
      static const std::string s_name( "AAPairDistance");
      return s_name;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a specified histogram file
    //! @param HISTOGRAM_FILENAME filename of the histogram to be used
    //! @param SCHEME scheme to be used
    AAPairDistance::AAPairDistance
    (
      const std::string &HISTOGRAM_FILENAME,
      const std::string &SCHEME
    ) :
      m_HistogramFileName( HISTOGRAM_FILENAME),
      m_Scheme( SCHEME),
      m_EnergyFunctionMap(),
      m_DistanceCutoff( 0.0)
    {
      // read the histogram file and store the energy functions
      ReadEnergyFunctionMap();
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new AAPairDistance object that is copied from this one
    AAPairDistance *AAPairDistance::Clone() const
    {
      return new AAPairDistance( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AAPairDistance::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate amino acid pairing potential for given amino acid pair
    //! @param AMINO_ACID_A first amino acid of interest
    //! @param AMINO_ACID_B second amino acid of interest
    //! @return amino acid pairing potential for given amino acid pair
    double AAPairDistance::operator()
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B
    ) const
    {
      // calculate the distance
      const double distance( biol::FirstSidechainAtomDistance( AMINO_ACID_A, AMINO_ACID_B));

      // return the score
      return operator()( AMINO_ACID_A, AMINO_ACID_B, distance);
    }

    //! @brief calculate amino acid pairing potential for given amino acid pair
    //! @param AMINO_ACID_A first amino acid of interest
    //! @param AMINO_ACID_B second amino acid of interest
    //! @param DISTANCE distance between the amino acid pair
    //! @return amino acid pairing potential for given amino acid pair
    double AAPairDistance::operator()
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B,
      const double DISTANCE
    ) const
    {
      // check that distance is defined
      if( !util::IsDefined( DISTANCE))
      {
        return 0.0;
      }

      // construct biol::AAType pair
      storage::Pair< biol::AAType, biol::AAType> type_pair( AMINO_ACID_A.GetType(), AMINO_ACID_B.GetType());

      // search the map for this aa type pair and store the iterator
      storage::Map
      <
        storage::Pair< biol::AAType, biol::AAType>,
        util::ShPtr< math::CubicSplineDamped>
      >::const_iterator
      itr_find
      (
        m_EnergyFunctionMap.Find( type_pair)
      );

      // if the itr is not valid
      if( itr_find == m_EnergyFunctionMap.End())
      {
        // return undefined
        return double( 0);
      }

      // now call the scoring function for the found energy function
      return itr_find->second->operator()( DISTANCE);
    }

    //! @brief calculate clash score for a protein model
    //! @param MODEL the protein model of interest
    //! @return amino acid pairing potential for given protein
    double AAPairDistance::operator()( const assemble::ProteinModel &MODEL) const
    {
      assemble::VoxelGridAA s_voxel_grid( m_DistanceCutoff);

      s_voxel_grid.SetObjects( MODEL.GetAminoAcids());
      double score( 0.0);

      storage::Vector< storage::Triplet< util::SiPtr< const biol::AABase>, util::SiPtr< const biol::AABase>, double> >
        res( s_voxel_grid.GetNeighbors( m_DistanceCutoff));
      for( auto itr( res.Begin()), itr_end( res.End()); itr != itr_end; ++itr)
      {
        score += operator()( *itr->First(), *itr->Second(), itr->Third());
      }
      return score;
    }

    //! @brief calculate clash score for SSEs
    //! @param SSE_A, SSE_B the SSEs to check for clashes
    //! @return clash score
    double AAPairDistance::operator()( const assemble::SSE &SSE_A, const assemble::SSE &SSE_B) const
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
          score += this->operator()( **itr_b, *pr.First(), pr.Second());
        }
      }
      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &AAPairDistance::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_HistogramFileName, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // read the histogram file and store the energy functions
      ReadEnergyFunctionMap();

      // end
      return ISTREAM;
    }

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @param INDENT indentation
    //! @return returns the output stream
    std::ostream &AAPairDistance::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_HistogramFileName, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
    //! @param AMINO_ACID_A first amino acid of interest
    //! @param AMINO_ACID_B second amino acid of interest
    //! @param OSTREAM the std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &
    AAPairDistance::WriteDetailedSchemeAndValues
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B,
      std::ostream &OSTREAM
    ) const
    {
      // calculate the CB Distance
      const double cb_distance( biol::FirstSidechainAtomDistance( AMINO_ACID_A, AMINO_ACID_B));

      // if distance is not defined return 0
      if( !util::IsDefined( cb_distance))
      {
        return OSTREAM;
      }

      // write Scheme
      OSTREAM << AMINO_ACID_A.GetSeqID() << '\t'
              << AMINO_ACID_A.GetType()->GetThreeLetterCode() << '\t'
              << AMINO_ACID_B.GetSeqID() << '\t'
              << AMINO_ACID_B.GetType()->GetThreeLetterCode() << '\t'
              << cb_distance << '\t'
              << operator()( AMINO_ACID_A, AMINO_ACID_B) << '\n';

      // end
      return OSTREAM;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AAPairDistance::GetSerializer() const
    {
      io::Serializer serializer;

      serializer.SetClassDescription( "Scores amino acid pair distances.");
      serializer.AddInitializer
      (
        "histogram filename",
        "path to file where the statistics and in consequence the energy potentials are read from",
        io::Serialization::GetAgent( &m_HistogramFileName),
        GetDefaultHistogramFilename()
      );
      serializer.AddInitializer
      (
        "distance cutoff",
        "distance cutoff above which score will always be 0",
        io::Serialization::GetAgent( &m_DistanceCutoff),
        "0.0"
      );

      return serializer;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief read map of amino acid pair energies based on distance from histogram files
    void AAPairDistance::ReadEnergyFunctionMap()
    {
      // reset distance cutoff to 0
      m_DistanceCutoff = 0.0;

      // store map of all aa pairs with their histogram
      storage::Map
      <
        storage::Pair< biol::AAType, biol::AAType>,
        math::Histogram
      > histogram_map;

      // read file with all histograms for each pair of aa types
      io::IFStream read;
      io::File::MustOpenIFStream( read, Score::AddHistogramPath( m_HistogramFileName));

      // initialize temporary strings to read
      std::string tmp_a, tmp_b;

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
        read >> histogram_map[ aa_type_pair];
      }

      // keep track of sum of normalized distributions of first aa with all other aas and second aa with all other aas
      storage::Map
      <
        storage::Pair< biol::AAType, biol::AAType>,
        math::Histogram
      > histogram_sum_map;

      for
      (
        storage::Map
        <
          storage::Pair< biol::AAType, biol::AAType>,
          math::Histogram
        >::const_iterator itr( histogram_map.Begin()), itr_end( histogram_map.End());
        itr != itr_end;
        ++itr
      )
      {
        // normalize the histogram before combining, to remove bias to certain aa types
        math::Histogram current_histogram( itr->second);
        current_histogram.Normalize();

        // sum all distributions for each amino acid pair to all other amino acid pairs - for each pair this should add up to a sum to 40 histogram combinations
        for
        (
          biol::AATypes::const_iterator
            aa_type_itr( biol::GetAATypes().Begin()), aa_type_itr_end( biol::GetAATypes().VAL.GetIterator());
          aa_type_itr <= aa_type_itr_end;
          ++aa_type_itr
        )
        {
          // sum for the first amino acid to all other
          if( *aa_type_itr <= itr->first.First())
          {
            BCL_Assert
            (
              ( histogram_sum_map[ storage::Pair< biol::AAType, biol::AAType>( *aa_type_itr, itr->first.First())]).Combine( current_histogram),
              "unable to combine histograms of different parameters"
            );
          }
          else
          {
            BCL_Assert
            (
              ( histogram_sum_map[ storage::Pair< biol::AAType, biol::AAType>( itr->first.First(), *aa_type_itr)]).Combine( current_histogram),
              "unable to combine histograms of different parameters"
            );
          }

          // sum for the second amino acid to all other
          if( *aa_type_itr <= itr->first.Second())
          {
            BCL_Assert
            (
              ( histogram_sum_map[ storage::Pair< biol::AAType, biol::AAType>( *aa_type_itr, itr->first.Second())]).Combine( current_histogram),
              "unable to combine histograms of different parameters"
            );
          }
          else
          {
            BCL_Assert
            (
              ( histogram_sum_map[ storage::Pair< biol::AAType, biol::AAType>( itr->first.Second(), *aa_type_itr)]).Combine( current_histogram),
              "unable to combine histograms of different parameters"
            );
          }
        }
      }

      // normalize the histogram sum
      for
      (
        storage::Map
        <
          storage::Pair< biol::AAType, biol::AAType>,
          math::Histogram
        >::iterator itr( histogram_sum_map.Begin()), itr_end( histogram_sum_map.End());
        itr != itr_end;
        ++itr
      )
      {
        itr->second.Normalize();
      }

      // write background distributions
      if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
      {
        util::SiPtrVector< const math::Histogram> histograms;
        storage::Vector< std::string> histograms_descriptors;

        for
        (
          storage::Map
          <
            storage::Pair< biol::AAType, biol::AAType>,
            math::Histogram
          >::const_iterator itr( histogram_sum_map.Begin()), itr_end( histogram_sum_map.End());
          itr != itr_end;
          ++itr
        )
        {
          // skip redundant pairs
          if( itr->first.First() > itr->first.Second())
          {
            continue;
          }
          histograms.PushBack( itr->second);
          histograms_descriptors.PushBack( itr->first.First()->GetThreeLetterCode() + " " + itr->first.Second()->GetThreeLetterCode());
        }

        // write splines to gnuplot file
        io::OFStream write;
        io::File::MustOpenOFStream( write, "aa_pair_distance_background.gnuplot");
        math::GnuplotHeatmap heatmap;
        heatmap.SetFromHistograms( histograms, false, false);
        heatmap.SetTicsY( histograms_descriptors, true, 1);
        heatmap.SetPixelAndRatio( 1080, 8000, util::GetUndefined< double>());
        heatmap.SetTitleAndLabel( "amino acid pair distance background", "distance [" + math::GnuplotHeatmap::s_AngstromSymbolGnuplot + "]", "amino acid pair type", "SUM(p)");
        heatmap.SetFont( "arialbd", 16);
        heatmap.SetFilename( "aa_pair_distance_background");
        heatmap.WriteScript( write);
        io::File::CloseClearFStream( write);
      }

      // generate energy distribution for all possible amino acid pairs
      for
      (
        biol::AATypes::const_iterator
          aa_type_itr1( biol::GetAATypes().Begin()), aa_type_itr_end( biol::GetAATypes().VAL.GetIterator());
        aa_type_itr1 <= aa_type_itr_end;
        ++aa_type_itr1
      )
      {
        for
        (
          biol::AATypes::const_iterator aa_type_itr2( aa_type_itr1);
          aa_type_itr2 <= aa_type_itr_end;
          ++aa_type_itr2
        )
        {
          const storage::Pair< biol::AAType, biol::AAType> current_aa_pair( *aa_type_itr1, *aa_type_itr2);

          // create spline for the current distribution
          util::ShPtr< math::CubicSplineDamped> current_spline
          (
            new math::CubicSplineDamped
            (
              EnergyDistribution::AAPairPotential( histogram_map[ current_aa_pair], histogram_sum_map[ current_aa_pair], true)
//              EnergyDistribution::AAPairPotential( histogram_map[ current_aa_pair], histogram_sum_map[ current_aa_pair], false)
            )
          );

          // store the spline in the map
          m_EnergyFunctionMap[ storage::Pair< biol::AAType, biol::AAType>( *aa_type_itr1, *aa_type_itr2)] = current_spline;

          // reverse the pair and insert it
          m_EnergyFunctionMap[ storage::Pair< biol::AAType, biol::AAType>( *aa_type_itr2, *aa_type_itr1)] = current_spline;
        }
      }
      m_DistanceCutoff = 8.0;

      // close the stream
      io::File::CloseClearFStream( read);
    }

    //! @brief set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @return ERROR_STREAM stream with which to write errors
    bool AAPairDistance::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERROR_STREAM
    )
    {
      ReadEnergyFunctionMap();
      return true;
    }

  } // namespace score
} // namespace bcl
