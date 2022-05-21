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
#include "score/bcl_score_loop_angle.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    //! @brief returns default file where loop angle histograms and in consequence the energy potentials are read from
    //! @return default file where the statistics and in consequence the energy potentials are read from
    const std::string &LoopAngle::GetDefaultHistogramFilename()
    {
      static const std::string s_default_histogram_filename( "loop_angle.histograms");
      return s_default_histogram_filename;
    }

    //! @brief returns default file where the loop angle table is stored
    //! @return default file where the loop angle table data is stored
    const std::string &LoopAngle::GetDefaultTableFilename()
    {
      static const std::string s_default_table_filename( "loop_angle.tbl");
      return s_default_table_filename;
    }

    //! @brief returns the maximum sequence distance between SSEs to be considered as consecutive
    //! @return maximum sequence distance between SSEs to be considered as consecutive
    size_t LoopAngle::GetDefaultMaxmimumSequenceDistance()
    {
      static const size_t s_max_sequence_distance( 20); // set according to the cos(angle) cutoff plots
      return s_max_sequence_distance;
    }

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &LoopAngle::GetDefaultScheme()
    {
      static const std::string s_default_scheme( "loop_angle");
      return s_default_scheme;
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return  the name of the object when used in a dynamic context
    const std::string &LoopAngle::GetAlias() const
    {
      static const std::string s_name( "LoopAngle");
      return s_name;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer LoopAngle::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "scores the loop angle between two SSEs");
      serializer.AddInitializer
      (
        "histogram filename",
        "file from which statistics and in consequence energy potentials are read",
        io::Serialization::GetAgent( &m_HistogramFileName)
      );
      serializer.AddInitializer
      (
        "max sequence distance",
        "maximum sequence distance between SSEs to be considered consecutive",
        io::Serialization::GetAgent( &m_MaxSequenceDistance)
      );

      return serializer;
    }

    //! @brief constructor from a specified histogram file and scheme
    //! @param HISTOGRAM_FILENAME filename of the histogram to be used
    //! @param MAX_SEQ_DISTANCE maximum sequence distance between SSEs to be considered as consecutive
    //! @param SCHEME scheme to be used
    LoopAngle::LoopAngle
    (
      const std::string &HISTOGRAM_FILENAME,
      const size_t &MAX_SEQ_DISTANCE,
      const std::string &SCHEME
    ) :
      m_HistogramFileName( HISTOGRAM_FILENAME),
      m_MaxSequenceDistance( MAX_SEQ_DISTANCE),
      m_Scheme( SCHEME),
      m_EnergyFunctionShortLoops(),
      m_EnergyFunctionLongLoops()
    {
      ReadEnergyVector();
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new Loop copied from this one
    LoopAngle *LoopAngle::Clone() const
    {
      return new LoopAngle( *this);
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LoopAngle::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the histogram filename
    //! @return the histogram filename
    const std::string &LoopAngle::GetHistogramFileName() const
    {
      return m_HistogramFileName;
    }

    //! @brief gets the maximum sequence distance for which the two SSEs are considered consecutive for this score
    //! @return the maximum sequence distance considered
    const size_t &LoopAngle::GetMaxSequenceDistance() const
    {
      return m_MaxSequenceDistance;
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &LoopAngle::GetScheme() const
    {
      return m_Scheme;
    }

    //! @brief get energy functions
    //! @return VectorND of short and long loop cubic spline
    storage::VectorND< 2, math::CubicSplineDamped> LoopAngle::GetEnergyFunctions() const
    {
      return storage::VectorND< 2, math::CubicSplineDamped>( m_EnergyFunctionShortLoops, m_EnergyFunctionLongLoops);
    }

    //! @brief get score type
    //! @return score type
    ProteinModel::Type LoopAngle::GetType() const
    {
      return ProteinModel::e_Structure;
    }

    //! @brief returns the angle between the given SSEs respective to the given position
    //! @param SSE_FIRST first SSE of interest
    //! @param SSE_SECOND first SSE of interest
    //! @param CENTER_OF_MASS the center of mass of the protein the SSEs are in
    //! @return angle between the given SSEs respective to the given position
    double LoopAngle::CalculateCosAngle
    (
      const assemble::SSE &SSE_FIRST,
      const assemble::SSE &SSE_SECOND,
      const linal::Vector3D &CENTER_OF_MASS
    ) const
    {
      // the angle is calculated between the end points of the loop
      const linal::Vector3D prev_end_of_z( SSE_FIRST.EndOfZ());
      const linal::Vector3D next_begin_of_z( SSE_SECOND.BeginOfZ());

      return linal::ProjAngleCosinus( CENTER_OF_MASS, prev_end_of_z, next_begin_of_z);
    }

    //! @brief the loop angle between the two given SSEs
    //! @param PROTEIN_MODEL the protein model to be scored
    //! @return score for the angle between consecutive SSEs
    double LoopAngle::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      std::stringstream dummy_stream;
      return ScoreLoops( PROTEIN_MODEL, dummy_stream, false);
    }

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &LoopAngle::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_HistogramFileName, ISTREAM);
      io::Serialize::Read( m_MaxSequenceDistance, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // read the histogram file and store the energy functions
      ReadEnergyVector();

      return ISTREAM;
    }

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @param INDENT indentation
    //! @return returns the output stream
    std::ostream &LoopAngle::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_HistogramFileName, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MaxSequenceDistance, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      return OSTREAM;
    }

    //! @brief write detailed scheme and values to OSTREAM
    //! @param MODEL protein models to be used to evaluate the function
    //! @param OSTREAM std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &LoopAngle::WriteDetailedSchemeAndValues
    (
      const assemble::ProteinModel &MODEL,
      std::ostream &OSTREAM
    ) const
    {
      ScoreLoops( MODEL, OSTREAM, true);
      return OSTREAM;
    }

    //! @brief set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @return ERROR_STREAM stream with which to write errors
     bool LoopAngle::ReadInitializeSuccessHook
     (
       const util::ObjectDataLabel &LABEL,
       std::ostream &ERROR_STREAM
     )
    {
      ReadEnergyVector();
      return true;
    }

    //! @brief read energy distribution for scoring pairs of AASequences
    void LoopAngle::ReadEnergyVector()
    {
      // read histogram
      io::IFStream read;
      io::File::MustOpenIFStream( read, Score::AddHistogramPath( m_HistogramFileName));
      math::Histogram loop_angle_short_loops_histogram, loop_angle_long_loops_histogram;
      read >> loop_angle_short_loops_histogram;
      read >> loop_angle_long_loops_histogram;
      io::File::CloseClearFStream( read);

      // create energy function from histogram
      m_EnergyFunctionShortLoops = EnergyDistribution::EnergyfunctionFromHistogram( loop_angle_short_loops_histogram);
      m_EnergyFunctionLongLoops = EnergyDistribution::EnergyfunctionFromHistogram( loop_angle_long_loops_histogram);
    }

    //! @brief helper function called by WriteDetailedSchemeAndValues and operator() so that the code remains in sync
    //! @param MODEL the protein model to be scored
    //! @param OSTREAM the output stream to write the detailed scheme to for this chain
    //! @param DO_WRITE set to true to actually write to the output stream; otherwise, nothing will be written
    //! @return the final score
    double LoopAngle::ScoreLoops( const assemble::ProteinModel &MODEL, std::ostream &OSTREAM, const bool &WRITE) const
    {
      // the center of mass is used as reference point for the angle calculation
      const linal::Vector3D center_of_mass( MODEL.GetCenterOfMass());

      // only select helices and strands from the model since the coils are the loops
      storage::Set< biol::SSType> sse_types( biol::GetSSTypes().GetHelixTypes());
      sse_types.Insert( biol::GetSSTypes().STRAND);

      // iterate over the SSEs and score the angle for consecutive ones
      const util::SiPtrVector< const assemble::SSE> sses( MODEL.GetSSEs( sse_types));
      double score( 0.0);
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator itr( sses.Begin()), itr_end( sses.End());
        itr != itr_end && itr + 1 != itr_end;
        ++itr
      )
      {
        // if the sequence distance between the current SSEs exceeds the provided maximum distance,
        // disregard this iteration
        const assemble::SSE &sse_a( **itr);
        const assemble::SSE &sse_b( **( itr + 1));
        const int seq_id_a( sse_a.GetLastAA()->GetSeqID());
        const int seq_id_b( sse_b.GetFirstAA()->GetSeqID());
        const size_t sequence_distance( seq_id_b - seq_id_a);

        // score the angle between the current SSE and it's successor
        const double cos_angle( CalculateCosAngle( sse_a, sse_b, center_of_mass));
        score += sequence_distance <= m_MaxSequenceDistance ?
          m_EnergyFunctionShortLoops( cos_angle) : m_EnergyFunctionLongLoops( cos_angle);

        if( WRITE)
        {
          OSTREAM << "Loop between " << sse_a.GetIdentification() << " and " << sse_b.GetIdentification()
            << " sequence_distance=" << sequence_distance << " cos(angle)="
            << cos_angle << " cumulative_score=" << score << "\n";
        }
      }

      // normalize the summed up score by the number of sses
      score /= sses.GetSize();
      if( WRITE)
      {
        OSTREAM << "Normalized score=" << score << " by number_of_sse=" << sses.GetSize();
      }

      return score;
    }
  } // namespace score
} // namespace bcl
