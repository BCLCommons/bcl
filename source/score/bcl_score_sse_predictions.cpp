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
#include "score/bcl_score_sse_predictions.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_running_average_sd.h"
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
    const util::SiPtr< const util::ObjectInterface> SSEPredictions::s_Instance
    (
      util::Enumerated< SSEPredictionInterface>::AddInstance( new SSEPredictions())
    );

    //! @brief returns default file where the statistics and in consequence the energy potentials are read from
    //! @return default file where the statistics and in consequence the energy potentials are read from
    const std::string &SSEPredictions::GetDefaultHistogramFilename()
    {
      // static string
      static const std::string s_default_histogram_filename( "sspred_sstype_nres_average_mean.bcl");

      // end
      return s_default_histogram_filename;
    }

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &SSEPredictions::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "ssepred");

      // end
      return s_default_scheme;
    }

    //! @brief multiple of the confidence interval for the z-score
    const double SSEPredictions::s_DefaultConfidenceThreshold( 0.5); // mean and higher probability sses should have a negative score

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SSEPredictions::SSEPredictions() :
      m_Method( sspred::GetMethods().e_Undefined),
      m_ConfidenceThreshold( s_DefaultConfidenceThreshold),
      m_Scheme( GetDefaultScheme())
    {
    }

    //! @brief constructor from a ssmethod and confidence threshold
    //! @param SS_METHOD ssmethod to use
    //! @param CONFIDENCE_THRESHOLD the threshold in units of z-score, above which the score becomes negative
    //! @param SCHEME scheme to be used
    SSEPredictions::SSEPredictions
    (
      const sspred::Method &SS_METHOD,
      const double CONFIDENCE_THRESHOLD,
      const std::string &SCHEME
    ) :
      m_Method( SS_METHOD),
      m_ConfidenceThreshold( CONFIDENCE_THRESHOLD),
      m_Scheme( SCHEME)
    {
      if( m_Method.IsDefined())
      {
        // read the histogram file and store the energy functions
        ReadEnergyVector();
      }
    }

    //! @brief Clone function
    //! @return pointer to new SSEPredictions
    SSEPredictions *SSEPredictions::Clone() const
    {
      return new SSEPredictions( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSEPredictions::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &SSEPredictions::GetScheme() const
    {
      return m_Scheme;
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &SSEPredictions::GetAlias() const
    {
      static const std::string s_name( "SSEPredictions");
      return s_name;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SSEPredictions::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Score the log likelihodd of each amino acids prediction");
      serializer.AddInitializer
      (
        "method",
        "sspred method use in SSE prediction evaluation",
        io::Serialization::GetAgent( &m_Method),
        sspred::GetMethods().e_Undefined
      );
      serializer.AddInitializer
      (
        "confidence threshold",
        "In units of z-score, at which z-score, the confidence is negative ?",
        io::Serialization::GetAgent( &m_ConfidenceThreshold),
        util::Format()( s_DefaultConfidenceThreshold)
      );
      return serializer;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that calculates the score for a given SSE
    //! @param SSE SSE of interest
    //! @param MEMBRANE membrane object
    //! @return score calculated for the given SSE
    storage::Pair< double, size_t> SSEPredictions::operator()
    (
      const assemble::SSE &SSE,
      const biol::Membrane &MEMBRANE
    ) const
    {
      double score( 0);

      const storage::Map< biol::SSType, math::PiecewiseFunction>::const_iterator
        itr( m_Potentials.Find( SSE.GetType()));

      if( itr == m_Potentials.End())
      {
        return storage::Pair< double, size_t>( score, 0);
      }

      // iterate over all amino acids in sequence
      for
      (
        biol::AASequence::const_iterator aa_itr( SSE.Begin()), aa_itr_end( SSE.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        score += ScoreAminoAcid( **aa_itr, SSE.GetType(), itr->second);
      }

      // end
      return storage::Pair< double, size_t>( score, SSE.GetSize());
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
    //! @param SSE SSE of interest
    //! @param MEMBRANE membrane object
    //! @param OSTREAM the std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &SSEPredictions::WriteDetailedSchemeAndValues
    (
      const assemble::SSE &SSE,
      const biol::Membrane &MEMBRANE,
      std::ostream &OSTREAM
    ) const
    {
      const storage::Map< biol::SSType, math::PiecewiseFunction>::const_iterator
        itr( m_Potentials.Find( SSE.GetType()));

      if( itr == m_Potentials.End())
      {
        return OSTREAM;
      }

      // iterate over all amino acids in sequence
      for
      (
        biol::AASequence::const_iterator aa_itr( SSE.Begin()), aa_itr_end( SSE.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        OSTREAM << ( *aa_itr)->GetIdentification() << '\t' << SSE.GetType()->GetOneLetterCode() << '\t' << ScoreAminoAcid( **aa_itr, SSE.GetType(), itr->second) << '\n';
      }

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEPredictions::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Method             , ISTREAM);
      io::Serialize::Read( m_ConfidenceThreshold, ISTREAM);
      io::Serialize::Read( m_Scheme             , ISTREAM);

      if( m_Method.IsDefined())
      {
        // read the histogram file and store the energy functions
        ReadEnergyVector();
      }

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SSEPredictions::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Method             , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ConfidenceThreshold, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme             , OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief score a single amino acid
    //! @param AMIO_ACID amino acid to score
    //! @param SS_TYPE the sstype amino acid prediction will be evaluated for
    //! @param PiecewiseFunction the PiecewiseFunction to use
    //! @return the score
    double SSEPredictions::ScoreAminoAcid( const biol::AABase &AMINO_ACID, const biol::SSType &SS_TYPE, const math::PiecewiseFunction &FUNC) const
    {
      // if Z score has sigma of zero, return score of zero
      auto pred( AMINO_ACID.GetData()->GetSSPrediction( m_Method));
      if( !pred.IsDefined())
      {
        return 0.0;
      }

      // get the raw predictions
      linal::Vector3D predictions( pred->GetThreeStatePrediction());
      // translate them into probabilities using the localPPV
      for( int i( 0); i < 3; ++i)
      {
        predictions( i) = m_Potentials.GetValue( biol::SSType( i))( predictions( i));
      }
      predictions.SetToSum( 1.0);

      // evaluate
      return -predictions( SS_TYPE);
    }

    //! @brief read energy distribution for scoring pairs of sses
    void SSEPredictions::ReadEnergyVector()
    {
      storage::Map< sspred::Method, storage::Map< biol::SSType, math::PiecewiseFunction> > funcs;

      io::IFStream read_str;
      io::File::MustOpenIFStream( read_str, Score::AddHistogramPath( GetDefaultHistogramFilename()));
      read_str >> funcs;
      io::File::CloseClearFStream( read_str);

      // search for method
      auto meth_itr( funcs.Find( m_Method));

      if( meth_itr == funcs.End())
      {
        BCL_MessageCrt
        (
          "there is no statistics for method: " + m_Method.GetName() + " in file: " + GetDefaultHistogramFilename() +
          " ==> no SSPred score for this method!!"
        );
        return;
      }
      m_Potentials = meth_itr->second;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool SSEPredictions::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      ReadEnergyVector();
      return true;
    }
  } // namespace score
} // namespace bcl
