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
#include "score/bcl_score_environment_predictions.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_atom.h"
#include "biol/bcl_biol_membrane.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> EnvironmentPredictions::s_Instance
    (
      GetObjectInstances().AddInstance( new EnvironmentPredictions())
    );

    //! @brief returns default file where the statistics and in consequence the energy potentials are read from
    //! @return default file where the statistics and in consequence the energy potentials are read from
    const std::string &EnvironmentPredictions::GetDefaultHistogramFilename()
    {
      // static string
      static const std::string s_default_histogram_filename( "sspred_env_mean.bcl");

      // end
      return s_default_histogram_filename;
    }

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &EnvironmentPredictions::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "ssepred_env");

      // end
      return s_default_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    EnvironmentPredictions::EnvironmentPredictions() :
      m_Method( sspred::GetMethods().e_Undefined),
      m_Scheme( GetDefaultScheme()),
      m_Potentials()
    {
    }

    //! @brief constructor from a ssmethod and confidence threshold
    //! @param SS_METHOD ssmethod to use
    //! @param CONFIDENCE_THRESHOLD the threshold in units of z-score, above which the score becomes negative
    //! @param SCHEME scheme to be used
    EnvironmentPredictions::EnvironmentPredictions
    (
      const sspred::Method &SS_METHOD,
      const std::string &SCHEME
    ) :
      m_Method( SS_METHOD),
      m_Scheme( SCHEME),
      m_Potentials()
    {
      if( m_Method.IsDefined())
      {
        // read the histogram file and store the energy functions
        ReadEnergyVector();
      }
    }

    //! @brief Clone function
    //! @return pointer to new EnvironmentPredictions
    EnvironmentPredictions *EnvironmentPredictions::Clone() const
    {
      return new EnvironmentPredictions( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &EnvironmentPredictions::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that calculates the score for a given SSE
    //! @param THIS_SSE SSE of interest
    //! @param MEMBRANE membrane object
    //! @return score calculated for the given SSE
    storage::Pair< double, size_t> EnvironmentPredictions::operator()
    (
      const assemble::SSE &THIS_SSE,
      const biol::Membrane &MEMBRANE
    ) const
    {
      // initialize the score
      double score( 0.0);

      // iterate over all amino acids in sequence
      for
      (
        biol::AASequence::const_iterator aa_itr( THIS_SSE.Begin()),
          aa_itr_end( THIS_SSE.End());
        aa_itr != aa_itr_end; ++aa_itr
      )
      {
        // get the environment type
        const biol::EnvironmentType current_env_type
        (
          MEMBRANE.IsDefined()
          ? MEMBRANE.DetermineEnvironmentType( ( *aa_itr)->GetAtom( biol::GetAtomTypes().CA).GetCoordinates())->GetReducedType()
          : biol::GetEnvironmentTypes().e_Solution
        );

        // find the potential
        auto itr( m_Potentials.Find( current_env_type));

        // add the score
        if( itr != m_Potentials.End())
        {
          score += ScoreAminoAcid( **aa_itr, current_env_type, itr->second);
        }
      }

      // end
      return storage::Pair< double, size_t>( score, THIS_SSE.GetSize());
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
    //! @param THIS_SSE SSE of interest
    //! @param MEMBRANE membrane object
    //! @param OSTREAM the std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &EnvironmentPredictions::WriteDetailedSchemeAndValues
    (
      const assemble::SSE &THIS_SSE,
      const biol::Membrane &MEMBRANE,
      std::ostream &OSTREAM
    ) const
    {
      // iterate over all amino acids in sequence
      for
      (
        biol::AASequence::const_iterator aa_itr( THIS_SSE.Begin()),
          aa_itr_end( THIS_SSE.End());
        aa_itr != aa_itr_end; ++aa_itr
      )
      {
        // get the environment type
        const biol::EnvironmentType current_env_type
        (
          MEMBRANE.IsDefined() ?
            MEMBRANE.DetermineEnvironmentType
            (
              ( *aa_itr)->GetAtom( biol::GetAtomTypes().CA).GetCoordinates()
            )->GetReducedType() :
            biol::GetEnvironmentTypes().e_Solution
        );

        // find the potential
        auto itr( m_Potentials.Find( current_env_type));

        // add the score
        if( itr != m_Potentials.End())
        {
          OSTREAM << ( *aa_itr)->GetIdentification() << '\t' << current_env_type->GetTwoLetterCode() << '\t'
                  << ScoreAminoAcid( **aa_itr, current_env_type, itr->second) << '\n';
        }
      }

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &EnvironmentPredictions::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Method             , ISTREAM);
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
    std::ostream &EnvironmentPredictions::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Method             , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme             , OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief score a single amino acid
    //! @param AMIO_ACID amino acid to score
    //! @param ENV_TYPE the environment type the amino acid prediction will be evaluated for
    //! @param ZSCORE the zscore to use
    //! @return the score
    double EnvironmentPredictions::ScoreAminoAcid
    (
      const biol::AABase &AMINO_ACID,
      const biol::EnvironmentType &ENV_TYPE,
      const math::PiecewiseFunction &ZSCORE
    ) const
    {
      auto pred( AMINO_ACID.GetData()->GetSSPrediction( m_Method));
      if( !pred.IsDefined())
      {
        return 0.0;
      }
      // get the raw predictions
      linal::Vector3D predictions( pred->GetThreeStateTMPrediction());
      // translate them into probabilities using the localPPV
      predictions( 0) = m_Potentials.GetValue( biol::GetEnvironmentTypes().e_MembraneCore)( predictions( 0));
      predictions( 1) = m_Potentials.GetValue( biol::GetEnvironmentTypes().e_Transition)( predictions( 1));
      predictions( 2) = m_Potentials.GetValue( biol::GetEnvironmentTypes().e_Solution)( predictions( 2));
      predictions.SetToSum( 1.0);

      // evaluate
      return -predictions( ENV_TYPE->GetReducedIndex());
    }

    //! @brief read energy distribution for scoring pairs of sses
    void EnvironmentPredictions::ReadEnergyVector()
    {
      storage::Map< sspred::Method, storage::Map< biol::EnvironmentType, math::PiecewiseFunction> > funcs;

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
          " ==> no SSPred Environment score for this method!!"
        );
        return;
      }
      m_Potentials = meth_itr->second;
    }

  } // namespace score
} // namespace bcl
