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
#include "score/bcl_score_phi_psi.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
#include "biol/bcl_biol_atom.h"
#include "biol/bcl_biol_membrane.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_bicubic_spline.h"
#include "math/bcl_math_gnuplot_heatmap.h"
#include "math/bcl_math_histogram_2d.h"
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
    const util::SiPtr< const util::ObjectInterface> PhiPsi::s_Instance
    (
      util::Enumerated
      <
        math::BinaryFunctionInterfaceSerializable< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> >
      >::AddInstance( new PhiPsi())
    );

    //! @brief returns default file where the statistics and in consequence the energy potentials are read from
    //! @return the default filename for the default phi_psi_angle statistics
    const std::string &PhiPsi::GetDefaultHistogramFilename()
    {
      // static string
      static const std::string s_default_histogram_filename( "phi_psi_angles_by_sstype.histogram2D");

      // end
      return s_default_histogram_filename;
    }

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &PhiPsi::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "phi_psi");

      // end
      return s_default_scheme;
    }

    //! @brief get the name of the object
    //! @return the name of the object
    const std::string &PhiPsi::GetAlias() const
    {
      static const std::string s_name( "PhiPsi");
      return s_name;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer PhiPsi::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "scores the agreement of phi psi values with expected probabilities for an SSE");
      serializer.AddInitializer
      (
        "histogram filename",
        "path to where the statistics and energy potentials are read from",
        io::Serialization::GetAgent( &m_HistogramFileName),
        GetDefaultHistogramFilename()
       );
      serializer.AddOptionalInitializer
      (
        "ss types",
        "types of SSEs whose phi psi angles will be scored independently of the ss type",
        io::Serialization::GetAgent( &m_SSTypes)
       );

      return serializer;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PhiPsi::PhiPsi() :
      m_Scheme( GetDefaultScheme()),
      m_HistogramFileName(),
      m_EnergyMapSoluble(),
      m_EnergyMapMembrane(),
      m_AATypeEnergyMap(),
      m_SSTypes()
    {
    }

    //! @brief constructor from a specified histogram file
    //! @param SCHEME scheme to be used
    //! @param HISTOGRAM_FILENAME filename of the histogram to be used
    //! @param SSTYPE_SET the types of sses whose phi psi angles will be score independently of the sstype
    PhiPsi::PhiPsi
    (
      const std::string &SCHEME,
      const std::string &HISTOGRAM_FILENAME,
      const storage::Set< biol::SSType> &SSTYPE_SET
    ) :
      m_Scheme( SCHEME),
      m_HistogramFileName( HISTOGRAM_FILENAME),
      m_EnergyMapSoluble(),
      m_EnergyMapMembrane(),
      m_AATypeEnergyMap(),
      m_SSTypes( SSTYPE_SET)
    {
      // read the histogram file and store the energy functions
      ReadEnergyVector();
    }

    //! @brief virtual copy constructor
    PhiPsi *PhiPsi::Clone() const
    {
      return new PhiPsi( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PhiPsi::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief Operator that scores the phi_psi bending of various SSEs in a protein
    //! @param THIS_SSE SSE to be scored
    //! @param MEMBRANE membrane object
    //! @return pair of score and number of scored amino acids
    storage::Pair< double, size_t> PhiPsi::operator()
    (
      const assemble::SSE &THIS_SSE,
      const biol::Membrane &MEMBRANE
    ) const
    {
      // ensure that the amino acid chain length is at least three
      if( THIS_SSE.GetSize() < 3)
      {
        return storage::Pair< double, size_t>( double( 0.0), 0);
      }

      // true if the phi and psi angles of the sse should be scored without taking into account the type of sse
      if( m_SSTypes.Find( THIS_SSE.GetType()) != m_SSTypes.End())
      {
        return ScorePhiPsiSSTypeIndependent( THIS_SSE);
      }

      // find the vector of histograms with corresponding sstype
      const storage::Map
      <
        biol::SSType, storage::Map< biol::AAType, math::BicubicSpline>
      >::const_iterator map_itr_sol( m_EnergyMapSoluble.Find( THIS_SSE.GetType())),
                        map_itr_mem( m_EnergyMapMembrane.Find( THIS_SSE.GetType()));

      // make sure it's found
      if( map_itr_sol == m_EnergyMapSoluble.End())
      {
        BCL_MessageCrt( "unexpected sstype found " + THIS_SSE.GetType().GetName());
        return storage::Pair< double, size_t>( double( 0.0), 0);
      }

      // end
      return ScoreSequencePhiPsi( THIS_SSE, map_itr_sol->second, map_itr_mem->second, MEMBRANE);
    }

    //! @brief scores an sse according to phi and psi without taking into account the SS Type
    //! @param THIS_SSE will be scored for agreement of phi and psi probabilities
    //! @return pair of double and size_t which is the score and number of scored amino acids, respectively
    storage::Pair< double, size_t> PhiPsi::ScorePhiPsiSSTypeIndependent( const assemble::SSE &THIS_SSE) const
    {
      return ScoreSequencePhiPsi( THIS_SSE, m_AATypeEnergyMap, m_AATypeEnergyMap, biol::Membrane::GetUndefinedMembrane());
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &PhiPsi::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Scheme           , ISTREAM);
      io::Serialize::Read( m_HistogramFileName, ISTREAM);
      io::Serialize::Read( m_SSTypes          , ISTREAM);

      // read the histogram files and store the energy functions
      ReadEnergyVector();

      // end
      return ISTREAM;
    }

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @param INDENT indentation
    //! @return returns the output stream
    std::ostream &PhiPsi::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Scheme           , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_HistogramFileName, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SSTypes          , OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read from this class
    //! @param ERROR_STREAM stream with which to write errors
    bool PhiPsi::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERROR_STREAM
      )
    {
      ReadEnergyVector();
      return true;
    }

    //! @brief read the energy distribution for scoring phi_psi angles
    void PhiPsi::ReadEnergyVector()
    {
      // open the histogram file into a stream
      m_EnergyMapSoluble.Reset();
      m_EnergyMapMembrane.Reset();
      m_AATypeEnergyMap.Reset();

      // initialize read
      io::IFStream read;
      io::File::MustOpenIFStream( read, Score::AddHistogramPath( m_HistogramFileName));

      util::SiPtrVector< storage::Map< biol::SSType, storage::Map< biol::AAType, math::BicubicSpline> > > soluble_membrane_maps;
      soluble_membrane_maps.PushBack( util::ToSiPtrNonConst( m_EnergyMapSoluble));
      soluble_membrane_maps.PushBack( util::ToSiPtrNonConst( m_EnergyMapMembrane));
      storage::Vector< std::string> suffices( storage::Vector< std::string>::Create( "soluble", "membrane"));
      auto itr_suffices( suffices.Begin());
      for
      (
        auto itr_env( soluble_membrane_maps.Begin()), itr_env_end( soluble_membrane_maps.End());
        itr_env != itr_env_end;
        ++itr_env, ++itr_suffices
      )
      {
        std::string env_type;
        read >> env_type;
        // store aatype_sstype_histograms before deriving energy functions
        storage::Map< biol::AAType, storage::Map< biol::SSType, math::Histogram2D> > aatype_sstype_histograms;

        math::Histogram2D background;

        // for each sstype
        for
        (
          biol::SSTypes::const_iterator
            sstype_itr( biol::GetSSTypes().Begin()), sstype_itr_end( biol::GetSSTypes().COIL.GetIterator() + 1);
          sstype_itr != sstype_itr_end; ++sstype_itr
        )
        {
          // initialize SSType to be read
          biol::SSType current_sstype;

          // read sstype
          read >> current_sstype;

          // make sure it's the same as the expected
          BCL_Assert
          (
            current_sstype == *sstype_itr,
            "Unexpected sstype read from file!" + sstype_itr->GetName() + " != " + current_sstype.GetName()
          );

          // for each amino acids type read histogram
          for
          (
            biol::AATypes::const_iterator
              aa_itr( biol::GetAATypes().Begin()),
              aa_itr_end( biol::GetAATypes().GetEnumIteratorFromIndex( biol::AATypes::s_NumberStandardAATypes));
            aa_itr != aa_itr_end;
            ++aa_itr
          )
          {
            // initialize aatype to be read
            biol::AAType current_aatype;

            // read one letter code that represents the aatype
            std::string tmp;
            read >> tmp;
            current_aatype = biol::GetAATypes().AATypeFromOneLetterCode( tmp[ 0]);

            // assert that the aatypes in the histogram file are in the same order
            BCL_Assert
            (
              *aa_itr == current_aatype,
              "unexpected aatype read from file! " + aa_itr->GetName() + " != " + current_aatype.GetName()
            );

            // read Histogram2D from stream
            math::Histogram2D &current_phi_psi_histogram( aatype_sstype_histograms[ *aa_itr][ current_sstype]);
            read >> current_phi_psi_histogram;

            // normalize
            current_phi_psi_histogram.Normalize();

            // add to background
            background.Combine( current_phi_psi_histogram);
          }
        }

        // normalize background
        background.Normalize();

        if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
        {
          io::OFStream write;
          io::File::MustOpenOFStream( write, "phi_psi_background_" + *itr_suffices + ".gnuplot");
          math::GnuplotHeatmap heatmap;
          heatmap.SetFromHistogram( background, false, false);
          heatmap.SetTitleAndLabel( "phi psi background " + *itr_suffices, "phi [" + math::Angle::s_DegreeSymbolGnuplot + "]", "psi [" + math::Angle::s_DegreeSymbolGnuplot + "]", "p");
          const linal::Vector< double> binning
          (
            linal::FillVector< double>
            (
              background.GetNumberOfBinsX() + 1,
              math::Angle::Degree( background.GetBoundariesX().First()),
              math::Angle::Degree( background.GetBinSizeXY().First())
            )
          );
          heatmap.SetTicsX( math::GnuplotHeatmap::TicsFromBinning( binning, 1, util::Format().FFP( 0).W( 3).Fill( ' ').R()), false, 1);
          heatmap.SetTicsY( math::GnuplotHeatmap::TicsFromBinning( binning, 1, util::Format().FFP( 0).W( 3).Fill( ' ').R()), false, 1);
          heatmap.SetFont( "arialbd", 16);
          heatmap.SetPixelAndRatio( 1080, 800, -1);
          heatmap.SetRotationXTics( 90.0);
          heatmap.SetFilename( "phi_psi_background");
          heatmap.WriteScript( write);
          io::File::CloseClearFStream( write);
        }

        // generate scores

        for
        (
          storage::Map< biol::AAType, storage::Map< biol::SSType, math::Histogram2D> >::const_iterator
            aa_itr( aatype_sstype_histograms.Begin()), aa_itr_end( aatype_sstype_histograms.End());
          aa_itr != aa_itr_end;
          ++aa_itr
        )
        {
          // sum of all sstypes, for sstype independent score
          math::Histogram2D sstype_sum;
          for
          (
            storage::Map< biol::SSType, math::Histogram2D>::const_iterator
              ss_itr( aa_itr->second.Begin()), ss_itr_end( aa_itr->second.End());
            ss_itr != ss_itr_end;
            ++ss_itr
          )
          {
            // add to sstype independent sum
            sstype_sum.Combine( ss_itr->second);

            ( **itr_env)[ ss_itr->first][ aa_itr->first] = EnergyDistribution::PhiPsiAnglePotential( ss_itr->second, background, false);
            //          m_EnergyMap[ ss_itr->first][ aa_itr->first] = EnergyDistribution::PhiPsiAnglePotential( ss_itr->second, background, false);
          }

          if( itr_suffices == suffices.Begin())
          {
            sstype_sum.Normalize();
            m_AATypeEnergyMap[ aa_itr->first] = EnergyDistribution::PhiPsiAnglePotential( sstype_sum, background, false);
            //        m_AATypeEnergyMap[ aa_itr->first] = EnergyDistribution::PhiPsiAnglePotential( sstype_sum, background, false);
          }
        }
      }

      // close and clear read stream
      io::File::CloseClearFStream( read);
    }

    //! @brief scores the phi psi of a single residue according to a given energy map
    //! @param AATYPE_ENERGY the energy used to score the phi psi of the residue
    //! @param AA_BASE the residue of interest
    //! @param AA_PREVIOUS the residue previous in sequence to the residue of interest
    //! @param AA_NEXT the residue following in sequence to the residue of interest
    //! @return double which is the score of the phi and psi angles for the residue of interest
    double PhiPsi::ScoreAAPhiPsi
    (
      const math::BicubicSpline &AATYPE_ENERGY,
      const biol::AABase &AA_BASE,
      const biol::AABase &AA_PREVIOUS,
      const biol::AABase &AA_NEXT
    )
    {
      // calculate phi psi for selected amino acid
      const storage::Pair< double, double> phi_psi_pair
      (
        AA_BASE.CalculatePhiPsi
        (
          AA_PREVIOUS.GetAtom( biol::GetAtomTypes().C),
          AA_NEXT.GetAtom( biol::GetAtomTypes().N)
        )
      );

      // check that angles are defined
      if( !util::IsDefined( phi_psi_pair.First()) || !util::IsDefined( phi_psi_pair.Second()))
      {
        return util::GetUndefinedDouble();
      }

      // get the score from the bicubicspline for this amino acid add this score to the sum
      const double score( AATYPE_ENERGY.F( linal::MakeVector( phi_psi_pair.Second(), phi_psi_pair.First())));

      // end
      return score;
    }

    //! @brief scores the phi and psi angles of an entire aa sequence according to a given energy map
    //! @param AA_SEQUENCE the amino acid sequence which will be scored
    //! @param AATYPE_ENERGY_MAP the energy map used to score the amino acid sequence
    //! @return pair of double and size_t which are the score and number of residues scored, respectively
    storage::Pair< double, size_t> PhiPsi::ScoreSequencePhiPsi
    (
      const biol::AASequence &AA_SEQUENCE,
      const storage::Map< biol::AAType, math::BicubicSpline> &AATYPE_ENERGY_MAP_SOLUBLE,
      const storage::Map< biol::AAType, math::BicubicSpline> &AATYPE_ENERGY_MAP_MEMBRANE,
      const biol::Membrane &MEMBRANE
    )
    {
      // initialize the sum
      double score_sum( 0.0);
      size_t scored_entities( 0);

      // iterate over amino acids
      for
      (
        biol::AASequence::const_iterator
          aa_itr( AA_SEQUENCE.Begin() + 1), aa_itr_end( AA_SEQUENCE.End() - 1);
        aa_itr != aa_itr_end; ++aa_itr
      )
      {
        const storage::Map< biol::AAType, math::BicubicSpline> &map
        (
          MEMBRANE.IsDefined() &&
          MEMBRANE.DetermineEnvironmentType( ( *aa_itr)->GetCA().GetCoordinates())->GetReducedType()
          == biol::GetEnvironmentTypes().e_MembraneCore
          ? AATYPE_ENERGY_MAP_MEMBRANE
          : AATYPE_ENERGY_MAP_SOLUBLE
        );

        // checks whether amino acid type is contained within the energy map.
        const storage::Map< biol::AAType, math::BicubicSpline>::const_iterator spline_itr
        (
          map.Find( ( *aa_itr)->GetType())
        );

        // make sure a bicubic spline exists for the selected amino acid type
        if( spline_itr == map.End())
        {
          continue;
        }

        // get the score from the bicubicspline for this amino acid
        const double current_score( ScoreAAPhiPsi( spline_itr->second, **aa_itr, **( aa_itr - 1), **( aa_itr + 1)));

        // true if "current_score" is defined
        if( util::IsDefined( current_score))
        {
          // add this score to the sum and increment "scored_entities"
          score_sum += current_score;
          ++scored_entities;
        }
      }

      // end
      return storage::Pair< double, size_t>( score_sum, scored_entities);
    }

  } // namespace score
} // namespace bcl

