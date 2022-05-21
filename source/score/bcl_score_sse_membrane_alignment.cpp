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
#include "score/bcl_score_sse_membrane_alignment.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
#include "biol/bcl_biol_membrane.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_cubic_spline.h"
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
    const std::string &SSEMembraneAlignment::GetDefaultHistogramFilename()
    {
      // static string
      static const std::string s_default_histogram_filename( "sse_membrane_alignment.histograms");

      // end
      return s_default_histogram_filename;
    }

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &SSEMembraneAlignment::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "ssealign");

      // end
      return s_default_scheme;

    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a specified histogram file
    //! @param HISTOGRAM_FILENAME filename of the histogram to be used
    //! @param SCHEME scheme to be used
    SSEMembraneAlignment::SSEMembraneAlignment
    (
      const std::string &HISTOGRAM_FILENAME,
      const std::string &SCHEME
    ) :
      m_HistogramFileName( HISTOGRAM_FILENAME),
      m_Scheme( SCHEME),
      m_EnergyFunctions()
    {
      // read the histogram file and store the energy functions
      ReadEnergyFunctions();
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new SSEMembraneAlignment object copied from this one
    SSEMembraneAlignment *SSEMembraneAlignment::Clone() const
    {
      return new SSEMembraneAlignment( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSEMembraneAlignment::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &SSEMembraneAlignment::GetScheme() const
    {
      return m_Scheme;
    }

    //! @brief access energy functions
    //! @return map for each sstype, that has a map for each environmenttype with a cubic spline
    const storage::Map< biol::SSType, storage::Map< biol::EnvironmentType, math::CubicSplineDamped> > &
    SSEMembraneAlignment::GetEnergyFunctions() const
    {
      return m_EnergyFunctions;
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &SSEMembraneAlignment::GetAlias() const
    {
      static const std::string s_name( "SSEMembraneAlignment");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SSEMembraneAlignment::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Scores the alignment of transmembrane SSEs with the membrane normal.");
      serializer.AddInitializer
      (
        "histogram filename",
        "path to the histogram file that is used to construct the potential function",
        io::Serialization::GetAgent( &m_HistogramFileName)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculates the projection angle of the given SSE to the membrane plane
    //! @param SSE_GEOMETRY sse geometry derived class of interest
    //! @param AXIS membrane axis
    //! @param MEMBRANE membrane object
    //! @return calculated angle
    double SSEMembraneAlignment::AngleToMembranePlane
    (
      const assemble::SSEGeometryInterface &SSE_GEOMETRY,
      const coord::Axis &AXIS,
      const biol::Membrane &MEMBRANE
    )
    {
      return math::Absolute( linal::ProjAngle( MEMBRANE.GetNormal(), SSE_GEOMETRY.GetAxis( AXIS)) - math::g_Pi / 2.0);
    }

    //! @brief calculate weight for the given strand, that it is not turned in membrane
    //! @param STRAND sse geometry derived class of interest
    //! @param MEMBRANE membrane object
    //! @return calculated weight
    double SSEMembraneAlignment::WeightXAxis
    (
      const assemble::SSEGeometryInterface &STRAND,
      const biol::Membrane &MEMBRANE
    )
    {
      if( STRAND.GetType() == biol::GetSSTypes().HELIX)
      {
        return double( 1);
      }
      else if( STRAND.GetType() == biol::GetSSTypes().STRAND)
      {
        return math::WeightBetweenZeroAndPi_ThreeSections( AngleToMembranePlane( STRAND, coord::GetAxes().e_X, MEMBRANE) * 2.0);
      }
      else
      {
        return double( 0);
      }

      return util::GetUndefined< double>();
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that calculates the score for a given SSE
    //! @param SSE SSE of interest
    //! @param MEMBRANE membrane object
    //! @return score calculated for the given SSE
    storage::Pair< double, size_t> SSEMembraneAlignment::operator()
    (
      const assemble::SSE &SSE,
      const biol::Membrane &MEMBRANE
    ) const
    {
      // make sure it's helix or strand
      if( !SSE.GetType()->IsStructured())
      {
        return storage::Pair< double, size_t>( double( 0), 0);
      }

      // initialize score and count
      storage::Pair< double, size_t> score_count( 0.0, 0);

      // iterate over all fragments
      for
      (
        util::ShPtrVector< assemble::SSEGeometryInterface>::const_iterator frag_itr( SSE.GetFragments().Begin()),
          frag_itr_end( SSE.GetFragments().End());
        frag_itr != frag_itr_end;
        ++frag_itr
      )
      {
        // determine environment type and weight in gap
        const storage::Pair< biol::EnvironmentType, double> environment_weight
        (
          MEMBRANE.DetermineEnvironmentTypeAndWeight( ( *frag_itr)->GetCenter())
        );

        BCL_MessageDbg
        (
          "ssemembranealignment: " +
          environment_weight.First().GetName() + " " + ( *frag_itr)->GetIdentification() +
          " origin " + util::Format()( ( *frag_itr)->GetCenter())
        );

        const double angle_to_membrane_plane( AngleToMembranePlane( **frag_itr, coord::GetAxes().e_Z, MEMBRANE));
        double energy( 0);

        if( !environment_weight.First()->IsGap())
        {
          const storage::Map< biol::SSType, storage::Map< biol::EnvironmentType, math::CubicSplineDamped> >::const_iterator
            sse_type_itr
            (
              m_EnergyFunctions.Find
              (
                ( *frag_itr)->GetType()
              )
            );
          if( sse_type_itr == m_EnergyFunctions.End())
          {
            BCL_MessageDbg( "sse type not found in map of energy functions");
            continue;
          }
          const storage::Map< biol::EnvironmentType, math::CubicSplineDamped>::const_iterator environment_type_itr
          (
            sse_type_itr->second.Find( environment_weight.First())
          );

          if( environment_type_itr == sse_type_itr->second.End())
          {
            BCL_MessageDbg
            (
              "environment type " + environment_weight.First()->GetName() + " not found for ss type: " +
              sse_type_itr->first->GetName()
            );
            continue;
          }

          // calculate the energy
          energy = environment_type_itr->second( angle_to_membrane_plane);
        }
        else
        {
          // if it is a gap type, energy for adjacent env type have to be determined
          const biol::EnvironmentType env_type_left( environment_weight.First().GetIndex() - 1);
          const biol::EnvironmentType env_type_right( environment_weight.First().GetIndex() + 1);

          const double energy_left
          (
            m_EnergyFunctions.Find
            (
              ( *frag_itr)->GetType()
            )->second.Find( env_type_left)->second( angle_to_membrane_plane)
          );
          const double energy_right
          (
            m_EnergyFunctions.Find
            (
              ( *frag_itr)->GetType()
            )->second.Find( env_type_right)->second( angle_to_membrane_plane)
          );

          // weight the two energies depending on how close the z-coordinate is to left or right
          energy += energy_left * environment_weight.Second();
          energy += energy_right * ( 1.0 - environment_weight.Second());
        }

        // weight of the x axes orientation (necessary for strands)
        if( ( *frag_itr)->GetType() == biol::GetSSTypes().STRAND)
        {
          energy *= WeightXAxis( **frag_itr, MEMBRANE);
        }

        score_count.First() += energy;
        ++score_count.Second();
      }

      // energy and number of scored entities
      return score_count;
    }

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    //! @return return code indicating success or failure
    bool SSEMembraneAlignment::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      ReadEnergyFunctions();

      return true;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &SSEMembraneAlignment::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_HistogramFileName, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // read the histogram file and store the energy functions
      ReadEnergyFunctions();

      // end
      return ISTREAM;
    }

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @return returns the output stream
    std::ostream &SSEMembraneAlignment::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_HistogramFileName, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
    //! @param SSE SSE of interest
    //! @param MEMBRANE membrane object
    //! @param OSTREAM the std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &
    SSEMembraneAlignment::WriteDetailedSchemeAndValues
    (
      const assemble::SSE &SSE,
      const biol::Membrane &MEMBRANE,
      std::ostream &OSTREAM
    ) const
    {
      //write sstype, Angel to membrane plane first and last amino acid and value to the STREAM
      OSTREAM << SSE.GetType() << '\t'
              << AngleToMembranePlane( SSE, coord::GetAxes().e_Z, MEMBRANE) << '\t'
              << SSE.GetFirstAA()->GetSeqID() << '\t'
              << SSE.GetFirstAA()->GetType()->GetThreeLetterCode() << '\t'
              << SSE.GetLastAA()->GetSeqID() << '\t'
              << SSE.GetLastAA()->GetType()->GetThreeLetterCode() << '\t'
              << operator()( SSE, MEMBRANE).First() << '\n';

      //end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief read energy distribution for scoring sse membrane alignment
    void
    SSEMembraneAlignment::ReadEnergyFunctions()
    {
      // reset the energy functions
      m_EnergyFunctions.Reset();

      //initialize sstype, and vector to store the histograms read from the file
      storage::Map< biol::SSType, storage::Map< biol::EnvironmentType, math::Histogram> > sse_membrane_alignment;

      // initialize read
      io::IFStream read;
      io::File::MustOpenIFStream( read, Score::AddHistogramPath( m_HistogramFileName));

      //read alignment histograms
      read >> sse_membrane_alignment;

      // reset stream
      io::File::CloseClearFStream( read);

      // iterate over two sstypes
      for
      (
        biol::SSTypes::const_iterator sstype_itr( biol::GetSSTypes().Begin()),
          sstype_itr_end( biol::GetSSTypes().COIL.GetIterator());
        sstype_itr != sstype_itr_end;
        ++sstype_itr
      )
      {
        // iterate over every environment and copy the counts
        for
        (
          biol::EnvironmentTypes::const_iterator env_itr( biol::GetEnvironmentTypes().GetReducedTypes().Begin()),
            env_itr_end( biol::GetEnvironmentTypes().GetReducedTypes().End());
          env_itr != env_itr_end; ++env_itr
        )
        {
          // calculate spline
          m_EnergyFunctions[ *sstype_itr][ *env_itr] =
            EnergyDistribution::AngleAlignmentPotential( sse_membrane_alignment[ *sstype_itr][ *env_itr]);

        } // env_itr
      } // sstype_itr
    }

  } // namespace score
} // namespace bcl
