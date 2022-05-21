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
#include "score/bcl_score_aa_neighborhood_exposure.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"
#include "biol/bcl_biol_aa_base.h"
#include "biol/bcl_biol_atom.h"
#include "biol/bcl_biol_membrane.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_histogram.h"
#include "score/bcl_score_energy_distribution.h"
#include "score/bcl_score_read_histograms.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AANeighborhoodExposure::s_Instance
    (
      util::Enumerated< AANeighborhoodInterface>::AddInstance( new AANeighborhoodExposure())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AANeighborhoodExposure::AANeighborhoodExposure() :
      m_AAExposure()
    {
    }

    //! @brief constructor from aa exposure function
    //! @param SP_AA_EXPOSURE ShPtr to AA exposure function to be used for aa exposure measure and scoring
    AANeighborhoodExposure::AANeighborhoodExposure
    (
      const assemble::AAExposureInterface &SP_AA_EXPOSURE
    ) :
      m_AAExposure( SP_AA_EXPOSURE)
    {
      ReadEnergyVector();
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new AANeighborhoodExposure copied from this one
    AANeighborhoodExposure *AANeighborhoodExposure::Clone() const
    {
      return new AANeighborhoodExposure( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AANeighborhoodExposure::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &AANeighborhoodExposure::GetScheme() const
    {
      // write the scheme of the aa exposure function being used
      return m_AAExposure->GetScheme();
    }

    //! @brief access to the minimal sequence separation
    //! @return minimal sequence separation for neighbors of that exposure score between amino acids in the same chain
    size_t AANeighborhoodExposure::GetMinimalSequenceSeparation() const
    {
      return m_AAExposure->GetMinimalSequenceSeparation();
    }

    //! @brief access to the distance cutoff
    //! @return distance cutoff above which the neighbor does not have influence on the score anymore
    double AANeighborhoodExposure::GetDistanceCutoff() const
    {
      return m_AAExposure->GetDistanceCutoff();
    }

    //! @brief access to the membrane potentials for given environment
    //! @param ENVIRONMENT environment for the potentials
    //! @return reference to membrane environment potentials
    const storage::Map< biol::AAType, math::CubicSplineDamped> &AANeighborhoodExposure::GetMembranePotentials
    (
      const biol::EnvironmentType &ENVIRONMENT
    ) const
    {
      return m_EnergyFunctionsMembrane.GetValue( ENVIRONMENT);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculates the sum of exposures of all amino acids for the given ProteinModel
    //! @param AA_NEIGHBOR_LIST neighbor list which's center amino acid is scored in the context
    //! @param MEMBRANE if it is defined, the score can be determined based on the membrane environment
    //! @return the exposure score for this neighbor list
    double AANeighborhoodExposure::operator()
    (
      const assemble::AANeighborList &AA_NEIGHBOR_LIST,
      const util::SiPtr< const biol::Membrane> &MEMBRANE
    ) const
    {
      const biol::AABase &current_aa( *AA_NEIGHBOR_LIST.GetCenterAminoAcid());

      // skip undefined aa types
      if( !current_aa.GetType().IsDefined() || !current_aa.GetType()->IsNaturalAminoAcid())
      {
        return 0.0;
      }

      // check that the current aminoacid has a defined coordinate
      if( !current_aa.GetFirstSidechainAtom().GetCoordinates().IsDefined())
      {
        return 0.0;
      }

      // calculate exposure for the current amino acid
      const double current_exposure( m_AAExposure->operator()( AA_NEIGHBOR_LIST));

      // score the exposure
      if( MEMBRANE.IsDefined())
      {
        const storage::Pair< biol::EnvironmentType, double> environment_weight
        (
          MEMBRANE->DetermineEnvironmentTypeAndWeight( current_aa.GetFirstSidechainAtom().GetCoordinates())
        );

        // for non-gap regions use the energy functions
        if( !environment_weight.First()->IsGap())
        {
          return m_EnergyFunctionsMembrane.Find( environment_weight.First())->second.Find( current_aa.GetType())->second( current_exposure);
        }
        // for gaps use a combination of bordering regions depending on the given weight

        // if it is a gap type, energy for adjacent env type have to be determined
        const biol::EnvironmentType env_type_left( environment_weight.First().GetIndex() - 1);
        const biol::EnvironmentType env_type_right( environment_weight.First().GetIndex() + 1);

        const double energy_left( m_EnergyFunctionsMembrane.Find( env_type_left)->second.Find( current_aa.GetType())->second( current_exposure));
        const double energy_right( m_EnergyFunctionsMembrane.Find( env_type_right)->second.Find( current_aa.GetType())->second( current_exposure));

        // weight the two energies depending on how close the z-coordinate is to left or right
        return   energy_left * environment_weight.Second()
               + energy_right * ( 1.0 - environment_weight.Second());
      }

      // soluble protein
      return m_EnergyFunctionsSoluble.Find( current_aa.GetType())->second( current_exposure);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @param INDENT indentation
    //! @return returns the output stream
    std::ostream &AANeighborhoodExposure::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_AAExposure, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &AANeighborhoodExposure::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_AAExposure, ISTREAM);

      ReadEnergyVector();

      // end
      return ISTREAM;
    }

    //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
    //! @param AA_NEIGHBOR_LIST neighbor list which's center amino acid is scored in the context
    //! @param OSTREAM the std::ostream to be written to
    //! @param MEMBRANE if it is defined, the score can be determined based on the membrane environment
    //! @return std::ostream which was written to
    std::ostream &
    AANeighborhoodExposure::WriteDetailedSchemeAndValues
    (
      const assemble::AANeighborList &AA_NEIGHBOR_LIST,
      const util::SiPtr< const biol::Membrane> &MEMBRANE,
      std::ostream &OSTREAM
    ) const
    {
      const biol::AABase &current_aa( *AA_NEIGHBOR_LIST.GetCenterAminoAcid());

      // skip undefined aa types
      if( !current_aa.GetType().IsDefined() || !current_aa.GetType()->IsNaturalAminoAcid())
      {
        return OSTREAM;
      }

      // check that the current aminoacid has a defined coordinate
      if( !current_aa.GetFirstSidechainAtom().GetCoordinates().IsDefined())
      {
        return OSTREAM;
      }

      // calculate exposure for the current amino acid
      const double current_exposure( m_AAExposure->operator()( AA_NEIGHBOR_LIST));

      // write seq id of current AA
      OSTREAM << current_aa.GetIdentification() << '\t';
      // write exposure
      OSTREAM << current_exposure << ": ";

      // score the exposure
      if( MEMBRANE.IsDefined())
      {
        const storage::Pair< biol::EnvironmentType, double> environment_weight
        (
          MEMBRANE->DetermineEnvironmentTypeAndWeight( current_aa.GetFirstSidechainAtom().GetCoordinates())
        );

        double current_score( 0.0);
        // for non-gap regions use the energy functions
        if( !environment_weight.First()->IsGap())
        {
          current_score = m_EnergyFunctionsMembrane.Find( environment_weight.First())->second.Find( current_aa.GetType())->second( current_exposure);
        }
        // for gaps use a combination of bordering regions depending on the given weight
        else
        {
          // if it is a gap type, energy for adjacent env type have to be determined
          const biol::EnvironmentType env_type_left( environment_weight.First().GetIndex() - 1);
          const biol::EnvironmentType env_type_right( environment_weight.First().GetIndex() + 1);

          const double energy_left( m_EnergyFunctionsMembrane.Find( env_type_left)->second.Find( current_aa.GetType())->second( current_exposure));
          const double energy_right( m_EnergyFunctionsMembrane.Find( env_type_right)->second.Find( current_aa.GetType())->second( current_exposure));

          // weight the two energies depending on how close the z-coordinate is to left or right
          current_score =   energy_left * environment_weight.Second()
                          + energy_right * ( 1.0 - environment_weight.Second());
        }

        // write according score
        OSTREAM << current_score << '\t'
                // write z-coordinate
                << current_aa.GetFirstSidechainAtom().GetCoordinates().Z() << '\t'
                // write environment
                << environment_weight.First().GetName() << '\t'
                // write weight
                << environment_weight.Second() << '\n';

        // end
        return OSTREAM;
      }

      // soluble protein
      return OSTREAM << m_EnergyFunctionsSoluble.Find( current_aa.GetType())->second( current_exposure) << '\n';
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AANeighborhoodExposure::GetSerializer() const
    {
      io::Serializer serializer;

      serializer.SetClassDescription( "Scores the exposure of an amino acid in the context of its neighbors.");
      serializer.AddInitializer
      (
        "aa exposure",
        "AAExposure function to be used for calculations",
        io::Serialization::GetAgent( &m_AAExposure)
      );

      return serializer;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief read the energy vector for the AminoAcid neighbor counts
    void
    AANeighborhoodExposure::ReadEnergyVector()
    {
      // initialize read
      io::IFStream read;
      const std::string filename( Score::AddHistogramPath( m_AAExposure->GetHistogramFileName()));
      io::File::MustOpenIFStream( read, filename);

      // read boundaries
      math::Range< double> threshold_low_high;
      read >> threshold_low_high;

      // read sequence separation
      size_t minimal_sequence_separation( 0);
      read >> minimal_sequence_separation;

      m_AAExposure->SetThresholdRange( threshold_low_high);
      m_AAExposure->SetMinimalSequenceSeparation( minimal_sequence_separation);

      // read membrane and environment dependent potential
      const storage::Vector< storage::Vector< math::Histogram> >
        membrane_aa_histograms( ReadHistograms::ReadMembraneDependentEnvironmentHistograms( read));

      // close and clear stream
      io::File::CloseClearFStream( read);

      // derive energy distribution for membrane aa environment potential and store it
      m_EnergyFunctionsMembrane = EnergyDistribution::AAMembraneEnvironmentPotential( membrane_aa_histograms);
      m_EnergyFunctionsSoluble = EnergyDistribution::AAEnvironmentPotential( membrane_aa_histograms( biol::GetEnvironmentTypes().e_Solution));
    }

  } // namespace score
} // namespace bcl
