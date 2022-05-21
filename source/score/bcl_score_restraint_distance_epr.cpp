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
#include "score/bcl_score_restraint_distance_epr.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container_generator_protein_model.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_histogram_2d.h"
#include "score/bcl_score_energy_distribution.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from a restraint set, histogram file, and scheme
    //! @param RESTRAINTS EPR restraints for scoring the protein model
    //! @param HISTOGRAM_FILENAME name of the histogram file used to interpolate the potential function
    //! @param SCHEME scheme of this score
    RestraintDistanceEPR::RestraintDistanceEPR
    (
      const util::ShPtrVector< restraint::AtomDistance> &RESTRAINTS,
      const std::string &HISTOGRAM_FILENAME,
      const std::string &SCHEME
    ) :
      m_Scheme( SCHEME),
      m_EnergyFunction( CreateEnergyFunction( HISTOGRAM_FILENAME)),
      m_Restraints( RESTRAINTS),
      m_NeighborCalculator
      (
        assemble::AANeighborListContainerGeneratorProteinModel::AANeighborListGenerator( 11.4, 0, true, false)
      ),
      m_ExposureCalculator()
    {
    }

    //! @brief returns a pointer to a new RestraintDistanceEPR
    //! @return pointer to a new RestraintDistanceEPR
    RestraintDistanceEPR *RestraintDistanceEPR::Clone() const
    {
      return new RestraintDistanceEPR( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &RestraintDistanceEPR::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the scheme of this score
    //! @return the scheme of this score
    const std::string &RestraintDistanceEPR::GetScheme() const
    {
      return m_Scheme;
    }

    //! @brief returns default file name where the statistics and in consequence the energy potentials are read from
    //! @return default file name where the statistics and in consequence the energy potentials are read from
    const std::string &RestraintDistanceEPR::GetDefaultHistogramFilename()
    {
      static const std::string s_default_histogram_filename( "sl-cb_distances_geom.histograms");
      return s_default_histogram_filename;
    }

    //! @brief returns the default scheme of this score
    //! @return the default scheme of this score
    const std::string &RestraintDistanceEPR::GetDefaultScheme()
    {
      static const std::string s_default_scheme( "epr_distance_geom");
      return s_default_scheme;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief scores the agreement of a given protein model with distance measurements obtained from an EPR experiment
    //! @param PROTEIN_MODEL protein model for which to compute the agreement
    //! @return normalized agreement score with -1 being the best and 0 being the worst agreement
    double RestraintDistanceEPR::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // compute the exposure of the residues in the protein model
      const assemble::AANeighborListContainer neighbor_lists( ( *m_NeighborCalculator)( PROTEIN_MODEL));

      // sum up the energy scores for the restraints
      double sum_score( 0.0);
      for
      (
        util::ShPtrVector< restraint::AtomDistance>::const_iterator it( m_Restraints.Begin()), it_end( m_Restraints.End());
        it != it_end;
        ++it
      )
      {
        // get the restraint data
        const restraint::AtomDistance &restraint( **it);
        const double exp_distance( restraint.GetDistance()->GetDistance());
        const restraint::DataPairwise &data( restraint.GetData());
        const double model_distance( data.EuclidianDistance( PROTEIN_MODEL));

        // get the coordinates of the Ca and Cb atoms of the spin labeling sites
        const biol::AABase &aa_first( *data.First()->LocateAA( PROTEIN_MODEL));
        const biol::AABase &aa_second( *data.Second()->LocateAA( PROTEIN_MODEL));
        const biol::Atom &ca_first( aa_first.GetAtom( biol::GetAtomTypes().CA));
        const linal::Vector3D &ca_first_coord( ca_first.GetCoordinates());
        const linal::Vector3D &cb_first_coord( data.First()->Locate( PROTEIN_MODEL));
        const biol::Atom &ca_second( aa_second.GetAtom( biol::GetAtomTypes().CA));
        const linal::Vector3D &ca_second_coord( ca_second.GetCoordinates());
        const linal::Vector3D &cb_second_coord( data.Second()->Locate( PROTEIN_MODEL));

        // compute the projection angles between the CaCb and CaCa vectors of the spin labeling sites
        const double proj_first_cos( linal::ProjAngleCosinus( ca_first_coord, cb_first_coord, ca_second_coord));
        const double proj_second_cos( linal::ProjAngleCosinus( ca_second_coord, cb_second_coord, ca_first_coord));

        // compute the neighbor counts of the spin labeling sites
        const double nc_first( m_ExposureCalculator( neighbor_lists.Find( aa_first)->second));
        const double nc_second( m_ExposureCalculator( neighbor_lists.Find( aa_second)->second));

        // compute the energy score and add it to the sum score for all restraints
        const double aggregate( ( ( nc_first * proj_first_cos) + ( nc_second * proj_second_cos)) / 8.0);
        const double dsl_dbb( exp_distance - model_distance);
        const double score( ( *m_EnergyFunction)( aggregate, dsl_dbb));
        sum_score += score;
      }

      return sum_score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads in members from stream
    //! @param ISTREAM stream to read members from
    //! @return returns the input stream
    std::istream &RestraintDistanceEPR::Read( std::istream &ISTREAM)
    {
      // read members from input stream
      io::Serialize::Read( m_Scheme, ISTREAM);
      io::Serialize::Read( m_EnergyFunction, ISTREAM);
      io::Serialize::Read( m_Restraints, ISTREAM);
      io::Serialize::Read( m_NeighborCalculator, ISTREAM);
      io::Serialize::Read( m_ExposureCalculator, ISTREAM);

      return ISTREAM;
    }

    //! @brief write members into a stream
    //! @param OSTREAM stream to write members into
    //! @INDENT number of indentations to use
    //! @return returns the output stream
    std::ostream &RestraintDistanceEPR::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members into output stream
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_EnergyFunction, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_Restraints, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_NeighborCalculator, OSTREAM, INDENT) << std::endl;
      io::Serialize::Write( m_ExposureCalculator, OSTREAM, INDENT);

      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief creates the energy function used to score the agreement of model with the EPR data from the given histogram
    //! @param HISTOGRAM_FILENAME histogram from which to create the energy function
    //! @return shared pointer to the energy function created from the given histogram
    util::ShPtr< math::BicubicSpline> RestraintDistanceEPR::CreateEnergyFunction( const std::string &HISTOGRAM_FILENAME)
    {
      // read in the histogram
      math::Histogram2D histogram;
      io::IFStream read;
      io::File::MustOpenIFStream( read, Score::AddHistogramPath( HISTOGRAM_FILENAME));
      read >> histogram;
      io::File::CloseClearFStream( read);
      histogram.NormalizeRows();

      // create the energy function
      const math::Range< double> distance_range( -14.0, 14.0);
      util::ShPtr< math::BicubicSpline> sp_energy_function
      (
        new math::BicubicSpline( EnergyDistribution::EPRDistance( histogram, distance_range))
      );

      return sp_energy_function;
    }

  } // namespace score
} // namespace bcl
