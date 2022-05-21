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
#include "assemble/bcl_assemble_aa_neighbor_vector.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_list.h"
#include "biol/bcl_biol_aa_base.h"
#include "biol/bcl_biol_atom.h"
#include "linal/bcl_linal_vector_3d_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! @brief returns the default lower and higher thresholds
    //! @return the default lower and higher thresholds
    const math::Range< double> &AANeighborVector::GetDefaultThresholdLowHigh()
    {
      // static range to hold values
      static const math::Range< double> s_default_threshold_low_high( 3.3, 11.1);

      // end
      return s_default_threshold_low_high;
    }

    //! @brief returns default file where the statistics and in consequence the energy potentials are read from
    //! @return default file where the statistics and in consequence the energy potentials are read from
    const std::string &AANeighborVector::GetDefaultHistogramFilename()
    {
      // static string
      static const std::string s_default_histogram_filename( "neighbor_vector_sasa_all_chains_membrane.histograms");

      // end
      return s_default_histogram_filename;
    }

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &AANeighborVector::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "aaneighvector");

      // end
      return s_default_scheme;

    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! @param SCHEME scheme to be used
    AANeighborVector::AANeighborVector( const std::string &SCHEME) :
      m_Scheme( SCHEME),
      m_ThresholdLowHigh( GetDefaultThresholdLowHigh()),
      m_MininalSequenceSeparation( s_DefaultMinimalSequenceSeparation)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new AANeighborVector object copied from this one
    AANeighborVector *AANeighborVector::Clone() const
    {
      return new AANeighborVector( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AANeighborVector::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &AANeighborVector::GetScheme() const
    {
      return m_Scheme;
    }

    //! @brief min and max exposure measure
    //! @return range in which exposure can be
    const math::Range< double> &AANeighborVector::GetRange() const
    {
      static const math::Range< double> s_range( 0, 1);
      return s_range;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculates the neighbor vector for a given amino acid within the specified sequence separation and threshold
    //! @param AA_NEIGHBOR_LIST AANeighborList to be used for a single amino acid
    //! @param THRESHOLD_LOW_HIGH lower and higher distance thresholds
    //! @param NORMALIZE_BY_NEIGHBOR_COUNT flag to normalize by neighbor count
    //! @return the neighbor vector for a given amino acid within the specified sequence separation and threshold
    linal::Vector3D
    AANeighborVector::NeighborVector
    (
      const AANeighborList &AA_NEIGHBOR_LIST,
      const math::Range< double> &THRESHOLD_LOW_HIGH,
      const bool NORMALIZE_BY_NEIGHBOR_COUNT
    )
    {
      // check that the given neighbor list has the proper parameters
      const double dist_cutoff( AA_NEIGHBOR_LIST.GetDistanceCutoff());
      const double thresh_high( THRESHOLD_LOW_HIGH.GetMax());
      BCL_Assert
      (
        dist_cutoff == thresh_high, "given neighbor list does not have the proper distance cutoff : " +
        util::Format()( dist_cutoff) + " != " + util::Format()( thresh_high)
      );

      // current first sidechain atom coordinate
      const linal::Vector3D &current_cb_coordinate( AA_NEIGHBOR_LIST.GetCenterAminoAcid()->GetFirstSidechainAtom().GetCoordinates());

      // start with the assumption that the current AA has no neighbors
      double neighbor_count( 0.0);

      // vector variable used to hold the sums of all vectors between the current AA and neighboring AAs
      linal::Vector3D vector_sum;

      //loop over all AAs to calculate env_radius
      for
      (
        AANeighborList::const_iterator
          dist_aa_itr( AA_NEIGHBOR_LIST.Begin()),
          dist_aa_itr_end( AA_NEIGHBOR_LIST.End());
        dist_aa_itr != dist_aa_itr_end;
        ++dist_aa_itr
      )
      {
        // determine the weight assigned to this neighbor
        const double neighbor_weight
        (
          coord::NeighborWeight( dist_aa_itr->Second(), THRESHOLD_LOW_HIGH.GetMin(), THRESHOLD_LOW_HIGH.GetMax())
        );

        neighbor_count += neighbor_weight;

        // normalize each neighbor vector such that each gets an initial weight of 1 (before the neighbor_weight)
        vector_sum += ( ( dist_aa_itr->First()->GetFirstSidechainAtom().GetCoordinates() - current_cb_coordinate).Normalize() * neighbor_weight);

      } // aa loop

      if( NORMALIZE_BY_NEIGHBOR_COUNT && ( neighbor_count > 0.0))
      {
        vector_sum /= neighbor_count;
      }

      // end
      return vector_sum;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate neighbor vector potential for a given amino acid and its AANeighborList
    //! @return neighbor vector potential for a given amino acid and its AANeighborList
    double AANeighborVector::operator()( const AANeighborList &AA_NEIGHBOR_LIST) const
    {
      // calculate RadiusEnvirmonment for current aminoacid
      return NeighborVector( AA_NEIGHBOR_LIST, m_ThresholdLowHigh, true).Norm();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &AANeighborVector::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Scheme          , ISTREAM);
      io::Serialize::Read( m_ThresholdLowHigh, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @return returns the output stream
    std::ostream &AANeighborVector::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Scheme          , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ThresholdLowHigh, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
