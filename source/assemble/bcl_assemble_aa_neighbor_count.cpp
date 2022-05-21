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
#include "assemble/bcl_assemble_aa_neighbor_count.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_list.h"
#include "biol/bcl_biol_aa_base.h"

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
    const math::Range< double> &AANeighborCount::GetDefaultThresholdLowHigh()
    {
      // static range to hold thresholds
      static const math::Range< double> s_default_threshold_low_high( 4.0, 11.4);

      // end
      return s_default_threshold_low_high;
    }

    //! @brief default minimal sequence separation
    const size_t AANeighborCount::GetDefaultMinimalSequenceSeparation()
    {
      // default minimal sequence separation
      static const size_t s_default_minimal_sequence_separation( 2);

      return s_default_minimal_sequence_separation;
    }

    //! @brief returns default file where the statistics and in consequence the energy potentials are read from
    //! @return default file where the statistics and in consequence the energy potentials are read from
    const std::string &AANeighborCount::GetDefaultHistogramFilename()
    {
      // static string
      static const std::string s_default_histogram_filename( "neighbor_count_sasa_all_chains_membrane.histograms");

      // end
      return s_default_histogram_filename;
    }

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &AANeighborCount::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "aaneigh");

      // end
      return s_default_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! @param SCHEME scheme to be used
    AANeighborCount::AANeighborCount( const std::string &SCHEME) :
      m_Scheme( SCHEME),
      m_ThresholdLowHigh( GetDefaultThresholdLowHigh()),
      m_MininalSequenceSeparation( GetDefaultMinimalSequenceSeparation())
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new AANeighborCount object copied from this one
    AANeighborCount *AANeighborCount::Clone() const
    {
      return new AANeighborCount( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AANeighborCount::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &AANeighborCount::GetScheme() const
    {
      return m_Scheme;
    }

    //! @brief min and max exposure measure
    //! @return range in which exposure can be
    const math::Range< double> &AANeighborCount::GetRange() const
    {
      static const math::Range< double> s_range( 0, 50);
      return s_range;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate neighbor potential for a given amino acid and its AANeighborList
    //! @return neighbor potential for a given amino acid and its AANeighborList
    double AANeighborCount::operator()( const AANeighborList &AA_NEIGHBOR_LIST) const
    {
      // check that the given neighbor list has the proper parameters
      BCL_Assert
      (
        AA_NEIGHBOR_LIST.GetDistanceCutoff() == m_ThresholdLowHigh.GetMax(),
        "given neighbor list does not have the proper distance cutoff or sequence separation: " +
        util::Format()( AA_NEIGHBOR_LIST.GetDistanceCutoff()) + " == " + util::Format()( m_ThresholdLowHigh.GetMax())
      );

      // start with the assumption that the current AA has no neighbors
      double neighbor_count( 0.0);

      //loop over all AAs to calculate env_radius
      for
      (
        AANeighborList::const_iterator
          aa_itr( AA_NEIGHBOR_LIST.Begin()), aa_itr_end( AA_NEIGHBOR_LIST.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        neighbor_count +=
          coord::NeighborWeight( aa_itr->Second(), m_ThresholdLowHigh.GetMin(), m_ThresholdLowHigh.GetMax());
      }

      // return number of neighbors for that AMINO_ACID
      return neighbor_count;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &AANeighborCount::Read( std::istream &ISTREAM)
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
    std::ostream &AANeighborCount::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ThresholdLowHigh, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
