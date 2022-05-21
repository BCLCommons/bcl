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
#include "contact/bcl_contact_statistics.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container_generator_protein_model.h"
#include "biol/bcl_biol_aa_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {

    //! @brief StatisticType as string
    //! @param STATISTIC_TYPE the StatisticType
    //! @return the string for the StatisticType
    const std::string &Statistics::GetStatisticTypeDescriptor( const StatisticType &STATISTIC_TYPE)
    {
      static const std::string s_descriptors[] =
      {
        "NumberContacts",
        "NumberContactsShort",
        "NumberContactsMid",
        "NumberContactsLong",
        "RatioContactsShort",
        "RatioContactsMid",
        "RatioContactsLong",
        GetStaticClassName< StatisticType>()
      };

      return s_descriptors[ STATISTIC_TYPE];
    }

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> Statistics::s_Instance
    (
      GetObjectInstances().AddInstance( new Statistics())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Statistics::Statistics() :
      m_SequenceSeparationRange(),
      m_UseRatio(),
      m_NeighborGenerator()
    {
    }

    //! @brief constructor from a statistics type to be calculated and cache boolean
    //! @param STATISTIC_TYPE StatisticType to be calculated
    //! @param CACHE whether a neighbor generator with cache should be used
    Statistics::Statistics
    (
      const StatisticType &STATISTIC_TYPE,
      const bool CACHE
    ) :
      m_SequenceSeparationRange(),
      m_UseRatio(),
      m_NeighborGenerator
      (
        assemble::AANeighborListContainerGeneratorProteinModel::AANeighborListGenerator
        (
          g_ContactCbDistanceCutoff, g_ContactMinSequenceSeparation, true, CACHE
        )
      )
    {
      // switch over the statistic type given
      switch( STATISTIC_TYPE)
      {
        case e_NumberContacts:
        {
          m_SequenceSeparationRange = GetDefaultSequenceSeparationRange();
          m_UseRatio = false;
          break;
        }
        case e_NumberContactsShort:
        {
          m_SequenceSeparationRange = GetDefaultSequenceSeparationShortRange();
          m_UseRatio = false;
          break;
        }
        case e_NumberContactsMid:
        {
          m_SequenceSeparationRange = GetDefaultSequenceSeparationMidRange();
          m_UseRatio = false;
          break;
        }
        case e_NumberContactsLong:
        {
          m_SequenceSeparationRange = GetDefaultSequenceSeparationLongRange();
          m_UseRatio = false;
          break;
        }
        case e_RatioContactsShort:
        {
          m_SequenceSeparationRange = GetDefaultSequenceSeparationShortRange();
          m_UseRatio = true;
          break;
        }
        case e_RatioContactsMid:
        {
          m_SequenceSeparationRange = GetDefaultSequenceSeparationMidRange();
          m_UseRatio = true;
          break;
        }
        case e_RatioContactsLong:
        {
          m_SequenceSeparationRange = GetDefaultSequenceSeparationLongRange();
          m_UseRatio = true;
          break;
        }
        case s_NumberStatisticType:
        {
          break;
        }
      }
    }

    //! @brief constructor from a sequence separation range, whether to calculate ratio and whether to use cache
    //! @param SEQUENCE_SEPARATION_RANGE Sequence separation range
    //! @param CALCULATE_RATIO boolean indicating whether just the counts or the ratio wrt to total number of contacts should be calculated
    //! @param CACHE whether a neighbor generator with cache should be used
    Statistics::Statistics
    (
      const math::Range< size_t> &SEQUENCE_SEPARATION_RANGE,
      const bool CALCULATE_RATIO,
      const bool CACHE
    ) :
      m_SequenceSeparationRange( SEQUENCE_SEPARATION_RANGE),
      m_UseRatio( CALCULATE_RATIO),
      m_NeighborGenerator
      (
        assemble::AANeighborListContainerGeneratorProteinModel::AANeighborListGenerator
        (
          g_ContactCbDistanceCutoff, g_ContactMinSequenceSeparation, true, CACHE
        )
      )
      {
      }

    //! @brief Clone function
    //! @return pointer to new Statistics
    Statistics *Statistics::Clone() const
    {
      return new Statistics( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Statistics::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  /////////////////
  // operations  //
  /////////////////

    //! @brief calculate the contact statistics value for the given AANeighborListContainer
    //! @param NEIGHBOR_CONTAINER AANeighborListContainer of interest
    //! @return calculated statistics value for the given AANeighborListContainer
    double Statistics::Calculate( const assemble::AANeighborListContainer &NEIGHBOR_CONTAINER) const
    {
      // make a copy of the container
      assemble::AANeighborListContainer container
      (
        NEIGHBOR_CONTAINER,
        g_ContactCbDistanceCutoff,
        g_ContactMinSequenceSeparation,
        true
      );

      // if ratio is asked
      if( m_UseRatio)
      {
        return GetContactsRatio( m_SequenceSeparationRange, container);
      }
      // otherwise just return counts
      else
      {
        return double( GetNumberContacts( m_SequenceSeparationRange, container) / 2);
      }
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate the contact statistics value for the given ProteinModel
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @return calculuated statistics value for the given model
    double Statistics::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // calculate the neighbor list for the given ProteinModel and calculate and return the statistic value
      return Calculate( m_NeighborGenerator->operator()( PROTEIN_MODEL));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Statistics::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SequenceSeparationRange, ISTREAM);
      io::Serialize::Read( m_UseRatio, ISTREAM);
      io::Serialize::Read( m_NeighborGenerator, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Statistics::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SequenceSeparationRange, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_UseRatio, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NeighborGenerator, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculate the number of contacts within the given sequence separation range for the given NeighborListContainer
    //! @param SEQUENCE_SEPARATION_RANGE Sequence separation range to be used
    //! @param NEIGHBOR_CONTAINER AANeighborListContainer of interest
    //! @return the number of contacts in the given NeighborContainer
    size_t Statistics::GetNumberContacts
    (
      const math::Range< size_t> &SEQUENCE_SEPARATION_RANGE,
      const assemble::AANeighborListContainer &NEIGHBOR_CONTAINER
    )
    {
      // initialize count
      size_t count( 0);

      // iterate over the neighbor list container
      for
      (
        assemble::AANeighborListContainer::const_iterator
          itr( NEIGHBOR_CONTAINER.Begin()), itr_end( NEIGHBOR_CONTAINER.End());
        itr != itr_end; ++itr
      )
      {
        // create a reference on the center amino acid
        const biol::AABase &center_aa( *itr->second.GetCenterAminoAcid());

        // iterate over the neighbors listed
        for
        (
          assemble::AANeighborList::const_iterator neigh_itr( itr->second.Begin()), neigh_itr_end( itr->second.End());
          neigh_itr != neigh_itr_end; ++neigh_itr
        )
        {
          // if it's a different chain or the sequence separation is within the range then skip over
          if
          (
            center_aa.GetChainID() != neigh_itr->First()->GetChainID() ||
            SEQUENCE_SEPARATION_RANGE.IsWithin( biol::SequenceSeparation( center_aa, *neigh_itr->First()))
          )
          {
            // then increment count
            ++count;
          }
        }
      }

      // end
      return count;
    }

    //! @brief calculate and return the ratio of all contacts that are within given sequence separation range
    //! @param SEQUENCE_SEPARATION_RANGE sequence separation range to be used for calculating ratio
    //! @param NEIGHBOR_CONTAINER AANeighborListContainer of interest
    //! @return the ratio of all contacts that are within given sequence separation range
    double Statistics::GetContactsRatio
    (
      const math::Range< size_t> &SEQUENCE_SEPARATION_RANGE,
      const assemble::AANeighborListContainer &NEIGHBOR_CONTAINER
    )
    {
      // calculate the counts with the sequence separation
      const double counts_range( GetNumberContacts( SEQUENCE_SEPARATION_RANGE, NEIGHBOR_CONTAINER));

      // divide by total number of counts and return
      return 100.0 * counts_range / NEIGHBOR_CONTAINER.GetNumberNeighbors();
    }

  } // namespace contact
} // namespace bcl
