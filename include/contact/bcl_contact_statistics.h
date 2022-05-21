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

#ifndef BCL_CONTACT_STATISTICS_H_
#define BCL_CONTACT_STATISTICS_H_

// include the namespace header
#include "bcl_contact.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_function_interface_serializable.h"
#include "math/bcl_math_range.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Statistics
    //! @brief helper class to calculate contact related statistics from a given ProteinModel
    //! @details This class can calculate counts or the ratio of contacts for a given ProteinModel. By default it
    //! considers amino acids with at least 6 residues sequence separation as possible contacts. Using the convenience
    //! constructor which takes the enum StatisticType, a variety of statistics can be calculated such as counts/ratios
    //! of short/mid/long range contacts
    //!
    //! @see @link example_contact_statistics.cpp @endlink
    //! @author karakam
    //! @date Sep 22, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Statistics :
      public math::FunctionInterfaceSerializable< assemble::ProteinModel, double>
    {

    public:

    ///////////
    // types //
    ///////////

      //! @brief type of statistics to be calculated
      enum StatisticType
      {
        e_NumberContacts,
        e_NumberContactsShort,
        e_NumberContactsMid,
        e_NumberContactsLong,
        e_RatioContactsShort,
        e_RatioContactsMid,
        e_RatioContactsLong,
        s_NumberStatisticType
      };

      //! @brief StatisticType as string
      //! @param STATISTIC_TYPE the StatisticType
      //! @return the string for the StatisticType
      static const std::string &GetStatisticTypeDescriptor( const StatisticType &STATISTIC_TYPE);

    private:

    //////////
    // data //
    //////////

      //! range of sequence separation to be used when identifying contacts
      math::Range< size_t> m_SequenceSeparationRange;

      //! boolean indicating whether just the counts or the ratio wrt to total number of contacts should be calculated
      bool m_UseRatio;

      //! pointer to the neighbor list container generator
      util::ShPtr< math::FunctionInterfaceSerializable< assemble::ProteinModel, assemble::AANeighborListContainer> > m_NeighborGenerator;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Statistics();

      //! @brief constructor from a statistics type to be calculated and cache boolean
      //! @param STATISTIC_TYPE StatisticType to be calculated
      //! @param CACHE whether a neighbor generator with cache should be used
      Statistics
      (
        const StatisticType &STATISTIC_TYPE,
        const bool CACHE
      );

      //! @brief constructor from a sequence separation range, whether to calculate ratio and whether to use cache
      //! @param SEQUENCE_SEPARATION_RANGE Sequence separation range
      //! @param CALCULATE_RATIO boolean indicating whether just the counts or the ratio wrt to total number of contacts should be calculated
      //! @param CACHE whether a neighbor generator with cache should be used
      Statistics
      (
        const math::Range< size_t> &SEQUENCE_SEPARATION_RANGE,
        const bool CALCULATE_RATIO,
        const bool CACHE
      );

      //! @brief Clone function
      //! @return pointer to new Statistics
      Statistics *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get sequence separation range used
      //! @return sequence separation range used
      const math::Range< size_t> &GetSequenceSeparationRange() const
      {
        return m_SequenceSeparationRange;
      }

      //! @brief return whether ratio is calculated
      //! @return whether ratio is calculated
      bool GetUseRatio() const
      {
        return m_UseRatio;
      }

    /////////////////
    // operations  //
    /////////////////

      //! @brief calculate the contact statistics value for the given AANeighborListContainer
      //! @param NEIGHBOR_CONTAINER AANeighborListContainer of interest
      //! @return calculated statistics value for the given AANeighborListContainer
      double Calculate( const assemble::AANeighborListContainer &NEIGHBOR_CONTAINER) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief calculate the contact statistics value for the given ProteinModel
      //! @param PROTEIN_MODEL ProteinModel of interest
      //! @return calculated statistics value for the given model
      double operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief calculate the number of contacts within the given sequence distance range for the given NeighborListContainer
      //! @param SEQUENCE_SEPARATION_RANGE Sequence separation range to be used
      //! @param NEIGHBOR_CONTAINER AANeighborListContainer of interest
      //! @return the number of contacts in the given NeighborContainer
      static size_t GetNumberContacts
      (
        const math::Range< size_t> &SEQUENCE_SEPARATION_RANGE,
        const assemble::AANeighborListContainer &NEIGHBOR_CONTAINER
      );

      //! @brief calculate and return the ratio of all contacts that are within given sequence separation range
      //! @param SEQUENCE_SEPARATION_RANGE sequence separation range to be used for calculating ratio
      //! @param NEIGHBOR_CONTAINER AANeighborListContainer of interest
      //! @return the ratio of all contacts that are within given sequence separation range
      static double GetContactsRatio
      (
        const math::Range< size_t> &SEQUENCE_SEPARATION_RANGE,
        const assemble::AANeighborListContainer &NEIGHBOR_CONTAINER
      );

    }; // class Statistics

  } // namespace contact
} // namespace bcl

#endif // BCL_CONTACT_STATISTICS_H_
