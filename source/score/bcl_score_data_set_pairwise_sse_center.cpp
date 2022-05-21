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
#include "assemble/bcl_assemble_sse_pool.h"
#include "restraint/bcl_restraint_data_set_pairwise.h"
#include "restraint/bcl_restraint_locator_coordinates_first_side_chain_atom.h"
#include "score/bcl_score_data_set_pairwise_sse_center.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> DataSetPairwiseSSECenter::s_Instance
    (
      GetObjectInstances().AddInstance( new DataSetPairwiseSSECenter())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &DataSetPairwiseSSECenter::GetDefaultScheme()
    {
      static const std::string s_scheme( "sse_center");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    DataSetPairwiseSSECenter::DataSetPairwiseSSECenter( const std::string &SCHEME) :
      m_Data(),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor taking pool
    //! @param SSE_POOL the pool the sse definitions will come from
    //! @param SCHEME optional scheme
    DataSetPairwiseSSECenter::DataSetPairwiseSSECenter( const assemble::SSEPool &SSE_POOL, const std::string &SCHEME) :
      m_Data(),
      m_Scheme( SCHEME)
    {
      m_Data = CalculateResiduePositionWeight( SSE_POOL);
    }

    //! @brief Clone function
    //! @return pointer to new DataSetPairwiseSSECenter
    DataSetPairwiseSSECenter *DataSetPairwiseSSECenter::Clone() const
    {
      return new DataSetPairwiseSSECenter( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DataSetPairwiseSSECenter::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &DataSetPairwiseSSECenter::GetScheme() const
    {
      return m_Scheme;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate the score of a data set
    //! @param DATA data set to be scored
    //! @return the score of the current data set
    double DataSetPairwiseSSECenter::operator()( const restraint::DataSetPairwise &DATA) const
    {
      if( DATA.IsEmpty())
      {
        return 0;
      }

      double score( 0);

      // iterate over the data pairs to score them
      for
      (
        restraint::DataSetPairwise::const_iterator itr( DATA.Begin()), itr_end( DATA.End());
        itr != itr_end;
        ++itr
      )
      {
        storage::Map
        <
          util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, double,
          assemble::LocatorAtomCoordinatesInterface::PtrResidueLessThan
        >::const_iterator resi_itr_a( m_Data.Find( itr->First()));

        storage::Map
        <
          util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, double,
          assemble::LocatorAtomCoordinatesInterface::PtrResidueLessThan
        >::const_iterator resi_itr_b( m_Data.Find( itr->Second()));

        // if either residue could not found in the map
        if( resi_itr_a == m_Data.End() || resi_itr_b == m_Data.End())
        {
          ++score;
          continue;
        }

        double weight_a( resi_itr_a->second);
        double weight_b( resi_itr_b->second);

        double current_score( weight_a + weight_b);

        score += current_score;
      }

      score /= double( DATA.GetSize());

      --score;

      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DataSetPairwiseSSECenter::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Data, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DataSetPairwiseSSECenter::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Data, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculates the weight of every residue in helical sses according to its position relative to the ends
    //! @param SSE_POOL the pool of sses that will be used to calculate the weights
    //! @return map that has for every residue in helical sses the weight associated according to its position relative
    //!         to the end of the sse
    storage::Map
    <
      util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, double,
      assemble::LocatorAtomCoordinatesInterface::PtrResidueLessThan
    > DataSetPairwiseSSECenter::CalculateResiduePositionWeight( const assemble::SSEPool &SSE_POOL)
    {
      // the sses from the pool
      storage::Set< util::SiPtr< const assemble::SSE>, assemble::SSELessThanNoOverlap> sses
      (
        SSE_POOL.GetRandomNonOverlappingSet()
      );

      // map to hold the residues and their weights
      storage::Map
      <
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, double,
        assemble::LocatorAtomCoordinatesInterface::PtrResidueLessThan
      > residue_weights;

      // iterate through the sses
      for
      (
        storage::Set< util::SiPtr< const assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
          sse_itr( sses.Begin()), sse_itr_end( sses.End());
          sse_itr != sse_itr_end;
          ++sse_itr
      )
      {
        // true if the sse is not helix, go to next sse
        if( ( *sse_itr)->GetType() != biol::GetSSTypes().HELIX)
        {
          continue;
        }

        // iterate over the residues of the sse to assign their weights
        for
        (
          biol::AASequence::const_iterator aa_itr( ( *sse_itr)->Begin()), aa_itr_end( ( *sse_itr)->End());
          aa_itr != aa_itr_end;
          ++aa_itr
        )
        {

          util::ShPtr< assemble::LocatorAtomCoordinatesInterface> current_resi
          (
            new restraint::LocatorCoordinatesFirstSideChainAtom( ( *aa_itr)->GetChainID(), ( *aa_itr)->GetSeqID())
          );

          const double weight( CalculateWeight( **sse_itr, **aa_itr));

          bool success
          (
            residue_weights.Insert
            (
              std::pair< util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, double>( current_resi, weight)
            ).second
          );

          BCL_Assert( success, "was not able to insert residue " + current_resi->GetIdentification());
        }
      }

      return residue_weights;
    }

    //! @brief calculates the weight of the given residue according to its position away from the end of the given sse
    //! @param SSE the sse whose ends will be used to calculate the weight
    //! @param AA_BASE the residue of interest whose weight will be calculated
    //! @return double which is the weight of the residue
    double DataSetPairwiseSSECenter::CalculateWeight( const assemble::SSE &SSE, const biol::AABase &AA_BASE)
    {
      double aa_resi_num( AA_BASE.GetSeqID());

      double sse_middle_resi( ( SSE.GetLastAA()->GetSeqID() - SSE.GetFirstAA()->GetSeqID()) / 2.0 + SSE.GetFirstAA()->GetSeqID());

      const double weight( math::Absolute( aa_resi_num - sse_middle_resi) / SSE.GetSize());

      BCL_MessageDbg
      (
        "resi is " + util::Format()( aa_resi_num) + " sse is " + SSE.GetIdentification() +
        " weight is " + util::Format()( weight) + " middle resi is " + util::Format()( sse_middle_resi)
        + " sse size is " + util::Format()( SSE.GetSize())
      );

      return weight;
    }

  } // namespace score

} // namespace bcl
