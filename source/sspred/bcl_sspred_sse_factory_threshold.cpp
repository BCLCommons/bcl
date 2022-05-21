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
#include "sspred/bcl_sspred_sse_factory_threshold.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_pool.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sspred
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> SSEFactoryThreshold::s_Instance
    (
      GetObjectInstances().AddInstance( new SSEFactoryThreshold())
    );

    //! @brief returns default map for sse thresholds
    //! @return default map for sse thresholds
    const storage::Map< biol::SSType, double> &SSEFactoryThreshold::GetDefaultThresholdsMap()
    {
      // initialize empty default map
      static storage::Map< biol::SSType, double> s_default_map;

      // set the default map values
      if( s_default_map.IsEmpty())
      {
        s_default_map[ biol::GetSSTypes().HELIX] = 0.5;
        s_default_map[ biol::GetSSTypes().STRAND] = 0.5;
      }

      // return default map
      return s_default_map;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SSEFactoryThreshold::SSEFactoryThreshold() :
      m_SSMethod( GetMethods().e_Undefined),
      m_SSEThresholds( GetDefaultThresholdsMap()),
      m_ExtendSSEs( false)
    {
    }

    //! @brief constructor from specified SSMethod
    //! @param SS_METHOD set of ssmethod to be used for generating SSEs from AASequence
    //! @param SSE_THRESHOLDS threshold values to classify ss predictions as helix or strand or not
    //! @param EXTEND_SSES whether or not to add systematically extended copies of sses in pool
    SSEFactoryThreshold::SSEFactoryThreshold
    (
      const Method &SS_METHOD,
      const storage::Map< biol::SSType, double> &SSE_THRESHOLDS,
      const bool EXTEND_SSES
    ) :
      m_SSMethod( SS_METHOD),
      m_SSEThresholds( SSE_THRESHOLDS),
      m_ExtendSSEs( EXTEND_SSES)
    {
    }

    //! @brief virtual copy constructor
    SSEFactoryThreshold *SSEFactoryThreshold::Clone() const
    {
      return new SSEFactoryThreshold( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSEFactoryThreshold::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that returns a set of SSEs for the given AASequence
    //! @brief SEQUENCE AASequence from which the SSEPool is going to be built
    //! @return SSEPool built from provided SEQUENCE
    assemble::SSEPool
    SSEFactoryThreshold::operator()( const biol::AASequence &SEQUENCE) const
    {
      // initialize an SSEPool
      assemble::SSEPool sse_pool;

      // iterate over methods
      for( biol::SSTypes::const_iterator ss_itr( biol::GetSSTypes().Begin()), ss_itr_end( biol::GetSSTypes().COIL.GetIterator()); ss_itr != ss_itr_end; ++ss_itr)
      {
        // find the threshold
        const storage::Map< biol::SSType, double>::const_iterator thresh_itr( m_SSEThresholds.Find( *ss_itr));

        // no threshold given - not considered
        if( thresh_itr == m_SSEThresholds.End())
        {
          continue;
        }

        // get the SSEPool created using this method and insert to sse_pool
        sse_pool.InsertElements( IdentifySSERegions( SEQUENCE, m_SSMethod, *ss_itr, thresh_itr->second));
      }

      // if extend sses flag is set
      if( m_ExtendSSEs)
      {
        // add extended versions of the sses to the pool
        sse_pool = SystematicallyExtendSSEs( sse_pool, SEQUENCE);
      }

      // return
      return sse_pool;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEFactoryThreshold::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SSMethod     , ISTREAM);
      io::Serialize::Read( m_SSEThresholds, ISTREAM);
      io::Serialize::Read( m_ExtendSSEs   , ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &SSEFactoryThreshold::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SSMethod     , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SSEThresholds, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ExtendSSEs   , OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief identifies the regions with most likely probability to be in SSE and returns them in a ssepool
    //! @param SEQUENCE AASequence from which SSEs will be identified
    //! @param METHOD the prediction method considered
    //! @param SS_TYPE SSType of interest
    //! @param THRESHOLD the prediction threshold, that the prediction needs to be above to be considered for that sstype
    //! @return SSEPool composed of SSEs generated from SEQUENCE
    assemble::SSEPool SSEFactoryThreshold::IdentifySSERegions
    (
      const biol::AASequence &SEQUENCE,
      const Method &METHOD,
      const biol::SSType &SS_TYPE,
      const double THRESHOLD
    )
    {
      // initialize SSEPool
      assemble::SSEPool sse_pool;

      // iterate over each residue in the sequence
      biol::AASequence::const_iterator aa_itr( SEQUENCE.Begin()), aa_itr_end( SEQUENCE.End());
      while( aa_itr != aa_itr_end)
      {
        // locate the end of the elongated SSERegion
        const biol::AASequence::const_iterator sse_end_itr
        (
          ElongateSSERegion( aa_itr, aa_itr_end, SS_TYPE, METHOD, THRESHOLD)
        );

        // could not find any residue
        if( aa_itr == sse_end_itr)
        {
          // move on to next aa
          ++aa_itr;
          continue;
        }

        // create the sse from the subsequence of SEQUENCE using iterators and ss_type
        util::ShPtr< assemble::SSE> this_sse
        (
          new assemble::SSE
          (
            biol::AASequence
            (
              util::ShPtrVector< biol::AABase>( aa_itr, sse_end_itr),
              SEQUENCE.GetChainID(),
              SEQUENCE.GetFastaHeader()
            ),
            SS_TYPE
          )
        );

        // set to origin and ideal coordinates
        this_sse->SetToIdealConformationAtOrigin();

        // insert this sse into the pool
        sse_pool.Insert( this_sse);

        // update the aa_itr to sse_end_itr
        aa_itr = sse_end_itr;
      }

      // return
      return sse_pool;
    }

    //! @brief elongates the SSE region of type SS_TYPE and beginning at AA_ITR
    //! @param AA_ITR iterator to the beginning of a new SSE region
    //! @param AA_ITR_END end iterator for the sequence
    //! @param SS_TYPE SSE type for which the window is going to be maximized
    //! @param METHOD the prediction method considered
    //! @param THRESHOLD the prediction threshold, that the prediction needs to be above to be extended
    //! @return iterator to the end of the region which SSE can be elongated
    biol::AASequence::const_iterator
    SSEFactoryThreshold::ElongateSSERegion
    (
      const biol::AASequence::const_iterator &AA_ITR,
      const biol::AASequence::const_iterator &AA_ITR_END,
      const biol::SSType &SS_TYPE,
      const Method &METHOD,
      const double THRESHOLD
    )
    {
      // put the end iterator one past the beginning iterator
      biol::AASequence::const_iterator aa_itr( AA_ITR);

      // loop until you hit end of sequence or find the end of this sse
      while( aa_itr != AA_ITR_END)
      {
        // prediction for that amino acid
        const util::SiPtr< const MethodInterface> sp_prediction( ( *aa_itr)->GetSSPrediction( METHOD));

        // is there a prediction for that amino acid
        if( !sp_prediction.IsDefined())
        {
          break;
        }

        // is the prediction above the threshold
        if( sp_prediction->GetThreeStatePrediction()( SS_TYPE) < THRESHOLD)
        {
          break;
        }

        // go to next aa
        ++aa_itr;
      }

      // return
      return aa_itr;
    }

    //! @brief systematically extends the SSEs in the given POOL and inserts them into POOL
    //! @param POOL SSEPool for which SSEs will be extended
    //! @param SEQUENCE AASequence on which the extensions will be based
    //! @return SSEPool which has SSEs in the original pool and the extended SSEs
    assemble::SSEPool SSEFactoryThreshold::SystematicallyExtendSSEs( assemble::SSEPool &POOL, const biol::AASequence &SEQUENCE) const
    {
      // make a copy of the pool
      assemble::SSEPool new_pool( POOL);

      // iterate over sses in the pool
      for
      (
          assemble::SSEPool::const_iterator sse_itr( POOL.Begin()), sse_itr_end( POOL.End());
        sse_itr != sse_itr_end; ++sse_itr
      )
      {
        // initialize empty set of sses to store the extended copies
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThan> extended_sses;

        // identify aa in front and behind the sse
        const util::ShPtr< biol::AABase> aa_in_front( SEQUENCE.GetAA( ( *sse_itr)->GetFirstAA()->GetSeqID() - 2));
        const util::ShPtr< biol::AABase> aa_in_back( SEQUENCE.GetAA( ( *sse_itr)->GetLastAA()->GetSeqID()));

        // make a copy of the sse
        util::ShPtr< assemble::SSE> front_extended_sse( ( *sse_itr)->Clone());
        // prepend one amino acid to the front
        front_extended_sse->Prepend( *aa_in_front);
        // insert this front extended aa into extended_sses
        extended_sses.Insert( front_extended_sse);

        // make a copy of the sse
        util::ShPtr< assemble::SSE> back_extended_sse( ( *sse_itr)->Clone());
        // append one amino acid ro the end of the sse
        back_extended_sse->Append( *aa_in_back);
        // insert this back extended sse into extended_sses
        extended_sses.Insert( back_extended_sse);

        // make a copy of the sse
        util::ShPtr< assemble::SSE> front_and_back_extended_sse( ( *sse_itr)->Clone());
        // prepend one amino acid to the front
        front_and_back_extended_sse->Prepend( *aa_in_front);
        // to the sse that already has one aa extended to the front, add one aa to the back
        front_and_back_extended_sse->Append( *aa_in_back);
        // insert this doubly extended sse into extended_sses
        extended_sses.Insert( front_and_back_extended_sse);

        // extend this sse and insert the resultant set of sses into SSEPool
        new_pool.InsertElements( extended_sses.Begin(), extended_sses.End());
      }

      // return new pool
      return new_pool;
    }

  } // namespace sspred
} // namespace bcl
