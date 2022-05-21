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
#include "assemble/bcl_assemble_protein_model_inverter.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_sequence_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProteinModelInverter::s_Instance
    (
      GetObjectInstances().AddInstance( new ProteinModelInverter())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from whether to use cache and a cache size
    //! @param USE_CACHE whether to use cache, by default it is set to false
    //! @param MODEL_CACHE_SIZE Maximum number of models to keep in cache, default is 2
    //! @param COIL_CACHE_SIZE Maximum number of coils to keep in cache, default is 30
    ProteinModelInverter::ProteinModelInverter
    (
      const bool USE_CACHE,
      const size_t MODEL_CACHE_SIZE,
      const size_t COIL_CACHE_SIZE
    ) :
      m_UseCache( USE_CACHE),
      m_ModelCacheSize( MODEL_CACHE_SIZE),
      m_CoilCacheSize( COIL_CACHE_SIZE),
      m_ExtendedSequences(),
      m_Coils(),
      m_InvertedModels()
    {
    }

    //! @brief Clone function
    //! @return pointer to new ProteinModelInverter
    ProteinModelInverter *ProteinModelInverter::Clone() const
    {
      return new ProteinModelInverter( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinModelInverter::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief generates the description string for a given chain
    //! @param CHAIN Chain of interest
    //! @return description string
    std::string ProteinModelInverter::GenerateDescription( const Chain &CHAIN)
    {
      // initialize string with the address of the sequence
      std::string description( util::Format()( CHAIN.GetSequence().GetPointer()));

      // static variable to hold valid SSTypes set
      static const storage::Set< biol::SSType> valid_types( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND);

      // get SSEs
      const util::SiPtrVector< const SSE> sses( CHAIN.GetSSEs( valid_types));

      // iterate over chains
      for
      (
        util::SiPtrVector< const SSE>::const_iterator sse_itr( sses.Begin()), sse_itr_end( sses.End());
        sse_itr != sse_itr_end; ++sse_itr
      )
      {
        // generate description for chain
        description += "_" + util::Format()( ( *sse_itr)->GetFirstAA()->GetSeqID()) +
                       "_" + util::Format()( ( *sse_itr)->GetLastAA()->GetSeqID());
      }

      // end
      return description;
    }

    //! @brief generates the description string for a given model
    //! @param MODEL Model of interest
    //! @return description string
    std::string ProteinModelInverter::GenerateDescription( const ProteinModel &MODEL)
    {
      // initialize string
      std::string description;

      // iterate over chains
      for
      (
        util::ShPtrVector< Chain>::const_iterator
          chain_itr( MODEL.GetChains().Begin()), chain_itr_end( MODEL.GetChains().End());
        chain_itr != chain_itr_end; ++chain_itr
      )
      {
        // generate description for chain
        description += GenerateDescription( **chain_itr);
      }

      // end
      return description;
    }

    //! @brief returns whether any extended loop with the given begin and end SeqIDs exist for the given AASequence
    //! @param SP_SEQUENCE ShPtr to AASequence of interest
    //! @param BEGIN_SEQ_ID seqid for the first amino acid
    //! @param END_SEQ_ID seqid for the last amino acid
    //! @return whether any extended loop with the given begin and end SeqIDs exist for the given AASequence
    bool ProteinModelInverter::DoesContainCoil
    (
      const util::ShPtr< biol::AASequence> &SP_SEQUENCE,
      const size_t BEGIN_SEQ_ID,
      const size_t END_SEQ_ID
    ) const
    {
      // search for the sequence
      storage::Map< util::ShPtr< biol::AASequence>, util::ShPtrList< SSE> >::const_iterator seq_itr
      (
        m_Coils.Find( SP_SEQUENCE)
      );

      // if not found
      if( seq_itr == m_Coils.End())
      {
        return false;
      }

      // if sequence is found, then search and return if the SSE is found
      return
        std::find_if
        (
          seq_itr->second.Begin(),
          seq_itr->second.End(),
          SSECompareByIdentity( BEGIN_SEQ_ID, END_SEQ_ID, SP_SEQUENCE->GetChainID(), biol::GetSSTypes().COIL)
        ) != seq_itr->second.End();
    }

    //! @brief reset the cached storage
    void ProteinModelInverter::Reset()
    {
      // reset the members
      m_ExtendedSequences.Reset();
      m_Coils.Reset();
      m_InvertedModels.Reset();
    }

  ////////////////
  // operations //
  ///////////////

    //! @brief return inverted model for a given ProteinModel
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @return inverted ProteinModel
    util::ShPtr< ProteinModel> ProteinModelInverter::GetInvertedModel( const ProteinModel &PROTEIN_MODEL) const
    {
      // if cache is not used
      if( !m_UseCache)
      {
        return ConstructInvertedModel( PROTEIN_MODEL);
      }
      // otherwise generate model with cache
      else
      {
        return GetInvertedModelFromCache( PROTEIN_MODEL);
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProteinModelInverter::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_UseCache, ISTREAM);
      io::Serialize::Read( m_ModelCacheSize, ISTREAM);
      io::Serialize::Read( m_CoilCacheSize, ISTREAM);
      // do not read other data members since they should be regenerated
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProteinModelInverter::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_UseCache, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ModelCacheSize, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_CoilCacheSize, OSTREAM, INDENT);
      // do not write other data members since they can be regenerated
      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief get a cached coil using begin seqid and end seqid
    //! @param SP_SEQUENCE ShPtr to AASequence of the chain inverted
    //! @param BEGIN_SEQ_ID seqid for the first amino acid
    //! @param END_SEQ_ID seqid for the last amino acid
    //! @return coil with specified begin and end seqids
    util::ShPtr< SSE> ProteinModelInverter::GetCoilFromCache
    (
      const util::ShPtr< biol::AASequence> &SP_SEQUENCE,
      const size_t BEGIN_SEQ_ID,
      const size_t END_SEQ_ID
    ) const
    {
      // create reference on the cached coils
      util::ShPtrList< SSE> &coils( m_Coils[ SP_SEQUENCE]);

      // search for this SSE
      storage::List< util::ShPtr< SSE> >::const_iterator coil_itr
      (
        std::find_if
        (
          coils.Begin(),
          coils.End(),
          SSECompareByIdentity( BEGIN_SEQ_ID, END_SEQ_ID, SP_SEQUENCE->GetChainID(), biol::GetSSTypes().COIL)
        )
      );

      // if found
      if( coil_itr != coils.End())
      {
        // then return it
        return *coil_itr;
      }

      // otherwise create the coil
      util::ShPtr< SSE> sp_new_coil
      (
        new SSE
        (
          m_ExtendedSequences[ SP_SEQUENCE]->SubSequence( BEGIN_SEQ_ID - 1, END_SEQ_ID - BEGIN_SEQ_ID + 1),
          biol::GetSSTypes().COIL
        )
      );

      // insert this SSE and return it
      coils.PushBack( sp_new_coil);
      return sp_new_coil;
    }

    //! @brief returns a chain that contains extended loops containing all amino acids that are not in the given chain
    //! @param CHAIN Chain to be inverted
    //! @return inverted chain
    util::ShPtr< Chain> ProteinModelInverter::GetInvertedChainFromCache( const Chain &CHAIN) const
    {
      // create reference on the sequence of the Chain
      const util::ShPtr< biol::AASequence> &sp_sequence( CHAIN.GetSequence());

      // check to see if an extended sequence exists for this chain
      storage::Map< util::ShPtr< biol::AASequence>, util::ShPtr< biol::AASequence> >::const_iterator extended_itr
      (
        m_ExtendedSequences.Find( CHAIN.GetSequence())
      );

      // if it does not exists
      if( extended_itr == m_ExtendedSequences.End())
      {
        // construct the extended sequence and store it
        extended_itr =
          m_ExtendedSequences.Insert( std::make_pair( sp_sequence, ConstructExtendedSequence( *sp_sequence))).first;
      }

      // create a reference on the extended sequence
      const util::ShPtr< biol::AASequence> &sp_extended_sequence( extended_itr->second);

      // get SSEs
      util::SiPtrVector< const SSE> sses( CHAIN.GetSSEs());

      // initialize a vector to store the inverted SSEs
      util::ShPtrVector< SSE> coils;

      // if no SSEs meaning an empty chain
      if( sses.IsEmpty())
      {
        // insert only the sequence that covers the full sequence
        coils.PushBack( util::ShPtr< SSE>( new SSE( *sp_extended_sequence, biol::GetSSTypes().COIL)));

        // construct a chain and return it
        return util::ShPtr< Chain>( new Chain( sp_sequence, coils));
      }

      // initialize iterators
      util::SiPtrVector< const SSE>::const_iterator sse_itr( sses.Begin());
      const util::SiPtrVector< const SSE>::const_iterator sse_itr_end( sses.End());

      // initialize previous seqid
      int prev_seqid( 0);

      // iterate over the SSEs
      for( ; sse_itr != sse_itr_end; ++sse_itr)
      {
        // store the begin and end id
        const int begin_seqid( ( *sse_itr)->GetFirstAA()->GetSeqID());

        // if there is any region between this and the previous sse
        if( begin_seqid > prev_seqid + 1)
        {
          // get the inverted SSE and insert it
          coils.PushBack( GetCoilFromCache( sp_sequence, prev_seqid + 1, begin_seqid - 1));
        }
        // update prev_seqid with the seqid of the last amino acid in this SSE
        prev_seqid = ( *sse_itr)->GetLastAA()->GetSeqID();

      }

      // if theres is any loop left after the last SSE
      if( prev_seqid < sp_sequence->GetLastAA()->GetSeqID())
      {
        // insert one more SSE that covers the last prev_seqid to the end of the sequence
        coils.PushBack( GetCoilFromCache( sp_sequence, prev_seqid + 1, sp_sequence->GetSize()));
      }

      // initialize a new chain and return it
      return util::ShPtr< Chain>( new Chain( sp_sequence, coils));
    }

    //! @brief returns an inverted model
    //! @param MODEL ProteinModel to be inverted
    //! @return ShPtr to inverted model
    util::ShPtr< ProteinModel> ProteinModelInverter::GetInvertedModelFromCache( const ProteinModel &MODEL) const
    {
      // generate the string description to be used as key
      const std::string key( GenerateDescription( MODEL));

      // construct a predicate that matches the key
      storage::PairEqualFirst< std::string> predicate( key);

      // search for this predicate to see if there is already a model cached with that description
      storage::List< storage::Pair< std::string, util::ShPtr< ProteinModel> > >::const_iterator
        model_itr( std::find_if( m_InvertedModels.Begin(), m_InvertedModels.End(), predicate));

      // if found
      if( model_itr != m_InvertedModels.End())
      {
        // then return it directly
        return model_itr->Second();
      }

      // otherwise construct a new model and description pair
      storage::Pair< std::string, util::ShPtr< ProteinModel> > sp_new_key_model
      (
        key, util::ShPtr< ProteinModel>( new ProteinModel())
      );

      // iterate over chains
      for
      (
        util::ShPtrVector< Chain>::const_iterator
          chain_itr( MODEL.GetChains().Begin()), chain_itr_end( MODEL.GetChains().End());
        chain_itr != chain_itr_end; ++chain_itr
      )
      {
        // construct and insert a new chain
        sp_new_key_model.Second()->Insert( GetInvertedChainFromCache( **chain_itr));

        // create non-const reference on the coils list for this sequence
        util::ShPtrList< SSE> &coils( m_Coils[ ( *chain_itr)->GetSequence()]);

        // while the size of the coils is larger than the coil cache size
        while( coils.GetSize() > m_CoilCacheSize)
        {
          // remove the first element
          coils.PopFront();
        }
      }

      // set the protein model data
      sp_new_key_model.Second()->SetProteinModelData( MODEL.GetProteinModelData());

      // insert the new model into cache
      m_InvertedModels.PushBack( sp_new_key_model);

      // if model cache size is reached
      if( m_InvertedModels.GetSize() > m_ModelCacheSize)
      {
        // remove first element
        m_InvertedModels.PopFront();
      }

      // return new model
      return sp_new_key_model.Second();
    }

    //! @brief static function to generate an extended chain from a given AASequence
    //! @param AA_SEQUENCE AASequence of interest
    //! @return ShPtr to AASequence in the extended chain conformation
    util::ShPtr< biol::AASequence> ProteinModelInverter::ConstructExtendedSequence( const biol::AASequence &AA_SEQUENCE)
    {
      // construct transformation matrix that will move the sequence to be parallel to xy plane and at z coordinate 500
      // this is required in order to make sure for membrane protein, the inverted model is not in the membrane
      math::TransformationMatrix3D transform( coord::GetAxes().e_X, math::g_Pi / 2.0);
      transform( linal::Vector3D( 0.0, 0.0, 500.0));

      // calculate from the chain id the z transform
      linal::Vector3D z_translation( 0.0, 0.0, ( AA_SEQUENCE.GetChainID() - 'A') * 25);
      transform( z_translation);

      // hardcopy the sequence
      util::ShPtr< biol::AASequence> sp_new_sequence( AA_SEQUENCE.HardCopy());

      // idealize the sequence
      biol::AASequenceFactory::IdealizeSequence( *sp_new_sequence, biol::GetSSTypes().STRAND);

      // transform sequence
      sp_new_sequence->Transform( transform);

      // return the sequence
      return sp_new_sequence;
    }

    //! @brief static function to generate an inverted chain for a given Chain
    //! @param CHAIN Chain of interest
    //! @return ShPtr to inverted Chain
    util::ShPtr< Chain> ProteinModelInverter::ConstructInvertedChain( const Chain &CHAIN)
    {
      // create reference on the sequence
      const util::ShPtr< biol::AASequence> &sp_sequence( CHAIN.GetSequence());

      // construct extended sequence
      util::ShPtr< biol::AASequence> sp_extended_sequence( ConstructExtendedSequence( *sp_sequence));

      // get SSEs from the Chain
      util::SiPtrVector< const SSE> sses( CHAIN.GetSSEs());

      // initialize a vector to store the inverted SSEs
      util::ShPtrVector< SSE> inverted_sses;

      // if no SSEs meaning an empty chain
      if( sses.IsEmpty())
      {
        // insert only the sequence that covers the full sequence
        inverted_sses.PushBack( util::ShPtr< SSE>( new SSE( *sp_extended_sequence, biol::GetSSTypes().COIL)));

        // construct a chain and return it
        return util::ShPtr< Chain>( new Chain( sp_sequence, inverted_sses));
      }

      // initialize iterators
      util::SiPtrVector< const SSE>::const_iterator sse_itr( sses.Begin());
      const util::SiPtrVector< const SSE>::const_iterator sse_itr_end( sses.End());

      // initialize previous seqid
      int prev_seqid( 0);

      // iterate over the SSEs
      for( ; sse_itr != sse_itr_end; ++sse_itr)
      {
        // store the begin and end id
        const int begin_seqid( ( *sse_itr)->GetFirstAA()->GetSeqID());

        // if there is any region between this and the previous sse
        if( begin_seqid > prev_seqid + 1)
        {
          // construct the SSE
          util::ShPtr< SSE> sp_new_sse
          (
            new SSE
            (
              sp_extended_sequence->SubSequence( prev_seqid, begin_seqid - prev_seqid - 1), biol::GetSSTypes().COIL
            )
          );

          // insert it into the model
          inverted_sses.PushBack( sp_new_sse);
        }
        // update prev_seqid with the seqid of the last amino acid in this SSE
        prev_seqid = ( *sse_itr)->GetLastAA()->GetSeqID();
      }

      // if there is any loop left after the last SSE
      if( prev_seqid < sp_sequence->GetLastAA()->GetSeqID())
      {
        // construct the SSE
        util::ShPtr< SSE> sp_new_sse
        (
          new SSE
          (
            sp_extended_sequence->SubSequence( prev_seqid, CHAIN.GetSequence()->GetSize() - prev_seqid),
            biol::GetSSTypes().COIL
          )
        );

        // insert it into the model
        inverted_sses.PushBack( sp_new_sse);
      }

      // initialize new chain
      util::ShPtr< Chain> sp_new_chain( new Chain( sp_sequence, inverted_sses));

      // return new chain
      return sp_new_chain;
    }

    //! @brief static function to generate an inverted model for a given ProteinModel
    //! @param MODEL ProteinMOdel of interest
    //! @return ShPtr to inverted ProteinMOdel
    util::ShPtr< ProteinModel> ProteinModelInverter::ConstructInvertedModel( const ProteinModel &MODEL)
    {
      // otherwise construct a new model and description pair
      util::ShPtr< ProteinModel> sp_new_model( new ProteinModel());

      // iterate over chains
      for
      (
        util::ShPtrVector< Chain>::const_iterator
          chain_itr( MODEL.GetChains().Begin()), chain_itr_end( MODEL.GetChains().End());
        chain_itr != chain_itr_end; ++chain_itr
      )
      {
        // construct inverted chain and insert it into the model
        sp_new_model->Insert( ConstructInvertedChain( **chain_itr));
      }

      // set the protein model data
      sp_new_model->SetProteinModelData( MODEL.GetProteinModelData());

      // end
      return sp_new_model;
    }

  } // namespace assemble
} // namespace bcl
