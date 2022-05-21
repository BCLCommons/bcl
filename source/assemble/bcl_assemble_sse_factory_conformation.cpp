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
#include "assemble/bcl_assemble_sse_factory_conformation.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_chain.h"
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "biol/bcl_biol_dssp.h"
#include "biol/bcl_biol_membrane.h"
#include "math/bcl_math_mutate_result.h"
#include "sspred/bcl_sspred_ci_phi_psi.h"
#include "storage/bcl_storage_object_nd_hash_map.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> SSEFactoryConformation::s_Instance
    (
      GetObjectInstances().AddInstance( new SSEFactoryConformation())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SSEFactoryConformation::SSEFactoryConformation() :
      m_MinSSELengths( 5, 3, 0)
    {
    }

    //! @brief constructor from provided minimal sse lengths
    //! @param MIN_SSE_LENGTHS minimum size of SSEs to be considered
    SSEFactoryConformation::SSEFactoryConformation( const storage::VectorND< 3, size_t> &MIN_SSE_LENGTHS) :
      m_MinSSELengths( MIN_SSE_LENGTHS)
    {
    }

    //! @brief virtual copy constructor
    SSEFactoryConformation *SSEFactoryConformation::Clone() const
    {
      return new SSEFactoryConformation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSEFactoryConformation::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that returns a set of SSEs for the given AASequence
    //! @brief SEQUENCE AASequence from which the SSEPool is going to be built
    //! @return SSEPool built from provided SEQUENCE
    SSEPool
    SSEFactoryConformation::operator()( const biol::AASequence &SEQUENCE) const
    {
      // initialize sse_pool
      SSEPool sse_pool;

      // collection of identified sses
      Domain identified_sses;

      // sequence needs to have at least s_MinimalSequenceLength amino acids
      if( SEQUENCE.GetSize() < s_MinimalSequenceLength)
      {
        return sse_pool;
      }

      ProteinModel model;
      util::ShPtr< Chain> tmp_chain( new Chain( util::CloneToShPtr( SEQUENCE)));
      tmp_chain->Insert( util::ShPtr< SSE>( new SSE( SEQUENCE, biol::GetSSTypes().COIL)));
      model.Insert( tmp_chain);

      const biol::DSSP dssp;
      math::MutateResult< ProteinModel> result( dssp( model));
      storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap> dssp_sses
      (
        result.GetArgument().IsDefined()
        ? result.GetArgument()->GetChains().FirstElement()->GetData()
        : storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap>()
      );

      // copy the sequence
      biol::AASequence seq_copy( SEQUENCE);

      sspred::CIPhiPsi().Calculate( seq_copy, util::SiPtr< const biol::Membrane>());

      util::ShPtr< SSE> current_sse;

      // find sstype for each amino acid and generate long sequences of a certain sstype
      for
      (
        biol::AASequence::const_iterator
          current_itr( SEQUENCE.Begin()),
          end_itr( SEQUENCE.End()), copy_itr( seq_copy.Begin());
        current_itr != end_itr;
        ++current_itr, ++copy_itr
      )
      {
        // calculate phi and psi and store it in map
        util::SiPtr< const sspred::MethodInterface> ciphipsi_pred
        (
          ( *copy_itr)->GetSSPrediction( sspred::GetMethods().e_CIPhiPsi)
        );
        const biol::SSType current_sstype( ciphipsi_pred.IsDefined() ? ciphipsi_pred->GetOneStateSSPrediction() : biol::GetSSTypes().COIL);

        // there was no sses started yet
        if( !current_sse.IsDefined())
        {
          current_sse =
            util::ShPtr< SSE>
            (
              new SSE( biol::AASequence( ( *current_itr)->GetAAClass(), 0, SEQUENCE.GetChainID()), current_sstype)
            );
        }

        // elongate existing sse of identical sstype
        if( current_sse->GetType() == current_sstype)
        {
          current_sse->PushBack( current_itr->HardCopy());
        }

        // different sstypes - insert sse into domain and create new sse with that aminoacid
        else
        {
          // if sse is not long enough, set it to type coil, otherwise calculate body
          if
          (
            (
              ( current_sse->GetType() == biol::GetSSTypes().HELIX) &&
              ( current_sse->GetSize() < m_MinSSELengths( biol::GetSSTypes().HELIX))
            )
            ||
            (
              ( current_sse->GetType() == biol::GetSSTypes().STRAND) &&
              ( current_sse->GetSize() < m_MinSSELengths( biol::GetSSTypes().STRAND))
            )
          )
          {
            current_sse->SetType( biol::GetSSTypes().COIL);
          }
          // else if sse is long enough
          else
          {
            current_sse->SetGeometry();
          }

          if( current_sse->GetType()->IsStructured())
          {
            // insert current sse
            current_sse->SetChainID( SEQUENCE.GetChainID());
            auto itr_dssp( dssp_sses.Find( current_sse));
            while( itr_dssp != dssp_sses.End())
            {
              dssp_sses.RemoveElement( itr_dssp);
              itr_dssp = dssp_sses.Find( current_sse);
            }
            identified_sses.Insert( current_sse);
          }

          // create new current sse for current aa
          current_sse = util::ShPtr< SSE>
          (
            new SSE( biol::AASequence( ( *current_itr)->GetAAClass(), 0, SEQUENCE.GetChainID()), current_sstype)
          );
          current_sse->PushBack( current_itr->HardCopy());
        }
      }

      // this is the same routine, when all aas were processed since the last sse is still not inserted in the sse list
      // if sse is not long enough, set it to type coil, otherwise calculate body
      if
      (
        (
          ( current_sse->GetType() == biol::GetSSTypes().HELIX) &&
          ( current_sse->GetSize() < m_MinSSELengths( biol::GetSSTypes().HELIX))
        )
        ||
        (
          ( current_sse->GetType() == biol::GetSSTypes().STRAND) &&
          ( current_sse->GetSize() < m_MinSSELengths( biol::GetSSTypes().STRAND))
        )
      )
      {
        current_sse->SetType( biol::GetSSTypes().COIL);
      }
      else
      {
        current_sse->SetGeometry();
      }

      if( current_sse->GetType()->IsStructured())
      {
        // insert current sse
        current_sse->SetChainID( SEQUENCE.GetChainID());
        auto itr_dssp( dssp_sses.Find( current_sse));
        while( itr_dssp != dssp_sses.End())
        {
          dssp_sses.RemoveElement( itr_dssp);
          itr_dssp = dssp_sses.Find( current_sse);
        }
        identified_sses.Insert( current_sse);
      }

      if( !dssp_sses.IsEmpty())
      {
        for( auto itr_dssp( dssp_sses.Begin()), itr_dssp_end( dssp_sses.End()); itr_dssp != itr_dssp_end; ++itr_dssp)
        {
          if
          (
            (
              ( ( *itr_dssp)->GetType() == biol::GetSSTypes().HELIX) &&
              ( ( *itr_dssp)->GetSize() >= m_MinSSELengths( biol::GetSSTypes().HELIX))
            )
            ||
            (
              ( ( *itr_dssp)->GetType() == biol::GetSSTypes().STRAND) &&
              ( ( *itr_dssp)->GetSize() >= m_MinSSELengths( biol::GetSSTypes().STRAND))
            )
          )
          {
            util::ShPtr< SSE> current_sse
            (
              new SSE
              (
                biol::AASequence( SEQUENCE.GetFirstAA()->GetAAClass(), 0, SEQUENCE.GetChainID()),
                ( *itr_dssp)->GetType()
              )
            );
            for
            (
              auto itr_sse( ( *itr_dssp)->GetData().Begin()), itr_sse_end( ( *itr_dssp)->GetData().End());
              itr_sse != itr_sse_end;
              ++itr_sse
            )
            {
              current_sse->PushBack( *itr_sse);
            }
            current_sse->SetGeometry();
            identified_sses.Insert( current_sse);
          }
        }
      }

      // insert the sses in the domain to the pool
      sse_pool.InsertElements( identified_sses.GetData().Begin(), identified_sses.GetData().End());

      // return set of sses
      return sse_pool;
    }

  ////////////////
  // operations //
  ////////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEFactoryConformation::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_MinSSELengths, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &SSEFactoryConformation::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_MinSSELengths, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
