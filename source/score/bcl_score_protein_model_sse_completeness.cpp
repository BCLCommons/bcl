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
#include "score/bcl_score_protein_model_sse_completeness.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "biol/bcl_biol_atom.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProteinModelSSECompleteness::s_Instance
    (
      util::Enumerated< ProteinModel>::AddInstance( new ProteinModelSSECompleteness( false))
    );

    const util::SiPtr< const util::ObjectInterface> ProteinModelSSECompleteness::s_AAInstance
    (
      util::Enumerated< ProteinModel>::AddInstance( new ProteinModelSSECompleteness( true))
    );
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief returns a pointer to a new ProteinModelSSECompleteness
    //! @return pointer to a new ProteinModelSSECompleteness
    ProteinModelSSECompleteness *ProteinModelSSECompleteness::Clone() const
    {
      return new ProteinModelSSECompleteness( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &ProteinModelSSECompleteness::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the scheme of this score
    //! @return the scheme of this score
    const std::string &ProteinModelSSECompleteness::GetScheme() const
    {
      static const std::string s_default_scheme( "sse_completeness"), s_aa_scheme( "sse_completeness_aas");
      return m_CountAAs ? s_aa_scheme : s_default_scheme;
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ProteinModelSSECompleteness::GetAlias() const
    {
      return GetScheme();
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProteinModelSSECompleteness::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        m_CountAAs
        ? "# of aas in sses that could still be added"
        : "# of sses that could still be added"
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief scores the completeness of the given protein model
    //! @detail the completeness is calculated by computing the ratio of residues with defined coordinates and the
    //! total number of residues in the sequence.
    //! @param PROTEIN_MODEL protein model for which to compute the completeness
    //! @return completeness of the given protein model with -1.0 being complete and 0.0 being empty
    double ProteinModelSSECompleteness::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // get the pool from the given ProteinModel and make sure it is valid
      const util::ShPtr< assemble::SSEPool> sp_pool
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Pool)
      );
      if( !sp_pool.IsDefined())
      {
        return 0.0;
      }

      // find all sse's not yet in the model
      util::SiPtrList< const assemble::SSE> non_overlapping
      (
        sp_pool->GetNonOverlappingSSEs( PROTEIN_MODEL)
      );

      size_t total_sse_sizes( 0);
      if( m_CountAAs)
      {
        for( auto itr( non_overlapping.Begin()), itr_end( non_overlapping.End()); itr != itr_end; ++itr)
        {
          total_sse_sizes += ( *itr)->GetSize();
        }
      }
      // Return # of missing aas
      return m_CountAAs ? total_sse_sizes : non_overlapping.GetSize();
    }

  } // namespace score
} // namespace bcl
