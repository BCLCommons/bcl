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
#include "fold/bcl_fold_mutate_protein_model_compress.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelCompress::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateProteinModelCompress())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelCompress::MutateProteinModelCompress() :
      m_CompressionRange( math::Range< double>( 0.97, 1.00)),
      m_Scheme( GetStaticClassName< MutateProteinModelCompress>())
    {
    }

    //! @brief constructor from a specified compression factor and a scheme
    //! @param COMPRESSION_FACTOR compression factor
    //! @param SCHEME Scheme to be used
    MutateProteinModelCompress::MutateProteinModelCompress
    (
      const double COMPRESSION_FACTOR,
      const std::string &SCHEME
    ) :
      m_CompressionRange( COMPRESSION_FACTOR, COMPRESSION_FACTOR),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor from a compression factor range
    //! @param COMPRESSION_FACTOR_RANGE compression factor range composed of an lower and an upper value, and a scheme
    //! @param SCHEME Scheme to be used
    MutateProteinModelCompress::MutateProteinModelCompress
    (
      const math::Range< double> &COMPRESSION_FACTOR_RANGE,
      const std::string &SCHEME
    ) :
      m_CompressionRange( COMPRESSION_FACTOR_RANGE),
      m_Scheme( SCHEME)
    {
    }

    //! @brief clone
    MutateProteinModelCompress *MutateProteinModelCompress::Clone() const
    {
      return new MutateProteinModelCompress( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateProteinModelCompress::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns compression range
    //! @return compression range
    const math::Range< double> MutateProteinModelCompress::GetCompressionRange() const
    {
      return m_CompressionRange;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
    //! @param PROTEIN_MODEL ProteinModel which will be mutated
    //! @return MutateResult ProteinModel after mutation
    math::MutateResult< assemble::ProteinModel> MutateProteinModelCompress::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // static empty result
      static util::ShPtr< assemble::ProteinModel> s_empty_model;

      // make sure the protein model at least has 2 SSEs
      if( PROTEIN_MODEL.GetSSEs().GetSize() < 2)
      {
        // end
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // determine the scaling factor
      // if it's just one value then use it, if it's a range then pick a random factor from the specified range
      const double compression_factor
      (
        m_CompressionRange.GetMin() == m_CompressionRange.GetMax() ?
          m_CompressionRange.GetMin() : random::GetGlobalRandom().Double( m_CompressionRange)
      );

      // copy the protein model
      util::ShPtr< assemble::ProteinModel> new_model( new assemble::ProteinModel());

      // calculate the center of mass
      const linal::Vector3D model_center( PROTEIN_MODEL.GetCenterOfSSEs());

      // iterate over all chains
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator chain_itr( PROTEIN_MODEL.GetChains().Begin()),
          chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end; ++chain_itr
      )
      {
        // initialize a new chain
        util::ShPtr< assemble::Chain> sp_chain( new assemble::Chain( ( *chain_itr)->GetSequence()));

        // iterate over all the SSEs in the chain
        for
        (
          storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
            sse_itr( ( *chain_itr)->GetData().Begin()), sse_itr_end( ( *chain_itr)->GetData().End());
          sse_itr != sse_itr_end; ++sse_itr
        )
        {
          // make a copy of the SSE
          util::ShPtr< assemble::SSE> new_sse( ( *sse_itr)->Clone());

          // get the translation axis
          const linal::Vector3D translation_axis( ( *sse_itr)->GetCenter() - model_center);

          // apply the translation along the calculated axis
          new_sse->Translate( translation_axis * ( compression_factor - double( 1.0)));

          // insert the SSE
          sp_chain->Insert( new_sse);
        }

        // insert this chain into the model
        new_model->Insert( sp_chain);
      }

      // end
      return math::MutateResult< assemble::ProteinModel>( new_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateProteinModelCompress::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_CompressionRange, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &MutateProteinModelCompress::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_CompressionRange, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
