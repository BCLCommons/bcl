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
#include "fold/bcl_fold_mutate_protein_model_sse_add_multiple.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain.h"
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
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelSSEAddMultiple::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateProteinModelSSEAddMultiple())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelSSEAddMultiple::MutateProteinModelSSEAddMultiple() :
      m_SSEPool(),
      m_Placement(),
      m_Scheme( GetStaticClassName< MutateProteinModelSSEAddMultiple>())
    {
    }

    //! @brief constructor from a SSEPool, a PlacementDomainInterface and a scheme
    //! @param SSE_POOL ShPtr to pool of SSEs to be used
    //! @param PLACEMENT placement for the SSEs in the pool into the model
    //! @param SCHEME Scheme to be used
    MutateProteinModelSSEAddMultiple::MutateProteinModelSSEAddMultiple
    (
      const util::ShPtr< assemble::SSEPool> &SSE_POOL,
      const PlacementDomainInterface &PLACEMENT,
      const std::string &SCHEME
    ) :
      m_SSEPool( SSE_POOL),
      m_Placement( PLACEMENT.Clone()),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor from a ShPtr to a SSEPool, a ShPtr to PlacementDomainInterface and a scheme
    //! @param SSE_POOL ShPtr to pool of SSEs to be used
    //! @param SP_PLACEMENT ShPtr to placement for the SSEs in the pool into the model
    //! @param SCHEME Scheme to be used
    MutateProteinModelSSEAddMultiple::MutateProteinModelSSEAddMultiple
    (
      const util::ShPtr< assemble::SSEPool> &SSE_POOL,
      const util::ShPtr< PlacementDomainInterface> &SP_PLACEMENT,
      const std::string &SCHEME
    ) :
      m_SSEPool( SSE_POOL),
      m_Placement( SP_PLACEMENT),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateProteinModelSSEAddMultiple
    MutateProteinModelSSEAddMultiple *MutateProteinModelSSEAddMultiple::Clone() const
    {
      return new MutateProteinModelSSEAddMultiple( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateProteinModelSSEAddMultiple::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator taking a protein model and returning a mutate object of protein model type
    //! @param PROTEIN_MODEL protein model of interest
    //! @return MutateResult with ProteinModel after the mutate
    math::MutateResult< assemble::ProteinModel> MutateProteinModelSSEAddMultiple::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // static empty model
      static util::ShPtr< assemble::ProteinModel> s_empty_model( new assemble::ProteinModel());

      // return empty model if no SSEs in pool
      if( m_SSEPool->GetSize() == 0)
      {
        BCL_MessageStd( "Skipping Adding Multiple SSEs: no SSEs in pool");

        // return empty result
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // initialize temp SSEs set
      const storage::Set< util::SiPtr< const assemble::SSE>, assemble::SSELessThanNoOverlap> sses
      (
        m_SSEPool->GetRandomNonOverlappingSet()
      );

      // initialize ShPtr Set
      storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap> sp_sses;

      // iterate through the SSEs
      for
      (
        storage::Set< util::SiPtr< const assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
          sse_itr( sses.Begin()), sse_itr_end( sses.End());
        sse_itr != sse_itr_end; ++sse_itr
      )
      {
        sp_sses.Insert( util::ShPtr< assemble::SSE>( ( *sse_itr)->Clone()));
      }

      // get the map of transformation matrices from the placement
      storage::Pair
      <
        storage::Map< util::SiPtr< const assemble::SSE>, math::TransformationMatrix3D>,
        bool
      > transformations
      (
        m_Placement->Place( assemble::Domain( sp_sses))
      );

      // if the placement was not successful
      if( !transformations.Second())
      {
        BCL_MessageStd
        (
          "The provided transformations were not successful, skipping the add multiple SSEs"
        );

        // return empty result
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // create a new protein model using the chain information
      util::ShPtr< assemble::ProteinModel> new_model( new assemble::ProteinModel( PROTEIN_MODEL.GetEmptyChains()));

      // iterate through the SSEs
      for
      (
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
          sse_itr( sp_sses.Begin()), sse_itr_end( sp_sses.End());
        sse_itr != sse_itr_end; ++sse_itr
      )
      {
        // create an iterator on the map
        storage::Map< util::SiPtr< const assemble::SSE>, math::TransformationMatrix3D>::iterator map_itr
        (
          transformations.First().Find( util::SiPtr< const assemble::SSE>( **sse_itr))
        );

        // if the iterator is valid
        if( map_itr != transformations.First().End())
        {
          // create temp SSE
          util::ShPtr< assemble::SSE> this_sse( ( *sse_itr)->Clone());

          // set this sse to origin
          this_sse->Transform( math::Inverse( this_sse->GetOrientation()));

          // apply the proper transformation
          this_sse->Transform( map_itr->second);

          // insert the SSE into the protein model
          new_model->Insert( this_sse);
        }
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
    std::istream &MutateProteinModelSSEAddMultiple::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SSEPool, ISTREAM);
      io::Serialize::Read( m_Placement, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &MutateProteinModelSSEAddMultiple::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SSEPool, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Placement, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
