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
#include "fold/bcl_fold_mutate_protein_model_sse_move.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "math/bcl_math_mutate_result.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelSSEMove::s_Instance
    (
      util::Enumerated< math::MutateInterface< assemble::ProteinModel> >::AddInstance( new MutateProteinModelSSEMove())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelSSEMove::MutateProteinModelSSEMove() :
      m_Locator(),
      m_Placement(),
      m_Scheme( GetStaticClassName< MutateProteinModelSSEMove>())
    {
    }

    //! @brief constructor from a locator and placement and a scheme
    //! @param LOCATOR ShPtr to LocatorInterface to be used
    //! @param PLACEMENT ShPtr to PlacementInterface to be used
    //! @param SCHEME Scheme to be used
    MutateProteinModelSSEMove::MutateProteinModelSSEMove
    (
      const find::LocatorCriteriaInterface
      <
        util::SiPtr< const assemble::SSE>, assemble::DomainInterface, assemble::DomainInterface
      > &LOCATOR,
      const PlacementInterface< assemble::SSE, assemble::ProteinModel> &PLACEMENT,
      const std::string &SCHEME
    ) :
      m_Locator( LOCATOR),
      m_Placement( PLACEMENT),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor from ShPtrs to a locator and placement and a scheme
    //! @param SP_LOCATOR ShPtr to LocatorInterface to be used
    //! @param SP_PLACEMENT ShPtr to PlacementInterface to be used
    //! @param SCHEME Scheme to be used
    MutateProteinModelSSEMove::MutateProteinModelSSEMove
    (
      const util::ShPtr
      <
        find::LocatorCriteriaInterface
        <
          util::SiPtr< const assemble::SSE>, assemble::DomainInterface, assemble::DomainInterface
        >
      > &SP_LOCATOR,
      const util::ShPtr< PlacementInterface< assemble::SSE, assemble::ProteinModel> > &SP_PLACEMENT,
      const std::string &SCHEME
    ) :
      m_Locator( *SP_LOCATOR),
      m_Placement( *SP_PLACEMENT),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateProteinModelSSEMove
    MutateProteinModelSSEMove *MutateProteinModelSSEMove::Clone() const
    {
      return new MutateProteinModelSSEMove( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateProteinModelSSEMove::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get scheme
    //! @return scheme
    const std::string &MutateProteinModelSSEMove::GetScheme() const
    {
      return m_Scheme;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that moves an SSE
    //! @param PROTEIN_MODEL protein model interest
    //! @return MutateResult with ProteinModel after the mutate
    math::MutateResult< assemble::ProteinModel> MutateProteinModelSSEMove::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // create the ShPtr to model
      util::ShPtr< assemble::ProteinModel> sp_model;

      // try to locate sse
      const util::SiPtr< const assemble::SSE> sp_located_sse( m_Locator->Locate( PROTEIN_MODEL, PROTEIN_MODEL));

      // no sse located
      if( !sp_located_sse.IsDefined())
      {
        return math::MutateResult< assemble::ProteinModel>( sp_model, *this);
      }

      // find placement for the sse
      const storage::Pair< math::TransformationMatrix3D, bool> place( m_Placement->Place( *sp_located_sse, PROTEIN_MODEL));

      // is there a placement defined
      if( !place.Second())
      {
        return math::MutateResult< assemble::ProteinModel>( sp_model, *this);
      }

      // create the transformation
      math::TransformationMatrix3D transform( sp_located_sse->GetOrientation());
      transform.Invert();
      transform( place.First());

      // create a copy of the sse to tranform
      util::ShPtr< assemble::SSE> sp_sse( sp_located_sse->Clone());

      // transform to new location
      sp_sse->Transform( transform);

      // create new protein model
      sp_model = util::ShPtr< assemble::ProteinModel>( PROTEIN_MODEL.Clone());

      // replace
      BCL_Assert( sp_model->Replace( sp_sse), "cannot replace sse: " + sp_sse->GetIdentification());

      // end
      return math::MutateResult< assemble::ProteinModel>( sp_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace fold
} // namespace bcl
