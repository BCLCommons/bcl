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
#include "fold/bcl_fold_mutate_protein_model_sse_remove.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "find/bcl_find_locator_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelSSERemove::s_Instance
    (
      util::Enumerated< math::MutateInterface< assemble::ProteinModel> >::AddInstance( new MutateProteinModelSSERemove())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelSSERemove::MutateProteinModelSSERemove() :
      m_Locator(),
      m_Scheme( GetStaticClassName< MutateProteinModelSSERemove>())
    {
    }

    //! @brief constructor from a locator and a scheme
    //! @param LOCATOR function that chooses the sse
    //! @param SCHEME Scheme to be used
    MutateProteinModelSSERemove::MutateProteinModelSSERemove
    (
      const util::ShPtr< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> > &LOCATOR,
      const std::string &SCHEME
    ) :
      m_Locator( *LOCATOR),
      m_Scheme( SCHEME)
    {
    }

    //! @brief clone
    MutateProteinModelSSERemove *MutateProteinModelSSERemove::Clone() const
    {
      return new MutateProteinModelSSERemove( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateProteinModelSSERemove::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this mutate
    //! @return the scheme for this mutate
    const std::string &MutateProteinModelSSERemove::GetScheme() const
    {
      return m_Scheme;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
    //! @param PROTEIN_MODEL protein model interest
    //! @return MutateResult with ProteinModel after the mutate
    math::MutateResult< assemble::ProteinModel> MutateProteinModelSSERemove::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // static empty model
      static util::ShPtr< assemble::ProteinModel> s_empty_model;

      // if there is only one sse in the protein model do not remove it and return
      if( PROTEIN_MODEL.GetNumberSSEs() == 1)
      {
        // warn user
        BCL_MessageVrb( "Only one SSE in protein model, skipping remove");

        // return empty result
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // choose a random sse and make copy
      const util::SiPtr< const assemble::SSE> located_sse( m_Locator->Locate( PROTEIN_MODEL));

      // if sse cannot be located
      if( !located_sse.IsDefined())
      {
        // warn user
        BCL_MessageVrb( "could not find sse in model to be removed: " + util::Format()( *m_Locator));

        // return empty result
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // report selected sse from pool to be added
      BCL_MessageVrb( "sse to be removed: " + located_sse->GetIdentification());

      // copy the protein model
      util::ShPtr< assemble::ProteinModel> new_model( PROTEIN_MODEL.Clone());

      // replace the sse with the mutated hardcopy of the same sse
      new_model->Remove( *located_sse);

      // end
      return math::MutateResult< assemble::ProteinModel>( new_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateProteinModelSSERemove::Read( std::istream &ISTREAM)
    {
      // read members
      // return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &MutateProteinModelSSERemove::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
