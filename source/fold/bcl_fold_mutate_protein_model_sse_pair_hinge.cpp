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
#include "fold/bcl_fold_mutate_protein_model_sse_pair_hinge.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "coord/bcl_coord_movable_eccentric.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelSSEPairHinge::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateProteinModelSSEPairHinge())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelSSEPairHinge::MutateProteinModelSSEPairHinge() :
      m_Collector(),
      m_Move(),
      m_MoveHinge( false),
      m_Scheme( GetStaticClassName< MutateProteinModelSSEPairHinge>())
    {
    }

    //! @brief constructor from a CollectorInterface, MoveInterface and a scheme
    //! @param COLLECTOR function that chooses the sses
    //! @param MOVE function that performs the move on the sse
    //! @param MOVE_HINGE whether the hinge should be also moved
    //! @param SCHEME Scheme to be used
    MutateProteinModelSSEPairHinge::MutateProteinModelSSEPairHinge
    (
      const find::CollectorInterface< storage::List< storage::VectorND< 2, util::SiPtr< const assemble::SSE> > >, assemble::DomainInterface> &COLLECTOR,
      const coord::MoveInterface &MOVE,
      const bool MOVE_HINGE,
      const std::string &SCHEME
    ) :
      m_Collector( COLLECTOR.Clone()),
      m_Move( MOVE.Clone()),
      m_MoveHinge( MOVE_HINGE),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor from ShPtrs to CollectorInterface, MoveInterface and a scheme
    //! @param SP_COLLECTOR ShPtr to function that chooses the sses
    //! @param SP_MOVE ShPtr to function that performs the move on the sse
    //! @param MOVE_HINGE whether the hinge should be also moved
    //! @param SCHEME Scheme to be used
    MutateProteinModelSSEPairHinge::MutateProteinModelSSEPairHinge
    (
      const util::ShPtr< find::CollectorInterface< storage::List< storage::VectorND< 2, util::SiPtr< const assemble::SSE> > >, assemble::DomainInterface> > &SP_COLLECTOR,
      const util::ShPtr< coord::MoveInterface> &SP_MOVE,
      const bool MOVE_HINGE,
      const std::string &SCHEME
    ) :
      m_Collector( *SP_COLLECTOR),
      m_Move( SP_MOVE),
      m_MoveHinge( MOVE_HINGE),
      m_Scheme( SCHEME)
    {
    }

    //! @brief clone
    MutateProteinModelSSEPairHinge *MutateProteinModelSSEPairHinge::Clone() const
    {
      return new MutateProteinModelSSEPairHinge( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateProteinModelSSEPairHinge::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
    //! @param PROTEIN_MODEL ProteinModel which will be mutated
    //! @return MutateResult ProteinModel after mutation
    math::MutateResult< assemble::ProteinModel> MutateProteinModelSSEPairHinge::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // static empty result
      static util::ShPtr< assemble::ProteinModel> s_empty_model;

      // collect possible sse pairs
      const storage::List< storage::VectorND< 2, util::SiPtr< const assemble::SSE> > > sse_pairs( m_Collector->Collect( PROTEIN_MODEL));

      // if no pair is available
      if( sse_pairs.IsEmpty())
      {
        // warn user
        BCL_MessageVrb( "unable to find a pair of sses");

        // return empty model
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // pick random sse pair
      const storage::VectorND< 2, util::SiPtr< const assemble::SSE> > &sse_pair
      (
        *random::GetGlobalRandom().Iterator( sse_pairs.Begin(), sse_pairs.End(), sse_pairs.GetSize())
      );

      // get a random number to decide which of the pair SSEs to assign as hinge
      const size_t index_for_hinge( random::GetGlobalRandom().Boolean());

      BCL_MessageVrb
      (
        "Hinge sse " + sse_pair( index_for_hinge)->GetIdentification() +
        "\nmutated_sse " + sse_pair( !index_for_hinge)->GetIdentification() +
        "with m_MoveHinge " + util::Format()( m_MoveHinge)
      );

      // create a new domain
      assemble::Domain this_domain;

      // initialize hinge to domain
      util::SiPtr< const coord::MovableInterface> hinge( this_domain);

      // insert the non-hinge sse into domain
      this_domain.Insert( util::ShPtr< assemble::SSE>( sse_pair( !index_for_hinge)->Clone()));

      // if m_MoveHinge boolean is set
      if( m_MoveHinge)
      {
        // insert the hinge into data
        this_domain.Insert( util::ShPtr< assemble::SSE>( sse_pair( index_for_hinge)->Clone()));
      }

      // constructor the MovableEccentric from domain and the hinge
      coord::MovableEccentric movable_eccentric( this_domain, *sse_pair( index_for_hinge));

      // move the MovableEccentric
      m_Move->Move( movable_eccentric);

      // copy the protein model
      util::ShPtr< assemble::ProteinModel> new_model( PROTEIN_MODEL.Clone());

      // replace the sses in the protein model contained in this domain with the moved copies
      new_model->Replace( this_domain.GetData());

      // end
      return math::MutateResult< assemble::ProteinModel>( new_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateProteinModelSSEPairHinge::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Collector, ISTREAM);
      io::Serialize::Read( m_Move, ISTREAM);
      io::Serialize::Read( m_MoveHinge, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &MutateProteinModelSSEPairHinge::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Collector, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Move, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MoveHinge, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
