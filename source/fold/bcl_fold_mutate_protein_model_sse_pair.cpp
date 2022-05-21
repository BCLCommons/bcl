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
#include "fold/bcl_fold_mutate_protein_model_sse_pair.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_geometry_packing.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelSSEPair::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateProteinModelSSEPair())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelSSEPair::MutateProteinModelSSEPair() :
      m_Collector(),
      m_TranslationRange(),
      m_RotationRange(),
      m_Scheme( GetStaticClassName< MutateProteinModelSSEPair>())
    {
    }

    //! @brief construct from Collector, max translation and rotation
    //! @param SP_COLLECTOR ShPtr to collector of SSE pairs
    //! @param MAX_TRANSLATION maximum translation allowed
    //! @param MAX_ROTATION maximum rotation allowed
    //! @param SCHEME Scheme to be used
    MutateProteinModelSSEPair::MutateProteinModelSSEPair
    (
      const util::ShPtr
      <
        find::CollectorInterface
        <
          storage::List< storage::VectorND< 2, util::SiPtr< const assemble::SSE> > >, assemble::DomainInterface
        >
      > &SP_COLLECTOR,
      const double MAX_TRANSLATION,
      const double MAX_ROTATION,
      const std::string &SCHEME
    ) :
      m_Collector( SP_COLLECTOR),
      m_TranslationRange( 0.0, MAX_TRANSLATION),
      m_RotationRange( 0.0, MAX_ROTATION),
      m_Scheme( SCHEME)
    {
    }

    //! @brief construct from Collector, min and max rotation and translation
    //! @param SP_COLLECTOR ShPtr to collector of SSE pairs
    //! @param TRANSLATION_RANGE range that specifies min and max translations
    //! @param ROTATION_RANGE range that specifies min and max rotations
    //! @param SCHEME Scheme to be used
    MutateProteinModelSSEPair::MutateProteinModelSSEPair
    (
      const util::ShPtr
      <
        find::CollectorInterface
        <
          storage::List< storage::VectorND< 2, util::SiPtr< const assemble::SSE> > >, assemble::DomainInterface
        >
      > &SP_COLLECTOR,
      const math::Range< double> &TRANSLATION_RANGE,
      const math::Range< double> &ROTATION_RANGE,
      const std::string &SCHEME
    ) :
      m_Collector( SP_COLLECTOR),
      m_TranslationRange( TRANSLATION_RANGE),
      m_RotationRange( ROTATION_RANGE),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateProteinModelSSEPair
    MutateProteinModelSSEPair *MutateProteinModelSSEPair::Clone() const
    {
      return new MutateProteinModelSSEPair( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateProteinModelSSEPair::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
    //! @param PROTEIN_MODEL protein model interest
    //! @return MutateResult with ProteinModel after the mutate
    math::MutateResult< assemble::ProteinModel> MutateProteinModelSSEPair::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // initialize empty model
      static util::ShPtr< assemble::ProteinModel> s_empty_model;

      // collect possible sse pairs
      const storage::List< storage::VectorND< 2, util::SiPtr< const assemble::SSE> > > sse_pairs
      (
        m_Collector->Collect( PROTEIN_MODEL)
      );

      // if no pair is available
      if( sse_pairs.IsEmpty())
      {
        // warn user
        BCL_MessageVrb( "unable to find a pair of sses");

        // return empty result
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // pick random sse pair
      const storage::VectorND< 2, util::SiPtr< const assemble::SSE> > &sse_pair
      (
        *random::GetGlobalRandom().Iterator( sse_pairs.Begin(), sse_pairs.End(), sse_pairs.GetSize())
      );

      // calculate the sse_packing
      const assemble::SSEGeometryPacking sse_pack( *sse_pair.First(), *sse_pair.Second());

      // pick either sse
      util::ShPtr< assemble::SSE> mutated_sse
      (
        sse_pair( !random::GetGlobalRandom().Boolean())->Clone()
      );

      BCL_MessageVrb( "mutated sse " + mutated_sse->GetIdentification());

      // generate TransformationMatrix3D and place sse in origin
      math::TransformationMatrix3D transform( -mutated_sse->GetCenter());

      // apply random translation around axis defined by the shortest connection between the two sses
      transform
      (
        math::RotationMatrix3D
        (
          sse_pack.GetShortestConnection().GetDirection(),
          random::GetGlobalRandom().Sign() *
          random::GetGlobalRandom().Double( m_RotationRange)
        )
      );

      // transform back to original position
      transform( mutated_sse->GetCenter());

      // translate along shortest connection
      transform
      (
        linal::Vector3D( linal::Vector3D( sse_pack.GetShortestConnection().GetDirection()).Normalize())
        * double( random::GetGlobalRandom().Sign())
        * random::GetGlobalRandom().Double( m_TranslationRange)
      );

      // apply transformation
      mutated_sse->Transform( transform);

      // make copy of proteinmodel
      util::ShPtr< assemble::ProteinModel> new_model( PROTEIN_MODEL.Clone());

      // replace mutated sse
      new_model->Replace( mutated_sse);

      // return
      return math::MutateResult< assemble::ProteinModel>( new_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateProteinModelSSEPair::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Collector, ISTREAM);
      io::Serialize::Read( m_TranslationRange, ISTREAM);
      io::Serialize::Read( m_RotationRange, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &MutateProteinModelSSEPair::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Collector, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_TranslationRange, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_RotationRange, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
