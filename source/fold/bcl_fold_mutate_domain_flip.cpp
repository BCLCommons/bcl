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
#include "fold/bcl_fold_mutate_domain_flip.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    // instantiate instance
    const util::SiPtr< const util::ObjectInterface> MutateDomainFlip::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateDomainFlip())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a single flip axis, flip internal and and a flip all boolean
    //! @param FLIP_AXIS Axis to apply the flip around
    //! @param FLIP_INTERNAL boolean to flip internally
    //! @param FLIP_ALL boolean to decide whether to flip all strands or only few
    //! @param USE_DIFFERENT_FLIP_AXES boolean to decide whether to use different flip axes for each flip
    MutateDomainFlip::MutateDomainFlip
    (
      const coord::Axis FLIP_AXIS,
      const bool FLIP_INTERNAL,
      const bool FLIP_ALL,
      const bool USE_DIFFERENT_FLIP_AXES
    ) :
      m_FlipAxes( FLIP_AXIS),
      m_FlipInternal( FLIP_INTERNAL),
      m_FlipAll( FLIP_ALL),
      m_UseDifferentFlipAxes( USE_DIFFERENT_FLIP_AXES)
    {
      // check the flip booleans
      BCL_Assert
      (
        m_FlipInternal || ( m_FlipAll && !m_FlipInternal), "Flip internal cannot be false if flip all was false!"
      );
    }

    //! @brief constructor from a set of flip axis, flip internal and and a flip all boolean
    //! @param FLIP_AXES Set of axes to apply the flip around, (randomly chosen at flip time)
    //! @param FLIP_INTERNAL boolean to flip internally
    //! @param FLIP_ALL boolean to decide whether to flip all strands or only few
    //! @param USE_DIFFERENT_FLIP_AXES boolean to decide whether to use different flip axes for each flip
    MutateDomainFlip::MutateDomainFlip
    (
      const storage::Set< coord::Axis> &FLIP_AXES,
      const bool FLIP_INTERNAL,
      const bool FLIP_ALL,
      const bool USE_DIFFERENT_FLIP_AXES
    ) :
      m_FlipAxes( FLIP_AXES),
      m_FlipInternal( FLIP_INTERNAL),
      m_FlipAll( FLIP_ALL),
      m_UseDifferentFlipAxes( USE_DIFFERENT_FLIP_AXES)
    {
      // make sure the given set of axes is not empty
      BCL_Assert( !m_FlipAxes.IsEmpty(), "The given set of flip axes is empty!");

      // check the flip booleans
      BCL_Assert
      (
        m_FlipInternal || ( m_FlipAll && !m_FlipInternal), "Flip internal cannot be false if flip all was false!"
      );
    }

    //! @brief Clone function
    //! @return pointer to new MutateDomainFlip
    MutateDomainFlip *MutateDomainFlip::Clone() const
    {
      return new MutateDomainFlip( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateDomainFlip::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return flip axes
    //! @return flip axes
    const storage::Set< coord::Axis> &MutateDomainFlip::GetFlipAxes() const
    {
      return m_FlipAxes;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that takes a Domain and return a mutated Domain
    //! @param THIS_DOMAIN Domain which will be mutated
    //! @return MutateResult with the mutated Domain
    math::MutateResult< assemble::Domain> MutateDomainFlip::operator()( const assemble::Domain &THIS_DOMAIN) const
    {
      // initialize an empty Domain ShPtr
      static util::ShPtr< assemble::Domain> s_empty_domain;

      // make sure the Domain has at least two SSEs
      if( THIS_DOMAIN.GetNumberSSEs() < 2)
      {
        // warn user and return
        BCL_MessageVrb( "The given domain has less than 2 SSEs, therefore skipping");
        return math::MutateResult< assemble::Domain>( s_empty_domain, *this);
      }

      // make a copy of the Domain
      util::ShPtr< assemble::Domain> new_domain( THIS_DOMAIN.Clone());

      // collect the list of SSEs to flip
      util::SiPtrVector< const assemble::SSE> sses_to_flip( new_domain->GetSSEs());

      // initialize transformation to be applied
      math::TransformationMatrix3D transformation;

      // determine the random flip axes
      coord::Axis flip_axis( *random::GetGlobalRandom().Iterator( m_FlipAxes.Begin(), m_FlipAxes.End(), m_FlipAxes.GetSize()));

      // if flipping externally
      if( !m_FlipInternal)
      {
        BCL_MessageDbg( "Flipping externally");

        // make a copy of the domain orientation
        math::TransformationMatrix3D domain_orientation( new_domain->GetOrientation());

        // build up the transformation matrix
        // first move the origin
        transformation = math::Inverse( domain_orientation);

        // apply the flip along the random axis
        transformation( flip_axis, math::g_Pi);

        // transform back
        transformation( domain_orientation);
      }
      // else if flipping internally
      else
      {
        // collect the list of SSEs to flip
        const size_t nr_sses( sses_to_flip.GetSize());
        BCL_MessageDbg( "Flipping internally ");

        // if flip all is not selected
        if( !m_FlipAll)
        {
          // determine how many SSEs to flip
          const size_t nr_flips( random::GetGlobalRandom().Random< size_t>( 1, nr_sses));
          const size_t nr_sses_unflipped( nr_sses - nr_flips);
          BCL_MessageDbg
          (
            "Flipping internally " + util::Format()( nr_flips) + " SSEs out of " + util::Format()( nr_sses)
          );

          // iterate number of SSEs to left unflipped times
          for( size_t ctr( 0); ctr < nr_sses_unflipped; ++ctr)
          {
            // remove a random element
            sses_to_flip.RemoveRandomElement();
          }
        }
      }

      // now iterate over the SSEs to flip
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
          sse_itr( sses_to_flip.Begin()), sse_itr_end( sses_to_flip.End());
        sse_itr != sse_itr_end; ++sse_itr
      )
      {
        // get the next strand
        const assemble::SSE &this_sse( **sse_itr);

        // if internal flip
        if( m_FlipInternal)
        {
          // if the boolean to use different axis is given
          if( m_UseDifferentFlipAxes)
          {
            // then get a new random axis
            flip_axis = *random::GetGlobalRandom().Iterator( m_FlipAxes.Begin(), m_FlipAxes.End(), m_FlipAxes.GetSize());
          }
          // update the transformation so that it's specific to this SSE
          // first step is to move to the origin
          transformation = math::Inverse( this_sse.GetOrientation());
          // rotate apply the flip
          transformation( flip_axis, math::g_Pi);
          // move back to the original location
          transformation( this_sse.GetOrientation());
        }

        BCL_MessageDbg
        (
          "Flipping " + this_sse.GetIdentification() + " around axis " + flip_axis.GetName()
        );

        // make a copy of this SSE and apply the transformation
        util::ShPtr< assemble::SSE> new_sse( this_sse.Clone());
        new_sse->Transform( transformation);

        // replace the corresponding SSE with the new_sse in the new_domain
        new_domain->Replace( new_sse);
      }

      // end
      return math::MutateResult< assemble::Domain>( new_domain, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateDomainFlip::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_FlipAxes, ISTREAM);
      io::Serialize::Read( m_FlipInternal, ISTREAM);
      io::Serialize::Read( m_FlipAll, ISTREAM);
      io::Serialize::Read( m_UseDifferentFlipAxes, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateDomainFlip::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // read members
      io::Serialize::Write( m_FlipAxes, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_FlipInternal, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_FlipAll, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_UseDifferentFlipAxes, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
