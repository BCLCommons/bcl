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
#include "fold/bcl_fold_mutate_protein_model_domain_add.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_fold_template.h"
#include "assemble/bcl_assemble_fold_template_handler.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    //! initialize deviation range
    const math::Range< double> MutateProteinModelDomainAdd::s_Deviation( 0.0, 5.0);

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelDomainAdd::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateProteinModelDomainAdd())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelDomainAdd::MutateProteinModelDomainAdd() :
      m_CollectorDomain(),
      m_FitToFoldTemplate( false),
      m_Scheme( GetStaticClassName< MutateProteinModelDomainAdd>())
    {
    }

    //! @brief construct from a domain collector and a scheme
    //! @param DOMAIN_COLLECTOR domain collector to be used to get a domain
    //! @param FIT_TO_FOLD_TEMPLATE bool whether to fit the domain to a fold template prior to adding
    //! @param SCHEME Scheme to be used
    MutateProteinModelDomainAdd::MutateProteinModelDomainAdd
    (
      const find::CollectorInterface< util::ShPtrVector< assemble::Domain>, assemble::DomainInterface> &DOMAIN_COLLECTOR,
      const bool FIT_TO_FOLD_TEMPLATE,
      const std::string &SCHEME
    ) :
      m_CollectorDomain( DOMAIN_COLLECTOR.Clone()),
      m_FitToFoldTemplate( FIT_TO_FOLD_TEMPLATE),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateProteinModelDomainAdd
    MutateProteinModelDomainAdd *MutateProteinModelDomainAdd::Clone() const
    {
      return new MutateProteinModelDomainAdd( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateProteinModelDomainAdd::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that takes a ProteinModel and return a mutated ProteinModel
    //! @param PROTEIN_MODEL ProteinModel which will be mutated
    //! @return MutateResult with the mutated ProteinModel
    math::MutateResult< assemble::ProteinModel> MutateProteinModelDomainAdd::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // static empty model
      static util::ShPtr< assemble::ProteinModel> s_empty_model;

      // get possible domains from the collector
      const util::ShPtrVector< assemble::Domain> domains( m_CollectorDomain->Collect( PROTEIN_MODEL));

      // if no domains were collected
      if( domains.IsEmpty())
      {
        // warn the user
        BCL_MessageStd( "No domains collected, skipping " + m_Scheme);

        // return empty result
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // get a random domain
      util::ShPtr< assemble::Domain> sp_domain
      (
        *random::GetGlobalRandom().Iterator( domains.Begin(), domains.End(), domains.GetSize())
      );

      // create a new protein model
      util::ShPtr< assemble::ProteinModel> new_model( PROTEIN_MODEL.Clone());

      // get the SSEs in the domain
      const util::SiPtrVector< const assemble::SSE> domain_sses( sp_domain->GetSSEs());
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator sse_itr( domain_sses.Begin()),
          sse_itr_end( domain_sses.End());
        sse_itr != sse_itr_end; ++sse_itr
      )
      {
        // remove the SSE from the protein model
        new_model->Remove( **sse_itr);
      }

      // if the domain should be fit into a fold template
      if( m_FitToFoldTemplate)
      {
        // get a random template
        const assemble::FoldTemplate fold_template
        (
          assemble::FoldTemplateHandler::GetRandomSubTemplate( sp_domain->GetSSEs())
        );

        // if a random template could not be found
        if( fold_template.GetGeometries().IsEmpty())
        {
          // warn the user
          BCL_MessageStd( "Unable to find suitable fold template, skipping " + m_Scheme);

          // return empty result
          return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
        }

        // get the domain from fitting the SSEs into the template
        sp_domain = util::ShPtr< assemble::Domain>( fold_template.FitSSEs( sp_domain->GetSSEs()).Clone());
      }

      // initialize transformation matrix with move from domain orientation to origin
      math::TransformationMatrix3D transform( math::Inverse( sp_domain->GetOrientation()));

      // create a range to choose the random rotation from
      const math::Range< double> angle_range( 0.0, 2.0 * math::g_Pi);

      // rotate the transformation matrix with a random rotation
      transform
      (
        math::RotationMatrix3D
        (
          random::GetGlobalRandom().Double( angle_range),
          random::GetGlobalRandom().Double( angle_range),
          random::GetGlobalRandom().Double( angle_range)
        )
      );

      // set the translation vector to be along the Z-axis of the rotation matrix
      const linal::Vector3D translation_vector( transform.GetRotation().GetAxis( coord::GetAxes().e_Z));

      // translate to the center of the protein model
      const linal::Vector3D protein_center( new_model->GetCenter());
      transform( protein_center);

      // move the domain to the center of the protein model
      sp_domain->Transform( transform);

      // create a line segment from the center pointing outwards
      const coord::LineSegment3D line_segment_forward
      (
        protein_center,
        protein_center + 1000.0 * translation_vector
      );

      // collect the coordinates in the protein model
      const util::SiPtrVector< const linal::Vector3D> coords_a
      (
        new_model->GetAtomCoordinates( storage::Set< biol::AtomType>( biol::GetAtomTypes().CA))
      );

      // get the distance to the edge of the protein model
      const double protein_edge_distance( GetMaxDistanceFromCenter( line_segment_forward, coords_a));

      // collect the coordinates in the domain
      const util::SiPtrVector< const linal::Vector3D> coords_b
      (
        sp_domain->GetAtomCoordinates( storage::Set< biol::AtomType>( biol::GetAtomTypes().CA))
      );

      // create a line segment from the center pointing outwards
      const coord::LineSegment3D line_segment_reverse
      (
        protein_center,
        protein_center - 1000.0 * translation_vector
      );

      // get the distance to the edge of the domain
      const double domain_edge_distance( GetMaxDistanceFromCenter( line_segment_reverse, coords_b));

      // get a random deviation for the translation
      const double random_deviation( random::GetGlobalRandom().Double( s_Deviation));

      // move the domain to the determined position
      sp_domain->Translate
      (
        ( protein_edge_distance + domain_edge_distance + random_deviation) * translation_vector
      );

      // iterate over the SSEs
      const util::SiPtrVector< const assemble::SSE> sses( sp_domain->GetSSEs());
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator sse_itr( sses.Begin()), sse_itr_end( sses.End());
        sse_itr != sse_itr_end; ++sse_itr
      )
      {
        // insert the sse into the protein model
        new_model->Insert( util::ShPtr< assemble::SSE>( ( *sse_itr)->Clone()));
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
    std::istream &MutateProteinModelDomainAdd::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_CollectorDomain, ISTREAM);
      io::Serialize::Read( m_FitToFoldTemplate, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateProteinModelDomainAdd::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_CollectorDomain, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_FitToFoldTemplate, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief takes a set of coordinates and returns the farthest distance from the beginning of the line segment
    //! @param LINE_SEGMENT line segment that contains the point and direction to be used
    //! @param COORDS coordinates to be checked
    //! @return the farthest distance from the beginning of the line segment
    double MutateProteinModelDomainAdd::GetMaxDistanceFromCenter
    (
      const coord::LineSegment3D &LINE_SEGMENT,
      const util::SiPtrVector< const linal::Vector3D> &COORDS
    )
    {
      // initialize the maximum distance from the center of the model
      double max_distance( 0.0);

      // iterate over the coordinates
      for
      (
        util::SiPtrVector< const linal::Vector3D>::const_iterator coord_itr( COORDS.Begin()),
          coord_itr_end( COORDS.End());
        coord_itr != coord_itr_end; ++coord_itr
      )
      {
        // if the coordinate is within 5 angstroms from the line segment
        if( coord::CalculateDistancePointFromLineSegment( LINE_SEGMENT, **coord_itr).First() < 5.0)
        {
          // update the max distance if this distance from the center is greater
          max_distance = std::max( max_distance, linal::Distance( LINE_SEGMENT.GetStartPoint(), **coord_itr));
        }
      }

      // end
      return max_distance;
    }

  } // namespace fold
} // namespace bcl
