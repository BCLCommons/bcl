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
#include "fold/bcl_fold_placement_domain_using_fold_template.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_fold_template.h"
#include "assemble/bcl_assemble_fold_template_handler.h"
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
    const util::SiPtr< const util::ObjectInterface> PlacementDomainUsingFoldTemplate::s_Instance
    (
      GetObjectInstances().AddInstance( new PlacementDomainUsingFoldTemplate())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor with body deviations
    //! @param BODY_DEVIATION positive values get a larger template, and negative values exclude SSEs
    PlacementDomainUsingFoldTemplate::PlacementDomainUsingFoldTemplate
    (
      const int &BODY_DEVIATION
    ) :
      m_BodyDeviation( BODY_DEVIATION)
    {
    }

    //! @brief Clone function
    //! @return pointer to new PlacementDomainUsingFoldTemplate
    PlacementDomainUsingFoldTemplate *PlacementDomainUsingFoldTemplate::Clone() const
    {
      return new PlacementDomainUsingFoldTemplate( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &PlacementDomainUsingFoldTemplate::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief determines transformation matrices for placing SSEs from a Domain into a protein model
    //! @param DOMAIN_TO_PLACE domain containing SSEs to be placed
    //! @param PROTEIN_MODEL model that domain will be placed into
    //! @return transformation matrices for placing SSEs from a Domain into a protein model
    storage::Pair< storage::Map< util::SiPtr< const assemble::SSE>, math::TransformationMatrix3D>, bool>
    PlacementDomainUsingFoldTemplate::Place
    (
      const assemble::DomainInterface &DOMAIN_TO_PLACE, const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      return Place( DOMAIN_TO_PLACE);
    }

    //! @brief determines transformation matrices for placing SSEs from a Domain into a protein model
    //! @param DOMAIN_TO_PLACE domain containing SSEs to be placed
    //! @return transformation matrices for placing SSEs from a Domain into a protein model
    storage::Pair< storage::Map< util::SiPtr< const assemble::SSE>, math::TransformationMatrix3D>, bool>
    PlacementDomainUsingFoldTemplate::Place
    (
      const assemble::DomainInterface &DOMAIN_TO_PLACE
    ) const
    {
      // initialize map
      storage::Map< util::SiPtr< const assemble::SSE>, math::TransformationMatrix3D> matrices;

      // copy SSEs
      util::SiPtrVector< const assemble::SSE> helices( DOMAIN_TO_PLACE.GetSSEs( biol::GetSSTypes().HELIX));
      util::SiPtrVector< const assemble::SSE> strands( DOMAIN_TO_PLACE.GetSSEs( biol::GetSSTypes().STRAND));

      // return empty vector if domain is empty
      if( helices.GetSize() == 0 && strands.GetSize() == 0)
      {
        BCL_MessageStd( "No SSEs in domain to place into a fold template");
        return storage::Pair< storage::Map< util::SiPtr< const assemble::SSE>, math::TransformationMatrix3D>, bool>
        (
          matrices, false
        );
      }

      // if sses should be excluded
      if( m_BodyDeviation < 0)
      {
        const size_t sses_to_exclude( -m_BodyDeviation);

        // warn the user and break if the number of sses to exclude is greater than or equal to the number of total sses
        if( sses_to_exclude >= helices.GetSize() + strands.GetSize())
        {
          BCL_MessageStd
          (
            "The number of SSEs to exclude from the fold template placement"
            " is larger than or equal to the total number of SSEs."
          );
          return storage::Pair< storage::Map< util::SiPtr< const assemble::SSE>, math::TransformationMatrix3D>, bool>
          (
            matrices, false
          );
        }

        // iterate through the SSEs to exclude
        for( size_t i( 0); i != sses_to_exclude; ++i)
        {
          // decide whether to remove helix or strand
          bool remove_helix( random::GetGlobalRandom().Boolean());

          // true if remove helix was chosen
          if( remove_helix)
          {
            // if there are no helical elements to remove, remove a strand instead
            if( helices.GetSize() == 0)
            {
              remove_helix = false;
            }
          }
          // otherwise a strand was chosen
          else
          {
            // if there are no strand elements to remove, remove a helix instead
            if( strands.GetSize() == 0)
            {
              remove_helix = true;
            }
          }

          // true if a helix is to be removed
          if( remove_helix)
          {
            // remove a random helical element
            helices.RemoveRandomElement();
            }
          // otherwise remove a strand
          else
          {
            // remove a random strand element
            strands.RemoveRandomElement();
          }
        }
      }

      // sort the SSEs by size
      std::sort( helices.Begin(), helices.End(), assemble::SSELessThanBySize());
      std::sort( strands.Begin(), strands.End(), assemble::SSELessThanBySize());

      // determine the number of helices and strands in the fold template
      storage::Pair< size_t, size_t> number_of_sses( helices.GetSize(), strands.GetSize());
      size_t added_helices( 0);
      size_t added_strands( 0);

      // if a larger template is requested
      if( m_BodyDeviation > 0)
      {
        // iterate through the sses to add
        for( int i( 0); i != m_BodyDeviation; ++i)
        {
          // decide whether to add a helix or strand
          const bool add_helix( random::GetGlobalRandom().Boolean());

          // increment proper sse
          add_helix ? ++added_helices : ++added_strands;
        }
      }

      // update the sse count
      number_of_sses.First() += added_helices;
      number_of_sses.Second() += added_strands;

      // store random fold template from the probability distribution
      assemble::FoldTemplate fold_template
      (
        assemble::FoldTemplateHandler::GetRandomTemplate( number_of_sses.First(), number_of_sses.Second())
      );
      BCL_MessageDbg( "Fold template generated from PDB: " + fold_template.GetPDBID());

      // calculate the topology
      fold_template.CalculateTopology();

      // if a larger template was requested
      if( m_BodyDeviation > 0)
      {
        // get a subdomain from the template
        fold_template =
          fold_template.GetSubDomain( number_of_sses.First() - added_helices, number_of_sses.Second() - added_strands);
      }

      // get the helical and strand geometries
      util::SiPtrVector< const assemble::SSEGeometryPhiPsi> helical_geometries( fold_template.GetHelicalGeometries());
      util::SiPtrVector< const assemble::SSEGeometryPhiPsi> strand_geometries( fold_template.GetStrandGeometries());
      // create iterators on the geometries in the template
      util::SiPtrVector< const assemble::SSEGeometryPhiPsi>::const_iterator helix_geometry_itr( helical_geometries.Begin());
      const util::SiPtrVector< const assemble::SSEGeometryPhiPsi>::const_iterator helix_geometry_itr_end( helical_geometries.End());
      util::SiPtrVector< const assemble::SSEGeometryPhiPsi>::const_iterator strand_geometry_itr( strand_geometries.Begin());
      const util::SiPtrVector< const assemble::SSEGeometryPhiPsi>::const_iterator strand_geometry_itr_end( strand_geometries.End());

      // iterate through the helices
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
          helix_itr( helices.Begin()), helix_itr_end( helices.End());
        helix_itr != helix_itr_end && helix_geometry_itr != helix_geometry_itr_end;
        ++helix_itr, ++helix_geometry_itr
      )
      {
        // create transformation matrix and move to new position
        math::TransformationMatrix3D transform( math::Inverse( ( *helix_itr)->GetOrientation()));
        transform( ( *helix_geometry_itr)->GetOrientation());

        //insert into map
        matrices[ *helix_itr] = transform;
      }

      // iterate through the strands
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
          strand_itr( strands.Begin()), strand_itr_end( strands.End());
        strand_itr != strand_itr_end && strand_geometry_itr != strand_geometry_itr_end;
        ++strand_itr, ++strand_geometry_itr
      )
      {
        // create transformation matrix and move to new position
        math::TransformationMatrix3D transform( math::Inverse( ( *strand_itr)->GetOrientation()));
        transform( ( *strand_geometry_itr)->GetOrientation());

        //insert into map
        matrices[ *strand_itr] = transform;
      }

      // return the transformation matrices
      return storage::Pair< storage::Map< util::SiPtr< const assemble::SSE>, math::TransformationMatrix3D>, bool>
      (
        matrices, true
      );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PlacementDomainUsingFoldTemplate::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_BodyDeviation, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PlacementDomainUsingFoldTemplate::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_BodyDeviation, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
