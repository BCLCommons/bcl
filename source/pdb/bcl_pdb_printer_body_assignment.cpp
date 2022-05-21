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
#include "pdb/bcl_pdb_printer_body_assignment.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "pdb/bcl_pdb_line.h"
#include "restraint/bcl_restraint_assignment.h"
#include "restraint/bcl_restraint_body.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace pdb
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> PrinterBodyAssignment::s_Instance
    (
      GetObjectInstances().AddInstance( new PrinterBodyAssignment())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PrinterBodyAssignment::PrinterBodyAssignment() :
      m_Restraints()
    {
    }

    //! @brief construct from restraint information
    //! @param RESTRAINTS holds the restraints (density map)
    PrinterBodyAssignment::PrinterBodyAssignment
    (
      const util::ShPtr
      <
        util::ShPtrVector< restraint::Body>
      > &RESTRAINTS
    ) :
      m_Restraints( RESTRAINTS)
    {
    }

    //! @brief Clone function
    //! @return pointer to new PrinterBodyAssignment
    PrinterBodyAssignment *PrinterBodyAssignment::Clone() const
    {
      return new PrinterBodyAssignment( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PrinterBodyAssignment::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes a protein model and returns PDB lines
    //! @param PROTEIN_MODEL protein model to print
    //! @return PDB lines
    util::ShPtrList< Line> PrinterBodyAssignment::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // initialize lines
      util::ShPtrList< Line> lines;

      // if no restraints
      if( !m_Restraints.IsDefined())
      {
        return lines;
      }

      BCL_MessageCrt
      (
        "writing PrinterProteinModelBodyAssignment::WriteBodySSEAssignmentInformation!!! "
      );

      // get the pool from the given ProteinModel and make sure it is valid
      const util::ShPtr< assemble::SSEPool> sp_pool( PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Pool));
      BCL_Assert( sp_pool.IsDefined(), "No pool stored for the given model");

      // create storage::Vector to Assignments "density_sse_assignment_list" to hold all of the density map assignments
      // This should only contain a single restraint::Assignment in the case where you only use a single density
      // map (don't forget the restraint::Assignment contains all of the information about a density map, i.e. not just a
      //  single density rod)
      storage::Vector
      <
        restraint::SSEAssignment
      > density_sse_assignment_list( GenerateAssignments( PROTEIN_MODEL));

      // create SiPtrVector "protein_model_sses" and initialize with the SSEs of "PROTEIN_MODEL"
      const util::SiPtrVector< const assemble::SSE> protein_model_sses( PROTEIN_MODEL.GetSSEs());

      // create an empty remark line
      util::ShPtr< Line> empty_remark_line( new Line( GetLineTypes().REMARK));
      empty_remark_line->Put( GetEntryTypes().REMARK_Number, util::Format()( size_t( s_RemarkNumber)));
      lines.PushBack( empty_remark_line);

      // iterate through all of the SSEs of the protein model in order to see which SSE corresponds with which body in
      // the assignment
      // Need to do this because the information about the SSE from the protein model is lost when it is put into the
      // assignment because the assignment just keeps the body information.
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
          sse_itr( protein_model_sses.Begin()), sse_itr_end( protein_model_sses.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // create linal::Vector3D "sse_center" and initialize with the origin of the SSE body
        const linal::Vector3D sse_center( ( *sse_itr)->GetCenter());

        // iterate through the vector of assignments (for one density map this for loop will only execute once)
        for
        (
          storage::Vector< restraint::SSEAssignment>::const_iterator
            assignment_itr( density_sse_assignment_list.Begin()), assignment_itr_end( density_sse_assignment_list.End());
          assignment_itr != assignment_itr_end;
          ++assignment_itr
        )
        {
          // iterate through the individual bodies of the assignment denoted by "assignment_itr"
          // iterate through the bodies of the sses of the protein model (stored in Group Collection of the Assignment)
          for
          (
            restraint::GroupCollection< size_t, assemble::SSE>::const_iterator
              assignment_sses_iter( assignment_itr->GetGroupCollection().Begin()),
              assignment_sses_iter_end( assignment_itr->GetGroupCollection().End());
            assignment_sses_iter != assignment_sses_iter_end;
            ++assignment_sses_iter
          )
          {
            // true if the center of the body behind "assignment_sses_iter"
            // (this is NOT the body of the density rod but the center of a body coming from an SSE in the protein model)
            // is equal to "sse_center"
            if( math::EqualWithinTolerance( ( *assignment_sses_iter->second.Begin())->GetCenter(), sse_center))
            {
              // get the index of the corresponding restraint body (i.e. density rod)
              const size_t density_rod_index( assignment_sses_iter->first);

              // compare the orientation of the density rod and the helix in the protein model
              // initialize bool that holds relative orientation of density rod and the helix in the protein model
              bool relative_orientation_rod_helix( true);

              // get vector along helical axis of sse from protein model
              const linal::Vector3D sse_axis( ( *sse_itr)->GetAxis( coord::GetAxes().e_Z));

              // get vector along helical axis of density rod
              const linal::Vector3D density_rod_axis
              (
                assignment_itr->GetRestraint()->operator()( density_rod_index)->GetAxis( coord::GetAxes().e_Z)
              );

              // check whether both vectors are parallel or antiparallel and store result in relative_orientation_rod_helix
              const double scalar_product( linal::ScalarProduct( sse_axis, density_rod_axis));
              BCL_MessageDbg( "the scalar product is: " + util::Format()( scalar_product));
              // check sign of scalar product (positive if parallel, negative if antiparallel)
              // parallel
              if( math::EqualWithinTolerance( scalar_product, 1.0))
              {
                relative_orientation_rod_helix = true;
              }
              // antiparallel
              else if( math::EqualWithinTolerance( scalar_product, -1.0))
              {
                relative_orientation_rod_helix = false;
              }
              // they should always be parallel or antiparallel
              else
              {
                BCL_MessageCrt( "magnitude of scalar product should always be 1, but it is: " + util::Format()( scalar_product));
                relative_orientation_rod_helix = false;
              }

              // find best match from pool for this particular sse
              // (this is necessary as sses can be changed by moves from their original pool state)
              const util::SiPtr< const assemble::SSE> best_match_from_pool( sp_pool->FindBestMatchFromPool( **sse_itr).First());

              std::string assignment_string
              (
                "map: " + util::Format()( assignment_itr - density_sse_assignment_list.Begin()) +
                " rod: " + util::Format()( density_rod_index)
              );

              if( best_match_from_pool.IsDefined())
              {
                // still iterate over pool to get the pool index of best_match_from_pool sse
                // create size_t "pool_sse_index_counter" and initialize with zero
                // this variable will be used to give an index of an SSE in the pool
                size_t pool_sse_index_counter( 0);

                // now iterate through the pool to see which element in the pool the best_match_from_pool sse corresponds
                // to; the element in the pool will be denoted by its index
                for
                (
                  storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThan>::const_iterator
                    pool_itr( sp_pool->Begin()), pool_itr_end( sp_pool->End());
                  pool_itr != pool_itr_end;
                  ++pool_itr, ++pool_sse_index_counter
                )
                {
                  if( *pool_itr == best_match_from_pool)
                  {
                    break;
                  }
                }

                // generate the assignment string
                assignment_string +=
                  " index: " + util::Format()( pool_sse_index_counter) +
                  " begin: " + util::Format()( best_match_from_pool->GetFirstAA()->GetSeqID()) +
                  " end: " + util::Format()( best_match_from_pool->GetLastAA()->GetSeqID()) +
                  " orientation: " + util::Format()( relative_orientation_rod_helix);
              }
              else
              {
                assignment_string += " no assignment found";
              }

              // add to pdb lines
              util::ShPtr< Line> new_line( empty_remark_line->Clone());
              new_line->Put( GetEntryTypes().REMARK_String, assignment_string);
              lines.PushBack( new_line);
            }
          }
        }
      }

      // only print out body sse agreement information if message level is sufficiently low
      if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Standard))
      {
        lines.Append( WriteBodySSEAgreementInformation( PROTEIN_MODEL));
      }

      // end
      return lines;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PrinterBodyAssignment::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PrinterBodyAssignment::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief write relative rotation and translation information between every SSE in the protein model and the body
    //! it is assigned to
    //! @param PROTEIN_MODEL protein model for which the translation and rotation agreement information is determined
    //! @return pdb lines which are written to
    util::ShPtrList< Line> PrinterBodyAssignment::WriteBodySSEAgreementInformation
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // initialize lines
      util::ShPtrList< Line> lines;

      // if no restraints
      if( !m_Restraints.IsDefined())
      {
        return lines;
      }

      // create an empty remark line
      util::ShPtr< Line> empty_remark_line( new Line( GetLineTypes().REMARK));
      empty_remark_line->Put( GetEntryTypes().REMARK_Number, util::Format()( size_t( s_RemarkNumber)));

      // create storage::Vector to Assignments "density_sse_assignment_list" to hold all of the density map assignments
      // This should only contain a single restraint::Assignment in the case where you only use a single density
      // map (don't forget the restraint::Assignment contains all of the information about a density map, i.e. not just a
      //  single density rod)
      storage::Vector
      <
        restraint::SSEAssignment
      > density_sse_assignment_list( GenerateAssignments( PROTEIN_MODEL));

      // iterate over all the assignments that are contained in the density_sse_assignment_list (in case there is more
      // than 1 density map)
      for
      (
        storage::Vector< restraint::SSEAssignment>::const_iterator
          assignment_itr( density_sse_assignment_list.Begin()), assignment_itr_end( density_sse_assignment_list.End());
        assignment_itr != assignment_itr_end; ++assignment_itr
      )
      {
        // iterate through the assignment denoted by assignment iterator
        for
        (
          restraint::GroupCollection< size_t, assemble::SSE>::const_iterator
            group_collection_itr( assignment_itr->GetGroupCollection().Begin()),
            group_collection_itr_end( assignment_itr->GetGroupCollection().End());
          group_collection_itr != group_collection_itr_end; ++group_collection_itr
        )
        {
          // get body of SSE in protein model
          const assemble::SSE &model_body( *group_collection_itr->second.FirstElement());

          // get body of restraint (density rod)
          const assemble::SSEGeometryInterface &native_body
          (
            *assignment_itr->GetRestraint()->operator()( group_collection_itr->first)
          );

          const math::TransformationMatrix3D difference
          (
            math::Inverse( model_body.GetOrientation())( native_body.GetOrientation())
          );

          // calculate the effective rotation angle between these two bodies
          const double rotation_angle( difference.GetRotation().EffectiveRotationAngle());
          const double translation
          (
            ( model_body.GetOrientation().GetTranslation() - native_body.GetOrientation().GetTranslation()).Norm()
          );

          // generate the assignment string
          const std::string assignment_string
          (
            "map: " + util::Format()( assignment_itr - density_sse_assignment_list.Begin()) +
            " rod: " + util::Format()( group_collection_itr->first) +
            " angle: " + util::Format()( math::Angle::Degree( rotation_angle)) +
            " translation: " + util::Format()( translation)
          );

          // add to pdb lines
          util::ShPtr< Line> new_line( empty_remark_line->Clone());
          new_line->Put( GetEntryTypes().REMARK_String, assignment_string);
          lines.PushBack( new_line);
        }
      }

      // end
      return lines;
    }

    //! @brief generate assignments from a protein model
    //! @param PROTEIN_MODEL protein model of interest
    //! @return assignments
    storage::Vector< restraint::SSEAssignment>
    PrinterBodyAssignment::GenerateAssignments( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // create storage::Vector to Assignments "density_sse_assignment_list" to hold all of the density map assignments
      // This should only contain a single restraint::Assignment in the case where you only use a single density
      // map (don't forget the restraint::Assignment contains all of the information about a density map, i.e. not just
      // a single density rod)
      storage::Vector
      <
        restraint::SSEAssignment
      > density_sse_assignment_list;

      // create SiPtrVector "protein_model_sses" and initialize with the SSEs of "PROTEIN_MODEL"
      const util::SiPtrVector< const assemble::SSE> protein_model_sses( PROTEIN_MODEL.GetSSEs());

      // iterate through "m_restraints" to build assignment objects
      // An assignment contains the information about the density map and which SSEs from the protein model go with which
      // density rods.
      // In the usual case, this loop will only execute a single time (i.e. only one density map is given)
      for
      (
        util::ShPtrVector< restraint::Body>::const_iterator itr( m_Restraints->Begin()), itr_end( m_Restraints->End());
        itr != itr_end;
        ++itr
      )
      {
        density_sse_assignment_list.PushBack( ( *itr)->GenerateAssignment( protein_model_sses));
      }

      return density_sse_assignment_list;
    }

  } // namespace pdb
} // namespace bcl
