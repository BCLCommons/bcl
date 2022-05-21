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
#include "restraint/bcl_restraint_handler_body.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "restraint/bcl_restraint_body.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> HandlerBody::s_Instance
    (
      GetObjectInstances().AddInstance( new HandlerBody())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    HandlerBody::HandlerBody() :
      m_DetermineOccupancy()
    {
    }

    //! @brief constructor from a ShPtr to a FunctionInterface defining how occupancy is determined
    //! @param OCCUPANCY ShPtr to a FunctionInterface defining how occupancy is determined
    HandlerBody::HandlerBody
    (
      const util::ShPtr< util::BinaryFunctionInterface< assemble::SSEGeometryInterface, assemble::SSE, bool> > &OCCUPANCY
    ) :
      m_DetermineOccupancy( OCCUPANCY)
    {
      BCL_Assert
      (
        m_DetermineOccupancy.IsDefined(),
        "HandlerBody::HandlerBody occupancy constructor : m_DetermineOccupancy is not defined"
      );
    }

    //! @brief virtual copy constructor
    HandlerBody *HandlerBody::Clone() const
    {
      return new HandlerBody( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &HandlerBody::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief CreateRestraintsBody is the function which creates the Restraints from an istream
    //! @param ISTREAM is the istream from which the restraints will be created
    //! @return returns a ShPtrVector of RestraintInterfaces
    util::ShPtrVector< Body>
    HandlerBody::CreateRestraintsBody( std::istream &ISTREAM) const
    {
      // create the density rods
      const util::ShPtr< util::ShPtrVector< assemble::SSEGeometryInterface> > restraint_bodies( GenerateDensityRods( ISTREAM));

      // put density rods into body restraint and return body restraint
      return util::ShPtrVector< Body>( 1, util::ShPtr< Body>( new Body( restraint_bodies, m_DetermineOccupancy)));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read restraint from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &HandlerBody::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_DetermineOccupancy, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write restraint to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &HandlerBody::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_DetermineOccupancy, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief CreateRestraintsBody is the function which creates the Restraints from an istream
    //! @param ISTREAM is the istream from which the restraints will be created
    //! @return returns a ShPtrVector of RestraintInterfaces
    util::ShPtr< util::ShPtrVector< assemble::SSEGeometryInterface> > HandlerBody::GenerateDensityRods( std::istream &ISTREAM) const
    {
      // read pdb file
      pdb::Handler pdb( ISTREAM);
      ISTREAM.clear();

      // create "factory" to create protein model with amino acids of type AABackBone
      pdb::Factory factory( biol::GetAAClasses().e_AABackBone);

      // create ProteinModel "protein_model" from "pdb"
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 9;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel protein_model( factory.ProteinModelFromPDB( pdb, ssetype_min_size));

      // get the secondary structure elements of "protein_model" which will be the body restraints
      util::SiPtrVector< const assemble::SSE> sses( protein_model.GetSSEs());

      BCL_MessageStd( "number of SSEs to create bodies from: " + util::Format()( sses.GetSize()));
      BCL_MessageStd( "number of chains: " + util::Format()( protein_model.GetNumberOfChains()));
      BCL_MessageStd( "number of sses: " + util::Format()( protein_model.GetNumberSSEs()));
      //BCL_MessageStd( "number of sses in chain a: " + util::Format()( protein_model.GetChain('A')->GetNumberSSEs()));
      //BCL_MessageStd( "number of sses in chain b: " + util::Format()( protein_model.GetChain('B')->GetNumberSSEs()));
      BCL_MessageStd( "number of aas in protein model: " + util::Format()( protein_model.GetNumberAAs()));

      // create a ShPtrVector of bodies which will be created from "bodies" and will be used to create the body
      // restraint
      util::ShPtr< util::ShPtrVector< assemble::SSEGeometryInterface> > restraint_bodies
      (
        new util::ShPtrVector< assemble::SSEGeometryInterface>()
      );

      // create coord::Bodies out of the SSEs of "sses" and insert the coord::Bodies into "bodies"
      for
      (
        util::SiPtrVector< const assemble::SSE>::iterator sse_itr( sses.Begin()), sse_itr_end( sses.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        util::ShPtr< assemble::SSE> current_body( ( *sse_itr)->Clone());

        BCL_MessageDbg
        (
          "sse aa size : " + util::Format()( current_body->GetSize()) + " and length : "
          + util::Format()( 2 * current_body->GetExtent( coord::GetAxes().e_Z))
        );

        // if the inserted body corresponded to a helix, then shorten its fragments
        if( current_body->GetType() == biol::GetSSTypes().HELIX)
        {
          BCL_MessageDbg
          (
            "restraint helix fragment size before change: " +
              util::Format()( current_body->GetGeometries().FirstElement()->GetExtent( coord::GetAxes().e_Z))
          );

          // reduce the size of the helix fragment so that end points fall into density rod
          current_body->SetFragmentGeometries( 5);

          BCL_MessageDbg
          (
            "restraint helix fragment size after change: " +
              util::Format()( current_body->GetGeometries().FirstElement()->GetExtent( coord::GetAxes().e_Z))
          );
        }

        // if the inserted body corresponded to a strand, then change its x and y extents and shorten its fragments
        else if( current_body->GetType() == biol::GetSSTypes().STRAND)
        {
          BCL_MessageDbg
          (
            "restraint strand fragment size before change: " +
              util::Format()( current_body->GetGeometries().FirstElement()->GetExtent( coord::GetAxes().e_Z))
          );

          // reduce the size of the strand fragment so that end points fall into density rod
          current_body->SetFragmentGeometries( 2);

          BCL_MessageDbg
          (
            "restraint strand fragment size after change: " +
            util::Format()( current_body->GetGeometries().FirstElement()->GetExtent( coord::GetAxes().e_Z))
          );

//             // set x and y to 1.6, leave z extent unchanged
//             // (this is necessary for the correct occupancy check, otherwise more than one body restraint could
//             // be occupied by the center of one sse, due to idealization))
//             bodies.LastElement().SetExtents
//             (
//               linal::Vector3D( 1.6, 1.6, bodies.LastElement().GetExtent( coord::GetAxes().e_Z))
//             );
        }
        restraint_bodies->PushBack( current_body);
      }

      BCL_MessageCrt
      (
        "number of body restraints created : " + util::Format()( restraint_bodies->GetSize())
      );

      return restraint_bodies;
    }

  } // namespace restraint
} // namespace bcl
