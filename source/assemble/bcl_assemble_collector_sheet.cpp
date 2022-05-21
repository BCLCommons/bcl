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
#include "assemble/bcl_assemble_collector_sheet.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_collector_topology_sheet.h"
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> CollectorSheet::s_Instance
    (
      GetObjectInstances().AddInstance( new CollectorSheet())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    CollectorSheet::CollectorSheet()
    {
    }

    //! @brief Clone function
    //! @return pointer to new Sheet
    CollectorSheet *CollectorSheet::Clone() const
    {
      return new CollectorSheet( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &CollectorSheet::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &CollectorSheet::GetAlias() const
    {
      static const std::string s_name( "CollectorSheet");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer CollectorSheet::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Collects beta-sheets in a protein model.");

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Collect returns all beta-sheets in the ProteinModel argument
    //! @param PROTEIN_MODEL ProteinModel from which SSEs will be collected
    //! @return returns ShPtrList of beta-sheets found in the given protein model
    util::ShPtrVector< Domain> CollectorSheet::Collect( const ProteinModel &PROTEIN_MODEL) const
    {
      // construct sheet vector
      util::ShPtrVector< Domain> sheet_vector;

      // initialize vector to hold strands
      util::SiPtrVector< const SSE> strand_vector( PROTEIN_MODEL.GetSSEs( biol::GetSSTypes().STRAND));

      // collect the SheetTopologies
      util::ShPtrVector< Topology> topology_vector( CollectorTopologySheet().Collect( strand_vector));

      // if there are no strands or no beta-sheet topologies were found
      if( topology_vector.IsEmpty())
      {
        return sheet_vector;
      }

      // iterate over each found topology
      for
      (
        util::ShPtrVector< Topology>::iterator
          topology_itr( topology_vector.Begin()), topology_itr_end( topology_vector.End());
        topology_itr != topology_itr_end; ++topology_itr
      )
      {
        // cast the order vector into a sse vector
        util::SiPtrVector< const SSE> sse_vector( ( *topology_itr)->GetElements());

        // make sure the cast was possible
        if( !sse_vector.IsDefined())
        {
          BCL_MessageCrt( "The cast of SSEGeometryInterfaces to SSEs failed");
          continue;
        }

        // construct the sse set
        storage::Set< util::ShPtr< SSE>, SSELessThanNoOverlap> sse_set;

        // iterate over SSE vector
        for
        (
          util::SiPtrVector< const SSE>::const_iterator sse_itr( sse_vector.Begin()), sse_itr_end( sse_vector.End());
          sse_itr != sse_itr_end; ++sse_itr
        )
        {
          // find the corresponding SSE
          util::ShPtr< SSE> sp_this_sse( PROTEIN_MODEL.FindSSE( **sse_itr));

          // make sure the SSE was found
          BCL_Assert( sp_this_sse.IsDefined(), "The SSE from Sheet topology was not found in the Protein Model");

          // insert into the set
          sse_set.Insert( sp_this_sse);
        }

        // construct the sheet and add it to the sheet vector
        util::ShPtr< Domain> sp_new_sheet( new Domain( sse_set, *topology_itr));

        // pushback into the sheet vector
        sheet_vector.PushBack( sp_new_sheet);
      }

      // end
      return sheet_vector;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &CollectorSheet::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &CollectorSheet::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
