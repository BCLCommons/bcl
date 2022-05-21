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
#include "score/bcl_score_protein_model_fragment_topology.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_collector_topology_combined.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_geometry_packing_pickers.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "assemble/bcl_assemble_sse_pool_agreement.h"
#include "assemble/bcl_assemble_topology.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
  //////////
  // data //
  //////////
    
    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ProteinModelFragmentTopology::s_Instance
    (
      GetObjectInstances().AddInstance( new ProteinModelFragmentTopology())
    );
   
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////
    
    //! @brief Clone function
    //! @return pointer to new ProteinModelFragmentTopology
    ProteinModelFragmentTopology *ProteinModelFragmentTopology::Clone() const
    {
      return new ProteinModelFragmentTopology( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinModelFragmentTopology::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate and return the percentage of recovered contacts
    //! @param TEMPLATE the known correct protein model
    //! @param MODEL the protein model to be evaluated against the TEMPLATE
    //! @return the percentage of recovered contacts
    double ProteinModelFragmentTopology::operator()
    ( 
      const assemble::ProteinModel &TEMPLATE, 
      const assemble::ProteinModel &MODEL
    ) const
    {
      // get sequences
      util::SiPtrVector< const biol::AASequence> template_seqs( TEMPLATE.GetSequences());
      util::SiPtrVector< const biol::AASequence> model_seqs( MODEL.GetSequences());

      size_t model_size( MODEL.GetNumberAAs());
      size_t template_size( TEMPLATE.GetNumberAAs());

      util::SiPtrVector< const biol::AASequence>::const_iterator template_seq_itr( template_seqs.Begin());
      for
      ( 
        util::SiPtrVector< const biol::AASequence>::const_iterator model_seq_itr( model_seqs.Begin()),
        model_seq_itr_end( model_seqs.End());
        model_seq_itr != model_seq_itr_end;
        model_seq_itr++, template_seq_itr++
      )
      {
        // make sure template and model have the same sequence
        const biol::AASequence model_seq( **model_seq_itr);
        const biol::AASequence template_seq( **template_seq_itr);

         BCL_Assert
        (
          model_seq == template_seq,
          "model and template must have the same sequence"
        );
      }

      //initialize matrices of fragment interactions for template and model
      linal::Matrix< double> model_matrix( MODEL.CalculateFragmentInteractions());
      linal::Matrix< double> template_matrix( TEMPLATE.CalculateFragmentInteractions());

      //calculate distance between the two matrices
      double temp_value( 0);
      double model_value( 0);

      // calculate euclidean distance between matices
      double norm( linal::Distance( model_matrix.AsVector(), template_matrix.AsVector()));

      return norm;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProteinModelFragmentTopology::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProteinModelFragmentTopology::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }
    
  /////////////////////
  // helper funtions //
  /////////////////////
    
  } // namespace score
} // namespace bcl
