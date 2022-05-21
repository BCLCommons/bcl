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
#include "fold/bcl_fold_mutate_protein_model_multiple_geometries.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_fold_template_handler.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_pool.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelMultipleGeometries::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateProteinModelMultipleGeometries())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelMultipleGeometries::MutateProteinModelMultipleGeometries() :
      m_FoldTemplate(),
      m_TemplateSizeProbabilities( 0.0, 1.0, 0.0),
      m_Scheme( GetStaticClassName< MutateProteinModelMultipleGeometries>())
    {
    }

    //! @brief constructor from a fold template and scheme
    //! @param FOLD_TEMPLATE fold template to be used
    //! @param SCHEME Scheme to be used
    MutateProteinModelMultipleGeometries::MutateProteinModelMultipleGeometries
    (
      const assemble::FoldTemplate &FOLD_TEMPLATE,
      const std::string &SCHEME
    ) :
      m_FoldTemplate( FOLD_TEMPLATE),
      m_TemplateSizeProbabilities( 0.0, 1.0, 0.0),
      m_Scheme( SCHEME)
    {
    }

    //! @brief constructor from a fold template handler and scheme
    //! @param PROBABILITIES probability of fitting into a small, equal, or large template, respectively
    //! @param SCHEME Scheme to be used
    MutateProteinModelMultipleGeometries::MutateProteinModelMultipleGeometries
    (
      const storage::VectorND< 3, double> &PROBABILITIES,
      const std::string &SCHEME
    ) :
      m_FoldTemplate(),
      m_TemplateSizeProbabilities( PROBABILITIES),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateProteinModelMultipleGeometries
    MutateProteinModelMultipleGeometries *MutateProteinModelMultipleGeometries::Clone() const
    {
      return new MutateProteinModelMultipleGeometries( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateProteinModelMultipleGeometries::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator taking a protein model and returning a mutate object of protein model type
    //! @param PROTEIN_MODEL protein model interest
    //! @return MutateResult with ProteinModel after the mutate
    math::MutateResult< assemble::ProteinModel> MutateProteinModelMultipleGeometries::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // initialize empty model
      static util::ShPtr< assemble::ProteinModel> s_empty_model;

      // get the pool from the given ProteinModel and make sure it is valid
      const util::ShPtr< assemble::SSEPool> sp_pool
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Pool)
      );
      BCL_Assert( sp_pool.IsDefined(), "No pool stored for the given model");

      // get non-overlapping SSEs from the pool
      storage::Set< util::SiPtr< const assemble::SSE>, assemble::SSELessThanNoOverlap> pool_no_overlap
      (
        sp_pool->GetRandomNonOverlappingSet()
      );

      // if no sses were found
      if( pool_no_overlap.IsEmpty())
      {
        // return an empty model
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      // initialize vector of helices and strands
      util::SiPtrVector< const assemble::SSE> combined_sses( pool_no_overlap.Begin(), pool_no_overlap.End());

      // copy the fold template member so that it can be changed if it is empty
      assemble::FoldTemplate fold_template( m_FoldTemplate);

      // check if the given fold template is empty
      if( fold_template.GetGeometries().IsEmpty())
      {
        // get a random number based on the probabilities (between 0 and the sum)
        double random_double
        (
          random::GetGlobalRandom().Double
          (
            math::Range< double>
            (
              0.0,
              m_TemplateSizeProbabilities.First() +
              m_TemplateSizeProbabilities.Second() +
              m_TemplateSizeProbabilities.Third()
            )
          )
        );

        // if fitting into a small template
        if( random_double <= m_TemplateSizeProbabilities.First())
        {
          // remove roughly 1/3 of the SSEs
          for( size_t i( 0), i_end( combined_sses.GetSize() / 3); i != i_end; ++i)
          {
            combined_sses.RemoveRandomElement();
          }
        }

        // if fitting into a small or equal template
        if( random_double <= m_TemplateSizeProbabilities.First() + m_TemplateSizeProbabilities.Second())
        {
          fold_template = assemble::FoldTemplateHandler::GetRandomTemplate( combined_sses);
        }
        // if fitting into a large template
        else
        {
          fold_template = assemble::FoldTemplateHandler::GetRandomSubTemplate( combined_sses);
        }
      }

      // if the fold template is empty, return an empty model
      if( fold_template.GetGeometries().IsEmpty())
      {
        return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
      }

      BCL_MessageVrb
      (
        "Using fold template from PDB ID " + fold_template.GetPDBID() + " for the mutate"
      );

      // fit the SSEs to the fold template which returns a domain
      const assemble::Domain fitted_domain( fold_template.FitSSEs( combined_sses));

      // get a new protein model
      util::ShPtr< assemble::ProteinModel> new_model( PROTEIN_MODEL.Clone());

      // iterate through the SSEs in the domain
      const util::SiPtrVector< const assemble::SSE> domain_sses( fitted_domain.GetSSEs());
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator sse_itr( domain_sses.Begin()),
          sse_itr_end( domain_sses.End());
        sse_itr != sse_itr_end; ++sse_itr
      )
      {
        // insert the SSE into the model
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
    std::istream &MutateProteinModelMultipleGeometries::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_FoldTemplate, ISTREAM);
      io::Serialize::Read( m_TemplateSizeProbabilities, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &MutateProteinModelMultipleGeometries::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_FoldTemplate, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_TemplateSizeProbabilities, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
