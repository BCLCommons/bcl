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
#include "contact/bcl_contact_recovery.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_list_container_generator_protein_model.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Recovery::s_Instance
    (
      GetObjectInstances().AddInstance( new Recovery())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a minimum sequence separation and cache boolean
    //! @param MIN_SEQUENCE_SEPARATION Minimum sequence separation
    //! @param CACHE boolean to decide whether a neighbor generator with cache should be used
    Recovery::Recovery
    (
      const size_t MIN_SEQUENCE_SEPARATION,
      const math::ContingencyMatrixMeasures::MeasureEnum &MEASURE
    ) :
      m_MinSequenceSeparation( MIN_SEQUENCE_SEPARATION),
      m_NeighborGenerator
      (
        assemble::AANeighborListContainerGeneratorProteinModel::AANeighborListGenerator
        (
          g_ContactCbDistanceCutoff,
          g_ContactMinSequenceSeparation,
          true,
          true
        )
      ),
      m_Measure( MEASURE),
      m_Normalization
      (
        util::EndsWith( math::ContingencyMatrixMeasures::GetMeasureName( m_Measure.GetMeasure()), "R")
        ? 100.0
        : 1.0
      )
      {
        BCL_Assert
        (
          MIN_SEQUENCE_SEPARATION >= g_ContactMinSequenceSeparation,
          "Given sequence separation " + util::Format()( m_MinSequenceSeparation) + " can't be smaller than " +
          util::Format()( g_ContactMinSequenceSeparation)
        );
      }

      //! @brief Clone function
      //! @return pointer to new Recovery
      Recovery *Recovery::Clone() const
      {
        return new Recovery( *this);
      }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Recovery::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate and return the percentage of recovered contacts
    //! @param TEMPLATE_MODEL Template model of interest
    //! @param PROTEIN_MODEL ProteinModel to be evaluated
    //! @return the percentage of recovered contacts
    double Recovery::operator()
    (
      const assemble::ProteinModel &TEMPLATE_MODEL,
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // generate the neighbor list containers for both models
      assemble::AANeighborListContainer neighbors_template( m_NeighborGenerator->operator()( TEMPLATE_MODEL));
      assemble::AANeighborListContainer neighbors_model( m_NeighborGenerator->operator()( PROTEIN_MODEL));

      // prune both by the min sequence separation
      if( m_MinSequenceSeparation > g_ContactMinSequenceSeparation)
      {
        neighbors_template.Prune( g_ContactCbDistanceCutoff, m_MinSequenceSeparation, true);
        neighbors_model.Prune( g_ContactCbDistanceCutoff, m_MinSequenceSeparation, true);
      }

      // calculate the contingency matrix
      const math::ContingencyMatrix contingency_matrix( CalculateContingencyMatrix( neighbors_template, neighbors_model));

      // if there were no positives found
      if( contingency_matrix.GetNumberActualPositives() == 0)
      {
        BCL_MessageStd( "No contacts found in the template model, returning 0");
        return double( 0.0);
      }

      return m_Normalization * m_Measure( contingency_matrix);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Recovery::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_MinSequenceSeparation, ISTREAM);
      io::Serialize::Read( m_Measure, ISTREAM);
      *this = Recovery( m_MinSequenceSeparation, m_Measure.GetMeasure());

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Recovery::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_MinSequenceSeparation, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Measure, OSTREAM, INDENT) << '\n';

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculate and return the contingency matrix for two neighbor list containers
    //! @param TEMPLATE_NEIGHBOR_CONTAINER AANeighborListContainer for a template model
    //! @param MODEL_NEIGHBOR_CONTAINER AANeighborListContainer for the model to be evaluated
    //! @return calculate and return the contingency matrix for two neighbor list containers
    math::ContingencyMatrix Recovery::CalculateContingencyMatrix
    (
      const assemble::AANeighborListContainer &TEMPLATE_NEIGHBOR_CONTAINER,
      const assemble::AANeighborListContainer &MODEL_NEIGHBOR_CONTAINER
    )
    {
      // store the number of contacts for template and model
      const size_t nr_contacts_template( TEMPLATE_NEIGHBOR_CONTAINER.GetNumberNeighbors());
      const size_t nr_contacts_model( MODEL_NEIGHBOR_CONTAINER.GetNumberNeighbors());
      const size_t nr_common_contacts( TEMPLATE_NEIGHBOR_CONTAINER.IntersectionSize( MODEL_NEIGHBOR_CONTAINER));
      const size_t min_seq_sep( TEMPLATE_NEIGHBOR_CONTAINER.GetMinimalSequenceSeparation());
      const size_t chain_size
      (
        std::max( std::max( MODEL_NEIGHBOR_CONTAINER.GetSize(), size_t( min_seq_sep + 1)), nr_contacts_model)
      );

      // initialize a contingency matrix from the intersect values
      const math::ContingencyMatrix matrix
      (
        nr_common_contacts,                        // true positives
        nr_contacts_model - nr_common_contacts,    // false positives
        nr_contacts_template - nr_common_contacts, // false negatives
        chain_size * ( chain_size - min_seq_sep - 1) - nr_contacts_model // true negatives
      );

      // end
      return matrix;
    }

  } // namespace contact
} // namespace bcl
