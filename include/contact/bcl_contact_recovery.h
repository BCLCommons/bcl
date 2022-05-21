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

#ifndef BCL_CONTACT_RECOVERY_H_
#define BCL_CONTACT_RECOVERY_H_

// include the namespace header
#include "bcl_contact.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "math/bcl_math_contingency_matrix.h"
#include "math/bcl_math_contingency_matrix_measures.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "util/bcl_util_sh_ptr.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Recovery
    //! @brief return the percentage of contacts recovered from one given model to the other one
    //! @details This class implements the contact recovery quality measure for comparing protein models. It measures
    //! the percentage of contacts from the first given protein model (ideally a native/template model) that were
    //! recovered in the second given model
    //!
    //! @see @link example_contact_recovery.cpp @endlink
    //! @author karakam
    //! @date Sep 22, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Recovery :
      public math::BinaryFunctionInterfaceSerializable< assemble::ProteinModel, assemble::ProteinModel, double>
    {

    private:

    //////////
    // data //
    //////////

      //! min sequence separation range
      size_t m_MinSequenceSeparation;

      //! pointer to the neighbor list generator
      util::ShPtr< math::FunctionInterfaceSerializable< assemble::ProteinModel, assemble::AANeighborListContainer> > m_NeighborGenerator;

      //! Measure to compute
      math::ContingencyMatrixMeasures m_Measure;

      //! normalization factor; 100 for rates, 1 otherwise
      double m_Normalization;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a minimum sequence separation and cache boolean
      //! @param MIN_SEQUENCE_SEPARATION Minimum sequence separation
      //! @param CACHE boolean to decide whether a neighbor generator with cache should be used
      Recovery
      (
        const size_t MIN_SEQUENCE_SEPARATION = g_ContactMinSequenceSeparation,
        const math::ContingencyMatrixMeasures::MeasureEnum &MEASURE = math::ContingencyMatrixMeasures::e_TruePositiveRate
      );

      //! @brief Clone function
      //! @return pointer to new Recovery
      Recovery *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief calculate and return the percentage of recovered contacts
      //! @param TEMPLATE_MODEL Template model of interest
      //! @param PROTEIN_MODEL ProteinModel to be evaluated
      //! @return the percentage of recovered contacts
      double operator()
      (
        const assemble::ProteinModel &TEMPLATE_MODEL,
        const assemble::ProteinModel &PROTEIN_MODEL
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief calculate and return the contingency matrix for two neighbor list containers
      //! @param TEMPLATE_NEIGHBOR_CONTAINER AANeighborListContainer for a template model
      //! @param MODEL_NEIGHBOR_CONTAINER AANeighborListContainer for the model to be evaluated
      //! @return calculate and return the contingency matrix for two neighbor list containers
      static math::ContingencyMatrix CalculateContingencyMatrix
      (
        const assemble::AANeighborListContainer &TEMPLATE_NEIGHBOR_CONTAINER,
        const assemble::AANeighborListContainer &MODEL_NEIGHBOR_CONTAINER
      );

    }; // class Recovery

  } // namespace contact
} // namespace bcl

#endif // BCL_CONTACT_RECOVERY_H_
