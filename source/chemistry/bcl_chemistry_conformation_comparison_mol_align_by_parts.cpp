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
#include "chemistry/bcl_chemistry_conformation_comparison_mol_align_by_parts.h"
#include "chemistry/bcl_chemistry_conformation_comparison_psi_flex_field.h"
// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_interface.h"
#include "chemistry/bcl_chemistry_priority_dihedral_angles.h"
#include "descriptor/bcl_descriptor_constants.h"
#include "io/bcl_io_serialization.h"
#include "model/bcl_model_interface_retrieve_from_file.h"
#include "model/bcl_model_kappa_nearest_neighbor.h"
#include "model/bcl_model_retrieve_interface.h"
#include "util/bcl_util_assert.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ConformationComparisonMolAlignByParts::s_Instance
    (
      util::Enumerated< ConformationComparisonInterface>::AddInstance( new ConformationComparisonMolAlignByParts())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor
    ConformationComparisonMolAlignByParts::ConformationComparisonMolAlignByParts() :
      ConformationComparisonPsiField(),
      m_MobileAtoms( "")
    {
    }

    //! virtual copy constructor
    ConformationComparisonMolAlignByParts *ConformationComparisonMolAlignByParts::Clone() const
    {
      return new ConformationComparisonMolAlignByParts( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ConformationComparisonMolAlignByParts::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &ConformationComparisonMolAlignByParts::GetAlias() const
    {
      static std::string s_name( "MolAlignByParts");
      return s_name;
    }

  /////////////////
  //  operations //
  /////////////////

    //! @brief make conformers of input molecules, randomly pair conformers
    //! @brief align and find best property rmsd between input molecules
    //! @param MOLECULE_A - first molecule from which conformational library will be generated
    //! @param MOLECULE_B - second molecule from which conformational library will be generated

    double ConformationComparisonMolAlignByParts::operator()
    (
      const ConformationInterface &MOLECULE_A,
      const ConformationInterface &MOLECULE_B
    ) const
    {
      // set the KeepIndices, which are those atom indices that will be included in scoring
      storage::Vector< size_t> keep_indices_a, keep_indices_b;
      this->GetNonMaskedAtoms( MOLECULE_A, MOLECULE_B, m_ExclusionIndicesA, m_ExclusionIndicesB, keep_indices_a, keep_indices_b);
      m_KeepIndicesA = keep_indices_a;
      m_KeepIndicesB = keep_indices_b;

      FragmentEnsemble ens_a, ens_b;

      return this->MolAlignByParts( ens_a, ens_b)( 0).Third();
    }

    storage::Vector< storage::Triplet< FragmentComplete, FragmentComplete, double> > ConformationComparisonMolAlignByParts::MolAlignByParts
    (
      const FragmentEnsemble &CONFORMERS_MOL_A,
      const FragmentEnsemble &CONFORMERS_MOL_B
    ) const
    {
      return storage::Vector< storage::Triplet< FragmentComplete, FragmentComplete, double> >();
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ConformationComparisonMolAlignByParts::GetSerializer() const
    {
      io::Serializer parameters( ConformationComparisonPsiField::GetSerializer());
      parameters.SetClassDescription( "Aligns and computes property RMSD for molecule conformer libraries by parts");

      parameters.AddInitializer
      (
        "sample_by_parts",
        "allow sampling of these atom indices during alignment",
        io::Serialization::GetAgent( &m_MobileAtoms),
        ""
      );

      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool ConformationComparisonMolAlignByParts::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      return ConformationComparisonPsiField::ReadInitializerSuccessHook( LABEL, ERR_STREAM);
    }
  } // namespace chemistry
} // namespace bcl
