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

#ifndef BCL_CHEMISTRY_FRAGMENT_ALIGN_TO_SCAFFOLD_H_
#define BCL_CHEMISTRY_FRAGMENT_ALIGN_TO_SCAFFOLD_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically
#include "descriptor/bcl_descriptor.fwd.hh"
#include "find/bcl_find.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_conformational_interface.h"
#include "bcl_chemistry_configurational_bond_type_data.h"
#include "bcl_chemistry_conformation_comparison_property_field_correlation.h"
#include "bcl_chemistry_conformation_graph_converter.h"
#include "bcl_chemistry_fragment_complete.h"
#include "bcl_chemistry_fragment_constitution_shared.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "find/bcl_find_pick_interface.h"
#include "math/bcl_math_mutate_interface.h"
#include "math/bcl_math_mutate_result.h"
#include "util/bcl_util_function_interface.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_si_ptr_list.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FragmentAlignToScaffold
    //! @brief Align molecules based on substructure.
    //!
    //! @see @link example_chemistry_fragment_align_to_scaffold.cpp @endlink
    //! @author brownbp1
    //! @date May 18, 2020
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentAlignToScaffold :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      // atom comparison type
      ConformationGraphConverter::AtomComparisonType m_AtomType;

      // bond comparison type
      ConfigurationalBondTypeData::Data m_BondType;

      // minimum isomorphism size for alignment to occur
      size_t m_MinIsoSize;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FragmentAlignToScaffold();

      //! @brief constructor with binding pocket to be created from MDL property
      //! @param ATOM_TYPE the atom comparison type to be used for substructure comparison
      //! @param BOND_TYPE the bond comparison type to be used for substructure comparison
      FragmentAlignToScaffold
      (
        const ConformationGraphConverter::AtomComparisonType &ATOM_TYPE,
        const ConfigurationalBondTypeData::Data &BOND_TYPE,
        const size_t &MIN_ISO_SIZE
      );

      //! @brief clone constructor
      FragmentAlignToScaffold *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief maximum common substructure alignment of small molecules
      //! @param TARGET_MOL the molecule to be aligned
      //! @param SCAFFOLD_MOL the molecule against which the target is aligned
      //! @return true if alignment occurs, false if the isomorphism size
      //! is below the minimum size allowed
      bool AlignToScaffold
      (
        FragmentComplete &TARGET_MOL,
        const FragmentComplete &SCAFFOLD_MOL,
        const storage::Vector< size_t> &TARGET_MOL_INDICES = storage::Vector< size_t>(),
        const storage::Vector< size_t> &SCAFFOLD_MOL_INDICES = storage::Vector< size_t>()
      ) const;

      //! @brief maximum common substructure alignment of small molecules with pose-dependent scoring
      //! @param TARGET_MOL the molecule to be aligned
      //! @param SCAFFOLD_MOL the molecule against which the target is aligned
      //! @param MDL the SDF file MDL property specifying the binding pocket filename
      //! @param BINDING_POCKET_FILENAME name of protein/protein binding pocket PDB file
      storage::Pair< bool, float> PoseSensitiveAlignToScaffold
      (
        FragmentComplete &TARGET_MOL,
        const FragmentComplete &SCAFFOLD_MOL,
        const descriptor::CheminfoProperty &SCORE,
        const std::string &MDL,
        const std::string &BINDING_POCKET_FILENAME
      );

      //! @brief maximum common substructure alignment of small molecule conformer ensembles
      //! @param TARGET_ENS the molecule ensemble to be aligned
      //! @param SCAFFOLD_MOL the molecule against which the targets are aligned
      //! @param COMPARER the metric to be used to compare alignments
      //! @return true if alignment occurs, false if the isomorphism size
      //! is below the minimum size allowed or the ensemble is empty
      storage::Vector< storage::Pair< bool, float> > AlignEnsembleToScaffold
      (
        FragmentEnsemble &TARGET_ENS,
        const FragmentComplete &SCAFFOLD_MOL,
        const util::Implementation< ConformationComparisonInterface> &COMPARER
          = util::Implementation< ConformationComparisonInterface>( "PropertyFieldDistance")
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief get a mutex
      static sched::Mutex &GetMutex();

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

    }; // class FragmentAlignToScaffold

  } // namespace chemistry
} // namespace bcl

#endif //BCL_CHEMISTRY_FRAGMENT_ALIGN_TO_SCAFFOLD_H_
