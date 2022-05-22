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

#ifndef BCL_CHEMISTRY_FRAGMENT_DERIVATIVE_MANAGER_H_
#define BCL_CHEMISTRY_FRAGMENT_DERIVATIVE_MANAGER_H_

// forward headers from bcl - sorted alphabetically

// headers from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_graph_converter.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_triplet.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_object_data_label.h"
#include "util/bcl_util_serializable_interface.h"
#include "util/bcl_util_sh_ptr.h"

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FragmentDerivativeManager
    //! @brief Stores, retrieves, and dumps data relating fragments during design.
    //!
    //! @see @link example_chemistry_fragment_derivative_manager.cpp @endlink
    //! @author brownbp1
    //! @date May 07, 2022
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentDerivativeManager :
      public util::SerializableInterface
    {

    //////////////////////
    // type definitions //
    //////////////////////

    public:

      //! rename for convenience and readability
      //! First subgraph should reference m_Molecules, second subgraph should reference m_FragmentDerivatives
      typedef storage::Pair< graph::Subgraph< size_t, size_t>, graph::Subgraph< size_t, size_t>> AtomMap;

    //////////
    // data //
    //////////

    private:

      //! Input molecules that serve as references for the fragment derivatives
      util::ShPtrVector< FragmentComplete> m_Molecules;

      //! Derivatives whose information relative to m_Molecules is stored
      util::ShPtrVector< FragmentComplete> m_FragmentDerivatives;

      //! Atom maps between reference and derivatives
      //! Outer vector indexed by m_Molecules, inner vector indexed by m_FragmentDerivatives
      storage::Vector< storage::Vector< AtomMap>> m_AtomMaps;

      //! One-hot vector mask that controls which fragment derivatives are used in analysis
      //! Outer vector indexed by m_Molecules, inner vector indexed by m_FragmentDerivatives
      storage::Vector< storage::Vector< size_t>> m_Mask;

      //! Features for m_Molecules
      storage::Vector< linal::Vector< float>> m_MoleculeFeatures;

      //! Features for m_FragmentDerivatives
      storage::Vector< linal::Vector< float>> m_FragmentDerivativesFeatures;

      //! Descriptors to use to compute features
      descriptor::CheminfoProperty m_Descriptors;

      //! Perform analysis on only uncommon atoms between reference and derivative
      bool m_Complement;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FragmentDerivativeManager();

      //! @brief constructor from FragmentDerivativeManager
      FragmentDerivativeManager( const FragmentDerivativeManager &MANAGER);

      //! @brief constructor from data
      FragmentDerivativeManager
      (
        const util::ShPtrVector< FragmentComplete> &MOLECULES,
        const util::ShPtrVector< FragmentComplete> &FRAGMENT_DERIVATIVES,
        const storage::Vector< storage::Vector< size_t>> &MASK,
            const descriptor::CheminfoProperty &DESCRIPTORS,
            const bool COMPLEMENT
          );

      //! @brief Clone function
      //! @return pointer to new FragmentDerivativeManager
      FragmentDerivativeManager *Clone() const
      {
        return new FragmentDerivativeManager( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get the name of this class
      //! @return the name of this class
      const std::string &GetAlias() const
      {
        static const std::string s_name( "FragmentDerivativeManager");
        return s_name;
      }

    ////////////////
    // operations //
    ////////////////

      void AddMolecule( const FragmentComplete &MOLECULE);

      void AddFragmentDerivative( const FragmentComplete &FRAGMENT_DERIVATIVE);

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief Set the members with LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook
      (
        const util::ObjectDataLabel &LABEL,
        std::ostream &ERR_STREAM
      );

      io::Serializer GetSerializer() const;

    }; // class FragmentDerivativeManager

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_FRAGMENT_DERIVATIVE_MANAGER_H_
