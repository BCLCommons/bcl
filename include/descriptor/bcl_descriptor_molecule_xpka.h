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

#ifndef BCL_DESCRIPTOR_MOLECULE_XPKA_H_
#define BCL_DESCRIPTOR_MOLECULE_XPKA_H_

// include the namespace header
#include "bcl_descriptor.h"
#include "bcl_descriptor_combine.h"
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "model/bcl_model_retrieve_interface.h"
#include "util/bcl_util_sh_ptr_list.h"

// include other forward headers - sorted alphabetically
#include "chemistry/bcl_chemistry.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_sequence.h"
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeXpKa
    //! // TODO
    //!
    //! @see @link example_descriptor_molecule_xp_ka.cpp @endlink
    //! @author decortja
    //! @date Mar 20, 2023
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeXpKa :
      public BaseSequence< chemistry::AtomConformationalInterface, float>
    {

    //////////
    // data //
    //////////

      //! true for acid, false for base
      bool m_AcidModel;

    public:

    //////////
    // data //
    //////////

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! @param  ACID_MODEL true if acid model to be used. Initialized to true.
      MoleculeXpKa( const bool &ACID_MODEL = true);

      //! @brief Clone function
      //! @return pointer to new MoleculeXpKa
      MoleculeXpKa *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get name of the current class
      //! @return name of the class
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @brief set return to 16 to explore all 16 buckets of the multitask classifier.
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const
      {
        return 16;
      }

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief calculate the descriptors
      //! @param STORAGE storage for the descriptor
      void Calculate( linal::VectorReference< float> &STORAGE);
//      void Calculate( storage::Vector< float> STORAGE);

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook
      (
        const util::ObjectDataLabel &LABEL,
        std::ostream &ERR_STREAM
      );

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class MoleculeXpKa

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_MOLECULE_XPKA_H_
