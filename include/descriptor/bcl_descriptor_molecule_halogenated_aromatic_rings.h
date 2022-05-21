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

#ifndef BCL_DESCRIPTOR_MOLECULE_HALOGENATED_AROMATIC_RINGS_H_
#define BCL_DESCRIPTOR_MOLECULE_HALOGENATED_AROMATIC_RINGS_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "chemistry/bcl_chemistry.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_sequence.h"
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeHalogenatedAromaticRings
    //! Determines the number of halogens on aromatic rings in a molecule
    //!
    //! @see @link example_descriptor_molecule_halogenated_aromatic_rings.cpp @endlink
    //! @author brownbp1
    //! @date Aug 23, 2021
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeHalogenatedAromaticRings :
      public BaseSequence< chemistry::AtomConformationalInterface, float>
    {

    public:

        //! Count types
        enum AromaticRingHalogens
        {
          e_Total, //! total number of aromatic ring halogens
          e_Max    //! max number of aromatic ring halogens on a single fragment
        };

    private:

        //! The halogen counting method we are using
        AromaticRingHalogens m_CountMethod;

        //! halogens that are allowed
        storage::Vector< chemistry::AtomType> m_AllowedHalogens;
        std::string m_AllowedHalogensString;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MoleculeHalogenatedAromaticRings();

      //! @brief constructor
      //! @brief COUNT_METHOD how to count halogens
      MoleculeHalogenatedAromaticRings
      (
        const AromaticRingHalogens &COUNT_METHOD
      );

      //! @brief constructor with halogen types
      //! @brief COUNT_METHOD how to count halogens
      //! @brief ALLOWED_HALOGENS the halogens included in counting
      MoleculeHalogenatedAromaticRings
      (
        const AromaticRingHalogens &COUNT_METHOD,
        const storage::Vector< chemistry::AtomType> &ALLOWED_HALOGENS
      );

      //! @brief Clone function
      //! @return pointer to new MoleculeHalogenatedAromaticRings
      MoleculeHalogenatedAromaticRings *Clone() const;

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
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const
      {
        return 1;
      }

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief calculate the descriptors
      //! @param STORAGE storage for the descriptor
      void Calculate( linal::VectorReference< float> &STORAGE);

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief calculate the max number of halogen atoms on a ring in a molecule
      //! @param MOLECULE the molecule of interest
      size_t CountMaxAromaticHalogensSingleFragment( const chemistry::FragmentComplete &MOLECULE) const;

      //! @brief calculate the total number of non-aromatic ring Cl, Br, I halogen atoms
      //! @param MOLECULE the molecule of interest
      size_t CountAromaticHalogens( const chemistry::FragmentComplete &MOLECULE) const;

      //! @brief set the halogens to be counted
      void SetAllowedHalogens( const storage::Vector< chemistry::AtomType> &ALLOWED_HALOGENS);

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERROR_STREAM the stream to write errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class MoleculeHalogenatedAromaticRings

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_MOLECULE_HALOGENATED_AROMATIC_RINGS_H_
