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

#ifndef BCL_DESCRIPTOR_AA_DSSP_INFO_H_
#define BCL_DESCRIPTOR_AA_DSSP_INFO_H_

// include the namespace header
#include "bcl_descriptor.h"

// other forward includes
#include "model/bcl_model.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"
#include "bcl_descriptor_type.h"
#include "biol/bcl_biol_aa_sequence.h"
#include "iterate/bcl_iterate_generic.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AADSSPInfo
    //! @brief AADSSPInfo gives the size of the in which this SSE resides, index of the SSE the given AA is in, the
    //! position (0-1.0) from N to C within the SSE, and the normalized distance from the center of the SSE
    //!
    //! @see @link example_descriptor_aa_dssp_info.cpp @endlink
    //! @author mendenjl
    //! @date May 12, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AADSSPInfo :
      public BaseElement< biol::AABase, float>
    {
    public:

    //////////
    // data //
    //////////

      //! custom DSSP info type
      enum DSSPInfoType
      {
        e_Accessibility,            //!< Residue water exposed surface area in A^2
        e_TotalEnergy,              //!< Total h-bond energy determined for this AA
        e_MaxHBondPotential,        //!< Best (lowest) energy for any neighbor
        e_MaxHBondNeighborOffset,   //!< Offset of the neighbor with max hbond energy
        e_IsInSheetCore,            //!< 1 if residue is h-bonded to another strand and at least one neighboring residue is too
        e_IsOnSheetEdge,            //!< 1 if the residue is labeled a strand but neither neighbor is bonded to a strand, though this residue is
        e_NonBondedStrandResidue,   //!< 1 if the residue is labeled a strand but this residue is not hbonded to a strand
        e_BondedStrandResidue,      //!< 1 if the residue is labeled a strand and this residue is hbonded to a strand
        // for any strand residue, exactly one of IsInSheetCore IsOnSheetEdge or AdjacentStrandResidue is true
        // BondedStrandResidue = (IsInSheetCore || IsOnSheetEdge)
        e_BondsParallel,            //!< True if the residue bonds to another strand oriented parallel to this strand
        e_BondsAntiparallel,        //!< True if the residue bonds to another strand oriented antiparallel to this strand
        e_BondsBridge,              //!< True if the residue bonds to another strand only at this residue
        // BondsParallel | BondsAntiparallel | BondsBridge = BondedStrandResidue
        s_NumberDSSPInfoTypes
      };

      //! @brief DSSPInfoType as string
      //! @param DSSP_INFO_TYPE the message level
      //! @return the DSSPInfoType as string
      static const std::string &GetDSSPInfoTypeString( const DSSPInfoType &DSSP_INFO_TYPE);

      //! @brief DSSPInfoType description
      //! @param DSSP_INFO_TYPE the description
      //! @return the DSSPInfoType's description
      static const std::string &GetDSSPInfoTypeDescription( const DSSPInfoType &DSSP_INFO_TYPE);

      //! DSSPInfoTypeEnum simplifies the usage of the DSSPInfoType enum of this class
      typedef util::WrapperEnum< DSSPInfoType, &GetDSSPInfoTypeString, s_NumberDSSPInfoTypes> DSSPInfoTypeEnum;

    private:

      //! Vector containing the mapping from all seq IDs to indices within the SSE info vectors
      storage::Map< int, size_t> m_SeqIDToIndexMap;

      //! Vector containing the dssp info requested
      linal::Vector< float> m_Info;

      //! sse type info desired
      DSSPInfoTypeEnum m_InfoTypeEnum;

      //! dssp file suffix
      std::string m_DsspExtension;

      //! single instance of the class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AADSSPInfo( const DSSPInfoTypeEnum &INFO_TYPE_ENUM);

      //! @brief virtual copy constructor
      AADSSPInfo *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const;

    private:

      //! @brief calculate the descriptors
      //! @param ELEMENT_A the element of interest
      //! @param STORAGE storage for the descriptor
      void Calculate
      (
        const iterate::Generic< const biol::AABase> &ELEMENT,
        linal::VectorReference< float> &STORAGE
      );

      //! @brief function to load files; should only be called the first time Calculate is called with a new sequence
      //! since the results are often in the cache
      void LoadFiles();

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      virtual void SetObjectHook();

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;
    }; // class AADSSPInfo

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_AA_DSSP_INFO_H_
