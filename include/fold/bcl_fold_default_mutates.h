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

#ifndef BCL_FOLD_DEFAULT_MUTATES_H_
#define BCL_FOLD_DEFAULT_MUTATES_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "sspred/bcl_sspred.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_fold_mutates.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DefaultMutates
    //! @brief mutates used in Default fold protocol
    //! @details Implements initialization of all the mutates used in the default protocol in addition to default
    //! mutate tree.
    //!
    //! @see @link example_fold_default_mutates.cpp @endlink
    //! @author karakam
    //! @date Nov 17, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DefaultMutates :
      public util::ObjectInterface
    {

    public:

    //////////
    // data //
    //////////

      // add sse
      Mutate e_AddSSENextToSSE;
      Mutate e_AddSSEShortLoop;
      Mutate e_AddStrandNextToSheet;

      // remove sse
      Mutate e_RemoveRandom;
      Mutate e_RemoveUnpairedStrand;

      // swap sses
      Mutate e_SwapSSEs;
      Mutate e_SwapSSEWithPool;
      Mutate e_SwapSSEWithPoolOverlap;

      // resize
      Mutate e_SSEResizeCTerm;
      Mutate e_SSEResizeNTerm;

      // split
      storage::Map< sspred::Method, Mutate> m_SSESplit;

      // move sse
      Mutate e_SSEMoveNext;
      Mutate e_SSEMoveShortLoop;
      Mutate e_SSEFurthestMoveNext;
      Mutate e_SSEBendRamachandran;
      Mutate e_SSEBendRandomSmall;
      Mutate e_SSEBendRandomLarge;
      Mutate e_SSETranslateSmall;
      Mutate e_SSETranslateXSmall;
      Mutate e_SSETranslateYSmall;
      Mutate e_SSETranslateZSmall;
      Mutate e_SSETranslateLarge;
      Mutate e_SSETranslateXLarge;
      Mutate e_SSETranslateYLarge;
      Mutate e_SSETranslateZLarge;
      Mutate e_SSERotateSmall;
      Mutate e_SSERotateXSmall;
      Mutate e_SSERotateYSmall;
      Mutate e_SSERotateZSmall;
      Mutate e_SSERotateLarge;
      Mutate e_SSERotateXLarge;
      Mutate e_SSERotateYLarge;
      Mutate e_SSERotateZLarge;
      Mutate e_SSETransformSmall;
      Mutate e_SSETransformLarge;

      // flip
      Mutate e_SSEFlipX;
      Mutate e_SSEFlipY;
      Mutate e_SSEFlipZ;

      // helix moves
      Mutate e_HelixMoveNext;
      Mutate e_HelixMoveShortLoop;
      Mutate e_HelixFurthestMoveNext;
      Mutate e_HelixTranslateXYSmall;
      Mutate e_HelixTranslateXYLarge;
      Mutate e_HelixTranslateZSmall;
      Mutate e_HelixTranslateZLarge;
      Mutate e_HelixRotateXYSmall;
      Mutate e_HelixRotateXYLarge;
      Mutate e_HelixRotateZSmall;
      Mutate e_HelixRotateZLarge;
      Mutate e_HelixTransformXYSmall;
      Mutate e_HelixTransformXYLarge;
      Mutate e_HelixTransformZSmall;
      Mutate e_HelixTransformZLarge;
      Mutate e_HelixFlipXY;
      Mutate e_HelixFlipZ;

      // strand moves
      Mutate e_StrandMoveNext;
      Mutate e_StrandFurthestMoveNext;
      Mutate e_StrandMoveSheet;
      Mutate e_StrandFurthestMoveSheet;
      Mutate e_StrandTranslateZSmall;
      Mutate e_StrandTranslateZLarge;
      Mutate e_StrandFlipX;
      Mutate e_StrandFlipY;
      Mutate e_StrandFlipZ;

      // sse pair
      Mutate e_SSEPairTranslateNoHingeSmall;
      Mutate e_SSEPairTranslateNoHingeLarge;
      Mutate e_SSEPairTranslateSmall;
      Mutate e_SSEPairTranslateLarge;
      Mutate e_SSEPairRotateSmall;
      Mutate e_SSEPairRotateLarge;
      Mutate e_SSEPairTransformSmall;
      Mutate e_SSEPairTransformLarge;

      // helix pair
      Mutate e_HelixPairRotateZSmallNoHinge;
      Mutate e_HelixPairRotateZSmallHinge;
      Mutate e_HelixPairRotateZLargeNoHinge;
      Mutate e_HelixPairRotateZLargeHinge;

      // helix domain
      Mutate e_HelixDomainShuffle;
      Mutate e_HelixDomainTranslateSmall;
      Mutate e_HelixDomainTranslateLarge;
      Mutate e_HelixDomainRotateSmall;
      Mutate e_HelixDomainRotateLarge;
      Mutate e_HelixDomainTransformSmall;
      Mutate e_HelixDomainTransformLarge;
      Mutate e_HelixDomainFlipExt;
      Mutate e_HelixDomainFlipInt;

      // sheet
      Mutate e_SheetRotateSmall;
      Mutate e_SheetRotateLarge;
      Mutate e_SheetTranslateSmall;
      Mutate e_SheetTranslateLarge;
      Mutate e_SheetTransformSmall;
      Mutate e_SheetTransformLarge;
      Mutate e_SheetPairStrands;
      Mutate e_SheetSwitchStrand;
      Mutate e_SheetFlipExt;
      Mutate e_SheetFlipInt;
      Mutate e_Sheet_flip_int_sub;
      Mutate e_Sheet_flip_int_sub_diff;
      Mutate e_Sheet_divide;
      Mutate e_Sheet_divide_sandwich;
      Mutate e_Sheet_twist_small;
      Mutate e_Sheet_twist_large;
      Mutate e_Sheet_shuffle;
      Mutate e_Sheet_cycle;
      Mutate e_Sheet_cycle_intact;
      Mutate e_Sheet_cycle_subset;
      Mutate e_Sheet_cycle_subset_intact;
      Mutate e_Sheet_register_fix;
      Mutate e_Sheet_register_shift;
      Mutate e_Sheet_register_shift_flip;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    private:

      //! @brief default constructor
      DefaultMutates();

    public:

      //! @brief Clone function
      //! @return pointer to new DefaultMutates
      DefaultMutates *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief access to single instance
      //! @return single instance of default mutates
      static DefaultMutates &GetInstance();

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief initialize the mutates and add them to Mutates enumerator
      void InitializeMutates();

      //! @brief modify the mutate tree used
      //! @param MUTATE_TREE MutateTree to be modified
      void ModifyMutateTree( MutateTree &MUTATE_TREE) const;

      //! @brief initializes the mutate tree for the given protein model
      //! @param TREE mutate tree that shall be optimized
      //! @param MODEL protein model for which the mutate tree shall be optimized
      void InitializeMutateTree
      (
        MutateTree &TREE, const assemble::ProteinModel &MODEL = assemble::ProteinModel()
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

    }; // class DefaultMutates

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_DEFAULT_MUTATES_H_
