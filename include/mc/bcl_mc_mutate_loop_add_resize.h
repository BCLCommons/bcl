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

#ifndef BCL_MC_MUTATE_LOOP_ADD_RESIZE_H_
#define BCL_MC_MUTATE_LOOP_ADD_RESIZE_H_

// include the namespace header
#include "bcl_mc.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_locator_loop.h"
#include "fold/bcl_fold_loop_library.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateLoopAddResize
    //! @brief Adds a loop to the given protein model
    //! @details Uses conformational hashing to add a loop to the given protein model.
    //!
    //! @see @link example_mc_mutate_loop_add_resize.cpp @endlink
    //! @author fischea
    //! @date Mar 8, 2016
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API MutateLoopAddResize :
      public math::MutateInterface< assemble::ProteinModel>
    {

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    private:

      //! locator for unbuilt loops in a given protein model
      fold::LocatorLoop m_LoopLocator;

      //! minimum lengths of helices and strands after resize
      storage::VectorND< 2, size_t> m_MinSizes;

      //! path to the loop template library
      std::string m_LoopLibraryFilename;

      //! library containing the loop templates
      util::ShPtr< fold::LoopLibrary> m_LoopLibrary;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief construct from members
      //! @param LOOP_LIBRARY_FILENAME path of the loop template library
      //! @param MIN_SIZES minimum lengths of helices and strands after resize
      MutateLoopAddResize
      (
        const std::string &LOOP_LIBRARY_FILENAME = GetDefaultLibraryPath(),
        const storage::VectorND< 2, size_t> MIN_SIZES = GetDefaultMinSizes()
      );

      //! @brief copy constructor
      //! @return pointer to a new MutateLoopAddResize
      MutateLoopAddResize *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

      //! @brief returns the default path of the loop template library
      //! @return the default path of the loop template library
      static const std::string &GetDefaultLibraryPath();

      //! @brief returns the default minimum lengths of helices and strands after resize
      //! @return the default minimum lengths of helices and strands after resize
      static storage::VectorND< 2, size_t> GetDefaultMinSizes();

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read of this class
      //! @param ERROR_STREAM stream with which to write errors
      //! @return return code indicating success or failure
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    ///////////////
    // operators //
    ///////////////

      //! @brief applies a mutation to the given protein model
      //! @detail this mutate randomly selects a missing loop in a protein model. the anchor SSEs of the selected loop are
      //! randomly resized and the resulting loop is constructed using a suitable conformation from a template library.
      //! @param MODEL protein model to which to apply the mutation
      //! @return result mutating the given protein model
      math::MutateResult< assemble::ProteinModel> operator()( const assemble::ProteinModel &MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief returns the anchor SSEs of a given loop
      //! @param MODEL protein model in which to find the anchor SSEs
      //! @param LOOP loop for which to return the anchor SSEs
      //! @return anchor SSEs of the given loop
      static storage::VectorND< 2, util::ShPtr< assemble::SSE> > GetAnchors
      (
        const assemble::ProteinModel &MODEL,
        const fold::LoopParameters &LOOP
      );

      //! @brief shrinks the given SSE by the given length
      //! @param SSE SSE to be shrunk
      //! @param LENGTH length by which the SSE will get shrunk
      //! @param DIRECTION terminus at which to shrink the SSE
      //! @return the shrunk SSE
      static util::ShPtr< assemble::SSE> ShrinkSSE
      (
        const assemble::SSE &SSE, size_t LENGTH, const biol::AASequenceFlexibility::SequenceDirection &DIRECTION
      );

      //! @brief randomly resizes the given SSEs
      //! @param SSES the SSEs that shall be resized
      //! @return resized SSEs
      static storage::VectorND< 2, util::ShPtr< assemble::SSE> > ResizeSSEs
      (
        const storage::VectorND< 2, util::ShPtr< assemble::SSE> > &SSES
      );

    }; // class MutateLoopAddResize

  } // namespace mc
} // namespace bcl

#endif // BCL_MC_MUTATE_LOOP_ADD_RESIZE_H_
