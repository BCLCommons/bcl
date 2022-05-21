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

#ifndef BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_ADD_H_
#define BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_ADD_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_fold_placement_interface.h"
#include "find/bcl_find_pick_criteria_interface.h"
#include "math/bcl_math_mutate_result.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateProteinModelSSEAdd
    //! @brief Mutate class to be used for adding SSEs to protein models
    //! @details This class uses the provided SSE locator and SSE placement classes to select and add that SSE in a specific
    //! way to the provided protein model.
    //!
    //! @see @link example_fold_mutate_protein_model_sse_add.cpp @endlink
    //! @author karakam
    //! @date Apr 9, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateProteinModelSSEAdd :
      public math::MutateInterface< assemble::ProteinModel>
    {

    private:

    //////////
    // data //
    //////////

      //! Pick Interface to be used to pick the SSE to be added
      util::Implementation
      <
        find::PickCriteriaInterface
        <
          util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::DomainInterface
        >
      >
      m_SSEPoolPicker;

      //! Placement to be used to place the SSE in the model
      util::Implementation< PlacementInterface< assemble::SSE, assemble::ProteinModel> > m_Placement;

      //! scheme
      std::string m_Scheme;

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
      MutateProteinModelSSEAdd();

      //! @brief constructor from a PickCriteriaInterface, PlacementInterface and a scheme
      //! @param PICKER picker to be used for picking sses from the pool
      //! @param PLACEMENT placement to place the picked sse in the model
      //! @param SCHEME Scheme to be used
      MutateProteinModelSSEAdd
      (
        const find::PickCriteriaInterface
        <
          util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::DomainInterface
        > &PICKER,
        const PlacementInterface< assemble::SSE, assemble::ProteinModel> &PLACEMENT,
        const std::string &SCHEME = GetStaticClassName< MutateProteinModelSSEAdd>()
      );

      //! @brief constructor from a PickCriteriaInterface, PlacementInterface add a scheme
      //! @param SP_PICKER ShPtr to picker to be used for picking sses from the pool
      //! @param SP_PLACEMENT ShPtr to placement to place the picked sse in the model
      //! @param SCHEME Scheme to be used
      MutateProteinModelSSEAdd
      (
        const util::ShPtr
        <
          find::PickCriteriaInterface
          <
            util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::DomainInterface
          >
        > &SP_PICKER,
        const util::ShPtr< PlacementInterface< assemble::SSE, assemble::ProteinModel> > &SP_PLACEMENT,
        const std::string &SCHEME = GetStaticClassName< MutateProteinModelSSEAdd>()
      );

      //! @brief Clone function
      //! @return pointer to new MutateProteinModelSSEAdd
      MutateProteinModelSSEAdd *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this mutate
      //! @return the scheme for this mutate
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

      //! @brief set the placement object
      //! @param SP_PLACEMENT the object to use to decide where to place the SSE
      void SetPlacement( const util::ShPtr< PlacementInterface< assemble::SSE, assemble::ProteinModel> > &SP_PLACEMENT)
      {
        m_Placement = *SP_PLACEMENT;
      }

      //! @brief get a short name for this class
      //! @return a short name for this class
      virtual const std::string &GetAlias() const
      {
        static const std::string s_alias( "MutateProteinModelSSEAdd");
        return s_alias;
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      virtual io::Serializer GetSerializer() const
      {
        io::Serializer serializer;
        serializer.SetClassDescription( "Adds an SSE to a protein model.");
        serializer.AddInitializer
        (
          "picker", "selects, which SSE will be added",
          io::Serialization::GetAgent( &m_SSEPoolPicker)
        );
        serializer.AddInitializer
        (
          "placement",
          "places the SSE in the protein model",
          io::Serialization::GetAgent( &m_Placement)
        );

        return serializer;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
      //! @param PROTEIN_MODEL protein model interest
      //! @return MutateResult with ProteinModel after the mutate
      math::MutateResult< assemble::ProteinModel> operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    }; // MutateProteinModelSSEAdd

  } // namespace fold
} // namespace bcl

#endif //BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_ADD_H_
