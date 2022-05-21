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

#ifndef BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_MOVE_H_
#define BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_MOVE_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_fold_placement_interface.h"
#include "find/bcl_find_locator_criteria_interface.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_mutate_interface.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateProteinModelSSEMove
    //! @brief This mutate picks an SSE from given model and moves it using a placement
    //! @details First a LocatorInterface is used to locate an SSE from the given model, then a PlacementInterface is
    //! used to determine a new position and then this placement is applied to the located SSE
    //!
    //! @see @link example_fold_mutate_protein_model_sse_move.cpp @endlink
    //! @author karakam, woetzen
    //! @date Aug 5, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateProteinModelSSEMove :
      public math::MutateInterface< assemble::ProteinModel>
    {

    private:

    //////////
    // data //
    //////////

      //! locator that determines which SSE will be moved
      util::Implementation
      <
        find::LocatorCriteriaInterface
        <
          util::SiPtr< const assemble::SSE>, assemble::DomainInterface, assemble::DomainInterface
        >
      > m_Locator;

      //! placement that determines the new position for the SSE
      util::Implementation< PlacementInterface< assemble::SSE, assemble::ProteinModel> > m_Placement;

      //! scheme
      std::string m_Scheme;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutateProteinModelSSEMove();

      //! @brief constructor from a locator and placement and a scheme
      //! @param LOCATOR ShPtr to LocatorInterface to be used
      //! @param PLACEMENT ShPtr to PlacementInterface to be used
      //! @param SCHEME Scheme to be used
      MutateProteinModelSSEMove
      (
        const find::LocatorCriteriaInterface
        <
          util::SiPtr< const assemble::SSE>, assemble::DomainInterface, assemble::DomainInterface
        > &LOCATOR,
        const PlacementInterface< assemble::SSE, assemble::ProteinModel> &PLACEMENT,
        const std::string &SCHEME = GetStaticClassName< MutateProteinModelSSEMove>()
      );

      //! @brief constructor from ShPtrs to a locator and placement and a scheme
      //! @param SP_LOCATOR ShPtr to LocatorInterface to be used
      //! @param SP_PLACEMENT ShPtr to PlacementInterface to be used
      //! @param SCHEME Scheme to be used
      MutateProteinModelSSEMove
      (
        const util::ShPtr
        <
          find::LocatorCriteriaInterface
          <
            util::SiPtr< const assemble::SSE>, assemble::DomainInterface, assemble::DomainInterface
          >
        > &SP_LOCATOR,
        const util::ShPtr< PlacementInterface< assemble::SSE, assemble::ProteinModel> > &SP_PLACEMENT,
        const std::string &SCHEME = GetStaticClassName< MutateProteinModelSSEMove>()
      );

      //! @brief Clone function
      //! @return pointer to new MutateProteinModelSSEMove
      MutateProteinModelSSEMove *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get scheme
      //! @return scheme
      const std::string &GetScheme() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const
      {
        static const std::string s_alias( "MutateProteinModelSSEMove");
        return s_alias;
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer serializer;
        serializer.SetClassDescription( "Moves an SSE in a protein model.");
        serializer.AddInitializer
        (
          "locator",
          "locates an SSE in a protein model",
          io::Serialization::GetAgent( &m_Locator)
        );
        serializer.AddInitializer
        (
          "placement",
          "places an SSE in a protein model",
          io::Serialization::GetAgent( &m_Placement)
        );

        return serializer;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that moves an SSE
      //! @param PROTEIN_MODEL protein model interest
      //! @return MutateResult with ProteinModel after the mutate
      math::MutateResult< assemble::ProteinModel> operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    }; // class MutateProteinModelSSEMove

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_PROTEIN_MODEL_SSE_MOVE_H_
