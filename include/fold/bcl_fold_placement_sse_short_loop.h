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

#ifndef BCL_FOLD_PLACEMENT_SSE_SHORT_LOOP_H_
#define BCL_FOLD_PLACEMENT_SSE_SHORT_LOOP_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_fold_placement_interface.h"
#include "find/bcl_find_locator_criteria_interface.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_sh_ptr.h"

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PlacementSSEShortLoop
    //! @brief class places an SSE into a protein model that contains another SSE connected by a short loop
    //! @details this class first identifies SSEs in the given protein model that has short loops to the given SSE
    //! that is not yet added to the model and picks one of them as the neighbor. Then it picks a location to place
    //! the given SSE next to the picked neighbor in a way that would not violate the loop closure by making sure
    //! the ends of the SSEs that are connected by loop are close to each other
    //!
    //! @see @link example_fold_placement_sse_short_loop.cpp @endlink
    //! @author karakam
    //! @date Apr 9, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PlacementSSEShortLoop :
      public PlacementInterface< assemble::SSE, assemble::ProteinModel>
    {

    //////////
    // data //
    //////////

      //! minimum loop distance
      static const double s_MinLoopDistance;

      //! maxiumum loop distance
      static const double s_MaxLoopDistance;

      //! maximum loop distance per residue
      static const double s_MaxDistancePerResidue;

      //! maximum loop length for defining short loops
      size_t m_MaxShortLoopLength;

      //! probability that determines how frequently the SSE will be added to the top of the SSE
      double m_AddToTopProbability;

      //! max hinge angle to be used when adding to top
      double m_MaxHingeAngleAddingToTop;

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
      PlacementSSEShortLoop();

      //! @brief constructor from a maximum short loop length
      //! @param MAX_SHORT_LOOP_LENGTH maximum number of residues between SSEs to be classified as short loop
      //! @param ADD_TO_TOP_PROBABILITY probability that determines how frequently the SSE will be added to the top of the SSE
      //! @param MAX_HINGE_ANGLE_ADDING_TO_TOP max hinge angle to be used when adding to top
      PlacementSSEShortLoop
      (
        const size_t MAX_SHORT_LOOP_LENGTH,
        const double ADD_TO_TOP_PROBABILITY = 0.25, // 1:3 ratio
        const double MAX_HINGE_ANGLE_ADDING_TO_TOP = ( math::g_Pi / 3.0) // 60 degrees
      );

      //! @brief Clone function
      //! @return pointer to new PlacementSSEShortLoop
      PlacementSSEShortLoop *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const
      {
        static const std::string s_alias( "PlacementSSEShortLoop");
        return s_alias;
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer serializer;
        serializer.SetClassDescription( "Places an SSE in a protein model.");
        serializer.AddInitializer
        (
          "max_short_loop_length",
          "maximum sequence length of a loop to be considered short",
          io::Serialization::GetAgent( &m_MaxShortLoopLength)
        );
        serializer.AddInitializer
        (
          "add_to_top_probability",
          "probability to add SSE to the top of the SSE",
          io::Serialization::GetAgent( &m_AddToTopProbability),
          "0.25"
        );
        serializer.AddInitializer
        (
          "max_hinge_angle",
          "maximum hinge angle when adding to top",
          io::Serialization::GetAgent( &m_MaxHingeAngleAddingToTop),
          "1.0472"
        );

        return serializer;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief generate placement for given SSE at a random orientation wrt to a located SSE in PROTEIN_MODEL
      //! @param SELECTED_SSE SiPtr to SSE to be placed
      //! @param PROTEIN_MODEL to which the SSE is going to be added
      //! @return the transformationmatrix3d to place the SSE in the PROTEIN_MODEL
      storage::Pair< math::TransformationMatrix3D, bool> Place
      (
        const assemble::SSE &SELECTED_SSE,
        const assemble::ProteinModel &PROTEIN_MODEL
      ) const;

      //! @brief generate placement for given SSE at a random orientation wrt to NEIGHBOR_SSE
      //! @param SELECTED_SSE SiPtr to SSE to be placed
      //! @param NEIGHBOR_SSE neighbor SSE
      //! @return the transformationmatrix3d to place the SSE next to NEIGHBOR_SSE
      storage::Pair< math::TransformationMatrix3D, bool> Place
      (
        const assemble::SSE &SELECTED_SSE,
        const assemble::SSE &NEIGHBOR_SSE
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    }; // class PlacementSSEShortLoop

  } // namespace fold
} // namespace bcl

#endif //BCL_FOLD_PLACEMENT_SSE_SHORT_LOOP_H_
