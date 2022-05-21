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

#ifndef BCL_FOLD_PLACEMENT_STRAND_NEXT_TO_SHEET_H_
#define BCL_FOLD_PLACEMENT_STRAND_NEXT_TO_SHEET_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_fold_placement_interface.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PlacementStrandNextToSheet
    //! @brief class allows finding the placement for a strand when adding to an existing sheet
    //! @details This class picks one of the beta-sheets in the given ProteinModel, and tries to find a correct placement
    //! that would add the given strand to the picked sheet.
    //! If there are only sheets of size 1, meaning individual strand or strands, then it creates the placements that
    //! put the strand next to strand in the sheet with a pre-given twist angle and pre-given distance ( adjusted
    //! by the deviation)
    //! If it can find a Sheet of 2 or more strands, it picks one the edges of the sheet and tries to mimic the
    //! transformation from the reference strand ( one before the edge strand in the sheet) to the edge strand. This
    //! mimicking includes the translation as well as the twist angle.
    //! In both cases, this placement also applies a 180 degree flip among X, Y, or Z axis of the new strand.
    //!
    //! @see @link example_fold_placement_strand_next_to_sheet.cpp @endlink
    //! @author karakam
    //! @date Mar 2, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PlacementStrandNextToSheet :
      public PlacementInterface< assemble::SSE, assemble::ProteinModel>
    {

    private:

    //////////
    // data //
    //////////

      //! probability of applying a flip the strand around a random axis
      double m_FlipProbability;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a flip probability and a distance deviation
      //! @param FLIP_PROBABILITY probability of applying a flip the strand around a random axis
      PlacementStrandNextToSheet( const double FLIP_PROBABILITY = 0.5);

      //! @brief Clone function
      //! @return pointer to new PlacementStrandNextToSheet
      PlacementStrandNextToSheet *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return flip probability
      //! @return flip probability
      double GetFlipProbability()
      {
        return m_FlipProbability;
      }

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief find the placement for the given strand SELECTED_STRAND with respect to a sheet in PROTEIN_MODEL
      //! @param SELECTED_STRAND SiPtr to strand to be placed
      //! @param PROTEIN_MODEL to which the SELECTED_STRAND is going to be added
      //! @return the transformationmatrix3d to place the SELECTED_STRAND in the PROTEIN_MODEL
      storage::Pair< math::TransformationMatrix3D, bool> Place
      (
        const assemble::SSE &SELECTED_STRAND,
        const assemble::ProteinModel &PROTEIN_MODEL
      ) const;

      //! @brief find the placement for the given strand SELECTED_STRAND with respect to a sheet in SHEET
      //! @param SELECTED_STRAND SiPtr to strand to be placed
      //! @param SHEET to which the SELECTED_STRAND is going to be added
      //! @return the transformationmatrix3d to place the SELECTED_STRAND in the SHEET
      storage::Pair< math::TransformationMatrix3D, bool> Place
      (
        const assemble::SSE &SELECTED_STRAND,
        const assemble::Domain &SHEET
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    }; // class PlacementStrandNextToSheet

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_PLACEMENT_STRAND_NEXT_TO_SHEET_H_ 
