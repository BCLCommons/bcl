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

#ifndef BCL_SCORE_AA_NEIGHBORHOOD_INTERFACE_H_
#define BCL_SCORE_AA_NEIGHBORHOOD_INTERFACE_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_binary_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AANeighborhoodInterface
    //! @brief score the neighborhood of an amino acid using its AANeighborlist
    //! @details The Neighborlist is generated using the information that is implemented with the distance cutoff and
    //!          minimal sequence separation.
    //!
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date Jun 26, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AANeighborhoodInterface :
      public math::BinaryFunctionInterfaceSerializable< assemble::AANeighborList, util::SiPtr< const biol::Membrane>, double>
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new AANeighborhoodInterface
      virtual AANeighborhoodInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief access to the minimal sequence separation
      //! @return minimal sequence separation for neighbors of that exposure score between amino acids in the same chain
      virtual size_t GetMinimalSequenceSeparation() const = 0;

      //! @brief access to the distance cutoff
      //! @return distance cutoff above which the neighbor does not have influence on the score anymore
      virtual double GetDistanceCutoff() const = 0;

    }; // class AANeighborhoodInterface

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_AA_NEIGHBORHOOD_INTERFACE_H_ 
