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

#ifndef BCL_SCORE_AA_PAIR_DISTANCE_INTERFACE_H_
#define BCL_SCORE_AA_PAIR_DISTANCE_INTERFACE_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_binary_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AAPairDistanceInterface
    //! @brief interface for scoring an aa pair and the distance between them
    //! @details scoring aa pairs is usually just of function of the distance. This interface exposes a function, where
    //!          one can just pass the two amino acids and their distance enabling the use of an aa pair from a neighbor
    //!          list without the necessity to recalculate the distance.
    //!
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date May 14, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AAPairDistanceInterface :
      public math::BinaryFunctionInterfaceSerializable< biol::AABase, biol::AABase, double>
    {

    public:

      //! @brief virtual copy constructor
      //! @return pointer to a new AANeighborhoodDistances object copied from this one
      virtual AAPairDistanceInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief access to the minimal sequence separation
      //! @return minimal sequence separation in sequence distance
      virtual size_t GetMinimalSequenceSeparation() const = 0;

      //! @brief access to the maximal distance that is considered
      //! @return the maximal distance between two amino acids that is scored
      virtual double GetDistanceCutoff() const = 0;

      //! @brief are amino acids of different chains considered
      //! @return return true, if amino acid pairs of different chains are considered
      virtual bool GetConsiderDifferentChain() const = 0;

    ///////////////
    // operators //
    ///////////////

      //! @brief calculate amino acid pairing potential for given amino acid pair
      //! @param AMINO_ACID_A first amino acid of interest
      //! @param AMINO_ACID_B second amino acid of interest
      //! @param DISTANCE distance between the amino acid pair
      //! @return amino acid pairing potential for given amino acid pair
      virtual double operator()
      (
        const biol::AABase &AMINO_ACID_A,
        const biol::AABase &AMINO_ACID_B,
        const double DISTANCE
      ) const = 0;

    }; // class AAPairDistanceInterface

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_AA_PAIR_DISTANCE_INTERFACE_H_
