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

#ifndef BCL_CHEMISTRY_CHIRALITY_H_
#define BCL_CHEMISTRY_CHIRALITY_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_wrapper_enum.h"

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @file bcl_chemistry_chirality.h
  //! @brief enumerates Chiralities of atoms
  //!
  //! @remarks example unnecessary
  //! @author mendenjl
  //! @date Feb 24, 2012
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace chemistry
  {

    enum Chirality
    {
      e_NonChiral,        //!< No chirality
      e_RChirality,       //!< rectus   - lowest priority substituent oriented away from viewer, remaining 3 substituents decrease in clockwise direction
      e_SChirality,       //!< sinister - lowest priority substituent oriented away from viewer, remaining 3 substituents decrease in counter-clockwise direction
      e_Chiral,           //!< Chiral, but unknown chirality
      e_UnknownChirality, //!< Any chirality or non-chiral
      e_UnknownRingChirality, //!< Potentially-chiral/entantiomeric/non-chiral center on a ring like in 1,3 dichlorohexane
      e_CisRingChirality,   //!< Only applies to saturated rings with exactly two asymmetric atoms, with substituents on same sides
      e_TransRingChirality, //!< Only applies to saturated rings with exactly two asymmetric atoms, with substituents on opposite sides
      s_NumberChiralities
    };

    //! @brief Chirality as string
    //! @param CHIRALITY the chirality desired
    //! @return the chirality as string
    BCL_API const std::string &GetChiralityName( const Chirality &CHIRALITY);

    //! ChiralityEnum simplifies the usage of the Chirality enum of this class
    typedef util::WrapperEnum< Chirality, &GetChiralityName, s_NumberChiralities> ChiralityEnum;

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_CHIRALITY_H_
