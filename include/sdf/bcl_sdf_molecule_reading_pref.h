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

#ifndef BCL_SDF_MOLECULE_READING_PREF_H_
#define BCL_SDF_MOLECULE_READING_PREF_H_

// include the namespace header
#include "bcl_sdf.h"

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_wrapper_enum.h"

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @file bcl_sdf_molecule_reading_pref.h
  //! @brief enumerates ways of handling hydrogens when loading fragments and other molecule reading preferences
  //!
  //! @remarks example unnecessary
  //! @author mendenjl
  //! @date May 16, 2016
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace sdf
  {

    enum HydrogenHandlingPref
    {
      e_SaturatePartial, //!< saturates some atoms with H
      e_Saturate,        //!< saturates molecule with H ignoring hydrogen count property
      e_Maintain,        //!< Neither add nor remove hydrogens
      e_Remove,          //!< Remove any hydrogens that are loaded in
      s_NumberHydrogenHandlingPrefs
    };

    //! @brief MoleculeReadingPref as string
    //! @param PREF the MoleculeReadingPref desired
    //! @return the MoleculeReadingPref as string
    BCL_API const std::string &GetHydrogensPrefName( const HydrogenHandlingPref &PREF);

    //! MoleculeReadingPrefEnum simplifies the usage of the MoleculeReadingPref enum of this class
    typedef util::WrapperEnum< HydrogenHandlingPref, &GetHydrogensPrefName, s_NumberHydrogenHandlingPrefs>
      HydrogenHandlingPrefEnum;

    //! @brief get access to a global flag defining whether hydrogens should be added
    //! @return access to a global flag defining whether hydrogens should be added
    BCL_API const util::ShPtr< command::FlagInterface> &GetAddHydrogensFlag();

    //! @brief get access to a global flag defining whether hydrogens should be removed
    //! @return access to a global flag defining whether hydrogens should be removed
    BCL_API const util::ShPtr< command::FlagInterface> &GetRemoveHydrogensFlag();

    //! @brief get access to the desired hydrogen handling based on the add/remove flags
    //! @return access to the desired hydrogen handling based on the add/remove flags
    BCL_API HydrogenHandlingPref GetCommandLineHydrogensPref();

    //! Neutralization metho
    enum NeutralizationPref
    {
      e_None,      //!< No effort will be made to neutralize the molecule
      e_BondOrder, //!< Bond order changes between formally charged atoms that result in neutral species will be employed
      e_pH,        //!< Neutralize lone, charged species by adding or removing H from them. Bond order changes are not employed
      e_BondOrderAndpH, //!< Allows for changes in bond order and changes in pH. Allows for changes in hydrogenation on
                        //! formerly neutral species if adjusting bond order afterwards would fix the charges
      e_BondOrderAndpHAromaticityLossOk, //!< Like BondOrderAndpH but allows aromaticity loss
      s_NumberNeutralizationPrefs,
      e_CmdLine //!< Special enum used as a default; means to use whatever was provided on the command line
    };

    //! @brief NeutralizationPref as string
    //! @param PREF the NeutralizationPref desired
    //! @return the NeutralizationPref as string
    BCL_API const std::string &GetNeutralizationPrefName( const NeutralizationPref &PREF);

    //! MoleculeReadingPrefEnum simplifies the usage of the MoleculeReadingPref enum of this class
    typedef util::WrapperEnum< NeutralizationPref, &GetNeutralizationPrefName, s_NumberNeutralizationPrefs>
      NeutralizationPrefEnum;

    //! @brief get access to a global flag defining whether charges should be neutralized wherever possible
    //! @return access to a global flag defining whether charges should be neutralized wherever possible
    BCL_API const util::ShPtr< command::FlagInterface> &GetNeutralizeChargesFlag();

    //! @brief get whether to try to neutralize charges on read in molecules
    //! @return whether to try to neutralize charges on read in molecules
    BCL_API NeutralizationPrefEnum GetCommandLineNeutralizationPref();

    //! @brief determine whether the preference allows alteration of the bond order
    //! @return true if the neutralization pref allows alteration of the bond order
    BCL_API bool GetNeutralizationPrefAllowsBondOrderChange( const NeutralizationPref &PREF);

    //! @brief determine whether the preference allows alteration of the pH
    //! @return true if the neutralization pref allows alteration of the pH
    BCL_API bool GetNeutralizationPrefAllowspHChange( const NeutralizationPref &PREF);

    //! @brief determine whether the preference allows alteration of the aromaticity
    //! @return true if the neutralization pref allows alteration of the aromaticity
    BCL_API bool GetNeutralizationPrefAllowsAromaticityChange( const NeutralizationPref &PREF);

    //! @brief get access to a global flag defining whether aromatic bonds are written to MDL instead of kekulized bonds
    //! @return access to a global flag defining aromatic bonds explicitly
    BCL_API const util::ShPtr< command::FlagInterface> &GetExplicitAromaticityFlag();

    //! @brief whether to explicitly write aromatic bonds instead of default kekulized bonds
    //! @return true if writing aromatic bonds explicitly
    BCL_API bool GetExplicitAromaticBondsPref();

    //! @brief add molecule reading and writing preferences to the command line flag
    //! @param CMD command to add the molecule reading  preference flags to
    BCL_API void AddMoleculeIOPrefFlags( command::Command &CMD);

    //! @brief add molecule reading preferences to the command line flag
    //! @param CMD command to add the molecule reading  preference flags to
    BCL_API void AddMoleculeReadingPrefFlags( command::Command &CMD);

    //! @brief add molecule writing preferences to the command line flag
    //! @param CMD command to add the molecule reading  preference flags to
    BCL_API void AddMoleculeWritingPrefFlags( command::Command &CMD);

  } // namespace sdf
} // namespace bcl

#endif // BCL_SDF_MOLECULE_READING_PREF_H_
