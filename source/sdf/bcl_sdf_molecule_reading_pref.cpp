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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "sdf/bcl_sdf_molecule_reading_pref.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace sdf
  {

    //! @brief MoleculeReadingPref as string
    //! @param PREF the MoleculeReadingPref desired
    //! @return the MoleculeReadingPref as string
    const std::string &GetPrefName( const HydrogenHandlingPref &PREF)
    {
      static const std::string s_names[] =
      {
        "SelectiveSaturate",
        "Saturate",
        "Maintain",
        "Remove",
        GetStaticClassName< HydrogenHandlingPref>()
      };

      return s_names[ PREF];
    }

    //! @brief get access to a global flag defining whether hydrogens should be added
    //! @return access to a global flag defining whether hydrogens should be added
    const util::ShPtr< command::FlagInterface> &GetAddHydrogensFlag()
    {
      static util::ShPtr< command::FlagInterface> s_h_add
        (
          new command::FlagStatic( "add_h", "add hydrogens to molecules when loaded")
        );

      return s_h_add;
    }

    //! @brief get access to a global flag defining whether hydrogens should be removed
    //! @return access to a global flag defining whether hydrogens should be removed
    const util::ShPtr< command::FlagInterface> &GetRemoveHydrogensFlag()
    {
      static util::ShPtr< command::FlagInterface> s_h_del
        (
          new command::FlagStatic( "remove_h", "remove hydrogens from molecules when loaded")
        );

      return s_h_del;
    }

    //! @brief get access to a global flag defining whether charges should be neutralized wherever possible
    //! @return access to a global flag defining whether charges should be neutralized wherever possible
    const util::ShPtr< command::FlagInterface> &GetNeutralizeChargesFlag()
    {
      static util::ShPtr< command::FlagInterface> s_neutralize
      (
        new command::FlagStatic
        (
          "neutralize",
          "neutralize charges; by default, if the flag is specified but no neutralization type is given, "
          "BondOrderAndpH will be used, otherwise, no neutralization is used. All neutralization algorithms preserve "
          "aromaticity except BondOrderAndpHAromaticityLossOk",
          command::Parameter
          (
            "method",
            "method used to neutralize charged atoms in the molecule",
            command::ParameterCheckSerializable( NeutralizationPrefEnum()),
            "BondOrderAndpH"
          )
        )
      );

      return s_neutralize;
    }

    //! @brief get access to a global flag defining whether hydrogens should be added
    //! @return access to a global flag defining whether hydrogens should be added
    const util::ShPtr< command::FlagInterface> &GetExplicitAromaticityFlag()
    {
      static util::ShPtr< command::FlagInterface> s_explicit_aro
        (
          new command::FlagStatic
          (
            "explicit_aromaticity",
            "write MDL bonds with aromatic bonds specified explicitly (as 4); "
            "alternatively, the default behavior during MDL writing is to "
            "kekulize aromatic rings (write alternating single (1) / double (2) "
            "bonds)"
            )
        );

      return s_explicit_aro;
    }

    //! @brief add molecule reading and writing preferences to the command line flag
    //! @param CMD command to add the molecule reading preference flags to
    void AddMoleculeIOPrefFlags( command::Command &CMD)
    {
      CMD.AddFlag( GetAddHydrogensFlag());
      CMD.AddFlag( GetRemoveHydrogensFlag());
      CMD.AddFlag( GetNeutralizeChargesFlag());
      CMD.AddFlag( GetExplicitAromaticityFlag());
    }

    //! @brief add molecule reading  preferences to the command line flag
    //! @param CMD command to add the molecule reading preference flags to
    void AddMoleculeReadingPrefFlags( command::Command &CMD)
    {
      CMD.AddFlag( GetAddHydrogensFlag());
      CMD.AddFlag( GetRemoveHydrogensFlag());
      CMD.AddFlag( GetNeutralizeChargesFlag());
    }

    //! @brief add molecule writing  preferences to the command line flag
    //! @param CMD command to add the molecule reading preference flags to
    void AddMoleculeWritingPrefFlags( command::Command &CMD)
    {
      CMD.AddFlag( GetExplicitAromaticityFlag());
    }

    //! @brief get access to the desired hydrogen handling based on the add/remove flags
    //! @return access to the desired hydrogen handling based on the add/remove flags
    HydrogenHandlingPref GetCommandLineHydrogensPref()
    {
      const bool add_h( GetAddHydrogensFlag()->GetFlag());
      const bool remove_h( GetRemoveHydrogensFlag()->GetFlag());
      BCL_Assert( !add_h || !remove_h, "Cannot add and remove hydrogens; choose one or the other");
      return add_h ? e_Saturate : remove_h ? e_Remove : e_Maintain;
    }

    //! @brief NeutralizationPref as string
    //! @param PREF the NeutralizationPref desired
    //! @return the NeutralizationPref as string
    const std::string &GetNeutralizationPrefName( const NeutralizationPref &PREF)
    {
      static const std::string s_names[] =
      {
        "None",
        "BondOrder",
        "pH",
        "BondOrderAndpH",
        "BondOrderAndpHAromaticityLossOk",
        GetStaticClassName< NeutralizationPref>()
      };
      return s_names[ PREF];
    }

    //! @brief determine whether the preference allows alteration of the bond order
    //! @return true if the neutralization pref allows alteration of the bond order
    bool GetNeutralizationPrefAllowsBondOrderChange( const NeutralizationPref &PREF)
    {
      static const bool s_bond_order_chng_allowed[] = { false, true, false, true, true, false};
      return s_bond_order_chng_allowed[ size_t( PREF)];
    }

    //! @brief determine whether the preference allows alteration of the pH
    //! @return true if the neutralization pref allows alteration of the pH
    bool GetNeutralizationPrefAllowspHChange( const NeutralizationPref &PREF)
    {
      static const bool s_chng_allowed[] = { false, false, true, true, true, false};
      return s_chng_allowed[ size_t( PREF)];
    }

    //! @brief determine whether the preference allows alteration of the pH
    //! @return true if the neutralization pref allows alteration of the pH
    bool GetNeutralizationPrefAllowsAromaticityChange( const NeutralizationPref &PREF)
    {
      static const bool s_chng_allowed[] = { false, false, false, false, true, false};
      return s_chng_allowed[ size_t( PREF)];
    }

    //! @brief whether to explicitly write aromatic bonds instead of default kekulized bonds
    //! @return true if writing aromatic bonds explicitly
    BCL_API bool GetExplicitAromaticBondsPref()
    {
      return GetExplicitAromaticityFlag()->GetFlag();
    }

    //! @brief get whether to try to neutralize charges on read in molecules
    //! @return whether to try to neutralize charges on read in molecules
    NeutralizationPrefEnum GetCommandLineNeutralizationPref()
    {
      if( !GetNeutralizeChargesFlag()->GetFlag())
      {
        return e_None;
      }
      NeutralizationPrefEnum type( GetNeutralizeChargesFlag()->GetFirstParameter()->GetValue());
      return type;
    }

  } // namespace sdf
} // namespace bcl
