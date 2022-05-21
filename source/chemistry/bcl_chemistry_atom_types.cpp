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
#include "chemistry/bcl_chemistry_atom_types.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_bond_lengths.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_enumerate.hpp"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    //! @brief construct all AtomTypes
    AtomTypes::AtomTypes() :
      m_MaxGasteigerAtomTypeIndex( 0),
      // neutral atoms
      H_S
      (
        AddEnum
        (
          "H_S",
          AtomTypeData
          (
            GetElementTypes().e_Hydrogen,               // Element type
            GetHybridOrbitalTypes().e_Unhybridized,     // Hybridization
            0,                                          // # hybridized sigma orbitals in bonds
            0,                                          // # non-binding hybrid orbitals
            storage::Set< AtomicOrbitalTypesEnum>( e_S),    // pi orbitals in bonds
            storage::Set< AtomicOrbitalTypesEnum>(),        // non-binding atomic orbitals
            13.60,
            0.75,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            0.387
          )
        )
      ),
      Li_S
      (
        AddEnum
        (
          "Li_S",
          AtomTypeData
          (
            GetElementTypes().e_Lithium,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_S),
            storage::Set< AtomicOrbitalTypesEnum>(),
            5.39,
            0.82,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Li_P
      (
        AddEnum
        (
          "Li_P",
          AtomTypeData
          (
            GetElementTypes().e_Lithium,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Px),
            storage::Set< AtomicOrbitalTypesEnum>(),
            3.54,
            0.56,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Be_SP
      (
        AddEnum
        (
          "Be_SP",
          AtomTypeData
          (
            GetElementTypes().e_Beryllium,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_S, e_Px),
            storage::Set< AtomicOrbitalTypesEnum>(),
            9.92,
            3.18,
            5.96,
            0.11,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Be_PP
      (
        AddEnum
        (
          "Be_PP",
          AtomTypeData
          (
            GetElementTypes().e_Beryllium,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py),
            storage::Set< AtomicOrbitalTypesEnum>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            6.11,
            0.76,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Be_DiDi
      (
        AddEnum
        (
          "Be_DiDi",
          AtomTypeData
          (
            GetElementTypes().e_Beryllium,
            GetHybridOrbitalTypes().e_SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            8.58,
            0.99,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Be_DiPi
      (
        AddEnum
        (
          "Be_DiPi",
          AtomTypeData
          (
            GetElementTypes().e_Beryllium,
            GetHybridOrbitalTypes().e_SP,
            1,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Py),
            storage::Set< AtomicOrbitalTypesEnum>(),
            8.02,
            0.92,
            6.04,
            0.43,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Be_TrTr
      (
        AddEnum
        (
          "Be_TrTr",
          AtomTypeData
          (
            GetElementTypes().e_Beryllium,
            GetHybridOrbitalTypes().e_SP2,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            7.61,
            0.59,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Be_TrPi
      (
        AddEnum
        (
          "Be_TrPi",
          AtomTypeData
          (
            GetElementTypes().e_Beryllium,
            GetHybridOrbitalTypes().e_SP2,
            1,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            7.38,
            0.63,
            6.06,
            0.54,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Be_TeTe
      (
        AddEnum
        (
          "Be_TeTe",
          AtomTypeData
          (
            GetElementTypes().e_Beryllium,
            GetHybridOrbitalTypes().e_SP3,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            7.18,
            0.51,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_SPP
      (
        AddEnum
        (
          "B_SPP",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_S, e_Px, e_Py),
            storage::Set< AtomicOrbitalTypesEnum>(),
            14.91,
            5.70,
            8.42,
            0.32,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_PPP
      (
        AddEnum
        (
          "B_PPP",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            8.40,
            3.46,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_DiDiPi
      (
        AddEnum
        (
          "B_DiDiPi",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().e_SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            12.55,
            2.12,
            8.23,
            0.44,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_DiPiPi
      (
        AddEnum
        (
          "B_DiPiPi",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().e_SP,
            1,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            11.66,
            2.56,
            8.41,
            1.89,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_TrTrTr
      (
        AddEnum
        (
          "B_TrTrTr",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().e_SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            11.29,
            1.38,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_TrTrPi
      (
        AddEnum
        (
          "B_TrTrPi",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().e_SP2,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            10.97,
            1.87,
            8.33,
            1.42,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_TeTeTe
      (
        AddEnum
        (
          "B_TeTeTe",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().e_SP3,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            10.43,
            1.53,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      C_SPPP
      (
        AddEnum
        (
          "C_SPPP",
          AtomTypeData
          (
            GetElementTypes().e_Carbon,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_S, e_Px, e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            21.01,
            8.91,
            11.27,
            0.34,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      C_DiDiPiPi
      (
        AddEnum
        (
          "C_DiDiPiPi",
          AtomTypeData
          (
            GetElementTypes().e_Carbon,
            GetHybridOrbitalTypes().e_SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            17.42,
            3.34,
            11.19,
            0.1,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            1.283
          )
        )
      ),
      C_TrTrTrPi
      (
        AddEnum
        (
          "C_TrTrTrPi",
          AtomTypeData
          (
            GetElementTypes().e_Carbon,
            GetHybridOrbitalTypes().e_SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            15.62,
            1.95,
            11.16,
            0.03,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            1.352 // = 1.352 if connected to >0 H's, 1.896 otherwise.  Assume for now that at least one bond is to an H
          )
        )
      ),
      C_TeTeTeTe
      (
        AddEnum
        (
          "C_TeTeTeTe",
          AtomTypeData
          (
            GetElementTypes().e_Carbon,
            GetHybridOrbitalTypes().e_SP3,
            4,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            14.61,
            1.34,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            1.061
          )
        )
      ),
      N_S2PPP
      (
        AddEnum
        (
          "N_S2PPP",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( e_S),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            13.94,
            0.84,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      N_SP2PP
      (
        AddEnum
        (
          "N_SP2PP",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_S, e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( e_Px),
            26.92,
            14.05,
            14.42,
            2.54,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      N_Di2DiPiPi
      (
        AddEnum
        (
          "N_Di2DiPiPi",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().e_SP,
            1,
            1,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            23.91,
            7.45,
            14.18,
            1.66,
            37.024,
            17.254,
            0.956
          )
        )
      ),
      N_DiDiPi2Pi
      (
        AddEnum
        (
          "N_DiDiPi2Pi",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().e_SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( e_Px),
            22.10,
            6.84,
            14.11,
            2.14,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            1.012 // N_Di2DiPiPi * N_TrTrTrPi2 / N_Tr2TrTrPi
          )
        )
      ),
      N_Tr2TrTrPi
      (
        AddEnum
        (
          "N_Tr2TrTrPi",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().e_SP2,
            2,
            1,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            20.60,
            5.14,
            14.12,
            1.78,
            34.645,
            15.107,
            1.030
          )
        )
      ),
      N_TrTrTrPi2
      (
        AddEnum
        (
          "N_TrTrTrPi2",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().e_SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>( e_Px),
            19.72,
            4.92,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            1.090
          )
        )
      ),
      N_Te2TeTeTe
      (
        AddEnum
        (
          "N_Te2TeTeTe",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().e_SP3,
            3,
            1,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            18.93,
            4.15,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            33.313,
            14.153,
            0.964
          )
        )
      ),
      O_S2P2PP
      (
        AddEnum
        (
          "O_S2P2PP",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py),
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_S, e_Pz),
            17.28,
            2.01,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      O_SP2P2P
      (
        AddEnum
        (
          "O_SP2P2P",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_S, e_Px),
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            36.07,
            18.44,
            18.53,
            3.40,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      O_Di2Di2PiPi
      (
        AddEnum
        (
          "O_Di2Di2PiPi",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().e_SP,
            0,
            2,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            17.28,
            2.01,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      O_Di2DiPi2Pi
      (
        AddEnum
        (
          "O_Di2DiPi2Pi",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().e_SP,
            1,
            1,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( e_Py),
            30.17,
            10.23,
            17.91,
            2.71,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            0.569 // similar to Tr2Tr2TrPi
          )
        )
      ),
      O_DiDiPi2Pi2
      (
        AddEnum
        (
          "O_DiDiPi2Pi2",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().e_SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            28.71,
            9.51,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            0.637 // similar to Te2Te2TeTe
          )
        )
      ),
      O_Tr2Tr2TrPi
      (
        AddEnum
        (
          "O_Tr2Tr2TrPi",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().e_SP2,
            1,
            2,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            26.65,
            7.49,
            17.70,
            2.47,
            42.534,
            20.154,
            0.569
          )
        )
      ),
      O_Tr2TrTrPi2
      (
        AddEnum
        (
          "O_Tr2TrTrPi2",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().e_SP2,
            2,
            1,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            26.14,
            7.32,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            0.274
          )
        )
      ),
      O_Te2Te2TeTe
      (
        AddEnum
        (
          "O_Te2Te2TeTe",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().e_SP3,
            2,
            2,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            24.39,
            6.11,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            40.358,
            18.708,
            0.637
          )
        )
      ),
      F_S2P2P2P
      (
        AddEnum
        (
          "F_S2P2P2P",
          AtomTypeData
          (
            GetElementTypes().e_Fluorine,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_S, e_Px, e_Py),
            20.86,
            3.50,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            26.368,
            13.378,
            0.296
          )
        )
      ),
      F_SP2P2P2
      (
        AddEnum
        (
          "F_SP2P2P2",
          AtomTypeData
          (
            GetElementTypes().e_Fluorine,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_S),
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py, e_Pz),
            38.24,
            24.37,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            26.368,
            13.378,
            0.296
          )
        )
      ),
      Na_S
      (
        AddEnum
        (
          "Na_S",
          AtomTypeData
          (
            GetElementTypes().e_Sodium,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_S),
            storage::Set< AtomicOrbitalTypesEnum>(),
            5.14,
            0.47,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Na_P
      (
        AddEnum
        (
          "Na_P",
          AtomTypeData
          (
            GetElementTypes().e_Sodium,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Px),
            storage::Set< AtomicOrbitalTypesEnum>(),
            3.04,
            0.09,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Mg_SP
      (
        AddEnum
        (
          "Mg_SP",
          AtomTypeData
          (
            GetElementTypes().e_Magnesium,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_S, e_Px),
            storage::Set< AtomicOrbitalTypesEnum>(),
            8.95,
            2.80,
            4.52,
            0.06,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Mg_PP
      (
        AddEnum
        (
          "Mg_PP",
          AtomTypeData
          (
            GetElementTypes().e_Magnesium,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py),
            storage::Set< AtomicOrbitalTypesEnum>(),
            5.65,
            0.01,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Mg_DiDi
      (
        AddEnum
        (
          "Mg_DiDi",
          AtomTypeData
          (
            GetElementTypes().e_Magnesium,
            GetHybridOrbitalTypes().e_SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            7.10,
            1.08,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Mg_DiPi
      (
        AddEnum
        (
          "Mg_DiPi",
          AtomTypeData
          (
            GetElementTypes().e_Magnesium,
            GetHybridOrbitalTypes().e_SP,
            1,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            7.30,
            0.78,
            5.09,
            0.03,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Mg_TrTr
      (
        AddEnum
        (
          "Mg_TrTr",
          AtomTypeData
          (
            GetElementTypes().e_Magnesium,
            GetHybridOrbitalTypes().e_SP2,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            6.54,
            0.52,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Mg_TrPi
      (
        AddEnum
        (
          "Mg_TrPi",
          AtomTypeData
          (
            GetElementTypes().e_Magnesium,
            GetHybridOrbitalTypes().e_SP2,
            1,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            6.75,
            0.38,
            5.27,
            0.02,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Mg_TeTe
      (
        AddEnum
        (
          "Mg_TeTe",
          AtomTypeData
          (
            GetElementTypes().e_Magnesium,
            GetHybridOrbitalTypes().e_SP3,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            6.28,
            0.32,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Al_SPP
      (
        AddEnum
        (
          "Al_SPP",
          AtomTypeData
          (
            GetElementTypes().e_Aluminum,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_S, e_Px, e_Py),
            storage::Set< AtomicOrbitalTypesEnum>(),
            12.27,
            4.92,
            6.47,
            1.37,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Al_PPP
      (
        AddEnum
        (
          "Al_PPP",
          AtomTypeData
          (
            GetElementTypes().e_Aluminum,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            6.50,
            4.89,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Al_DiDiPi
      (
        AddEnum
        (
          "Al_DiDiPi",
          AtomTypeData
          (
            GetElementTypes().e_Aluminum,
            GetHybridOrbitalTypes().e_SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            9.91,
            2.61,
            6.36,
            1.45,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Al_DiPiPi
      (
        AddEnum
        (
          "Al_DiPiPi",
          AtomTypeData
          (
            GetElementTypes().e_Aluminum,
            GetHybridOrbitalTypes().e_SP,
            1,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            9.39,
            3.66,
            6.49,
            3.13,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Al_TrTrTr
      (
        AddEnum
        (
          "Al_TrTrTr",
          AtomTypeData
          (
            GetElementTypes().e_Aluminum,
            GetHybridOrbitalTypes().e_SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            8.83,
            2.11,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Al_TrTrPi
      (
        AddEnum
        (
          "Al_TrTrPi",
          AtomTypeData
          (
            GetElementTypes().e_Aluminum,
            GetHybridOrbitalTypes().e_SP2,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            8.65,
            2.94,
            6.43,
            2.58,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Al_TeTeTe
      (
        AddEnum
        (
          "Al_TeTeTe",
          AtomTypeData
          (
            GetElementTypes().e_Aluminum,
            GetHybridOrbitalTypes().e_SP3,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            8.17,
            2.58,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Si_SPPP
      (
        AddEnum
        (
          "Si_SPPP",
          AtomTypeData
          (
            GetElementTypes().e_Silicon,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_S, e_Px, e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            17.31,
            6.94,
            9.19,
            2.82,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Si_DiDiPiPi
      (
        AddEnum
        (
          "Si_DiDiPiPi",
          AtomTypeData
          (
            GetElementTypes().e_Silicon,
            GetHybridOrbitalTypes().e_SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            14.06,
            4.07,
            9.18,
            2.20,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Si_TrTrTrPi
      (
        AddEnum
        (
          "Si_TrTrTrPi",
          AtomTypeData
          (
            GetElementTypes().e_Silicon,
            GetHybridOrbitalTypes().e_SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            12.61,
            3.20,
            9.17,
            2.00,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Si_TeTeTeTe
      (
        AddEnum
        (
          "Si_TeTeTeTe",
          AtomTypeData
          (
            GetElementTypes().e_Silicon,
            GetHybridOrbitalTypes().e_SP3,
            4,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            11.82,
            2.78,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      P_S2PPP
      (
        AddEnum
        (
          "P_S2PPP",
          AtomTypeData
          (
            GetElementTypes().e_Phosphorus,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( e_S),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            10.73,
            1.42,
            31.172,
            18.612,
            util::GetUndefined< double>()
          )
        )
      ),
      P_SP2PP
      (
        AddEnum
        (
          "P_SP2PP",
          AtomTypeData
          (
            GetElementTypes().e_Phosphorus,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_S, e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( e_Px),
            20.20,
            8.48,
            12.49,
            1.98,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      P_Di2DiPiPi
      (
        AddEnum
        (
          "P_Di2DiPiPi",
          AtomTypeData
          (
            GetElementTypes().e_Phosphorus,
            GetHybridOrbitalTypes().e_SP,
            1,
            1,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            17.53,
            4.95,
            11.61,
            1.68,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            1.525 // P_Te2TeTeTe * N_Di2DiPiPi/N_Te2TeTeTe
          )
        )
      ),
      P_DiDiPi2Pi
      (
        AddEnum
        (
          "P_DiDiPi2Pi",
          AtomTypeData
          (
            GetElementTypes().e_Phosphorus,
            GetHybridOrbitalTypes().e_SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( e_Px),
            16.78,
            4.77,
            11.89,
            2.02,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      P_Tr2TrTrPi
      (
        AddEnum
        (
          "P_Tr2TrTrPi",
          AtomTypeData
          (
            GetElementTypes().e_Phosphorus,
            GetHybridOrbitalTypes().e_SP2,
            2,
            1,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            15.59,
            3.74,
            11.64,
            1.80,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            1.643 // P_Te2TeTeTe * N_Tr2TrTrPi/N_Te2TeTeTe
          )
        )
      ),
      P_TrTrTrPi2
      (
        AddEnum
        (
          "P_TrTrTrPi2",
          AtomTypeData
          (
            GetElementTypes().e_Phosphorus,
            GetHybridOrbitalTypes().e_SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>( e_Px),
            15.18,
            3.76,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            1.739 // P_Te2TeTeTe * N_TrTrTrPi2/N_Te2TeTeTe
          )
        )
      ),
      P_Te2TeTeTe
      (
        AddEnum
        (
          "P_Te2TeTeTe",
          AtomTypeData
          (
            GetElementTypes().e_Phosphorus,
            GetHybridOrbitalTypes().e_SP3,
            3,
            1,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            14.57,
            3.24,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            24.041,
            12.095,
            1.538 // Miller 1990 contains contains a typo; the lewis-diagram is clearly of P_Te2TeTeTe, but they call it P_TeTeTeTe
                  // Presumably they meant Te2 since TeTeTeTe would be charged
          )
        )
      ),
      S_S2P2PP
      (
        AddEnum
        (
          "S_S2P2PP",
          AtomTypeData
          (
            GetElementTypes().e_Sulfur,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py),
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_S, e_Pz),
            12.39,
            2.38,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            22.977,
            11.053,
            util::GetUndefined< double>()
          )
        )
      ),
      S_SP2P2P
      (
        AddEnum
        (
          "S_SP2P2P",
          AtomTypeData
          (
            GetElementTypes().e_Sulfur,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_S, e_Px),
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            20.08,
            11.54,
            13.32,
            3.50,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      S_Di2Di2PiPi
      (
        AddEnum
        (
          "S_Di2Di2PiPi",
          AtomTypeData
          (
            GetElementTypes().e_Sulfur,
            GetHybridOrbitalTypes().e_SP,
            0,
            2,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            12.39,
            2.38,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      S_Di2DiPi2Pi
      (
        AddEnum
        (
          "S_Di2DiPi2Pi",
          AtomTypeData
          (
            GetElementTypes().e_Sulfur,
            GetHybridOrbitalTypes().e_SP,
            1,
            1,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( e_Py),
            17.78,
            6.96,
            12.86,
            2.94,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            3.729 // similar to Tr2Tr2TrPi
          )
        )
      ),
      S_DiDiPi2Pi2
      (
        AddEnum
        (
          "S_DiDiPi2Pi2",
          AtomTypeData
          (
            GetElementTypes().e_Sulfur,
            GetHybridOrbitalTypes().e_SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            17.42,
            6.80,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      S_Tr2Tr2TrPi
      (
        AddEnum
        (
          "S_Tr2Tr2TrPi",
          AtomTypeData
          (
            GetElementTypes().e_Sulfur,
            GetHybridOrbitalTypes().e_SP2,
            1,
            2,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            16.33,
            5.43,
            12.70,
            2.76,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            3.729
          )
        )
      ),
      S_Tr2TrTrPi2
      (
        AddEnum
        (
          "S_Tr2TrTrPi2",
          AtomTypeData
          (
            GetElementTypes().e_Sulfur,
            GetHybridOrbitalTypes().e_SP2,
            2,
            1,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            16.27,
            5.49,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            2.700
          )
        )
      ),
      S_Te2Te2TeTe
      (
        AddEnum
        (
          "S_Te2Te2TeTe",
          AtomTypeData
          (
            GetElementTypes().e_Sulfur,
            GetHybridOrbitalTypes().e_SP3,
            2,
            2,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            15.50,
            4.77,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            27.728,
            13.638,
            3.000
          )
        )
      ),
      Cl_S2P2P2P
      (
        AddEnum
        (
          "Cl_S2P2P2P",
          AtomTypeData
          (
            GetElementTypes().e_Chlorine,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_S, e_Px, e_Py),
            15.03,
            3.73,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            26.368,
            13.378,
            2.315
          )
        )
      ),
      Cl_SP2P2P2
      (
        AddEnum
        (
          "Cl_SP2P2P2",
          AtomTypeData
          (
            GetElementTypes().e_Chlorine,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_S),
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py, e_Pz),
            24.02,
            14.45,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            26.368,
            13.378,
            2.315
          )
        )
      ),
      K_S
      (
        AddEnum
        (
          "K_S",
          AtomTypeData
          (
            GetElementTypes().e_Potassium,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_S),
            storage::Set< AtomicOrbitalTypesEnum>(),
            4.341,
            1.95,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      K_P
      (
        AddEnum
        (
          "K_P",
          AtomTypeData
          (
            GetElementTypes().e_Potassium,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Px),
            storage::Set< AtomicOrbitalTypesEnum>(),
            2.7,
            1.195,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      // cations
      Be_S
      (
        AddEnum
        (
          "Be_S",
          AtomTypeData
          (
            GetElementTypes().e_Beryllium,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_S),
            storage::Set< AtomicOrbitalTypesEnum>(),
            18.21,
            9.32,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Be_P
      (
        AddEnum
        (
          "Be_P",
          AtomTypeData
          (
            GetElementTypes().e_Beryllium,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Px),
            storage::Set< AtomicOrbitalTypesEnum>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            14.25,
            5.32,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_SP
      (
        AddEnum
        (
          "B_SP",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_S, e_Px),
            storage::Set< AtomicOrbitalTypesEnum>(),
            25.40,
            14.05,
            19.40,
            7.38,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_PP
      (
        AddEnum
        (
          "B_PP",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py),
            storage::Set< AtomicOrbitalTypesEnum>(),
            18.91,
            7.37,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_DiDi
      (
        AddEnum
        (
          "B_DiDi",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().e_SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            23.48,
            9.64,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_DiPi
      (
        AddEnum
        (
          "B_DiPi",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().e_SP,
            1,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Py),
            storage::Set< AtomicOrbitalTypesEnum>(),
            22.16,
            8.94,
            19.16,
            7.37,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_TrTr
      (
        AddEnum
        (
          "B_TrTr",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().e_SP2,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            21.72,
            8.33,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_TrPi
      (
        AddEnum
        (
          "B_TrPi",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().e_SP2,
            1,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            21.08,
            8.02,
            19.08,
            7.37,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_TeTe
      (
        AddEnum
        (
          "B_TeTe",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().e_SP3,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            20.93,
            7.88,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      C_SPP
      (
        AddEnum
        (
          "C_SPP",
          AtomTypeData
          (
            GetElementTypes().e_Carbon,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_S, e_Px, e_Py),
            storage::Set< AtomicOrbitalTypesEnum>(),
            33.03,
            19.42,
            23.93,
            9.91,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      C_PPP
      (
        AddEnum
        (
          "C_PPP",
          AtomTypeData
          (
            GetElementTypes().e_Carbon,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            23.29,
            11.65,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      C_DiDiPi
      (
        AddEnum
        (
          "C_DiDiPi",
          AtomTypeData
          (
            GetElementTypes().e_Carbon,
            GetHybridOrbitalTypes().e_SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Py),
            storage::Set< AtomicOrbitalTypesEnum>(),
            29.85,
            13.29,
            23.86,
            9.83,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      C_DiPiPi
      (
        AddEnum
        (
          "C_DiPiPi",
          AtomTypeData
          (
            GetElementTypes().e_Carbon,
            GetHybridOrbitalTypes().e_SP,
            1,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            28.16,
            12.96,
            23.61,
            10.78,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      C_TrTrTr
      (
        AddEnum
        (
          "C_TrTrTr",
          AtomTypeData
          (
            GetElementTypes().e_Carbon,
            GetHybridOrbitalTypes().e_SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            28.14,
            11.83,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      C_TrTrPi
      (
        AddEnum
        (
          "C_TrTrPi",
          AtomTypeData
          (
            GetElementTypes().e_Carbon,
            GetHybridOrbitalTypes().e_SP2,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            27.36,
            11.91,
            23.68,
            10.45,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      C_TeTeTe
      (
        AddEnum
        (
          "C_TeTeTe",
          AtomTypeData
          (
            GetElementTypes().e_Carbon,
            GetHybridOrbitalTypes().e_SP3,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            26.71,
            11.37,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      N_SPPP
      (
        AddEnum
        (
          "N_SPPP",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_S, e_Px, e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            41.84,
            25.59,
            28.69,
            12.48,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      N_DiDiPiPi
      (
        AddEnum
        (
          "N_DiDiPiPi",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().e_SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            37.00,
            17.24,
            28.70,
            12.06,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      N_TrTrTrPi
      (
        AddEnum
        (
          "N_TrTrTrPi",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().e_SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            34.62,
            15.09,
            28.71,
            11.96,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      N_TeTeTeTe
      (
        AddEnum
        (
          "N_TeTeTeTe",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().e_SP3,
            4,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            33.29,
            14.14,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      O_S2PPP
      (
        AddEnum
        (
          "O_S2PPP",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( e_S),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            34.15,
            14.61,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      O_SP2PP
      (
        AddEnum
        (
          "O_SP2PP",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_S, e_Px, e_Py),
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            51.41,
            32.29,
            34.22,
            15.86,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      O_Di2DiPiPi
      (
        AddEnum
        (
          "O_Di2DiPiPi",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().e_SP,
            1,
            1,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            46.80,
            23.45,
            34.19,
            15.24,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      O_DiDiPi2Pi
      (
        AddEnum
        (
          "O_DiDiPi2Pi",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().e_SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( e_Py),
            44.56,
            22.34,
            33.95,
            15.53,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      O_Tr2TrTrPi
      (
        AddEnum
        (
          "O_Tr2TrTrPi",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().e_SP2,
            2,
            1,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            42.49,
            20.15,
            34.08,
            15.30,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      O_TrTrTrPi2
      (
        AddEnum
        (
          "O_TrTrTrPi2",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().e_SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            41.39,
            19.64,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      O_Te2TeTeTe
      (
        AddEnum
        (
          "O_Te2TeTeTe",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().e_SP3,
            3,
            1,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            40.31,
            18.70,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Mg_S
      (
        AddEnum
        (
          "Mg_S",
          AtomTypeData
          (
            GetElementTypes().e_Magnesium,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_S),
            storage::Set< AtomicOrbitalTypesEnum>(),
            15.03,
            7.64,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Mg_P
      (
        AddEnum
        (
          "Mg_P",
          AtomTypeData
          (
            GetElementTypes().e_Magnesium,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Px),
            storage::Set< AtomicOrbitalTypesEnum>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            10.60,
            4.67,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Al_SP
      (
        AddEnum
        (
          "Al_SP",
          AtomTypeData
          (
            GetElementTypes().e_Aluminum,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_S, e_Px),
            storage::Set< AtomicOrbitalTypesEnum>(),
            20.15,
            11.32,
            13.48,
            5.99,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Al_PP
      (
        AddEnum
        (
          "Al_PP",
          AtomTypeData
          (
            GetElementTypes().e_Aluminum,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py),
            storage::Set< AtomicOrbitalTypesEnum>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            14.34,
            6.03,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Al_DiDi
      (
        AddEnum
        (
          "Al_DiDi",
          AtomTypeData
          (
            GetElementTypes().e_Aluminum,
            GetHybridOrbitalTypes().e_SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            17.47,
            8.00,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Al_DiPi
      (
        AddEnum
        (
          "Al_DiPi",
          AtomTypeData
          (
            GetElementTypes().e_Aluminum,
            GetHybridOrbitalTypes().e_SP,
            1,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Py),
            storage::Set< AtomicOrbitalTypesEnum>(),
            17.25,
            7.59,
            13.92,
            6.00,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Al_TrTr
      (
        AddEnum
        (
          "Al_TrTr",
          AtomTypeData
          (
            GetElementTypes().e_Aluminum,
            GetHybridOrbitalTypes().e_SP2,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            16.28,
            7.01,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Al_TrPi
      (
        AddEnum
        (
          "Al_TrPi",
          AtomTypeData
          (
            GetElementTypes().e_Aluminum,
            GetHybridOrbitalTypes().e_SP2,
            1,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            16.28,
            6.74,
            14.06,
            5.92,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Al_TeTe
      (
        AddEnum
        (
          "Al_TeTe",
          AtomTypeData
          (
            GetElementTypes().e_Aluminum,
            GetHybridOrbitalTypes().e_SP3,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            15.75,
            6.64,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Si_SPP
      (
        AddEnum
        (
          "Si_SPP",
          AtomTypeData
          (
            GetElementTypes().e_Silicon,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_S, e_Px, e_Py),
            storage::Set< AtomicOrbitalTypesEnum>(),
            24.68,
            14.93,
            16.56,
            8.61,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Si_PPP
      (
        AddEnum
        (
          "Si_PPP",
          AtomTypeData
          (
            GetElementTypes().e_Silicon,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            16.56,
            11.42,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Si_DiDiPi
      (
        AddEnum
        (
          "Si_DiDiPi",
          AtomTypeData
          (
            GetElementTypes().e_Silicon,
            GetHybridOrbitalTypes().e_SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Py),
            storage::Set< AtomicOrbitalTypesEnum>(),
            21.43,
            10.95,
            16.50,
            8.60,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Si_DiPiPi
      (
        AddEnum
        (
          "Si_DiPiPi",
          AtomTypeData
          (
            GetElementTypes().e_Silicon,
            GetHybridOrbitalTypes().e_SP,
            1,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            20.62,
            11.56,
            16.55,
            10.02,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Si_TrTrTr
      (
        AddEnum
        (
          "Si_TrTrTr",
          AtomTypeData
          (
            GetElementTypes().e_Silicon,
            GetHybridOrbitalTypes().e_SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            19.96,
            9.99,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Si_TrTrPi
      (
        AddEnum
        (
          "Si_TrTrPi",
          AtomTypeData
          (
            GetElementTypes().e_Silicon,
            GetHybridOrbitalTypes().e_SP2,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            19.62,
            10.57,
            16.53,
            9.54,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Si_TeTeTe
      (
        AddEnum
        (
          "Si_TeTeTe",
          AtomTypeData
          (
            GetElementTypes().e_Silicon,
            GetHybridOrbitalTypes().e_SP3,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            18.97,
            10.08,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      P_SPPP
      (
        AddEnum
        (
          "P_SPPP",
          AtomTypeData
          (
            GetElementTypes().e_Phosphorus,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_S, e_Px, e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            31.24,
            18.61,
            20.72,
            11.55,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      P_DiDiPiPi
      (
        AddEnum
        (
          "P_DiDiPiPi",
          AtomTypeData
          (
            GetElementTypes().e_Phosphorus,
            GetHybridOrbitalTypes().e_SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            27.01,
            14.05,
            20.69,
            10.96,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      P_TrTrTrPi
      (
        AddEnum
        (
          "P_TrTrTrPi",
          AtomTypeData
          (
            GetElementTypes().e_Phosphorus,
            GetHybridOrbitalTypes().e_SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            25.14,
            12.72,
            20.68,
            10.76,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      P_TeTeTeTe
      (
        AddEnum
        (
          "P_TeTeTeTe",
          AtomTypeData
          (
            GetElementTypes().e_Phosphorus,
            GetHybridOrbitalTypes().e_SP3,
            4,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            24.10,
            12.09,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      S_S2PPP
      (
        AddEnum
        (
          "S_S2PPP",
          AtomTypeData
          (
            GetElementTypes().e_Sulfur,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( e_S),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            22.91,
            11.05,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      S_SP2PP
      (
        AddEnum
        (
          "S_SP2PP",
          AtomTypeData
          (
            GetElementTypes().e_Sulfur,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_S, e_Px, e_Py),
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            35.18,
            21.13,
            24.49,
            11.98,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      S_Di2DiPiPi
      (
        AddEnum
        (
          "S_Di2DiPiPi",
          AtomTypeData
          (
            GetElementTypes().e_Sulfur,
            GetHybridOrbitalTypes().e_SP,
            1,
            1,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            31.57,
            16.09,
            23.70,
            11.51,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      S_DiDiPi2Pi
      (
        AddEnum
        (
          "S_DiDiPi2Pi",
          AtomTypeData
          (
            GetElementTypes().e_Sulfur,
            GetHybridOrbitalTypes().e_SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( e_Py),
            30.61,
            15.78,
            24.00,
            11.92,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      S_Tr2TrTrPi
      (
        AddEnum
        (
          "S_Tr2TrTrPi",
          AtomTypeData
          (
            GetElementTypes().e_Sulfur,
            GetHybridOrbitalTypes().e_SP2,
            2,
            1,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            28.99,
            14.38,
            23.74,
            11.65,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      S_TrTrTrPi2
      (
        AddEnum
        (
          "S_TrTrTrPi2",
          AtomTypeData
          (
            GetElementTypes().e_Sulfur,
            GetHybridOrbitalTypes().e_SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            28.51,
            14.33,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      S_Te2TeTeTe
      (
        AddEnum
        (
          "S_Te2TeTeTe",
          AtomTypeData
          (
            GetElementTypes().e_Sulfur,
            GetHybridOrbitalTypes().e_SP3,
            3,
            1,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            27.65,
            13.64,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      O_Te2Te2Te2Te
      (
        AddEnum
        (
          "O_Te2Te2Te2Te",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().e_SP3,
            1,
            3,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            6.11,
            0.0,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      P_TeTeTeTePi
      (
        AddEnum
        (
          "P_TeTeTeTePi",
          AtomTypeData
          (
            GetElementTypes().e_Phosphorus,
            GetHybridOrbitalTypes().e_SP3,
            4,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            17.704,
            5.694,
            5.385,
            -0.015,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            1.523 // P_Te2TeTeTe * N_TrTrTrPi2/N_Tr2TrTrPi * N_Te2TeTeTe / N_Tr2TrTrPi
          )
        )
      ),
      S_TeTeTeTePiPi
      (
        AddEnum
        (
          "S_TeTeTeTePiPi",
          AtomTypeData
          (
            GetElementTypes().e_Sulfur,
            GetHybridOrbitalTypes().e_SP3,
            4,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            20.59,
            6.69,
            5.39,
            -2.85,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Br_SP2P2P2
      (
        AddEnum
        (
          "Br_SP2P2P2",
          AtomTypeData
          (
            GetElementTypes().e_Bromine,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_S),
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py, e_Pz),
            22.081,
            14.315,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            26.368,
            13.378,
            3.013
          )
        )
      ),
      Br_S2P2P2P
      (
        AddEnum
        (
          "Br_S2P2P2P",
          AtomTypeData
          (
            GetElementTypes().e_Bromine,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_S, e_Px, e_Py),
            13.108,
            3.516,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            26.368,
            13.378,
            3.013
          )
        )
      ),
      I_SP2P2P2
      (
        AddEnum
        (
          "I_SP2P2P2",
          AtomTypeData
          (
            GetElementTypes().e_Iodine,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_S),
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py, e_Pz),
            18.01,
            13.23,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            26.368,
            13.378,
            5.415
          )
        )
      ),
      I_S2P2P2P
      (
        AddEnum
        (
          "I_S2P2P2P",
          AtomTypeData
          (
            GetElementTypes().e_Iodine,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_S, e_Px, e_Py),
            12.677,
            3.375,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            26.368,
            13.378,
            5.415
          )
        )
      ),
      Se_Te2Te2TeTe
      (
        AddEnum
        (
          "Se_Te2Te2TeTe",
          AtomTypeData
          (
            GetElementTypes().e_Selenium,
            GetHybridOrbitalTypes().e_SP3,
            2,
            2,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            20.908,
            10.469,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      S_Te2TeTeTePi
      (
        AddEnum
        (
          "S_Te2TeTeTePi",
          AtomTypeData
          (
            GetElementTypes().e_Sulfur,
            GetHybridOrbitalTypes().e_SP3,
            3,
            1,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            18.136,
            5.708,
            2.283,
            -4.393,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      N_TrTrTrPi2Pi
      (
        AddEnum
        (
          "N_TrTrTrPi2Pi",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().e_SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( e_Py),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Sn_TeTeTeTe
      (
        AddEnum
        (
          "Sn_TeTeTeTe",
          AtomTypeData
          (
            GetElementTypes().e_Tin,
            GetHybridOrbitalTypes().e_SP3,
            4,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            10.4,
            5.39,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Ge_TeTeTeTe
      (
        AddEnum
        (
          "Ge_TeTeTeTe",
          AtomTypeData
          (
            GetElementTypes().e_Germanium,
            GetHybridOrbitalTypes().e_SP3,
            4,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            11.48,
            4.66,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_TeTeTeTe
      (
        AddEnum
        (
          "B_TeTeTeTe",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().e_SP3,
            4,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            1.53,
            0.0,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_TrTrTrPi
      (
        AddEnum
        (
          "B_TrTrTrPi",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().e_SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            1.87,
            0.0,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Cl_S2P2P2P2
      (
        AddEnum
        (
          "Cl_S2P2P2P2",
          AtomTypeData
          (
            GetElementTypes().e_Chlorine,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_S, e_Px, e_Py, e_Pz),
            14.45,
            0.0,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Se_Di2DiPi2Pi
      (
        AddEnum
        (
          "Se_Di2DiPi2Pi",
          AtomTypeData
          (
            GetElementTypes().e_Selenium,
            GetHybridOrbitalTypes().e_SP,
            1,
            1,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( e_Py),
            17.29,
            6.44,
            13.06,
            2.28,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Te_Te2Te2TeTe
      (
        AddEnum
        (
          "Te_Te2Te2TeTe",
          AtomTypeData
          (
            GetElementTypes().e_Tellurium,
            GetHybridOrbitalTypes().e_SP3,
            2,
            2,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            15.11,
            4.20,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      I_S2P2P2P2
      (
        AddEnum
        (
          "I_S2P2P2P2",
          AtomTypeData
          (
            GetElementTypes().e_Iodine,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_S, e_Px, e_Py, e_Pz),
            13.38,
            0.0,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            26.368,
            13.378,
            5.573
          )
        )
      ),
      As_Te2TeTeTe
      (
        AddEnum
        (
          "As_Te2TeTeTe",
          AtomTypeData
          (
            GetElementTypes().e_Arsenic,
            GetHybridOrbitalTypes().e_SP3,
            3,
            1,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            12.80,
            3.81,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      N_TrTrTrPiPi
      (
        AddEnum
        (
          "N_TrTrTrPiPi",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().e_SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      P_TrTrTrPiPi
      (
        AddEnum
        (
          "P_TrTrTrPiPi",
          AtomTypeData
          (
            GetElementTypes().e_Phosphorus,
            GetHybridOrbitalTypes().e_SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      N_TeTeTeTePi
      (
        AddEnum
        (
          "N_TeTeTeTePi",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().e_SP3,
            4,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      N_DiDiPi2Pi2
      (
        AddEnum
        (
          "N_DiDiPi2Pi2",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().e_SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      N_Di2DiPi2Pi
      (
        AddEnum
        (
          "N_Di2DiPi2Pi",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().e_SP,
            1,
            1,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( e_Py),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      N_Tr2TrTrPi2
      (
        AddEnum
        (
          "N_Tr2TrTrPi2",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().e_SP2,
            2,
            1,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      N_Te2Te2TeTe
      (
        AddEnum
        (
          "N_Te2Te2TeTe",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().e_SP3,
            2,
            2,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      S_Te2Te2Te2Te
      (
        AddEnum
        (
          "S_Te2Te2Te2Te",
          AtomTypeData
          (
            GetElementTypes().e_Sulfur,
            GetHybridOrbitalTypes().e_SP3,
            1,
            3,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      S_OhOhOhOhOhOh
      (
        AddEnum
        (
          "S_OhOhOhOhOhOh",
          AtomTypeData
          (
            GetElementTypes().e_Sulfur,
            GetHybridOrbitalTypes().e_SP3,
            6,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            15.91,
            3.12,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      S_TeTeTeTePi
      (
        AddEnum
        (
          "S_TeTeTeTePi",
          AtomTypeData
          (
            GetElementTypes().e_Sulfur,
            GetHybridOrbitalTypes().e_SP3,
            4,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      C_Di2DiPiPi
      (
        AddEnum
        (
          "C_Di2DiPiPi",
          AtomTypeData
          (
            GetElementTypes().e_Carbon,
            GetHybridOrbitalTypes().e_SP,
            1,
            1,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      H_
      (
        AddEnum
        (
          "H_",
          AtomTypeData
          (
            GetElementTypes().e_Hydrogen,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            13.60 * 2.0,
            13.60,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Li_
      (
        AddEnum
        (
          "Li_",
          AtomTypeData
          (
            GetElementTypes().e_Lithium,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            5.39 * 2.0,
            5.39,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Na_
      (
        AddEnum
        (
          "Na_",
          AtomTypeData
          (
            GetElementTypes().e_Sodium,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            5.14 * 2.0,
            5.14,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      K_
      (
        AddEnum
        (
          "K_",
          AtomTypeData
          (
            GetElementTypes().e_Potassium,
            GetHybridOrbitalTypes().e_Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            4.341 * 2.0,
            4.341,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      )
    {
      m_MaxGasteigerAtomTypeIndex = GetEnumCount();
      DetermineElectronegativityRatios();
      CalculateElectronegativityValues();
      ConnectAtomTypesToBaseTypes();
      BondLengths::Initialize( GetEnumIteratorFromIndex( 0), GetEnumIteratorFromIndex( e_Undefined.GetIndex()));
      SetLinearBonds();
    }

  /////////////////
  // data access //
  /////////////////

    //! Map from element type & charge to corresponding base atom type
    std::map< std::pair< ElementType, short>, AtomType> &GetBasicAtomTypesMap()
    {
      static std::map< std::pair< ElementType, short>, AtomType> s_basic_atom_types;
      return s_basic_atom_types;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AtomTypes::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief # of gasteiger atom types. These types precede any pseudo-atomtypes introduced during standardization
    //! @return # of gasteiger atom types.
    size_t AtomTypes::GetNumberGasteigerTypes()
    {
      return GetAtomTypes().m_MaxGasteigerAtomTypeIndex;
    }

    //! @brief get an atom type for an atom for which only the element type and charge is known
    //! @param ELEMENT_TYPE the element type
    //! @param CHARGE the expected charge
    const AtomType &AtomTypes::GetAtomType( const ElementType &ELEMENT_TYPE, const short &CHARGE)
    {
      std::pair< ElementType, short> key( ELEMENT_TYPE, CHARGE);

      // always call GetEnums first, otherwise map may not have been populated yet
      GetEnums();

      std::map< std::pair< ElementType, short>, AtomType>::iterator itr( GetBasicAtomTypesMap().find( key));
      if( itr == GetBasicAtomTypesMap().end())
      {
        // create the name for the base type
        const std::string base_type_name( ELEMENT_TYPE.GetName() + "_" + util::Format()( CHARGE));

        // the element type does not exist with this charge in the map yet, so add it
        itr = GetBasicAtomTypesMap().insert
              (
                std::make_pair( key, GetEnums().AddEnum( base_type_name, AtomTypeData( ELEMENT_TYPE, CHARGE)))
              ).first;
        // link the new type up with to itself as the base type, since it is already a base type
        itr->second->m_BaseType = &itr->second;
      }

      return itr->second;
    }

    //! @brief get an atom type for an atom, even if the atom type is a non-gasteiger type
    //! @param ATOM_TYPE name of the atom type
    const AtomType &AtomTypes::GetAtomType( const std::string &ATOM_TYPE)
    {
      // try creating the atom type directly
      const AtomType &atom_type( GetAtomTypes().GetEnumFromName( ATOM_TYPE));
      if( atom_type.IsDefined())
      {
        // atom type existed return the desired atom type
        return atom_type;
      }
      else if( ATOM_TYPE.empty() || ATOM_TYPE == EnumDataType::GetUndefinedEnumName())
      {
        // user requested the undefined atom type
        return atom_type;
      }

      // atom type does not appear to exist yet, try to decompose it into an element type and charge
      const size_t split_pos( ATOM_TYPE.find( '_'));

      // if the delimiter was not present or was at the end of the string, return undefined
      if( split_pos >= ATOM_TYPE.size() - 1)
      {
        return atom_type;
      }

      // get the charge string
      const std::string charge_string( ATOM_TYPE.substr( split_pos + 1));

      // if the string is non-numeric, return the undefined atom type
      if( !util::IsNumerical( charge_string))
      {
        return atom_type;
      }

      // convert the charge string into the charge
      const short charge( util::ConvertStringToNumericalValue< short>( charge_string));

      // otherwise, get the element type
      ElementType element_type( ATOM_TYPE.substr( 0, split_pos));

      // create the element type
      return AtomTypes::GetAtomType( element_type, charge);
    }

    //! @brief sets up m_BaseAtomType for all atom types
    //! Only called in AtomTypes(); once
    void AtomTypes::ConnectAtomTypesToBaseTypes()
    {
      // connect the undefined type to itself
      e_Undefined->m_BaseType = &e_Undefined;

      // create a set of all element type - charge pairs in the atom types
      std::set< std::pair< ElementType, short> > keys;
      for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        keys.insert( std::make_pair( ( *itr)->GetElementType(), ( *itr)->GetFormalCharge()));
      }

      std::map< std::pair< ElementType, short>, AtomType> &base_atom_types_map( GetBasicAtomTypesMap());

      // add all the base types to the atoms types
      for
      (
        std::set< std::pair< ElementType, short> >::const_iterator itr( keys.begin()), itr_end( keys.end());
        itr != itr_end;
        ++itr
      )
      {
        // create the name for the base type
        const std::string base_type_name( itr->first.GetName() + "_" + util::Format()( itr->second));

        base_atom_types_map.insert
        (
          std::make_pair( *itr, AddEnum( base_type_name, AtomTypeData( itr->first, itr->second)))
        );
      }

      // connect all the base types to themselves
      for
      (
        std::map< std::pair< ElementType, short>, AtomType>::iterator
          itr( base_atom_types_map.begin()), itr_end( base_atom_types_map.end());
        itr != itr_end;
        ++itr
      )
      {
        itr->second->m_BaseType = &itr->second;
      }

      // connect all the gasteiger types to the corresponding base type
      for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        if( !( *itr)->IsGasteigerAtomType())
        {
          // bases of non-gasteiger types were set up in the prior for loop
          continue;
        }

        AtomType atom_type( *itr);
        atom_type->m_BaseType
          = &base_atom_types_map[ std::make_pair( atom_type->GetElementType(), atom_type->GetFormalCharge())];
      }
    }

    //! @brief Determine ratios between ionization potential and electronegativity for neutral atoms and cations
    void AtomTypes::DetermineElectronegativityRatios() const
    {
      // track the sum of the defined ground state electronegativity for each of the difference types, as well as the
      // charged state electronegativity (only when both charged state and neutral state are defined)
      double sum_neutral_sigma_ionization_potentials( 0.0);
      double sum_neutral_pi_ionization_potentials( 0.0);
      double sum_neutral_sigma_electron_affinities( 0.0);
      double sum_neutral_pi_electron_affinities( 0.0);
      double sum_cationic_sigma_ionization_potentials( 0.0);
      double sum_cationic_pi_ionization_potentials( 0.0);
      double num_neutral_sigma_ionization_potentials( 0.0);
      double num_neutral_pi_ionization_potentials( 0.0);
      double num_neutral_sigma_electron_affinities( 0.0);
      double num_neutral_pi_electron_affinities( 0.0);
      double num_cationic_sigma_ionization_potentials( 0.0);
      double num_cationic_pi_ionization_potentials( 0.0);

      for
      (
        storage::Vector< AtomType>::const_iterator itr_type( this->Begin()), itr_type_end( this->End());
        itr_type != itr_type_end;
        ++itr_type
      )
      {
        const AtomTypeData &data( **itr_type);
        const int charge( data.GetFormalCharge());

        if( charge < 0 || charge > 1)
        {
          continue;
        }

        const double s_ip( data.GetAtomTypeProperty( AtomTypeData::e_SigmaValenceStateIonizationPotential));
        const double s_ea( data.GetAtomTypeProperty( AtomTypeData::e_SigmaValenceStateElectronAffinity));
        const double p_ip( data.GetAtomTypeProperty( AtomTypeData::e_PiValenceStateIonizationPotential));
        const double p_ea( data.GetAtomTypeProperty( AtomTypeData::e_PiValenceStateElectronAffinity));

        if( charge == 1)
        {
          if( util::IsDefined( s_ip))
          {
            sum_cationic_sigma_ionization_potentials += s_ip;
            num_cationic_sigma_ionization_potentials += 1.0;
          }

          if( util::IsDefined( p_ip))
          {
            sum_cationic_pi_ionization_potentials += p_ip;
            num_cationic_pi_ionization_potentials += 1.0;
          }
        }
        else if( charge == 0)
        {
          if( util::IsDefined( s_ip) && util::IsDefined( s_ea))
          {
            sum_neutral_sigma_electron_affinities += s_ea;
            num_neutral_sigma_electron_affinities += 1.0;
            sum_neutral_sigma_ionization_potentials += s_ip;
            num_neutral_sigma_ionization_potentials += 1.0;
          }

          if( util::IsDefined( p_ip) && util::IsDefined( p_ea))
          {
            sum_neutral_pi_electron_affinities += p_ea;
            num_neutral_pi_electron_affinities += 1.0;
            sum_neutral_pi_ionization_potentials += p_ip;
            num_neutral_pi_ionization_potentials += 1.0;
          }
        }
      }

      const double ave_neutral_sigma_ionization_potentials( sum_neutral_sigma_ionization_potentials / num_neutral_sigma_ionization_potentials);
      const double ave_neutral_pi_ionization_potentials( sum_neutral_pi_ionization_potentials     / num_neutral_pi_ionization_potentials);
      const double ave_neutral_sigma_electron_affinities( sum_neutral_sigma_electron_affinities    / num_neutral_sigma_electron_affinities);
      const double ave_neutral_pi_electron_affinities( sum_neutral_pi_electron_affinities       / num_neutral_pi_electron_affinities);
      const double ave_cationic_sigma_ionization_potentials( sum_cationic_sigma_ionization_potentials / num_cationic_sigma_ionization_potentials);
      const double ave_cationic_pi_ionization_potentials( sum_cationic_pi_ionization_potentials    / num_cationic_pi_ionization_potentials);

      AtomTypeData::GetAverageNeutralSigmaIPToAnionIPRatio()
        = ave_neutral_sigma_ionization_potentials / ave_neutral_sigma_electron_affinities;
      AtomTypeData::GetAverageNeutralPiIPToAnionIPRatio()
        = ave_neutral_pi_ionization_potentials / ave_neutral_pi_electron_affinities;
      AtomTypeData::GetAverageCationSigmaIPToNeutralIPRatio()
        = ave_cationic_sigma_ionization_potentials / ave_neutral_sigma_ionization_potentials;
      AtomTypeData::GetAverageCationPiIPToNeutralIPRatio()
        = ave_cationic_pi_ionization_potentials / ave_neutral_pi_ionization_potentials;

    } // DetermineElectronegativityRatios

    //! @brief sets up m_ChargeToSigmaEN and m_ChargeToPiEN for all atom type data
    //! Only called in AtomTypes(); once
    void AtomTypes::CalculateElectronegativityValues()
    {
      // create a map with empty lists for each type difference
      std::map< AtomTypeData::TypeDifference, storage::List< AtomType> > undefined_map;
      undefined_map[ AtomTypeData::e_NumberBondingSOrbitals] = storage::List< AtomType>();
      undefined_map[ AtomTypeData::e_NumberBondingPOrbitals] = storage::List< AtomType>();
      undefined_map[ AtomTypeData::e_NumberLonePairOrbitals] = storage::List< AtomType>();

      // get a map from atomic number to a corresponding list of atom types
      storage::Map< size_t, storage::List< AtomType> > atomic_number_map;

      // iterate through atom types and store atom types by their atomic number
      for( AtomTypes::const_iterator itr_at( Begin()), itr_at_end( End()); itr_at != itr_at_end; ++itr_at)
      {
        atomic_number_map[ ( *itr_at)->GetElementType()->GetAtomicNumber()].PushBack( *itr_at);
      }

      // next, walk over each atomic number
      for
      (
        storage::Map< size_t, storage::List< AtomType> >::const_iterator
          itr_element( atomic_number_map.Begin()), itr_element_end( atomic_number_map.End());
        itr_element != itr_element_end;
        ++itr_element
      )
      {
        // make maps from hybrid orbital type, to a list of atom types with that hybridization
        // make separate maps for neutral atoms and ions
        storage::Map< HybridOrbitalType, storage::List< AtomType> > orbital_type_to_ions_map;
        storage::List< AtomType> neutrals_list;

        // load the maps with atom types for this element and hybridization state
        for
        (
          storage::List< AtomType>::const_iterator
            itr_type( itr_element->second.Begin()), itr_type_end( itr_element->second.End());
          itr_type != itr_type_end;
          ++itr_type
        )
        {
          const int charge( ( *itr_type)->GetFormalCharge());
          if( charge == 0)
          {
            neutrals_list.PushBack( *itr_type);
          }
          else
          {
            orbital_type_to_ions_map[ ( *itr_type)->GetHybridOrbitalType()].PushBack( *itr_type);
          }
        }

        // for each neutral atom type of this element, determine the relationship of each atom type to the others in the list
        for
        (
          storage::List< AtomType>::iterator
            itr_type( neutrals_list.Begin()), itr_type_end( neutrals_list.End());
          itr_type != itr_type_end;
          ++itr_type
        )
        {
          // get references to the type and data for easier access
          AtomType &type( *itr_type);
          AtomTypeData &data( *type);

          // a map from each neutral atom to the atom type ions that differs only by TypeDifference
          std::map< AtomTypeData::TypeDifference, storage::List< AtomType> >
            atom_types_to_neighbor_types( undefined_map);

          // make a list containing the atom types that may represent a sigma/pi ionized or ground state
          const storage::List< AtomType> &types_to_consider
          (
            orbital_type_to_ions_map[ data.GetHybridOrbitalType()] // neutral atoms consider ions
          );

          // look for compatible types in the list
          for
          (
            storage::List< AtomType>::const_iterator
              itr_cat( types_to_consider.Begin()), itr_cat_end( types_to_consider.End());
            itr_cat != itr_cat_end;
            ++itr_cat
          )
          {
            const AtomTypeData::TypeDifference diff( ( *itr_type)->DifferenceFrom( **itr_cat));

            if( diff == AtomTypeData::e_Other || diff == AtomTypeData::e_None)
            {
              // skip the type since it is more than one electron away from being the same as data
              continue;
            }

            storage::List< AtomType> &similar_types( atom_types_to_neighbor_types[ diff]);

            if( similar_types.IsEmpty())
            {
              similar_types.PushBack( *itr_cat);
            }
            else if( similar_types.FirstElement()->GetFormalCharge() > ( *itr_cat)->GetFormalCharge())
            {
              similar_types.PushBack( *itr_cat);
            }
            else if( ( *itr_cat)->GetFormalCharge() < similar_types.FirstElement()->GetFormalCharge())
            {
              similar_types.PushFront( *itr_cat);
            }
            else
            {
              // check whether the electronegativity was defined
              const double electronegativity( ( *itr_cat)->GetElectronegativity( diff));
              if( !util::IsDefined( electronegativity) || electronegativity == 0.0)
              {
                // add it to the end of the list because it cannot be used to determine the electronegativity parameters
                similar_types.PushBack( *itr_cat);
              }
              else
              {
                // get the electronegativity of the first element
                const double electronegativity_first( similar_types.FirstElement()->GetElectronegativity( diff));
                if( !util::IsDefined( electronegativity_first) || electronegativity_first == 0.0)
                {
                  // add it to the front of the list because this type is better than the existing first type in the list
                  // for determination of parameters
                  similar_types.PushFront( *itr_cat);
                }
                else
                {
                  BCL_MessageCrt
                  (
                    "AtomTypes indistiguishable types ERROR: " + type.GetName() + " <"
                    + AtomTypeData::GetTypeDifferenceName( diff) + "> " + itr_cat->GetName() + " OR "
                    + similar_types.FirstElement().GetName() + " OR " + similar_types.LastElement().GetName()
                  );
                }
              }
            }
          }

          // calculate the sigma en = f(charge)
          for
          (
            std::map< AtomTypeData::TypeDifference, storage::List< AtomType> >::iterator
              itr_diff( atom_types_to_neighbor_types.begin()), itr_diff_end( atom_types_to_neighbor_types.end());
            itr_diff != itr_diff_end;
            ++itr_diff
          )
          {
            const AtomTypeData::TypeDifference &diff( itr_diff->first);
            if
            (
              itr_diff->second.IsEmpty()
              || !itr_diff->second.FirstElement().IsDefined()
              || itr_diff->second.FirstElement()->GetElectronegativity( diff) == 0.0
            )
            {
              // set up the electronegativity using the average values electronegativity
              data.SetSimilarOrbitalElectronegativity( diff);
            }
            else
            {
              AtomType &similar_type( itr_diff->second.FirstElement());
              AtomTypeData &similar_data( *similar_type);

              // get the electronegativity from the cation, if one was found
              data.SetSimilarOrbitalElectronegativity( diff, similar_data.GetElectronegativity( diff));
            }
            for
            (
              storage::List< AtomType>::iterator
                itr_list( itr_diff->second.Begin()), itr_list_end( itr_diff->second.End());
              itr_list != itr_list_end;
              ++itr_list
            )
            {
              ( *itr_list)->SetSimilarOrbitalElectronegativity( diff, data);
            }
          }
        }
      } // for each element
    } // CalculateElectronegativityValues

    //! @brief set the atom types that are linear, namely all DiDiPiPi and O_DiDiPi2Pi2
    void AtomTypes::SetLinearBonds()
    {
      for( iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        if( ( *itr)->GetNumberBonds() == size_t( 2) && ( *itr)->GetNumberElectronsInBonds() >= size_t( 4))
        {
          ( *itr)->SetFormsOnlyLinearBonds();
        }
        else if( ( *itr)->GetNumberBonds() == size_t( 1) || *itr == O_DiDiPi2Pi2)
        {
          // single bonds are linear by definition
          // O_DiDiPi2Pi2 only comes up when connected to Si-O-Si
          ( *itr)->SetFormsOnlyLinearBonds();
        }
      }
    }

    const AtomTypes &GetAtomTypes()
    {
      return AtomTypes::GetEnums();
    }

    namespace
    {
      // trivial static instantiation
      //     const bool s_is_defined( !GetAtomTypes().e_Undefined.IsDefined());
    }
  } // namespace chemistry

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< chemistry::AtomTypeData, chemistry::AtomTypes>;

  } // namespace util
} // namespace bcl
