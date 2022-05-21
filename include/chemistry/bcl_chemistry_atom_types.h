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

#ifndef BCL_CHEMISTRY_ATOM_TYPES_H_
#define BCL_CHEMISTRY_ATOM_TYPES_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_type_data.h"
#include "util/bcl_util_enumerate.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AtomTypes
    //! @brief enumerator class for AtomTypes that separate atoms by their properties
    //! @details This is the enumerator class to be used for accessing Atom type information in Chemistry related classes.
    //! @details The values for sigma and pi electronegativity are taken from Bermann, Hinze, Structure & Binding, 1987
    //! @details available at http://www.springerlink.com/content/w300520690302287/fulltext.pdf
    //!
    //! @see @link example_chemistry_atom_types.cpp @endlink
    //! @author mueller, mendenjl
    //! @date 08/26/2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AtomTypes :
      public util::Enumerate< AtomTypeData, AtomTypes>
    {
      friend class util::Enumerate< AtomTypeData, AtomTypes>;

    //////////
    // data //
    //////////

    private:

      size_t m_MaxGasteigerAtomTypeIndex; //!< Index of last gasteiger atom type

    public:

      // declare all atom types
      // neutral
      const AtomType H_S;
      const AtomType Li_S;
      const AtomType Li_P;
      const AtomType Be_SP;
      const AtomType Be_PP;
      const AtomType Be_DiDi;
      const AtomType Be_DiPi;
      const AtomType Be_TrTr;
      const AtomType Be_TrPi;
      const AtomType Be_TeTe;
      const AtomType B_SPP;
      const AtomType B_PPP;
      const AtomType B_DiDiPi;
      const AtomType B_DiPiPi;
      const AtomType B_TrTrTr;
      const AtomType B_TrTrPi;
      const AtomType B_TeTeTe;
      const AtomType C_SPPP;
      const AtomType C_DiDiPiPi;
      const AtomType C_TrTrTrPi;
      const AtomType C_TeTeTeTe;
      const AtomType N_S2PPP;
      const AtomType N_SP2PP;
      const AtomType N_Di2DiPiPi;
      const AtomType N_DiDiPi2Pi;
      const AtomType N_Tr2TrTrPi;
      const AtomType N_TrTrTrPi2;
      const AtomType N_Te2TeTeTe;
      const AtomType O_S2P2PP;
      const AtomType O_SP2P2P;
      const AtomType O_Di2Di2PiPi;
      const AtomType O_Di2DiPi2Pi;
      const AtomType O_DiDiPi2Pi2;
      const AtomType O_Tr2Tr2TrPi;
      const AtomType O_Tr2TrTrPi2;
      const AtomType O_Te2Te2TeTe;
      const AtomType F_S2P2P2P;
      const AtomType F_SP2P2P2;
      const AtomType Na_S;
      const AtomType Na_P;
      const AtomType Mg_SP;
      const AtomType Mg_PP;
      const AtomType Mg_DiDi;
      const AtomType Mg_DiPi;
      const AtomType Mg_TrTr;
      const AtomType Mg_TrPi;
      const AtomType Mg_TeTe;
      const AtomType Al_SPP;
      const AtomType Al_PPP;
      const AtomType Al_DiDiPi;
      const AtomType Al_DiPiPi;
      const AtomType Al_TrTrTr;
      const AtomType Al_TrTrPi;
      const AtomType Al_TeTeTe;
      const AtomType Si_SPPP;
      const AtomType Si_DiDiPiPi;
      const AtomType Si_TrTrTrPi;
      const AtomType Si_TeTeTeTe;
      const AtomType P_S2PPP;
      const AtomType P_SP2PP;
      const AtomType P_Di2DiPiPi;
      const AtomType P_DiDiPi2Pi;
      const AtomType P_Tr2TrTrPi;
      const AtomType P_TrTrTrPi2;
      const AtomType P_Te2TeTeTe;
      const AtomType S_S2P2PP;
      const AtomType S_SP2P2P;
      const AtomType S_Di2Di2PiPi;
      const AtomType S_Di2DiPi2Pi;
      const AtomType S_DiDiPi2Pi2;
      const AtomType S_Tr2Tr2TrPi;
      const AtomType S_Tr2TrTrPi2;
      const AtomType S_Te2Te2TeTe;
      const AtomType Cl_S2P2P2P;
      const AtomType Cl_SP2P2P2;
      const AtomType K_S;
      const AtomType K_P;

      // cation
      const AtomType Be_S;
      const AtomType Be_P;
      const AtomType B_SP;
      const AtomType B_PP;
      const AtomType B_DiDi;
      const AtomType B_DiPi;
      const AtomType B_TrTr;
      const AtomType B_TrPi;
      const AtomType B_TeTe;
      const AtomType C_SPP;
      const AtomType C_PPP;
      const AtomType C_DiDiPi;
      const AtomType C_DiPiPi;
      const AtomType C_TrTrTr;
      const AtomType C_TrTrPi;
      const AtomType C_TeTeTe;
      const AtomType N_SPPP;
      const AtomType N_DiDiPiPi;
      const AtomType N_TrTrTrPi;
      const AtomType N_TeTeTeTe;
      const AtomType O_S2PPP;
      const AtomType O_SP2PP;
      const AtomType O_Di2DiPiPi;
      const AtomType O_DiDiPi2Pi;
      const AtomType O_Tr2TrTrPi;
      const AtomType O_TrTrTrPi2;
      const AtomType O_Te2TeTeTe;
      const AtomType Mg_S;
      const AtomType Mg_P;
      const AtomType Al_SP;
      const AtomType Al_PP;
      const AtomType Al_DiDi;
      const AtomType Al_DiPi;
      const AtomType Al_TrTr;
      const AtomType Al_TrPi;
      const AtomType Al_TeTe;
      const AtomType Si_SPP;
      const AtomType Si_PPP;
      const AtomType Si_DiDiPi;
      const AtomType Si_DiPiPi;
      const AtomType Si_TrTrTr;
      const AtomType Si_TrTrPi;
      const AtomType Si_TeTeTe;
      const AtomType P_SPPP;
      const AtomType P_DiDiPiPi;
      const AtomType P_TrTrTrPi;
      const AtomType P_TeTeTeTe;
      const AtomType S_S2PPP;
      const AtomType S_SP2PP;
      const AtomType S_Di2DiPiPi;
      const AtomType S_DiDiPi2Pi;
      const AtomType S_Tr2TrTrPi;
      const AtomType S_TrTrTrPi2;
      const AtomType S_Te2TeTeTe;

      // new atom types
      //! IMPORTANT: New atom types must be added to the *end* of this list. Adding them anywhere else will break the rotamer
      //!            libraries. The only exception is for atoms that have 0 bonds, which can be at the very end of this list
      //! TODO: Encode the current list of atom types in the rotamer library files rather than depending on their ordering
      //!       and then map the types in the given rotamer library file back to the current types
      const AtomType O_Te2Te2Te2Te;
      const AtomType P_TeTeTeTePi;
      const AtomType S_TeTeTeTePiPi;
      const AtomType Br_SP2P2P2;
      const AtomType Br_S2P2P2P;
      const AtomType I_SP2P2P2;
      const AtomType I_S2P2P2P;
      const AtomType Se_Te2Te2TeTe;
      const AtomType S_Te2TeTeTePi;
      const AtomType N_TrTrTrPi2Pi;
      const AtomType Sn_TeTeTeTe;
      const AtomType Ge_TeTeTeTe;
      const AtomType B_TeTeTeTe;
      const AtomType B_TrTrTrPi;
      const AtomType Cl_S2P2P2P2;
      const AtomType Se_Di2DiPi2Pi;
      const AtomType Te_Te2Te2TeTe;
      const AtomType I_S2P2P2P2;
      const AtomType As_Te2TeTeTe;
      const AtomType N_TrTrTrPiPi;
      const AtomType P_TrTrTrPiPi;
      const AtomType N_TeTeTeTePi;
      const AtomType N_DiDiPi2Pi2;
      const AtomType N_Di2DiPi2Pi;
      const AtomType N_Tr2TrTrPi2;
      const AtomType N_Te2Te2TeTe;
      const AtomType S_Te2Te2Te2Te;
      const AtomType S_OhOhOhOhOhOh; // octahedral sulfur
      const AtomType S_TeTeTeTePi;
      const AtomType C_Di2DiPiPi; // carbanion with triple bond; :C#

      // cations of group 1 elements (no hybrid orbitals)
      const AtomType H_;
      const AtomType Li_;
      const AtomType Na_;
      const AtomType K_;

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct all AtomTypes
      AtomTypes();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get an atom type for an atom for which only the element type and charge is known
      //! @param ELEMENT_TYPE the element type
      //! @param CHARGE the expected charge
      static const AtomType &GetAtomType( const ElementType &ELEMENT_TYPE, const short &CHARGE);

      //! @brief get an atom type for an atom, even if the atom type is a non-gasteiger type
      //! @param NAME name of the atom type
      static const AtomType &GetAtomType( const std::string &ATOM_TYPE);

      //! @brief # of gasteiger atom types. These types precede any pseudo-atomtypes introduced during standardization
      //! @return # of gasteiger atom types.
      static size_t GetNumberGasteigerTypes();

    private:

      //! @brief sets up m_BaseAtomType for all atom types
      //! Only called in AtomTypes(); once
      void ConnectAtomTypesToBaseTypes();

      //! @brief Determine ratios between ionization potential and electronegativity for neutral atoms and cations
      void DetermineElectronegativityRatios() const;

      //! @brief sets up m_ChargeToSigmaEN and m_ChargeToPiEN for all atom type data
      //! Only called in AtomTypes(); once
      void CalculateElectronegativityValues();

      //! @brief set the atom types that are linear, namely all DiDiPiPi and O_DiDiPi2Pi2
      void SetLinearBonds();

    }; // class AtomTypes

    BCL_API
    const AtomTypes &GetAtomTypes();

  } // namespace chemistry

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< chemistry::AtomTypeData, chemistry::AtomTypes>;

  } // namespace util
} // namespace bcl

#endif // BCL_CHEMISTRY_ATOM_TYPES_H_

