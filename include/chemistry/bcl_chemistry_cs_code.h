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

#ifndef BCL_CHEMISTRY_CS_CODE_H_
#define BCL_CHEMISTRY_CS_CODE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically
#include "descriptor/bcl_descriptor.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_conformational_interface.h"
#include "bcl_chemistry_atom_environment.h"
#include "bcl_chemistry_bond_constitutional.h"
#include "bcl_chemistry_conformation_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CsCode
    //! @brief This class contains atom constitution information
    //! @details  stores the atom type, bond type of bonds that atom has and valence bonds of the atom
    //!
    //! @see @link example_chemistry_cs_code.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date Mar 16, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CsCode :
      public util::ObjectInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief virtual copy constructor
      CsCode *Clone() const;

      //! @brief returns class name
      //! the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief compares two elements of the atom environment
      //! @return true if left-hand side has higher bond type or atomic number (for same bond type)
      static bool CompareAtomsAndBonds
      (
        const storage::Pair< ConfigurationalBondType, util::SiPtr< const AtomConformationalInterface> > &LHS,
        const storage::Pair< ConfigurationalBondType, util::SiPtr< const AtomConformationalInterface> > &RHS
      );

      //! @brief generate the code vector for each element of the atom environment
      //! @return the code vector for each element of the atom environment
      static storage::Vector< float> GenerateCodeVector
      (
        const ConformationInterface &MOLECULE,
        const AtomConformationalInterface &ATOM,
        const ConfigurationalBondType &BOND_TYPE,
        const size_t NR_HYDROGENS,
        const float  RING_CLOSURE,
        storage::Vector< descriptor::CheminfoProperty> &PROPERTIES
      );

      //! @brief parse the miscellaneous property string containing the carbon chemical shift (NMRshiftDB format)
      //! @return map between the atom index and the according chemical shift
      static storage::Map< size_t, float> ParseChemicalShifts( const ConformationInterface &MOLECULE);

      //! @brief implement a limited depth search from the atom of interest through the environment
      static void LimitedDepthFirstSearch
      (
        const ConformationInterface &MOLECULE,
        const AtomConformationalInterface &ATOM,
        const AtomEnvironment &ATOM_ENVIRONMENT,
        const size_t DEPTH,
        storage::Vector< float> &CODE,
        storage::Vector< descriptor::CheminfoProperty> &PROPERTIES,
        const bool &ONLY_FOLLOW_CONJUGATED_BONDS = false
      );

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief set bonds constitution of atom
      //! @param BONDS the bonds that the atom connects
      void SetBonds( const storage::Vector< BondConstitutional> &BONDS);

      //! @brief copy the atom and bonds using the specified difference in pointer address
      //! @param ATOM the atom to copy
      //! @param DIFF difference in begin of vector that atoms are located in
      void Copy( const CsCode &ATOM, const ptrdiff_t &DIFF);

    }; // class CsCode

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_CS_CODE_H_
