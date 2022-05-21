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

#ifndef BCL_CHEMISTRY_BOND_LENGTHS_H_
#define BCL_CHEMISTRY_BOND_LENGTHS_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_type_data.h"
#include "storage/bcl_storage_vector.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BondLengths
    //! @brief calculates and stores bond length information
    //! Determines bond lengths of atoms based on element type & bond type (Single/Double/Triple/Aromatic)
    //! by using the covalent radius around an atom of that type
    //!
    //! @see @link example_chemistry_bond_lengths.cpp @endlink
    //! @author mendenjl
    //! @date Jul 06, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API BondLengths :
      public util::ObjectInterface
    {
    /////////////
    // friends //
    /////////////

      friend class AtomTypes; //! calls initialize function, below, with enumerated members

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief convert bond order or aromatic into the corresponding radius
      //! @param BOND_ORDER_OR_AROMATIC bond type in notation: 1=single, 2=double, 3=triple, 4=aromatic
      static AtomTypeData::Properties BondOrderProperty( const size_t &BOND_ORDER_OR_AROMATIC);

      //! @brief set radius information for the atom types of a given phenotype
      //! @param ATOM_TYPES atom types for which to set the given property
      //! @param PROPERTY the property to set
      //! @param VALUE the value to set the property to
      static void SetPropertyInfo
      (
        storage::Vector< AtomType> &ATOM_TYPES,
        const AtomTypeData::Properties &PROPERTY,
        const double &VALUE
      );

      //! @brief initialization of bond lengths
      //! @param BEGIN iterator to beginning of atom types
      //! @param END iterator to end of atom types
      static void Initialize( std::vector< AtomType>::iterator BEGIN, std::vector< AtomType>::iterator END);

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief virtual copy constructor
      BondLengths *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get name of the current class
      //! @return name of the class
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief get the average bond length of two atom types joined by a bond of specified type
      //! @param ATOM_TYPE_A either atom type in the bond
      //! @param BOND_ORDER_OR_AROMATIC bond type in notation: 1=single, 2=double, 3=triple, 4=aromatic
      //! @param ATOM_TYPE_B the other atom type in the bond
      //! @return return the estimated bond length, in angstroms
      static double GetBondLength
      (
        const AtomType &ATOM_TYPE_A,
        const size_t &BOND_ORDER_OR_AROMATIC,
        const AtomType &ATOM_TYPE_B
      );

      //! @brief Get the covalent radius of an atom about a specific bond type
      //! @param ATOM_TYPE the atom type for which to estimate the covalent radius
      //! @param BOND_ORDER_OR_AROMATIC bond type in notation: 1=single, 2=double, 3=triple, 4=aromatic
      //! @return return the estimated covalent radius of the atom, in angstroms
      static double GetCovalentRadius( const AtomType &ATOM_TYPE, const size_t &BOND_ORDER_OR_AROMATIC);

      //! @brief Get the average covalent radius for a particular bond of an atom
      //! @param ELEMENT_TYPE the element type of the atom
      //! @param NUMBER_BONDS the total number of bonds for the atom
      //! @param NUMBER_ELECTRONS_IN_BONDS the total # of electrons in bonds around the atom
      //! @param BOND_ORDER_OR_AROMATIC bond type in sdf notation (1=single, 2=double, 3=triple, 4=aromatic)
      //! @return return the estimated covalent radius of the atom, in angstroms
      static double GetCovalentRadius
      (
        const ElementType &ELEMENT_TYPE,
        const size_t &NUMBER_BONDS,
        const size_t &NUMBER_ELECTRONS_IN_BONDS,
        const size_t &BOND_ORDER_OR_AROMATIC
      );

      //! @brief Get the average covalent radius of all the bonds about an atom
      //! @param ATOM the atom for which to estimate the covalent radius
      //! @return return the estimated covalent radius of the atom, in angstroms
      static double GetAverageCovalentRadius( const AtomConstitutionalInterface &ATOM);

      //! @brief Get the average covalent radius of all the bonds about an atom
      //! @param ATOM the atom for which to estimate the covalent radius
      //! @return return the estimated covalent radius of the atom, in angstroms
      static double GetAverageCovalentRadius( const AtomConfigurationalInterface &ATOM);

      //! @brief Get the average covalent radius of all the bonds about an atom
      //! @param ATOM the atom for which to estimate the covalent radius
      //! @return return the estimated covalent radius of the atom, in angstroms
      static double GetAverageCovalentRadius( const AtomConformationalInterface &ATOM);

      //! @brief solve for VDW radii
      //! @param SEP_MAP map containing strings with two atom types each and a set of distances between those types
      //! @param MIN_VALUES minimum values for the vdw-radii (any values smaller than 1 will be set up to 1)
      //! @return Map containing VDW radii computed for each atom type (assuming lp_solve is installed)
      static storage::Map< std::string, double> ComputeVdwRadii
      (
        storage::Map< std::string, storage::Vector< double> > &RADII,
        const storage::Map< std::string, double> &MIN_VALUES
      );

    ////////////////
    // operations //
    ////////////////

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief Estimate the covalent radius of an atom
      //! @param ATOM_TYPE the atom type
      //! @param SINGLE_BONDS the number of single bonds for this atom
      //! @param DOUBLE_BONDS the number of double bonds for this atom
      //! @param TRIPLE_BONDS the number of triple bonds for this atom
      //! @param AROMATIC_BONDS the number of aromatic bonds for this atom
      //! @return return the estimated covalent radius of the atom, in angstroms
      static double GetAveCovalentRadius
      (
        const AtomType &ATOM_TYPE,
        const size_t &SINGLE_BONDS,
        const size_t &DOUBLE_BONDS,
        const size_t &TRIPLE_BONDS,
        const size_t &AROMATIC_BONDS
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
      //! @param OSTREAM output stream
      //! @param INDENT number of indentations
      //! @return ostream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class BondLengths

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_BOND_LENGTHS_H_
