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
#include <chemistry/bcl_chemistry_atoms_complete_standardizer.h>
#include <chemistry/bcl_chemistry_bond_isometry_handler.h>
#include <chemistry/bcl_chemistry_stereocenters_handler.h>
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "chemistry/bcl_chemistry_rdkit_mol_utils.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief virtual copy constructor
    RdkitMolUtils *RdkitMolUtils::Clone() const
    {
      return new RdkitMolUtils( *this);
    }

    //! @brief returns class name
    //! the class name as const ref std::string
    const std::string &RdkitMolUtils::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }


  ////////////////
  // operations //
  ////////////////


  ///////////////
  // operators //
  ///////////////


    //! @brief converts an RDKit ROMol into a BCL FragmentComplete
    //! @param RDKIT_MOL the RDKit ROMol to be converted
    //! @return the FragmentComplete constructed from the RDKit ROMol
    std::shared_ptr< FragmentComplete> RdkitMolUtils::RDKitROMolToFragmentComplete
    (
      const ::RDKit::ROMol &RDKIT_MOL,
      const std::string &MOL_ID
    )
    {
      ::RDKit::RWMol editable_mol( RDKIT_MOL);
      return RDKitRWMolToFragmentComplete( editable_mol, MOL_ID);
    }

    std::shared_ptr< FragmentComplete> RdkitMolUtils::RDKitRWMolToFragmentComplete
    (
      const ::RDKit::RWMol &RDKIT_MOL,
      const std::string &MOL_ID
    )
    {
      // initialize data with which to construct our AtomVector
      storage::Vector< sdf::AtomInfo> bcl_atoms;
      storage::Vector< sdf::BondInfo> bcl_bonds;

      // retrieve coordinates of all atoms
      const std::vector<RDGeom::Point3D> RDKIT_MOL_positions
      (
        RDKIT_MOL.getConformer().getPositions().size() ?
            RDKIT_MOL.getConformer().getPositions() :
            std::vector<RDGeom::Point3D>( RDKIT_MOL.getNumAtoms(), RDGeom::Point3D( 0.0, 0.0, 0.0))
      );

      // retrieve ring information for atoms and bonds in our rdkit molecule
      std::shared_ptr< ::RDKit::RingInfo> rdkit_ringinfo( RDKIT_MOL.getRingInfo());

      // atoms
      size_t atom_index( 0);
      for
      (
          auto atom_itr( RDKIT_MOL.beginAtoms()), atom_itr_end( RDKIT_MOL.endAtoms());
          atom_itr != atom_itr_end;
          ++atom_itr, ++atom_index
      )
      {
        // determine atom type
        const int formal_charge( ( *atom_itr)->getFormalCharge());
        const std::string element( (*atom_itr)->getSymbol());
        const ElementType element_type( GetElementTypes().ElementTypeLookup( element));
        const AtomType atom_type( AtomTypes::GetAtomType( element_type, formal_charge));

        // retrieve coordinates of current atom
        const RDGeom::Point3D &position( RDKIT_MOL_positions[ atom_index]);
        const linal::Vector3D atom_coords(position.x, position.y, position.z);

        // save the AtomInfo object that will correspond to the current RDKit atom
        bcl_atoms.PushBack
        (
          sdf::AtomInfo( atom_type, e_UnknownChirality, atom_coords, true)
        );
      }

      // bonds
      size_t bond_index( 0);
      for
      (
          auto bond_itr( RDKIT_MOL.beginBonds()), bond_itr_end( RDKIT_MOL.endBonds());
          bond_itr != bond_itr_end;
          ++bond_itr, ++bond_index
      )
      {
        // the two atoms involved in the bond
        const size_t begin_atom_i( ( *bond_itr)->getBeginAtomIdx());
        const size_t end_atom_i( ( *bond_itr)->getEndAtomIdx());

        // check if both atoms are in a ring
        const bool ring_bond
        (
          !rdkit_ringinfo->atomMembers( begin_atom_i).empty() &&
          !rdkit_ringinfo->atomMembers( end_atom_i).empty() ?
              true : false
        );

        // conjugation and aromaticity
        const bool conjugated( (*bond_itr)->getIsConjugated());
        const bool aromatic( (*bond_itr)->getIsAromatic());

        // bond order as a simple numeric value
        const double bond_type( ( **bond_itr).getBondTypeAsDouble());

        // put it all together to determine the bond type
        chemistry::ConfigurationalBondType bcl_bond_type;

        // priority to aromaticity
        if( aromatic)
        {
          // standard bond orders
          if(  bond_type == double( 1.0))
          {
            bcl_bond_type = GetConfigurationalBondTypes().e_AromaticSingleBond;
          }
          else if( bond_type == double( 2.0))
          {
            bcl_bond_type = GetConfigurationalBondTypes().e_AromaticDoubleBond;
          }
          else if( bond_type == double( 3.0))
          {
            bcl_bond_type = GetConfigurationalBondTypes().e_AromaticTripleBond;
          }
          // catch all for remaining bond orders
          else
          {
            bcl_bond_type = GetConfigurationalBondTypes().e_AromaticBond;
          }
        }
        // second priority to conjugation, including ring designation
        else if( conjugated)
        {
          // standard bond orders
          if(  bond_type == double( 1.0))
          {
            ring_bond ?
                bcl_bond_type = GetConfigurationalBondTypes().e_ConjugatedSingleBondInRing :
                bcl_bond_type = GetConfigurationalBondTypes().e_ConjugatedSingleBond;
          }
          else if( bond_type == double( 2.0))
          {
            ring_bond ?
                bcl_bond_type = GetConfigurationalBondTypes().e_ConjugatedDoubleBondInRing :
                bcl_bond_type = GetConfigurationalBondTypes().e_ConjugatedDoubleBond;
          }
          else if( bond_type == double( 3.0))
          {
            ring_bond ?
                bcl_bond_type = GetConfigurationalBondTypes().e_ConjugatedTripleBondInRing :
                bcl_bond_type = GetConfigurationalBondTypes().e_ConjugatedTripleBond;
          }
          // catch all for remaining bond orders
          else
          {
            ring_bond ?
                bcl_bond_type = GetConfigurationalBondTypes().e_ConjugatedBondInRing :
                bcl_bond_type = GetConfigurationalBondTypes().e_ConjugatedBond;
          }
        }
        // non-aromatic, non-conjugated bonds
        else if( bond_type == double( 1.0))
        {
          ring_bond ?
              bcl_bond_type = GetConfigurationalBondTypes().e_NonConjugatedSingleBond :
              bcl_bond_type = GetConfigurationalBondTypes().e_NonConjugatedSingleBondInRing;
        }
        // undefined bond types otherwise
        else
        {
          bcl_bond_type = GetConfigurationalBondTypes().e_Undefined;
        }

        // note - no amide bond explicitly defined in RDKit

        // save the bond object for these two atoms
        bcl_bonds.PushBack
        (
          sdf::BondInfo
          (
            begin_atom_i < end_atom_i ? begin_atom_i : end_atom_i,
            begin_atom_i > end_atom_i ? begin_atom_i : end_atom_i,
            bcl_bond_type
          )
        );
      }

      // create our atom vector
      // build our complete atomvector from the atom and bond info vectors
      AtomVector< AtomComplete> atoms( bcl_atoms, bcl_bonds);

      // determine atom and bond types (fix all of our errors, e.g., amide bonds)
      AtomsCompleteStandardizer standardizer( atoms, MOL_ID, true); // force recalculation

      // add bond stereocenter/isometry info then build full molecule
      BondIsometryHandler::AddIsometryInformation( atoms, true);
      StereocentersHandler::AddChiralityFromConformation( atoms);
      return std::make_shared< FragmentComplete>( FragmentComplete( atoms, MOL_ID));
    }


    //! @brief converts a BCL FragmentComplete into an RDKit ROMol
    //! @param MOL the BCL FragmentComplete to be converted
    //! @param MOL_ID the molecule identifier on the new molecule
    //! @return the RDKit ROMol constructed from the FragmentComplete
    std::shared_ptr< ::RDKit::ROMol> RdkitMolUtils::FragmentCompleteToRDKitROMol
    (
      const FragmentComplete &MOL,
      const std::string &MOL_ID
    )
    {










      return std::shared_ptr< ::RDKit::ROMol>();
    }

    //! @brief converts a BCL FragmentComplete into an RDKit RWMol
    //! @param MOL the BCL FragmentComplete to be converted
    //! @param MOL_ID the molecule identifier on the new molecule
    //! @return the RDKit RWMol constructed from the FragmentComplete
    std::shared_ptr< ::RDKit::RWMol> RdkitMolUtils::FragmentCompleteToRDKitRWMol
    (
      const FragmentComplete &MOL,
      const std::string &MOL_ID
    )
    {
      std::shared_ptr< ::RDKit::ROMol> rdkit_mol( FragmentCompleteToRDKitROMol( MOL, MOL_ID));
      return std::make_shared< ::RDKit::RWMol>( ::RDKit::RWMol( *rdkit_mol));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RdkitMolUtils::Read( std::istream &ISTREAM)
    {
      // read member
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &RdkitMolUtils::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace chemistry
} // namespace bcl

