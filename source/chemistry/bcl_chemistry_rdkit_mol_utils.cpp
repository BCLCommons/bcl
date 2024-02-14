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
#include "chemistry/bcl_chemistry_rdkit_mol_utils.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_complete.h"
#include "chemistry/bcl_chemistry_atom_vector.h"
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_bond_isometry_handler.h"
#include "chemistry/bcl_chemistry_stereocenters_handler.h"

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
      const std::string &MOL_ID,
      const bool ALLOW_UNDEFINED_ATOMS
    )
    {
      ::RDKit::RWMol editable_mol( RDKIT_MOL);
      return RDKitRWMolToFragmentComplete( editable_mol, MOL_ID, ALLOW_UNDEFINED_ATOMS);
    }

    std::shared_ptr< FragmentComplete> RdkitMolUtils::RDKitRWMolToFragmentComplete
    (
      const ::RDKit::RWMol &RDKIT_MOL,
      const std::string &MOL_ID,
      const bool ALLOW_UNDEFINED_ATOMS
    )
    {
      // initialize data with which to construct our AtomVector
      storage::Vector< sdf::AtomInfo> bcl_atoms;
      storage::Vector< sdf::BondInfo> bcl_bonds;

      // retrieve coordinates of all atoms
      if( !RDKIT_MOL.getNumConformers())
      {
        BCL_MessageStd
        (
          "[WARNING] RdkitMolUtils::RDKitRWMolToFragmentComplete "
          "no 3D conformer stored on RDKit molecule. Setting all atom coordinates to origin."
        );
      }
      const std::vector< RDGeom::Point3D> rdkit_mol_positions
      (
        RDKIT_MOL.getNumConformers() ?
          RDKIT_MOL.getConformer().getPositions() :
          std::vector< RDGeom::Point3D>( RDKIT_MOL.getNumAtoms(), RDGeom::Point3D( 0.0, 0.0, 0.0))
      );

      // retrieve ring information for atoms and bonds in our rdkit molecule
      ::RDKit::RingInfo *const rdkit_ringinfo( RDKIT_MOL.getRingInfo());

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
        const std::string element( ( *atom_itr)->getSymbol());
        const ElementType element_type( GetElementTypes().ElementTypeLookup( element));
        const AtomType atom_type( AtomTypes::GetAtomType( element_type, formal_charge));

        // alert the user to undefined atom types and return null if not allowed
        if( atom_type.GetName() == GetAtomTypes().e_Undefined.GetName())
        {
          BCL_MessageVrb
          (
            "Undefined atom type conversion at index "
            + util::Format()( atom_index)
            + " with element symbol " + util::Format()( element)
            + " during creation of BCL molecule from RDKit molecule."
          );
          if( !ALLOW_UNDEFINED_ATOMS)
          {
            return std::shared_ptr< FragmentComplete>( nullptr);
          }
        }

        // retrieve coordinates of current atom
        const RDGeom::Point3D &position( rdkit_mol_positions[ atom_index]);
        const linal::Vector3D atom_coords( position.x, position.y, position.z);

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
        const bool conjugated( ( *bond_itr)->getIsConjugated());
        const bool aromatic( ( *bond_itr)->getIsAromatic());

        // bond order as a simple numeric value
        const double bond_type( ( **bond_itr).getBondTypeAsDouble());

        // put it all together to determine the bond type
        ConfigurationalBondType bcl_bond_type;

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
          if( bond_type == double( 1.0))
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

        // note - no amide bond explicitly defined in RDKit (special case of conjugated single bond)

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
      standardizer.SetConjugationOfBondTypes( atoms);

      // add bond stereocenter/isometry info then build full molecule
      BondIsometryHandler::AddIsometryInformation( atoms, true);
      StereocentersHandler::AddChiralityFromConformation( atoms);
      return std::make_shared< FragmentComplete>( FragmentComplete( atoms, MOL_ID));
    }

    //! @brief converts a BCL FragmentComplete into an RDKit RWMol
    //! @param MOL the BCL FragmentComplete to be converted
    //! @param UNDEFINED_SUBSTITUTE if the BCL atom is undefined (X), replace with this
    //! element type in the RDKit molecule; default is to leave undefined, which
    //! will skip this atom when building the new molecule
    //! @param SANITIZE perform RDKit-style standardization of molecule after an initial
    //! molecule is created from the BCL FragmentComplete
    //! @param SKIP_COORDS do not try to assign 3D coordinates from BCL object as a 3D conformer
    //! on the RDKit molecule; this means that we are effectively creating a molecule topology
    //! object from our BCL molecule; useful when just converting from FragmentComplete to SMILES
    //! @return the RDKit RWMol constructed from the FragmentComplete
    std::shared_ptr< ::RDKit::RWMol> RdkitMolUtils::FragmentCompleteToRDKitRWMol
    (
      const FragmentComplete &MOL,
      const ElementType &UNDEFINED_SUBSTITUTE,
      const bool SANITIZE,
      const bool SKIP_COORDS
    )
    {
      // initialize an empty RDKit molecule
      std::shared_ptr< ::RDKit::RWMol> rdkit_mol( new ::RDKit::RWMol());

      // RDKit atom from atomic number
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atoms( MOL.GetAtomsIterator());
        itr_atoms.NotAtEnd();
        ++itr_atoms
      )
      {
        // handle undefined elements in our BCL molecule by choosing to either skip them or replace them
        ElementType element_type( itr_atoms->GetElementType());
        if( itr_atoms->GetElementType() == GetElementTypes().e_Undefined)
        {
          // by default we skip undefined element types
          if( UNDEFINED_SUBSTITUTE == GetElementTypes().e_Undefined)
          {
            continue;
          }
          // but we can also substitute it for a dummy atom
          else
          {
            element_type = UNDEFINED_SUBSTITUTE;
          }
        }

        // RDKit requires that this is a raw pointer (non-owning for memory safety)
        ::RDKit::Atom *rdkit_atom_p( new ::RDKit::Atom( element_type->GetAtomicNumber()));
        rdkit_atom_p->setFormalCharge( itr_atoms->GetAtomType()->GetFormalCharge());
        rdkit_mol->addAtom
        (
          rdkit_atom_p, // raw pointer to atom
          false,        // do not set new atom to be the active atom
          true          // transfer ownership (very important that this is true with the owning pointer)
        );
      }

      // bonds
      const storage::Vector< sdf::BondInfo> bonds( MOL.GetAtomVector().GetBondInfo());
      for
      (
        auto bond_itr( bonds.Begin()), bond_itr_end( bonds.End());
        bond_itr != bond_itr_end;
        ++bond_itr
      )
      {
        // prepare to set bond type
        ::RDKit::Bond *rdkit_bond_p( new ::RDKit::Bond());

        // set the atoms involved in the bond
        rdkit_bond_p->setBeginAtomIdx( bond_itr->GetAtomIndexLow());
        rdkit_bond_p->setEndAtomIdx( bond_itr->GetAtomIndexHigh());

        // aromatic bonds
        if
        (
            bond_itr->GetConfigurationalBondType() == GetConfigurationalBondTypes().e_AromaticBond ||
            bond_itr->GetConfigurationalBondType() == GetConfigurationalBondTypes().e_AromaticSingleBond ||
            bond_itr->GetConfigurationalBondType() == GetConfigurationalBondTypes().e_AromaticDoubleBond ||
            bond_itr->GetConfigurationalBondType() == GetConfigurationalBondTypes().e_AromaticTripleBond
        )
        {
          rdkit_bond_p->setBondType( ::RDKit::Bond::BondType::AROMATIC);
          rdkit_bond_p->setIsAromatic( true);
          rdkit_bond_p->setIsConjugated( true);
        }

        // single bonds that contribute to conjugated systems
        else if
        (
            bond_itr->GetConfigurationalBondType() == GetConfigurationalBondTypes().e_AmideSingleBond ||
            bond_itr->GetConfigurationalBondType() == GetConfigurationalBondTypes().e_ConjugatedSingleBond ||
            bond_itr->GetConfigurationalBondType() == GetConfigurationalBondTypes().e_ConjugatedSingleBondInRing
        )
        {
          rdkit_bond_p->setBondType( ::RDKit::Bond::BondType::SINGLE);
          rdkit_bond_p->setIsAromatic( false);
          rdkit_bond_p->setIsConjugated( true);
        }

        // single bonds that are not part of a conjugated system
        else if
        (
            bond_itr->GetConfigurationalBondType() == GetConfigurationalBondTypes().e_NonConjugatedSingleBond ||
            bond_itr->GetConfigurationalBondType() == GetConfigurationalBondTypes().e_NonConjugatedSingleBondInRing
        )
        {
          rdkit_bond_p->setBondType( ::RDKit::Bond::BondType::SINGLE);
          rdkit_bond_p->setIsAromatic( false);
          rdkit_bond_p->setIsConjugated( false);
        }

        // double bonds
        else if
        (
            bond_itr->GetConfigurationalBondType() == GetConfigurationalBondTypes().e_ConjugatedDoubleBond ||
            bond_itr->GetConfigurationalBondType() == GetConfigurationalBondTypes().e_ConjugatedDoubleBond_E ||
            bond_itr->GetConfigurationalBondType() == GetConfigurationalBondTypes().e_ConjugatedDoubleBond_X ||
            bond_itr->GetConfigurationalBondType() == GetConfigurationalBondTypes().e_ConjugatedDoubleBond_Z ||
            bond_itr->GetConfigurationalBondType() == GetConfigurationalBondTypes().e_ConjugatedDoubleBondInRing ||
            bond_itr->GetConfigurationalBondType() == GetConfigurationalBondTypes().e_ConjugatedDoubleBondInRing_E ||
            bond_itr->GetConfigurationalBondType() == GetConfigurationalBondTypes().e_ConjugatedDoubleBondInRing_X ||
            bond_itr->GetConfigurationalBondType() == GetConfigurationalBondTypes().e_ConjugatedDoubleBondInRing_Z
        )
        {
          rdkit_bond_p->setBondType( ::RDKit::Bond::BondType::DOUBLE);
          rdkit_bond_p->setIsAromatic( false);
          rdkit_bond_p->setIsConjugated( true);
        }

        // triple bonds
        else if
        (
            bond_itr->GetConfigurationalBondType() == GetConfigurationalBondTypes().e_ConjugatedTripleBond ||
            bond_itr->GetConfigurationalBondType() == GetConfigurationalBondTypes().e_ConjugatedTripleBondInRing
        )
        {
          rdkit_bond_p->setBondType( ::RDKit::Bond::BondType::TRIPLE);
          rdkit_bond_p->setIsAromatic( false);
          rdkit_bond_p->setIsConjugated( true);
        }

        // give unspecified bonds to undefined and other bond types
        else
        {
          rdkit_bond_p->setBondType( ::RDKit::Bond::BondType::UNSPECIFIED);
        }

        // once the rdkit bond is setup, add it to the rdkit molecule
        rdkit_mol->addBond
        (
          rdkit_bond_p, // raw pointer to atom
          true          // transfer ownership (very important that this is true with the owning pointer)
        );
      }

      // rdkit standardization and return
      if( SANITIZE)
      {
        ::RDKit::MolOps::sanitizeMol( *rdkit_mol);
      }

      // currently the rdkit_mol is just a topology; this will add coordinates
      if( !SKIP_COORDS)
      {
        AddBCLConformersToRDKitRWMol( rdkit_mol, FragmentEnsemble( storage::List< FragmentComplete>( 1, MOL)));
      }
      return rdkit_mol;
    }

    //! @brief add coordinates from BCL conformers to RDKit molecule
    //! @details BCL and RDKit are organized differently with respect to
    //! storage of atomic coordinates. In the BCL, coordinates are a component
    //! of the AtomInfo object, which means each atom can independently have
    //! associated coordinate information. Each molecule contains both the
    //! topological and coordinate data, and thus a conformer ensemble is represented
    //! through a list of molecules (e.g., FragmentEnsemble). In contrast, RDKit
    //! molecule objects, such as RWMol, contain topology and property information,
    //! but no coordinate information. Conformers are stored on RDKit molecules as
    //! attributes, which means conformers need to be added to a molecule (rather than
    //! represented as individual molecules). This function will modify an RDKit molecule
    //! directly to add conformers based on an ensemble of BCL conformers.
    //! @param RDKIT_MOL the RDKit molecule to which conformers will be added
    //! @param MOL_ENS the BCL conformers that will be added to the RDKit molecule
    //! @param MAP_ATOMS if true, convert RDKit molecule to BCL molecule and perform
    //! substructure comparison to map atoms; default false assumes atoms are in the same
    //! order already, such as if the RDKit molecule was derived from the BCL molecule, or
    //! vice versa.
    //! @param ATOM_MAP user-supplied map between RDKit and BCL molecule atoms; key is the
    //! BCL atom index and value is the RDKit atom index
    void RdkitMolUtils::AddBCLConformersToRDKitRWMol
    (
      std::shared_ptr< ::RDKit::RWMol> RDKIT_MOL,
      const FragmentEnsemble &MOL_ENS,
      const bool MAP_ATOMS,
      const std::map< size_t, size_t> ATOM_MAP
    )
    {
      // verify that the molecules are the same-ish (dirty check; verifying size is the same)
      if( RDKIT_MOL->getNumAtoms( false) != MOL_ENS.GetMolecules().FirstElement().GetSize())
      {
        BCL_MessageStd
        (
          " Number of atoms in RDKit molecule (" + util::Format()( RDKIT_MOL->getNumAtoms( false)) + ")"
          " does not match number of atoms in BCL conformer (" + util::Format()( MOL_ENS.GetMolecules().FirstElement().GetSize()) + ");"
          " cannot add conformers to RDKit molecule."
        );
        return;
      }

      // substructure comparison to get matching atoms
      if( MAP_ATOMS)
      {
        // TODO: MCS return a SubGraph< size_t, size_t>
      }

      // add a conformer to RDKit molecule for each BCL conformer
      for
      (
          auto conf_itr( MOL_ENS.Begin()), conf_itr_end( MOL_ENS.End());
          conf_itr != conf_itr_end;
          ++conf_itr
      )
      {
        // initialize RDKit conformer with coordinates at origin for each atom
        ::RDKit::Conformer *rdkit_conformer_p( new ::RDKit::Conformer( RDKIT_MOL->getNumAtoms( false)));

        // assign an atom coordinate for each atom in BCL conformer
        for( size_t coord_index( 0), n_coords( conf_itr->GetAtomCoordinates().GetSize()); coord_index < n_coords; ++coord_index)
        {
          // default assumes that the atoms are in the same order between the two molecule types
          size_t rdkit_atom_id( coord_index);

          // optional - map relating atom sequences in molecules
          if( !ATOM_MAP.empty())
          {
            rdkit_atom_id = ATOM_MAP.find( coord_index)->second;
          }

          // set coordinate
          ::RDGeom::Point3D position
          (
            conf_itr->GetAtomCoordinates()( coord_index)->X(),
            conf_itr->GetAtomCoordinates()( coord_index)->Y(),
            conf_itr->GetAtomCoordinates()( coord_index)->Z()
          );
          rdkit_conformer_p->setAtomPos( rdkit_atom_id, position);
        }
        // add the conformer to our RDKit molecule
        RDKIT_MOL->addConformer( rdkit_conformer_p, true);
      }
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

