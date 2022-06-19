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
#include "chemistry/bcl_chemistry_fragment_evolve_base.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_bond_isometry_handler.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "chemistry/bcl_chemistry_sample_conformations.h"
#include "chemistry/bcl_chemistry_stereocenters_handler.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  /////////////////
  // data access //
  /////////////////

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // helper functions //
  //////////////////////

    namespace
    {
      // helper struct to keept track of bad-bond info
      struct BadBond
      {
        // Data
        ElementType m_Element1;
        ElementType m_Element2;
        size_t m_BondOrder; // 0 is any, 1, 2, 3 do what you expect, 4 is aromatic
        bool m_InRing;

        //! @brief default constructor
        BadBond() :
          m_Element1(),
          m_Element2(),
          m_BondOrder(),
          m_InRing()
        {
        }

        //! @brief full constructor
        BadBond
        (
          const ElementType &ELEMENT1,
          const ElementType &ELEMENT2,
          const int &BOND_ORDER,
          const bool &IN_RING
        ) :
          m_Element1( ELEMENT1),
          m_Element2( ELEMENT2),
          m_BondOrder( BOND_ORDER),
          m_InRing( IN_RING)
        {
        }
      };
    }

    //! @brief returns the number of heteroatom-heteroatom (same element) non-ring bonds in a molecule
    //! @param MOLECULE the molecule to inspect
    //! @return the number of bad bonds in MOLECULE
    size_t FragmentEvolveBase::NumberBadBonds( const FragmentComplete &MOLECULE)
    {
      static std::vector< BadBond> bad_bond_info;
      if( bad_bond_info.empty())
      {
        // To oxygen
        bad_bond_info.push_back( BadBond( GetElementTypes().e_Oxygen, GetElementTypes().e_Oxygen, 0, false));
        bad_bond_info.push_back( BadBond( GetElementTypes().e_Oxygen, GetElementTypes().e_Nitrogen, 1, false));

        // To sulfur
        bad_bond_info.push_back( BadBond( GetElementTypes().e_Sulfur, GetElementTypes().e_Sulfur, 1, true));
        bad_bond_info.push_back( BadBond( GetElementTypes().e_Sulfur, GetElementTypes().e_Sulfur, 1, false));
        bad_bond_info.push_back( BadBond( GetElementTypes().e_Sulfur, GetElementTypes().e_Sulfur, 1, false));

        // To nitrogen
        bad_bond_info.push_back( BadBond( GetElementTypes().e_Nitrogen, GetElementTypes().e_Nitrogen, 2, false));
        bad_bond_info.push_back( BadBond( GetElementTypes().e_Nitrogen, GetElementTypes().e_Nitrogen, 1, false));

        // To halogens
        bad_bond_info.push_back( BadBond( GetElementTypes().e_Nitrogen, GetElementTypes().e_Fluorine, 0, false));
        bad_bond_info.push_back( BadBond( GetElementTypes().e_Nitrogen, GetElementTypes().e_Chlorine, 0, false));
        bad_bond_info.push_back( BadBond( GetElementTypes().e_Nitrogen, GetElementTypes().e_Bromine, 0, false));
        bad_bond_info.push_back( BadBond( GetElementTypes().e_Nitrogen, GetElementTypes().e_Iodine, 0, false));

        bad_bond_info.push_back( BadBond( GetElementTypes().e_Oxygen, GetElementTypes().e_Fluorine, 0, false));
        bad_bond_info.push_back( BadBond( GetElementTypes().e_Oxygen, GetElementTypes().e_Chlorine, 0, false));
        bad_bond_info.push_back( BadBond( GetElementTypes().e_Oxygen, GetElementTypes().e_Bromine, 0, false));
        bad_bond_info.push_back( BadBond( GetElementTypes().e_Oxygen, GetElementTypes().e_Iodine, 0, false));

        bad_bond_info.push_back( BadBond( GetElementTypes().e_Sulfur, GetElementTypes().e_Fluorine, 0, false));
        bad_bond_info.push_back( BadBond( GetElementTypes().e_Sulfur, GetElementTypes().e_Chlorine, 0, false));
        bad_bond_info.push_back( BadBond( GetElementTypes().e_Sulfur, GetElementTypes().e_Bromine, 0, false));
        bad_bond_info.push_back( BadBond( GetElementTypes().e_Sulfur, GetElementTypes().e_Iodine, 0, false));
      }

      size_t bad_bonds( 0);
      iterate::Generic< const AtomConformationalInterface> itr_atom( MOLECULE.GetAtomsIterator());
      for( ; itr_atom.NotAtEnd(); ++itr_atom)
      {
        const ElementType &element( itr_atom->GetElementType());
        if
        (
          element == GetElementTypes().e_Carbon
          || element == GetElementTypes().e_Hydrogen
        )
        {
          continue;
        }
        const storage::Vector< BondConformational> &atom_bonds( itr_atom->GetBonds());
        for( size_t i( 0); i < atom_bonds.GetSize(); ++i)
        {
          const ConfigurationalBondType &bond_type( atom_bonds( i).GetBondType());
          const ElementType &neighbor_element( atom_bonds( i).GetTargetAtom().GetElementType());

          // Iterate through each bad bond and see if the bond matches
          for( size_t b = 0; b < bad_bond_info.size(); ++b)
          {
            const BadBond &bad_bond( bad_bond_info[ b]);
            if
            (
              ( ( element == bad_bond.m_Element1 && neighbor_element == bad_bond.m_Element2)
                || ( element == bad_bond.m_Element2 && neighbor_element == bad_bond.m_Element1))
              && bond_type->IsBondInRing() == bad_bond.m_InRing
            )
            {
              if
              (
                bad_bond.m_BondOrder == 0
                || ( bond_type->GetConjugation() == ConstitutionalBondTypeData::e_Aromatic && bad_bond.m_BondOrder == 4)
                || ( bond_type->GetNumberOfElectrons() == bad_bond.m_BondOrder * 2)
              )
              {
                ++bad_bonds;
              }
            }
          }
        }
      }
      return bad_bonds;
    }

    //! @brief generates the 3d conformation from the 2d conformation by doing a system call to corina
    //! @param MOLECULE the small molecule for which 3d coordinates will be generated
    //! @return the 3d generated small molecule or empty molecule if corina was unavailable
    util::ShPtr< FragmentComplete> FragmentEvolveBase::GetCorina3DCoordinates( const ConformationInterface &MOLECULE)
    {

#if !defined(__MINGW32__) && !defined(__WIN32__) && !defined(_WIN32)

      BCL_MessageVrb( "GetCorina3DCoordinates");

      if( MOLECULE.GetNumberAtoms() == 0)
      {
        BCL_MessageVrb( "GetCorina3DCoordinates: empty molecule as input");
        return util::ShPtr< FragmentComplete>();
      }

      time_t cur_time;
      time( &cur_time);
      std::string file_basename( util::Format()( &MOLECULE) + util::Format()( cur_time));

      io::DirectoryEntry filename2d( "/tmp/" + file_basename + "_gen2D.sdf");
      io::DirectoryEntry filename3d( "/tmp/" + file_basename + "_gen3D.sdf");

      io::OFStream out;
      if( !io::File::TryOpenOFStream( out, filename2d.GetFullName()))
      {
        BCL_MessageCrt( "unable to open " + filename2d.GetFullName() + " for writing");
        return util::ShPtr< FragmentComplete>();
      }
      MOLECULE.WriteMDL( out);
      io::File::CloseClearFStream( out);
      // generate meaningful coordinates for the ensemble of generated constitutions using corina until somebody implements code to do that into the bcl!!
      const std::string command( "corina -dwh " + filename2d.GetFullName() + " " + filename3d.GetFullName());
      const int error( system( command.c_str()));

      if( error != 0)
      {
        BCL_MessageCrt( "unable to execute command: " + command + " with error: " + util::Format()( error));
      }

      io::IFStream input;
      io::File::MustOpenIFStream( input, filename3d.GetFullName());
      util::ShPtr< FragmentComplete> fragment
      (
        new FragmentComplete( sdf::FragmentFactory::MakeFragment( input, sdf::e_Saturate))
      );
      io::File::CloseClearFStream( input);
      filename2d.Remove();
      filename3d.Remove();
      if( fragment.IsDefined() && !fragment->HasBadGeometry() && !fragment->HasNonGasteigerAtomTypes())
      {
        return fragment;
      }

      return util::ShPtr< FragmentComplete>();

#else // !defined(__MINGW32__) && !defined(__WIN32__) && !defined(_WIN32)
      return util::ShPtr< FragmentComplete>();
#endif

    }

    // Create 3D conformer of input molecule using BCL::Conf
    util::ShPtr< FragmentComplete> FragmentEvolveBase::MakeBCLConformer( const FragmentComplete &MOLECULE)
    {
      // Construct a SampleConformers object designed to select a single reasonably good 3D conformer
      static RotamerLibraryFile rotamer_library_file;
      static SampleConformations sample_confs
      (
        rotamer_library_file,
        "",
        0.25,  // tolerance
        1,     // number of conformations
        10, // number of iterations
        false,  // change chirality
        0.0,   // random dihedral change weight
        false,  // generate 3d
        0.10   // clash tolerance, auto-adjusted if necessary due to clashes
      );

      // Sample conformers
      FragmentEnsemble mol_ens( sample_confs( MOLECULE).First());

      //If we made conformers then return the first mol in the ensemble, otherwise return undefined mol
      util::ShPtr< FragmentComplete> molecule_3d;
      if( !mol_ens.IsEmpty())
      {
        util::ShPtr< FragmentComplete> new_mol( new FragmentComplete( *mol_ens.Begin()));
        molecule_3d = new_mol;
      }
      return molecule_3d;
    }

    //! @brief Finalizes a molecule by running it through the atom standardizer and getting a 3D conformation
    //! @param MOLECULE the molecule to finalize
    //! @param CORINA generate a 3D conformer with corina (requires system call to external program)
    //! @return a new FragmentComplete that has been cleaned
    util::ShPtr< FragmentComplete> FragmentEvolveBase::FinalizeMolecule( const FragmentComplete &MOLECULE, const bool CORINA)
    {

      if( !MOLECULE.GetNumberAtoms())
      {
        return util::ShPtr< FragmentComplete>();
      }

      AtomVector< AtomComplete> atom_vector( MOLECULE.GetAtomVector());
      const std::string mol_name;
      AtomsCompleteStandardizer atoms_standardizer( atom_vector, mol_name, true);

      BondIsometryHandler::AddIsometryInformation( atom_vector, true);
      StereocentersHandler::AddChiralityFromConformation( atom_vector);

      // Generate molecule with no defined bond connectivity
      FragmentComplete new_fragment( atom_vector, "");

      // Return the 3D conformer
      util::ShPtr< FragmentComplete> molecule_3D;

      // Generate 3D conformer
      CORINA ?
      molecule_3D = GetCorina3DCoordinates( new_fragment) :
      molecule_3D = MakeBCLConformer( new_fragment);

      if( molecule_3D.IsDefined())
      {
        if( molecule_3D->HasNonGasteigerAtomTypes())
        {
          BCL_MessageVrb( "FinalizeMolecule: Molecule has non-gasteiger atom types, cannot use it");
          return util::ShPtr< FragmentComplete>();
        }
        molecule_3D->SetName( MOLECULE.GetName());
        molecule_3D->GetStoredPropertiesNonConst() = MOLECULE.GetStoredProperties();
      }
      else
      {
        BCL_MessageVrb( "FinalizeMolecule: could not finalize molecule");
      }

      return molecule_3D;
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace chemistry
} // namespace bcl
