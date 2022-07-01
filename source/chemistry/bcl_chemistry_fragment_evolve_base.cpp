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
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
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

    //! @brief evaluates a molecule's topology to see if there are druglikeness violations
    //! @param MOLECULE the molecule to inspect
    //! @return false if the molecule fails any of the druglikeness checks in MoleculeDruglike; true otherwise
    bool FragmentEvolveBase::IsConstitutionDruglike( const FragmentComplete &MOLECULE)
    {
      return descriptor::GetCheminfoProperties().calc_IsMolDruglike->SumOverObject( MOLECULE)( 0);
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

    //! @brief determines what fragments would result from breaking a bond in a graph
    //! @param MOLECULE_GRAPH the graph that will have its bond broken
    //! @param FROM one vertex that makes up the bond to break
    //! @param TO the other vertex
    //! @return a list of vectors of indices which correspond to connected components of the graph
    storage::List< storage::Vector< size_t> > FragmentEvolveBase::CollectFragmentsFromBondBreakage
    (
      graph::ConstGraph< size_t, size_t> &MOLECULE_GRAPH,
      const size_t &FROM,
      const size_t &TO
    )
    {
      if( FROM >= MOLECULE_GRAPH.GetSize() || TO >= MOLECULE_GRAPH.GetSize() || FROM == TO)
      {
        return storage::List< storage::Vector< size_t> >();
      }

      // Save the bond info
      size_t bond_info( MOLECULE_GRAPH.GetEdgeData( FROM, TO));

      // Break the bond
      MOLECULE_GRAPH.RemoveEdge( FROM, TO);

      // Get the pieces of the graph
      storage::List< storage::Vector< size_t> > components( graph::Connectivity::GetComponents( MOLECULE_GRAPH));

      // Restore the bond
      MOLECULE_GRAPH.AddEdge( FROM, TO, bond_info);

      return components;
    }

    //! @brief determines what fragments would result from breaking a bond in a graph
    //! @param MOLECULE the molecule that will have a bond broken
    //! @param MOLECULE_GRAPH the graph MOLECULE
    //! @return a list of vectors of indices which correspond to connected components of the graph
    storage::List< storage::Vector< size_t> > FragmentEvolveBase::FragmentsFromRandomBondBreakage
    (
      const FragmentComplete &MOLECULE,
      graph::ConstGraph< size_t, size_t> &MOLECULE_GRAPH,
      const size_t &EDGE_TYPE
    )
    {
      // Make sure everything matches
      if( MOLECULE_GRAPH.GetSize() == 0 || MOLECULE_GRAPH.GetSize() != MOLECULE.GetNumberAtoms())
      {
        return storage::List< storage::Vector< size_t> >();
      }

      // Get a list of bonds of the molecule
      storage::Vector< sdf::BondInfo> bonds( MOLECULE.GetBondInfo());

      // Determine which ones can be broken; don't break ring bonds
      storage::Vector< size_t> available_bonds;
      available_bonds.AllocateMemory( bonds.GetSize());

      for( size_t pos( 0), end( bonds.GetSize()); pos < end; ++pos)
      {
        // Check to make sure the edge isn't in a ring, and it matches the edge type
        if
        (
          !bonds( pos).GetConstitutionalBondType()->IsBondInRing()
          && MOLECULE_GRAPH.GetEdgeData( bonds( pos).GetAtomIndexLow(), bonds( pos).GetAtomIndexHigh()) == EDGE_TYPE
        )
        {
          available_bonds.PushBack( pos);
        }
      }

      if( !available_bonds.GetSize())
      {
        return storage::List< storage::Vector< size_t> >();
      }

      size_t which_bond( random::GetGlobalRandom().Random< size_t>( available_bonds.GetSize() - 1));

      return CollectFragmentsFromBondBreakage
          (
            MOLECULE_GRAPH,
            bonds( available_bonds( which_bond)).GetAtomIndexLow(),
            bonds( available_bonds( which_bond)).GetAtomIndexHigh()
          );
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace chemistry
} // namespace bcl
