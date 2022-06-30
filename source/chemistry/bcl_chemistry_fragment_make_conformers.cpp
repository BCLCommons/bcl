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
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_classes.h"
#include "biol/bcl_biol_aa_side_chain_factory.h"
#include "chemistry/bcl_chemistry_aa_fragment_complete.h"
#include "chemistry/bcl_chemistry_conformation_interface.h"
#include "chemistry/bcl_chemistry_fragment_evolve_base.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "pdb/bcl_pdb_factory.h"
#include "quality/bcl_quality_rmsd.h"
#include "random/bcl_random_distribution_interface.h"
#include "sdf/bcl_sdf_fragment_factory.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_bond_isometry_handler.h"
#include "chemistry/bcl_chemistry_bond_lengths.h"
#include "chemistry/bcl_chemistry_conformation_comparison_psi_field.h"
#include "chemistry/bcl_chemistry_fragment_make_conformers.h"
#include "chemistry/bcl_chemistry_fragment_split_largest_component.h"
#include "chemistry/bcl_chemistry_fragment_split_rings_with_unsaturated_substituents.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "chemistry/bcl_chemistry_sample_conformations.h"
#include "chemistry/bcl_chemistry_stereocenters_handler.h"
#include "chemistry/bcl_chemistry_voxel_grid_atom.h"
#include "descriptor/bcl_descriptor_molecule_druglike.h"
#include "find/bcl_find_collector_interface.h"
#include "random/bcl_random_uniform_distribution.h"
// external includes - sorted alphabetically

#undef AddAtom
#undef ATOMS

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    // initialize static
    sched::Mutex &FragmentMakeConformers::GetMutex()
    {
      static sched::Mutex s_Mutex;
      return s_Mutex;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FragmentMakeConformers::FragmentMakeConformers() :
        m_RotamerLibrary( RotamerLibraryFile()),
        m_Comparer( "SymmetryRMSD"),
        m_ComparerTolerance( float( 0.25)),
        m_Cluster( true),
        m_NumberIterations( size_t( 2000)),
        m_ClashTolerance( float( 0.1)),
        m_Generate3D( true),
        m_Corina( false),
        m_Local( false)
    {
    }

    //! @brief constructor
    FragmentMakeConformers::FragmentMakeConformers
    (
      const RotamerLibraryFile &ROT_LIB,
      const std::string &COMPARER,
      const float &COMPARER_TOLERANCE,
      const bool &CLUSTER,
      const size_t &N_ITERATIONS,
      const float &CLASH_TOLERANCE,
      const bool &GENERATE_3D,
      const bool &CORINA,
      const bool &LOCAL
    ) :
        m_RotamerLibrary( ROT_LIB),
        m_Comparer( COMPARER),
        m_ComparerTolerance( COMPARER_TOLERANCE),
        //m_Cluster( m_Cluster),
        m_Cluster( CLUSTER),
        m_NumberIterations( N_ITERATIONS),
        m_ClashTolerance( CLASH_TOLERANCE),
        m_Generate3D( GENERATE_3D),
        m_Corina( CORINA),
        m_Local( LOCAL)
    {
    }

    //! @brief constructor for quickly making single BCL 3D conformer
     FragmentMakeConformers::FragmentMakeConformers
     (
       const RotamerLibraryFile &ROT_LIB,
       const size_t &N_ITERATIONS,
       const float &CLASH_TOLERANCE,
       const bool &GENERATE_3D,
       const bool &CORINA
     ) :
         m_RotamerLibrary( ROT_LIB),
         m_Comparer( ""),
         m_ComparerTolerance( 0.0),
         m_Cluster( false),
         m_NumberIterations( N_ITERATIONS),
         m_ClashTolerance( CLASH_TOLERANCE),
         m_Generate3D( GENERATE_3D),
         m_Corina( CORINA),
         m_Local( false)
     {
     }

    //! @brief clone constructor
    FragmentMakeConformers *FragmentMakeConformers::Clone() const
    {
      return new FragmentMakeConformers( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FragmentMakeConformers::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

  ////////////////
  // operations //
  ////////////////

    //! @brief system call to corina to generate 3d conformers
    //! @param MOLECULE the small molecule for which 3d coordinates will be generated
    //! @return the 3d generated small molecule or empty molecule if corina was unavailable
    util::ShPtr< FragmentComplete> FragmentMakeConformers::GetCorina3DCoordinates
    (
      const ConformationInterface &MOLECULE
    )
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

    //! @brief generate a 3D conformer of a molecule
    //! @param MOLECULE the molecule for which to generate a conformer
    //! @return a pointer to a new 3D molecule
    util::ShPtr< FragmentComplete> FragmentMakeConformers::MakeConformer
    (
      const FragmentComplete &MOLECULE
    )
    {
      // Construct a SampleConformers object designed to select a single reasonably okay-ish 3D conformer
      static SampleConformations sample_confs_single
      (
        m_RotamerLibrary,       // rotamer library file
        "",                     // conformation comparer type
        0.0,                    // conformational comparer tolerance
        1,                      // number of conformations
        m_NumberIterations,     // number of iterations
        false,                  // no change chirality
        0.0,                    // random dihedral change weight
        m_Generate3D,           // generate 3d?
        m_ClashTolerance,       // clash tolerance
        false                   // no clustering
      );

      // Sample conformers
      FragmentEnsemble mol_ens;

      // make corina conformer if desired
      FragmentComplete molecule;
      if( m_Corina)
      {
        util::ShPtr< FragmentComplete> molecule_ptr( GetCorina3DCoordinates( MOLECULE));
        if( molecule_ptr.IsDefined())
        {
          molecule = *molecule_ptr;
        }
        else
        {
          return util::ShPtr< FragmentComplete>();
        }
      }
      else
      {
        mol_ens = sample_confs_single( MOLECULE).First();
        if( mol_ens.GetSize())
        {
          molecule = mol_ens.GetMolecules().FirstElement();
        }
      }

      // make sure we pass the sample indices to new molecule
      molecule.GetStoredPropertiesNonConst().SetMDLProperty( "SampleByParts", MOLECULE.GetMDLProperty( "SampleByParts"));

      // if we did not use corina and our ensemble is empty, or we do not have a conformer, then return empty pointer
      if
      (
          ( mol_ens.IsEmpty() && !m_Corina) ||
          !molecule.GetSize()
      )
      {
        return util::ShPtr< FragmentComplete>();
      }
      // if we used corina and our molecule is there, then we are good to go
      else if( molecule.GetSize() && m_Corina)
      {
        return util::ShPtr< FragmentComplete>( new FragmentComplete( molecule));
      }
      //If we made conformers then return the first mol in the ensemble
      else if( ( !mol_ens.IsEmpty() && !m_Corina) && molecule.GetSize())
      {
        return util::ShPtr< FragmentComplete>( new FragmentComplete( molecule));
      }

      // return empty pointer
      return util::ShPtr< FragmentComplete>();
    }

    //! @brief generate an ensemble of 3D conformations of a molecule
    //! @param MOLECULE the molecule for which to generate conformers
    //! @return a pointer to a new 3D molecule
    util::ShPtr< FragmentEnsemble> FragmentMakeConformers::MakeConformers
    (
      const FragmentComplete &MOLECULE
    )
    {
      // Make an optimal conformer ensemble
      static SampleConformations sample_confs
      (
        m_RotamerLibrary,                            // rotamer library file
        m_Comparer,                                  // conformation comparer type
        m_ComparerTolerance,                         // conformational comparer tolerance
        std::min( size_t( 250), m_NumberIterations), // number of conformations
        m_NumberIterations,                          // number of iterations
        false,                                       // no change chirality
        0.0,                                         // random dihedral change weight
        m_Generate3D,                                // generate 3d?
        m_ClashTolerance,                            // clash tolerance
        m_Cluster                                    // cluster?
      );

      // restrict to local if specified
      if( m_Local)
      {
        // only sample bond angles/lengths
        sample_confs.SetSamplingPreferences( false, false, true, false);
      }

      // Sample conformers
      FragmentEnsemble mol_ens( sample_confs( MOLECULE).First());

      // Return any valid conformations
      if( !mol_ens.IsEmpty())
      {
        return util::ShPtr< FragmentEnsemble>( new FragmentEnsemble( mol_ens));
      }
      return util::ShPtr< FragmentEnsemble>();
    }

    //! @brief generate a 3D conformer of a molecule
    //! @param MOLECULE the molecule for which to generate a conformer
    //! @return a pointer to a new 3D molecule
    util::ShPtr< FragmentComplete> FragmentMakeConformers::MakeConformerBCLIdeal
    (
      const FragmentComplete &MOLECULE
    )
    {
      // Construct a SampleConformers object designed to select a single reasonably okay-ish 3D conformer
      static SampleConformations sample_confs_single_ideal
      (
        m_RotamerLibrary,                            // rotamer library file
        m_Comparer,                                  // conformation comparer type
        m_ComparerTolerance,                         // conformational comparer tolerance
        std::min( size_t( 250), m_NumberIterations), // number of conformations
        m_NumberIterations,                          // number of iterations
        false,                                       // no change chirality
        0.0,                                         // random dihedral change weight
        m_Generate3D,                                // generate 3d?
        m_ClashTolerance,                            // clash tolerance
        m_Cluster                                    // cluster?
      );

      // we want to sample the local space deeply
      sample_confs_single_ideal.SetSamplingPreferences( false, false, true, false);

      // Sample conformers
      FragmentComplete native_mol( MOLECULE);
      native_mol.SaturateWithH();
      FragmentEnsemble mol_ens( sample_confs_single_ideal( native_mol).First());

      // Find nearest BCL conf to input/native
      FragmentComplete molecule;
      if( mol_ens.GetSize())
      {
        // Make comparer and output matrix
        static util::Implementation< ConformationComparisonInterface> rmsd_calculator( m_Comparer);

        // make input molecule into ensemble
        native_mol.RemoveH();

        // Compare the conformers
        for
        (
            auto mol_itr( mol_ens.Begin()), mol_itr_end( mol_ens.End());
            mol_itr != mol_itr_end;
            ++mol_itr
        )
        {
          mol_itr->RemoveH();
          std::string rmsd( util::Format()( ( *rmsd_calculator)( *mol_itr, native_mol)));
          mol_itr->GetStoredPropertiesNonConst().SetMDLProperty( m_Comparer + "_to_input", rmsd);
          mol_itr->SaturateWithH();
        }

        // Get the lowest RMSD to input
        mol_ens.Sort( m_Comparer + "_to_input");

        molecule = mol_ens.GetMolecules().FirstElement();
      }
      else
      {
        return util::ShPtr< FragmentComplete>();
      }

      // make sure we pass the sample indices to new molecule
      molecule.GetStoredPropertiesNonConst().SetMDLProperty( "SampleByParts", MOLECULE.GetMDLProperty( "SampleByParts"));

      //If we made conformers then return the first mol in the ensemble, otherwise return undefined mol
      if( !mol_ens.IsEmpty() && molecule.GetSize())
      {
        return util::ShPtr< FragmentComplete>( new FragmentComplete( molecule));
      }
      return util::ShPtr< FragmentComplete>();
    }

    //! @brief generate a 3D conformer of a molecule
    //! @param MOLECULE the molecule for which to generate a conformer
    //! @param ATOM_INDICES minimize geometric center distance of these atoms to target
    //! @param POSITIONS geometric center of these coordinates is the target
    //! @return a pointer to a new 3D molecule with the minimum target distance
    util::ShPtr< FragmentComplete> FragmentMakeConformers::MakeMinDistConformer
    (
      const FragmentComplete &MOLECULE,
      const storage::Vector< size_t> &ATOM_INDICES,
      const storage::Vector< linal::Vector3D> &POSITIONS
    )
    {
      // get centroid position
      size_t p_sz( POSITIONS.GetSize());
      double x( 0.0), y( 0.0), z( 0.0), sq_distance( 0.0);
      for( size_t p( 0); p < p_sz; ++p)
      {
        x += POSITIONS( p).X();
        y += POSITIONS( p).Y();
        z += POSITIONS( p).Z();
      }
      linal::Vector3D centroid( x / double( p_sz), y / double( p_sz), z / double( p_sz));

      // Make an optimal conformer ensemble
      static SampleConformations sample_confs
      (
        m_RotamerLibrary,                            // rotamer library file
        m_Comparer,                                  // conformation comparer type
        m_ComparerTolerance,                         // conformational comparer tolerance
        std::min( size_t( 250), m_NumberIterations), // number of conformations
        m_NumberIterations,                          // number of iterations
        false,                                       // no change chirality
        0.0,                                         // random dihedral change weight
        m_Generate3D,                                // generate 3d?
        m_ClashTolerance,                            // clash tolerance
        m_Cluster                                    // cluster?
      );

      // restrict to local if specified
      if( m_Local)
      {
        // only sample bond angles/lengths
        sample_confs.SetSamplingPreferences( false, false, true, false);
      }

      // Sample conformers
//      BCL_Debug( MOLECULE.GetMDLProperty( "SampleByParts"));
      FragmentComplete mol( MOLECULE);
      mol.SaturateWithH();
//      BCL_Debug( mol.GetMDLProperty( "SampleByParts"));
      mol.GetStoredPropertiesNonConst().SetMDLProperty( "SampleByParts", MOLECULE.GetMDLProperty( "SampleByParts"));
//      BCL_Debug( mol.GetMDLProperty( "SampleByParts"));
      FragmentEnsemble mol_ens( sample_confs( mol).First());
//      BCL_Debug( centroid);
      if( mol_ens.GetSize())
      {
        // compute mean square distance to centroid for each conformer
        size_t a_sz( ATOM_INDICES.GetSize());
        for
        (
            auto conf_itr( mol_ens.Begin()), conf_itr_end( mol_ens.End());
            conf_itr != conf_itr_end;
            ++conf_itr
        )
        {
          math::RunningAverage< double> ave_sq_dist;
          for( size_t a( 0); a < a_sz; ++a)
          {
            sq_distance +=
                linal::SquareDistance
                (
                  conf_itr->GetAtomVector()( ATOM_INDICES( a)).GetPosition(),
                  centroid
                );
          }
          sq_distance /= double( a_sz);
          conf_itr->GetStoredPropertiesNonConst().SetMDLProperty( "SqrDistTarget", linal::Vector< double>( size_t( 1), sq_distance));
        }
        // Get the lowest distance to target
        mol_ens.Sort( "SqrDistTarget");
//        BCL_MessageStd( "SqrDistance Shortest: " + util::Format()( mol_ens.GetMolecules().FirstElement().GetMDLProperty( "SqrDistTarget")));
//        BCL_MessageStd( "SqrDistance Longest: " + util::Format()( mol_ens.GetMolecules().LastElement().GetMDLProperty( "SqrDistTarget")));
        FragmentComplete molecule( mol_ens.GetMolecules().FirstElement());
        molecule.GetStoredPropertiesNonConst().SetMDLProperty( "SqrDistTarget", molecule.GetMDLProperty( "SqrDistTarget"));
//        BCL_Debug( molecule.GetMDLProperty( "SqrDistTarget"));

        // make sure we pass the sample indices to new molecule
        molecule.GetStoredPropertiesNonConst().SetMDLProperty( "SampleByParts", MOLECULE.GetMDLProperty( "SampleByParts"));

        // done
        if( molecule.GetSize())
        {
          return util::ShPtr< FragmentComplete>( new FragmentComplete( molecule));
        }
        return util::ShPtr< FragmentComplete>();
      }
      else
      {
        return util::ShPtr< FragmentComplete>();
      }
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FragmentMakeConformers::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &FragmentMakeConformers::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
