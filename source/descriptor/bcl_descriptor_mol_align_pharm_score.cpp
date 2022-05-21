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
#include "descriptor/bcl_descriptor_mol_align_pharm_score.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "chemistry/bcl_chemistry_aa_fragment_complete.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "chemistry/bcl_chemistry_sample_conformations.h"
#include "command/bcl_command_command_state.h"
#include "descriptor/bcl_descriptor_constants.h"
#include "descriptor/bcl_descriptor_molecule_3da_smooth_sign_occlusion_code.h"
#include "io/bcl_io_ifstream.h"
#include "iterate/bcl_iterate_reflecting.h"
#include "math/bcl_math_running_min_max.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> MolAlignPharmScore::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new MolAlignPharmScore()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MolAlignPharmScore::MolAlignPharmScore() :
        m_AlignmentType( "flexible"),
        m_MaxAtomDistance( 1.0)
    {
      if( !command::CommandState::IsInStaticInitialization())
      {
        m_Properties.PushBack( CheminfoProperty( Constants< chemistry::AtomConformationalInterface, float>( 1.0)));
        m_PropertyWeights = linal::Vector< float>( size_t( 1), 1.0);
        m_Aligner.AssertRead
        (
          util::ObjectDataLabel
          (
            "(rigid mol b = true,"
            "number rigid trajectories=3,"
            "number flexible trajectories=3,"
            "fraction filtered initially=0.75,"
            "fraction filtered iteratively=0.50,"
            "iterations=100,"
            "filter iterations=50,"
            "refinement iterations=50,"
            "align to scaffold = false,"
            "conformer pairs = 10,"
            "number outputs = 1,"
            "sample conformers=SampleConformations("
            "conformation_comparer=RMSD,"
            "generate_3D=0,tolerance=0.50,rotamer_library = cod,"
            "max_iterations=100,max_conformations=10,"
            "cluster=true,"
            "clash_tolerance=0.4,"
            "change_chirality=0))"
          )
        );

        BCL_Assert
        (
          ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
          "Failed to create " + GetClassIdentifier()
        );
      }
    }

    //! @brief virtual copy constructor
    //! @return pointer to new MolAlignPharmScore
    MolAlignPharmScore *MolAlignPharmScore::Clone() const
    {
      return new MolAlignPharmScore( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &MolAlignPharmScore::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &MolAlignPharmScore::GetAlias() const
    {
      static const std::string s_name( "MolAlignPharmScore");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void MolAlignPharmScore::Calculate( linal::VectorReference< float> &STORAGE)
    {
      // Get our molecule
      util::SiPtr< const chemistry::ConformationInterface> current_mol( this->GetCurrentObject());
      chemistry::FragmentComplete molecule( *current_mol);

      // Perform alignment
      if( m_AlignmentType == "rigid")
      {
        m_Aligner.Prepare( molecule);
        storage::Pair< chemistry::FragmentComplete, double> alignment
        (
          m_Aligner.FieldOptimizeOrientation
          (
            molecule,
            m_Scaffold,
            100,
            100,
            false,
            m_MaxAtomDistance
          )
        );

        // create voxel grid of the newly aligned input molecule
        chemistry::VoxelGridAtom voxel_grid_molecule( m_MaxAtomDistance);
        voxel_grid_molecule.SetObjects( util::SiPtrVector< const chemistry::AtomConformationalInterface>( alignment.First().GetAtomsIterator().Begin(), alignment.First().GetAtomsIterator().End()));

        // Get neighbors of our input molecule to our scaffold molecule
        auto neighbors( m_VoxelGridScaffold.GetNeighborsIn( voxel_grid_molecule, m_MaxAtomDistance));

        // Get the number of unique atoms in the scaffold that have at least one neighbor in the input molecule
        storage::Set< size_t> scaffold_atoms_with_nbrs;
        for
        (
            auto neighbor_itr( neighbors.Begin()), neighbor_itr_end( neighbors.End());
            neighbor_itr != neighbor_itr_end;
            ++neighbor_itr
        )
        {
          scaffold_atoms_with_nbrs.Insert( m_Scaffold.GetAtomVector().GetAtomIndex( *( neighbor_itr->First())));
        }

        // Get neighbors of our input molecule to our scaffold molecule
        neighbors = voxel_grid_molecule.GetNeighborsIn( m_VoxelGridScaffold, m_MaxAtomDistance);
        BCL_Debug( neighbors.GetSize());

        // Get the number of unique atoms in the scaffold that have at least one neighbor in the input molecule
        storage::Set< size_t> molecule_atoms_with_nbrs;
        for
        (
            auto neighbor_itr( neighbors.Begin()), neighbor_itr_end( neighbors.End());
            neighbor_itr != neighbor_itr_end;
            ++neighbor_itr
        )
        {
          molecule_atoms_with_nbrs.Insert( alignment.First().GetAtomVector().GetAtomIndex( *( neighbor_itr->First())));
        }

        // assign result values
        STORAGE(0) = alignment.Second(); // RMSDX
        STORAGE(1) = float(scaffold_atoms_with_nbrs.GetSize()) / float(m_Scaffold.GetSize());
        STORAGE(2) = float(molecule_atoms_with_nbrs.GetSize()) / float(alignment.First().GetSize());
      }
      else if( m_AlignmentType == "flexible")
      {
        m_Aligner.Prepare( molecule);
        // make conformers
        auto molecule_confs(m_SampleConformations( molecule));
        auto scaffold_confs(m_SampleConformations( m_Scaffold));

        storage::Vector< storage::Triplet< chemistry::FragmentComplete, chemistry::FragmentComplete, double> > alignment
        (
          m_Aligner.FieldOptimizeOrientationFlex
          (
            molecule_confs.First(),
            scaffold_confs.First(),
            10,
            100,
            100,
            10,
            10,
            10,
            10,
            0.75,
            0.50,
            false
          )
        );

        // create voxel grid of the newly aligned input molecule
        chemistry::VoxelGridAtom voxel_grid_molecule( m_MaxAtomDistance);
        voxel_grid_molecule.SetObjects( util::SiPtrVector< const chemistry::AtomConformationalInterface>( alignment( 0).First().GetAtomsIterator().Begin(), alignment( 0).First().GetAtomsIterator().End()));

        chemistry::VoxelGridAtom voxel_grid_scaffold( m_MaxAtomDistance);
        voxel_grid_scaffold.SetObjects( util::SiPtrVector< const chemistry::AtomConformationalInterface>( alignment( 0).Second().GetAtomsIterator().Begin(), alignment( 0).Second().GetAtomsIterator().End()));

        // Get neighbors of our input molecule to our scaffold molecule
        auto neighbors( voxel_grid_scaffold.GetNeighborsIn( voxel_grid_molecule, m_MaxAtomDistance));
        BCL_Debug(voxel_grid_scaffold.GetNumberItems());
        BCL_Debug(voxel_grid_molecule.GetNumberItems());

        // Get the number of unique atoms in the scaffold that have at least one neighbor in the input molecule
        storage::Set< size_t> scaffold_atoms_with_nbrs;
        for
        (
            auto neighbor_itr( neighbors.Begin()), neighbor_itr_end( neighbors.End());
            neighbor_itr != neighbor_itr_end;
            ++neighbor_itr
        )
        {
          scaffold_atoms_with_nbrs.Insert( alignment( 0).Second().GetAtomVector().GetAtomIndex( *( neighbor_itr->First())));
        }

        // Get the number of unique atoms in the scaffold that have at least one neighbor in the input molecule
        storage::Set< size_t> molecule_atoms_with_nbrs;
        for
        (
            auto neighbor_itr( neighbors.Begin()), neighbor_itr_end( neighbors.End());
            neighbor_itr != neighbor_itr_end;
            ++neighbor_itr
        )
        {
          molecule_atoms_with_nbrs.Insert( alignment( 0).First().GetAtomVector().GetAtomIndex( *( neighbor_itr->Second())));
        }

        // assign result values
        BCL_Debug( scaffold_atoms_with_nbrs.GetSize());
        BCL_Debug( molecule_atoms_with_nbrs.GetSize());
        STORAGE(0) = alignment( 0).Third(); // RMSDX
        STORAGE(1) = float(scaffold_atoms_with_nbrs.GetSize()) / float(m_Scaffold.GetSize());
        STORAGE(2) = float(molecule_atoms_with_nbrs.GetSize()) / float(alignment( 0).First().GetSize());
      }
      else
      {
        BCL_MessageStd( "Alignment type specified, '" + util::Format()( m_AlignmentType) + "', is not valid. Please specify 'flexible' or 'rigid'.");
      }

    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MolAlignPharmScore::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Performs quick property-based molecular alignment with BCL::MolAlign and returns"
        "the RMSDX score as well as atomic overlap metrics."
      );
      parameters.AddInitializer
      (
        "scaffold",
        "the molecule to which the inputs are aligned",
        io::Serialization::GetAgentInputFilename( &m_ScaffoldFilename)
      );
      parameters.AddInitializer
      (
        "alignment type",
        "the type of alignment to perform; options are 'flexible' or 'rigid'",
        io::Serialization::GetAgent( &m_AlignmentType),
        "flexible"
      );
      parameters.AddInitializer
      (
        "max atom distance",
        "distance for computing MolAlign mutually matched atoms and VoxelGrid neighbors",
        io::Serialization::GetAgent( &m_MaxAtomDistance),
        "1.0"
      );
      parameters.AddInitializer
      (
        "properties",
        "atom properties to consider, use multiply(Constant(X),property y) for weighting",
        io::Serialization::GetAgent( &m_Properties),
        util::ObjectDataLabel
        (
          "Combine(Atom_Identity)"
        )
      );
      parameters.AddInitializer
      (
        "alignment options",
        "options for the flexible alignment procedure",
        io::Serialization::GetAgent( &m_Aligner),
        "(rigid mol b = true,"
        "number rigid trajectories=3,"
        "number flexible trajectories=3,"
        "fraction filtered initially=0.75,"
        "fraction filtered iteratively=0.50,"
        "iterations=100,"
        "filter iterations=50,"
        "refinement iterations=50,"
        "align to scaffold = false,"
        "conformer pairs = 10,"
        "number outputs = 1,"
        "sample conformers=SampleConformations("
        "conformation_comparer=RMSD,"
        "generate_3D=0,tolerance=0.50,rotamer_library = cod,"
        "max_iterations=100,max_conformations=10,"
        "cluster=true,"
        "clash_tolerance=0.4,"
        "change_chirality=0))"
      );
      parameters.AddInitializer
      (
        "property weights",
        "Weighting to give the properties. Typically these are the inverse standard deviations of the properties over a representative set of molecules, which can be computed w/ "
        "molecule:Properties -statistics",
        io::Serialization::GetAgent( &m_PropertyWeights),
        util::ObjectDataLabel( "(1.0)")
      );
      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool MolAlignPharmScore::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      // Check for scaffold small molecule
      if( m_ScaffoldFilename.size() > 0)
      {
        BCL_MessageVrb( "Scaffold filename: " + util::Format()( m_ScaffoldFilename));

        // Read in scaffold small molecule, remove hydrogen atoms for construction of voxel grid
        io::IFStream input_sdf;
        io::File::MustOpenIFStream( input_sdf, m_ScaffoldFilename);
        chemistry::FragmentEnsemble ensemble( input_sdf, sdf::e_Saturate);
        io::File::CloseClearFStream( input_sdf);

        // If they read in a file then it must have molecules in it
        BCL_Assert( ensemble.GetSize(), "Scaffold molecule file contains no molecules! Exiting...");

        // If they have more than one molecule in the ensemble, just use the first
        if( ensemble.GetSize() > 1)
        {
          BCL_MessageVrb( "Scaffold SDF contains multiple molecules. Just using the first molecule")
        }

        // Get scaffold
        m_Scaffold = *ensemble.Begin();
      }

      // Create voxel grid of our scaffold molecule
      m_VoxelGridScaffold = chemistry::VoxelGridAtom( m_MaxAtomDistance);
      m_VoxelGridScaffold.SetObjects( util::SiPtrVector< const chemistry::AtomConformationalInterface>( m_Scaffold.GetAtomsIterator(), m_Scaffold.GetAtomsIterator().End()));

      // Set up sample conformations object
      m_SampleConformations = chemistry::SampleConformations
          (
            chemistry::RotamerLibraryFile(),
            util::Implementation< chemistry::ConformationComparisonInterface>().GetAlias(),
            0.25, // tolerance
            10,     // number of conformations
            100, // number of iterations
            false,  // change chirality
            0.0,   // random dihedral change weight
            false,  // generate 3d
            0.10   // clash tolerance, auto-adjusted if necessary due to clashes
          );

      // Set properties for alignment
      m_Aligner.SetProperties( m_Properties);
      m_Aligner.SetPropertyWeights( m_PropertyWeights);
      m_Aligner.Prepare( m_Scaffold);

      return true;
    }

  } // namespace descriptor
} // namespace bcl
