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
#include "descriptor/bcl_descriptor_molecule_3da_pair_real_space_asymmetry.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_side_chain_factory.h"
#include "chemistry/bcl_chemistry_aa_fragment_complete.h"
#include "chemistry/bcl_chemistry_voxel_grid_atom.h"
#include "command/bcl_command_command_state.h"
#include "descriptor/bcl_descriptor_molecule_3da_smooth.h"
#include "descriptor/bcl_descriptor_molecule_3da_smooth_sign_occlusion_code.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_operations.h"
#include "pdb/bcl_pdb_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> Molecule3DAPairRealSpaceAsymmetry::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new Molecule3DAPairRealSpaceAsymmetry()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Molecule3DAPairRealSpaceAsymmetry::Molecule3DAPairRealSpaceAsymmetry() :
          m_NumberSteps( 14),
          m_StepSize( 0.50),
          m_Temperature( 5.0),
          m_Smooth( true),
          m_Interpolate( true),
          m_InternalDescriptorSize( 1)
    {
      if( !command::CommandState::IsInStaticInitialization())
      {

        BCL_Assert
        (
          ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
          "Failed to create " + GetClassIdentifier()
        );
      }
    }

    //! @brief constructor from number of steps, and mapped atom property
    Molecule3DAPairRealSpaceAsymmetry::Molecule3DAPairRealSpaceAsymmetry
    (
      const CheminfoProperty &ATOM_PROPERTY_A,
      const CheminfoProperty &ATOM_PROPERTY_B,
      const size_t NUMBER_STEPS,
      const float STEP_SIZE,
      const float TEMPERATURE,
      const bool SMOOTH,
      const bool INTERPOLATE
    ) :
      m_AtomPropertyA( ATOM_PROPERTY_A),
      m_AtomPropertyB( ATOM_PROPERTY_B),
      m_NumberSteps( NUMBER_STEPS),
      m_StepSize( STEP_SIZE),
      m_Temperature( TEMPERATURE),
      m_Smooth( SMOOTH),
      m_Interpolate( INTERPOLATE),
      m_InternalDescriptorSize( 1)
    {
      BCL_Assert
      (
        ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
        "Failed to create " + GetClassIdentifier()
      );
    }

    //! @brief virtual copy constructor
    //! @return pointer to new Molecule3DAPairRealSpaceAsymmetry
    Molecule3DAPairRealSpaceAsymmetry *Molecule3DAPairRealSpaceAsymmetry::Clone() const
    {
      return new Molecule3DAPairRealSpaceAsymmetry( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &Molecule3DAPairRealSpaceAsymmetry::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &Molecule3DAPairRealSpaceAsymmetry::GetAlias() const
    {
      static const std::string s_name( "3DAPairRealSpaceAsymmetry");
      return s_name;
    }

    //! @brief get number of steps of code
    //! @return number of steps in 2da code
    size_t Molecule3DAPairRealSpaceAsymmetry::GetNumberSteps() const
    {
      return m_NumberSteps;
    }

    //! @brief get step size of code
    //! @return step size of 3DA code
    float Molecule3DAPairRealSpaceAsymmetry::GetStepSize() const
    {
      return m_StepSize;
    }

    //! @brief get temperature of code
    //! @return const float  temperature of 3DA code
    const float &Molecule3DAPairRealSpaceAsymmetry::GetTemperature() const
    {
      return m_Temperature;
    }

    //! @brief get atom property A of code
    //! @return atom property A mapped in 3da code
    const CheminfoProperty &Molecule3DAPairRealSpaceAsymmetry::GetAtomPropertyA() const
    {
      return m_AtomPropertyA;
    }
    //! @brief get atom property B of code
    //! @return atom property B mapped in 3da code
    const CheminfoProperty &Molecule3DAPairRealSpaceAsymmetry::GetAtomPropertyB() const
    {
      return m_AtomPropertyB;
    }

    // intiialize static member data
    storage::Map< std::string, storage::Pair< chemistry::FragmentComplete, storage::Map< util::ObjectDataLabel, linal::Vector< float> > > >
        Molecule3DAPairRealSpaceAsymmetry::s_Pockets =
            storage::Map< std::string, storage::Pair< chemistry::FragmentComplete, storage::Map< util::ObjectDataLabel, linal::Vector< float> > > >();
    sched::Mutex Molecule3DAPairRealSpaceAsymmetry::s_Mutex = sched::Mutex();

  ////////////////
  // operations //
  ////////////////

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > Molecule3DAPairRealSpaceAsymmetry::GetInternalDescriptors()
    {
      return
          iterate::Generic< Base< chemistry::AtomConformationalInterface, float> >
      (
        &m_AtomPropertyA,
        &m_AtomPropertyA + 1
      );
    }

    //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
    //! @return the cache preference, assuming this feature has its normal dimension setting
    CachePreference Molecule3DAPairRealSpaceAsymmetry::GetNormalCachePreference() const
    {
      return e_IgnoreCache;
    }

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void Molecule3DAPairRealSpaceAsymmetry::Calculate( linal::VectorReference< float> &STORAGE)
    {
      m_DiscreteCode = 0.0;

      // Make molecule
      util::SiPtr< const chemistry::ConformationInterface> current_mol( this->GetCurrentObject());
      const chemistry::ConformationInterface &molecule( *current_mol);

      // Get ready to collect indices of atoms contributing to interactions
      storage::Vector< size_t> mol_a_indices( molecule.GetSize(), size_t( 0));
      storage::Vector< size_t> mol_b_indices( m_MolB.GetSize(), size_t( 0));

      // Compute descriptor on molecule
      linal::Vector< float> molecule_property( m_AtomPropertyA->CollectValuesOnEachElementOfObject( molecule));

      //Setup voxel grid for molecule
      double neighbor_distance( m_NumberSteps * m_StepSize);
      chemistry::VoxelGridAtom vg_p( neighbor_distance);
      vg_p.SetObjects( util::SiPtrVector< const chemistry::AtomConformationalInterface>( m_MolB.GetAtomsIterator(), m_MolB.GetAtomsIterator().End()));

      // compute products
      size_t mol_a_atom_index( 0);
      for
      (
          Iterator< chemistry::AtomConformationalInterface> itr_m( GetCurrentObject()->GetIterator());
          itr_m.NotAtEnd();
          ++itr_m, ++mol_a_atom_index
      )
      {
        const chemistry::AtomConformationalInterface &atom_m( *itr_m( 0));
        float prop_m( m_AtomPropertyA->operator ()( itr_m)( 0));

        // Find neighbors of atom_m in pocket
        auto neighbors_p( vg_p.GetNeighbors( atom_m.GetPosition(), neighbor_distance));
        for
        (
            auto itr_p( neighbors_p.Begin()), itr_p_end( neighbors_p.End());
            itr_p != itr_p_end;
            ++itr_p
        )
        {
          // Get pocket atom index
          size_t atom_index_vg( m_MolB.GetAtomIndex( *itr_p->First()));
          float prop_p( m_PocketProperties( atom_index_vg));

          // compute bin distance between atom pair
          size_t bin_distance( 4 * size_t( itr_p->Second() / m_StepSize));
          if( bin_distance / 4 < m_NumberSteps)
          {
            // store molecule a atom index for output
            if( m_GetMolAAtomIndices.size() > 0)
            {
              mol_a_indices( mol_a_atom_index) = size_t( 1);
            }

            // store molecule b atom index for output
            if( m_GetMolBAtomIndices.size() > 0)
            {
              mol_b_indices( m_MolB.GetAtomVector().GetAtomIndex( *itr_p->First())) = size_t( 1);
            }
          }
          // Assign product to bin
          Accumulate( itr_p->Second(), prop_m, prop_p);
        }
      }
      // apply smoothing to the discrete bin sums
      Smooth( STORAGE);
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief create a gaussian-smoothed signal from m_DiscreteCode and store it in STORAGE
    //! @param STORAGE storage for the gaussian-smoothed signal
    void Molecule3DAPairRealSpaceAsymmetry::Smooth( linal::VectorReference< float> &STORAGE) const
    {
      if( !m_Smooth)
      {
        STORAGE.CopyValues( m_DiscreteCode);
        return;
      }

      // generate the actual 3DA smooth code by applying the smoothing kernel over the 3DA
      int smoothing_vector_base_position( 0);
      for
      (
        linal::Vector< float>::iterator itr_3DA_smooth( STORAGE.Begin()),
        itr_3DA_smooth_end( STORAGE.End());
        itr_3DA_smooth != itr_3DA_smooth_end;
        ++itr_3DA_smooth, --smoothing_vector_base_position
      )
      {
        float neg_neg( 0), pos_pos( 0), neg_pos( 0), pos_neg( 0);
        int smoothing_vector_position( smoothing_vector_base_position);
        for
        (
          linal::Vector< float>::const_iterator itr_3DA_rough( m_DiscreteCode.Begin()),
          itr_3DA_rough_end( m_DiscreteCode.End());
          itr_3DA_rough != itr_3DA_rough_end;
          ++itr_3DA_rough, ++smoothing_vector_position
        )
        {
          const float exponential_factor( m_SmoothingCoefficients( std::abs( smoothing_vector_position)));
          neg_neg += *itr_3DA_rough * exponential_factor;
          pos_pos += *++itr_3DA_rough * exponential_factor;
          neg_pos += *++itr_3DA_rough * exponential_factor;
          pos_neg += *++itr_3DA_rough * exponential_factor;
        }
        *itr_3DA_smooth = neg_neg;
        *++itr_3DA_smooth = pos_pos;
        *++itr_3DA_smooth = neg_pos;
        *++itr_3DA_smooth = pos_neg;
      }
    }

    //! @brief add an observed distance/property value to m_DiscreteCode
    //! @param DISTANCE actual distance of the two atoms
    //! @param PROP_A property from atom A
    //! @param PROP_B property from atom B
    void Molecule3DAPairRealSpaceAsymmetry::Accumulate
    (
      const float &DISTANCE,
      const float &PROP_A,
      const float &PROP_B
    )
    {
      // compute product between atom properties
      float product = PROP_A * PROP_B;
      size_t sign_id( util::GetUndefinedSize_t());

      // check negative product cases (-/+, +/-)
      if( product < 0)
      {
        // invert sign on product since all bins will be positive
        product = -product;

        // case where the mol a is negative and mol b is positive (-/+)
        if( PROP_A < 0)
        {
          sign_id = 2;
        }
        // case where the mol a is positive and mol b is negative (+/-)
        else
        {
          sign_id = 3;
        }
      }
      // check positive product cases (-/-, +/+)
      else
      {
        // case where mol a and mol b are negative (-/-)
        if( PROP_A < 0)
        {
          sign_id = 0;
        }
        // case where mol a and mol b are positive (+/+)
        else
        {
          sign_id = 1;
        }
      }

      // no need to do anything if product is 0
      if( product == 0.0)
      {
        return;
      }

      // make sure we can place everything and I didn't fuck up
      BCL_Assert( sign_id != util::GetUndefinedSize_t(), "Cannot find sign bin for descriptor value");

//      // this was necessary for the 3DAPairConvol, but it has not
//      // been evaluated in the context of a simple pairwise 3DA like this
//      if( m_Sqrt)
//      {
//        product = math::Sqrt( product);
//      }

      if( !m_Interpolate && !m_Smooth)
      {
        size_t bin_distance( 4 * size_t( DISTANCE / m_StepSize));
        if( bin_distance / 4 < m_NumberSteps)
        {
          m_DiscreteCode( bin_distance + sign_id) += product;
        }
      }
      // calculate exponential function for the RDF equation
      else if( !m_Interpolate)
      {
        const size_t closest_step( float( DISTANCE + 0.5 * m_StepSize) / m_StepSize);
        if( closest_step < m_NumberSteps)
        {
          m_DiscreteCode( 4 * closest_step + sign_id) += product;
        }
      }
      else
      {
        // assign atom distance to the closest 3DA bin
        size_t closest_step_low( std::min( size_t( m_NumberSteps - 1), size_t( DISTANCE / m_StepSize)));
        size_t closest_step_high( closest_step_low + 1);

        const float dist_lower_step( m_StepSize * closest_step_low - DISTANCE);
        const float dist_higher_step( dist_lower_step + m_StepSize);

        float low_exp_factor
        (
          m_Smooth
          ? exp( -m_Temperature * math::Sqr( m_StepSize * closest_step_low - DISTANCE))
          : std::max( dist_higher_step, float( 0.0)) / m_StepSize
        );
        float high_exp_factor
        (
          m_Smooth
          ? exp( -m_Temperature * math::Sqr( m_StepSize * closest_step_high - DISTANCE))
          : ( dist_higher_step >= float( 0.0) ? -dist_lower_step / m_StepSize : 0.0)
        );
        if( closest_step_high == m_NumberSteps)
        {
          // generate rough 3DA for boundry distances
          m_DiscreteCode( 4 * closest_step_low + sign_id) += product * low_exp_factor;
        }
        else
        {
          // normalize rough 3DA code
          product /= ( low_exp_factor + high_exp_factor);
          // generate rough 3DA code for all other distances
          m_DiscreteCode( 4 * closest_step_low + sign_id) += product * low_exp_factor;
          m_DiscreteCode( 4 * closest_step_high + sign_id) += product * high_exp_factor;
        }
      }
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer Molecule3DAPairRealSpaceAsymmetry::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Relates the autocorrelations of two independent molecules"
      );

      parameters.AddInitializer
      (
        "property_a",
        "small molecule property",
        io::Serialization::GetAgent( &m_AtomPropertyA)
      );
      parameters.AddInitializer
      (
        "property_b",
        "protein pocket property",
        io::Serialization::GetAgent( &m_AtomPropertyB)
      );
      parameters.AddInitializer
      (
        "step size",
        "size of each step in angstroms",
        io::Serialization::GetAgentWithRange( &m_StepSize, 0.01, 100.0),
        "0.50"
      );
      parameters.AddInitializer
      (
        "steps",
        "# of steps/bins (each of size = step size) used in the function",
        io::Serialization::GetAgentWithRange( &m_NumberSteps, 1, 1000000),
        "14"
      );
      parameters.AddInitializer
      (
        "temperature",
        "increasing temperature spreads autocorrelation across more distant bins",
        io::Serialization::GetAgentWithRange( &m_Temperature, 0.0, 1000.0),
        "5.0"
      );
      parameters.AddInitializer
      (
        "gaussian",
        "whether to apply gaussian smoothing to the final curve. "
        "If set to false, temperature is ignored, interpolation is linear, and no gaussian smoothing is performed",
        io::Serialization::GetAgent( &m_Smooth),
        "True"
      );
      parameters.AddInitializer
      (
        "interpolate",
        "whether to interpolate values to the two nearest points; if false, all weight will be applied to the nearest bin",
        io::Serialization::GetAgent( &m_Interpolate),
        "True"
      );
      parameters.AddInitializer
      (
        "misc property id",
        "misc property name of corresponding protein binding pocket for current molecule",
        io::Serialization::GetAgent( &m_MiscPropertyString),
        ""
      );
      parameters.AddInitializer
      (
        "reference_conformer",
        "the molecule (e.g. protein binding pocket) whose autocorrelation is being compared against the input molecule descriptor set",
        io::Serialization::GetAgent( &m_MolBFilename),
        ""
      );
      parameters.AddInitializer
      (
        "mol_a_atom_indices_filename",
        "get the atom indices from molecule A contributing to the 3da interaction and output to this file; "
        "recommended to only use this for discrete properties, e.g. "
        "Partial(3DAPairRSAsym050(Multiply(Atom_HydrophobicTernary,Atom_Polarizability),Atom_IsInAromaticRingTernary),indices(23))",
        io::Serialization::GetAgent( &m_GetMolAAtomIndices),
        ""
      );
      parameters.AddInitializer
      (
        "mol_b_atom_indices_filename",
        "get the atom indices from molecule B contributing to the 3da interaction and output to this file; "
        "recommended to only use this for discrete properties, e.g. "
        "Partial(3DAPairRSAsym050(Multiply(Atom_HydrophobicTernary,Atom_Polarizability),Atom_IsInAromaticRingTernary),indices(23))",
        io::Serialization::GetAgent( &m_GetMolBAtomIndices),
        ""
      );
      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool Molecule3DAPairRealSpaceAsymmetry::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {

      // Ensure bin size/resolution specified
      if( m_NumberSteps == 0 || m_StepSize == 0.0)
      {
        ERR_STREAM << "m_NumberSteps equals zero - 3DA code will be empty!";
        return false;
      }

      // Check that a property is provided for the small molecule
        if( m_AtomPropertyA.IsDefined() && m_AtomPropertyA->GetNormalSizeOfFeatures() != 1)
        {
          ERR_STREAM
             << "Expected a property that returned 1 properties per atom, but property returns "
             << m_AtomPropertyA->GetNormalSizeOfFeatures()
             << " values per atom ( property was "
             << m_AtomPropertyA->GetString()
             << ")";
          return false;
        }

        // Check that a property is provided for the protein pocket
        if( m_AtomPropertyB.IsDefined() && m_AtomPropertyB->GetNormalSizeOfFeatures() != 1)
        {
          ERR_STREAM
             << "Expected a property that returned 1 properties per atom, but property returns "
             << m_AtomPropertyB->GetNormalSizeOfFeatures()
             << " values per atom ( property was "
             << m_AtomPropertyB->GetString()
             << ")";
          return false;
        }

      // Initialize our raw value bins and get smoothing coefficients
      m_DiscreteCode = linal::Vector< float>( GetNormalSizeOfFeatures(), 0.0);
      m_SmoothingCoefficients = Molecule3DASmooth::GetSmoothingCoefficientVector( m_NumberSteps, m_Temperature, m_StepSize);

      // Try reading in protein binding pocket from a filename specified in code object file
      if( !m_MolBFilename.empty())
      {
        // map filename to pocket file
        BCL_MessageVrb( "Using reference conformer in CodeObjectFile to generate pocket 3DA");

        // create a temporary inner map
        storage::Map< util::ObjectDataLabel, linal::Vector< float> > temp_map;

        // If we haven't seen this pocket at all yet, we need to make the pocket and compute some properties for it
        if( !s_Pockets.Has( m_MolBFilename))
        {
          // Read in PDB
          BCL_MessageVrb( "Caching pocket: " + m_MolBFilename);
          const pdb::Factory factory( biol::GetAAClasses().e_AAComplete);
          assemble::ProteinModel protein_model( factory.ProteinModelFromPDBFilename( m_MolBFilename, pdb::Factory::GetSSETypeMinSizes( 0, 0, 0)));

          //instantiate AASideChainFactory no hydrogens include backbone atoms for sidechain placement superimposition
          biol::AASideChainFactory side_chain_factory( false, true);

          // add side chains to model and set it to model
          protein_model = *side_chain_factory.ProteinModelWithSideChains( protein_model);
          chemistry::AAFragmentComplete aa_fragment( protein_model.GetAminoAcids(), true);
          aa_fragment.SaturateWithH();
          aa_fragment.RemoveAtomsUndefinedPos();

          // get final model
          m_MolB = aa_fragment;

          // compute properties
          BCL_MessageVrb( "Saving properties of base type: " + util::Format()( m_AtomPropertyB.GetLabel()));
          m_PocketProperties = m_AtomPropertyB->CollectValuesOnEachElementOfObject( m_MolB);

          // save everything to our map of maps
          temp_map[ m_AtomPropertyB.GetLabel()] = m_PocketProperties;
          s_Pockets[ m_MolBFilename] = std::make_pair( m_MolB, temp_map);
        }
        // if we have seen this pocket but not this property, then associate property with pocket
        else if
        (
            s_Pockets.Has( m_MolBFilename) &&
            !s_Pockets.Find( m_MolBFilename)->second.Second().Has( m_AtomPropertyB.GetLabel())
        )
        {
          // get pocket from map
          BCL_MessageVrb( "Retrieving pocket from cache: " + m_MolBFilename);
          m_MolB = s_Pockets.Find( m_MolBFilename)->second.First();

          // compute new properties to associate with the pocket
          BCL_MessageVrb( "Saving properties of base type: " + util::Format()( m_AtomPropertyB.GetLabel()));
          m_PocketProperties = m_AtomPropertyB->CollectValuesOnEachElementOfObject( m_MolB);

          // save everything to our map of maps
          temp_map[ m_AtomPropertyB.GetLabel()] = m_PocketProperties;
          s_Pockets[ m_MolBFilename] = std::make_pair( m_MolB, temp_map);
        }
        // If we have seen this pocket/property pair before, then get everything from the map
        else
        {
          BCL_MessageVrb( "Retrieving pocket from cache: " + m_MolBFilename);
          m_MolB = s_Pockets.Find( m_MolBFilename)->second.First();
          BCL_MessageVrb( "Retrieving properties of base type: " + util::Format()( m_AtomPropertyB.GetLabel()));
          m_PocketProperties = s_Pockets.Find( m_MolBFilename)->second.Second()[ m_AtomPropertyB.GetLabel()];
        }

        // Require that the protein has atoms and properties
        BCL_Assert( m_MolB.GetNumberAtoms() > 0, "The protein indicated by the misc property id label contains no atoms");
        BCL_Assert( m_PocketProperties.GetSize() > 0, "The protein property vector is empty");
      } // end trying to read MolB from the object file
      return true;
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void Molecule3DAPairRealSpaceAsymmetry::SetObjectHook()
    {
      // Get filename from MDL property
      util::SiPtr< const chemistry::ConformationInterface> current_mol( this->GetCurrentObject());
      const chemistry::ConformationInterface &molecule( *current_mol);
      std::string pocket_name( util::Strip( molecule.GetMDLProperty( m_MiscPropertyString), " \t\n\r"));

      // If the pocket is not specified in the CodeObjectFile
      // NOR is specified as an MDL property on the input SDF
      // AND we have seen 1 pocket already (with the associated properties in our CodeObjectFile),
      // then assume we are using the same pocket for every calculation, and we can
      // just find the correct property
      if( m_MolBFilename.empty())
      {
        s_Mutex.Lock();
        if( s_Pockets.GetSize() == 1 && pocket_name.empty())
        {
          BCL_MessageVrb( "Using pre-computed coordinates and properties");
          pocket_name = *( s_Pockets.GetKeysAsVector().Begin());
          m_MolB = s_Pockets.Find( pocket_name)->second.First();
          m_PocketProperties = m_AtomPropertyB->CollectValuesOnEachElementOfObject( m_MolB);
        }
        // If the pocket is not specified in the CodeObjectFile
        // AND is specified as an MDL property on the input SDF
        else if( !pocket_name.empty())
        {
          // map filename to pocket file
          BCL_MessageVrb( "Using MiscProperty pocket ID label to generate pocket 3DA");

          // create a temporary inner map
          storage::Map< util::ObjectDataLabel, linal::Vector< float> > temp_map;

          // If we haven't seen this pocket at all yet, we need to make the pocket and compute some properties for it
          if( !s_Pockets.Has( pocket_name))
          {
            // Read in PDB
            BCL_MessageVrb( "Caching pocket: " + pocket_name);
            const pdb::Factory factory( biol::GetAAClasses().e_AAComplete);
            assemble::ProteinModel protein_model( factory.ProteinModelFromPDBFilename( pocket_name, pdb::Factory::GetSSETypeMinSizes( 0, 0, 0)));

            //instantiate AASideChainFactory no hydrogens include backbone atoms for sidechain placement superimposition
            biol::AASideChainFactory side_chain_factory( false, true);

            // add side chains to model and set it to model
            protein_model = *side_chain_factory.ProteinModelWithSideChains( protein_model);
            chemistry::AAFragmentComplete aa_fragment( protein_model.GetAminoAcids(), true);
            aa_fragment.SaturateWithH();
            aa_fragment.RemoveAtomsUndefinedPos();

            // get final model
            m_MolB = aa_fragment;

            // compute properties
            BCL_MessageVrb( "Saving properties of base type: " + util::Format()( m_AtomPropertyB.GetLabel()));
            m_PocketProperties = m_AtomPropertyB->CollectValuesOnEachElementOfObject( m_MolB);

            // save everything to our map of maps
            temp_map[ m_AtomPropertyB.GetLabel()] = m_PocketProperties;
            s_Pockets[ pocket_name] = std::make_pair( m_MolB, temp_map);
          }
          // if we have seen this pocket but not this property, then associate property with pocket
          else if
          (
              s_Pockets.Has( pocket_name) &&
              !s_Pockets.Find( pocket_name)->second.Second().Has( m_AtomPropertyB.GetLabel())
          )
          {
            // get pocket from map
            BCL_MessageVrb( "Retrieving pocket from cache: " + pocket_name);
            m_MolB = s_Pockets.Find( pocket_name)->second.First();

            // compute new properties to associate with the pocket
            BCL_MessageVrb( "Saving properties of base type: " + util::Format()( m_AtomPropertyB.GetLabel()));
            m_PocketProperties = m_AtomPropertyB->CollectValuesOnEachElementOfObject( m_MolB);

            // save everything to our map of maps
            temp_map[ m_AtomPropertyB.GetLabel()] = m_PocketProperties;
            s_Pockets[ pocket_name] = std::make_pair( m_MolB, temp_map);
          }
          // If we have seen this pocket/property pair before, then get everything from the map
          else
          {
            BCL_MessageVrb( "Retrieving pocket from cache: " + pocket_name);
            m_MolB = s_Pockets.Find( pocket_name)->second.First();
            BCL_MessageVrb( "Retrieving properties of base type: " + util::Format()( m_AtomPropertyB.GetLabel()));
            m_PocketProperties = s_Pockets.Find( pocket_name)->second.Second()[ m_AtomPropertyB.GetLabel()];
          }

          // Require that the protein has atoms and properties
          BCL_MessageVrb( "Number of atoms in receptor: " + util::Format()( s_Pockets[ pocket_name].First().GetNumberAtoms()));
          BCL_Assert( s_Pockets[ pocket_name].First().GetNumberAtoms() > size_t( 0), "The protein indicated by the misc property id label contains no atoms");
          BCL_Assert( m_PocketProperties.GetSize() > size_t( 0), "The protein property vector is empty");
        }
        s_Mutex.Unlock();
      }
//        // Require that the protein has atoms and properties
//      BCL_MessageVrb("Final number of atoms in receptor: " + util::Format()( s_Pockets[ pocket_name].First().GetNumberAtoms()));
//      BCL_Assert( s_Pockets[pocket_name].First().GetNumberAtoms() > size_t( 0), "The protein indicated by the misc property id label contains no atoms");
//      BCL_Assert( m_PocketProperties.GetSize() > size_t( 0), "The protein property vector is empty");
    }

  } // namespace descriptor
} // namespace bcl
