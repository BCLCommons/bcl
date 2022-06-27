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
#include "math/bcl_math_angle.h"
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "descriptor/bcl_descriptor_molecule_inter_hbond_code.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_side_chain_factory.h"
#include "chemistry/bcl_chemistry_aa_fragment_complete.h"
#include "chemistry/bcl_chemistry_voxel_grid_atom.h"
#include "command/bcl_command_command_state.h"
#include "descriptor/bcl_descriptor_molecule_3da_smooth.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ifstream.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_operations.h"
#include "pdb/bcl_pdb_factory.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> MoleculeInterHBondCode::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new MoleculeInterHBondCode()
      )
    );

    // intiialize static member data
    storage::Map< std::string, storage::Pair< chemistry::FragmentComplete, storage::Map< util::ObjectDataLabel, linal::Vector< float> > > >
        MoleculeInterHBondCode::s_Pockets =
            storage::Map< std::string, storage::Pair< chemistry::FragmentComplete, storage::Map< util::ObjectDataLabel, linal::Vector< float> > > >();
    sched::Mutex MoleculeInterHBondCode::s_Mutex = sched::Mutex();

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoleculeInterHBondCode::MoleculeInterHBondCode() :
      m_AtomHBD( CheminfoProperty()),
      m_AtomHBA( CheminfoProperty()),
      m_NumberSteps( 14),
      m_StepSize( 0.50),
      m_InternalDescriptorSize( 1)
    {
      if( !command::CommandState::IsInStaticInitialization())
      {
        this->TryRead( util::ObjectDataLabel(), util::GetLogger());
        BCL_Assert
        (
          ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
          "Failed to create " + GetClassIdentifier()
        );
      }
    }

    //! @brief constructor from number of steps, and mapped atom property
    MoleculeInterHBondCode::MoleculeInterHBondCode
    (
      const CheminfoProperty &ATOM_HBD,
      const CheminfoProperty &ATOM_HBA,
      const size_t NUMBER_STEPS,
      const float STEP_SIZE
    ) :
      m_AtomHBD( ATOM_HBD),
      m_AtomHBA( ATOM_HBA),
      m_NumberSteps( NUMBER_STEPS),
      m_StepSize( STEP_SIZE),
      m_InternalDescriptorSize( 1)
    {
      BCL_Assert
      (
        ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
        "Failed to create " + GetClassIdentifier()
      );
    }

    //! @brief virtual copy constructor
    //! @return pointer to new MoleculeInterHBondCode
    MoleculeInterHBondCode *MoleculeInterHBondCode::Clone() const
    {
      return new MoleculeInterHBondCode( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &MoleculeInterHBondCode::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &MoleculeInterHBondCode::GetAlias() const
    {
      static const std::string s_name( "3DInterHBondCode");
      return s_name;
    }

    //! @brief get number of steps of code
    //! @return number of steps in 2da code
    size_t MoleculeInterHBondCode::GetNumberSteps() const
    {
      return m_NumberSteps;
    }

    //! @brief get step size of code
    //! @return step size of 3DA code
    float MoleculeInterHBondCode::GetStepSize() const
    {
      return m_StepSize;
    }

    //! @brief get angle size of code
    //! @return angle size of 3DA code
    size_t MoleculeInterHBondCode::GetAngleSize() const
    {
      return m_AngleSize;
    }

    //! @brief get atom property of code
    //! @return atom property mapped in 3da code
    const CheminfoProperty &MoleculeInterHBondCode::GetAtomProperty() const
    {
      return m_AtomHBD;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > MoleculeInterHBondCode::GetInternalDescriptors()
    {
      return
        iterate::Generic< Base< chemistry::AtomConformationalInterface, float> >
        (
          &m_AtomHBD,
          &m_AtomHBD + 1
        );
    }

    //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
    //! @return the cache preference, assuming this feature has its normal dimension setting
    CachePreference MoleculeInterHBondCode::GetNormalCachePreference() const
    {
      return e_IgnoreCache;
    }

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void MoleculeInterHBondCode::Calculate( linal::VectorReference< float> &STORAGE)
    {
      m_DiscreteCode = 0.0;

      // Make molecule
      util::SiPtr< const chemistry::ConformationInterface> current_mol( this->GetCurrentObject());
      const chemistry::ConformationInterface &molecule( *current_mol);
      chemistry::FragmentComplete mol( molecule);

      linal::Vector< float> molecule_hbd( m_AtomHBD->CollectValuesOnEachElementOfObject( molecule));
      linal::Vector< float> molecule_hba( m_AtomHBA->CollectValuesOnEachElementOfObject( molecule));

      // Get ready to collect indices of atoms contributing to interactions
      storage::Vector< size_t> mol_a_indices( molecule.GetSize(), size_t( 0));
      storage::Vector< size_t> mol_b_indices( m_MolB.GetSize(), size_t( 0));

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
        // compute property at molecule atom
        const chemistry::AtomConformationalInterface &atom_m( *itr_m( 0));

        // compute HBD property
        bool mol_hbd( true);
        float prop_m( m_AtomHBD->operator ()( itr_m)( 0));

        // if not a HBD atom, then compute HBA property
        if( !prop_m)
        {
          mol_hbd = false;
          prop_m = m_AtomHBA->operator ()( itr_m)( 0);
        }

        // if not a HBA atom either, then skip it and keep going
        if( !prop_m)
        {
          continue;
        }

        // Find neighbors of atom_m in pocket
        auto neighbors_p( vg_p.GetNeighbors( atom_m.GetPosition(), neighbor_distance));
        for
        (
            auto itr_p( neighbors_p.Begin()), itr_p_end( neighbors_p.End());
            itr_p != itr_p_end;
            ++itr_p
        )
        {
          // Get pocket atom index and property at that atom
          size_t atom_index_vg( m_MolB.GetAtomIndex( *itr_p->First()));

          // If atom_m is a HBD, then compute HBA for protein
          float prop_p;
          if( mol_hbd)
          {
            prop_p = m_PocketHBA( atom_index_vg);
          }
          // otherwise if we made it this far, compute the HBD property
          else
          {
            prop_p = m_PocketHBD( atom_index_vg);
          }

          // if the pocket atom is not a HBD or HBA, then continue
          if( !prop_p)
          {
            continue;
          }

          // make sure the bin distance does not exceed the number of bins we have allotted
          size_t bin_distance( size_t( itr_p->Second() / m_StepSize));
          if( bin_distance < m_NumberSteps)
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

            // Get angle formed by segments HBD-H and HBA-H
            linal::Vector3D mol_pos( mol.GetAtomVector()( mol_a_atom_index).GetPosition());
            linal::Vector3D pocket_pos( itr_p->First()->GetPosition());
            linal::Vector3D proton_pos;

            // find the proton on either the molecule or pocket atom
            if( mol_hbd)
            {
              for
              (
                  auto bond_itr( mol.GetAtomVector()( mol_a_atom_index).GetBonds().Begin()),
                  bond_itr_end( mol.GetAtomVector()( mol_a_atom_index).GetBonds().End());
                  bond_itr != bond_itr_end;
                  ++bond_itr
              )
              {
                if( bond_itr->GetTargetAtom().GetElementType() == chemistry::GetElementTypes().e_Hydrogen)
                {
                  proton_pos = bond_itr->GetTargetAtom().GetPosition();
                }
              }
            }
            else
            {
              for
              (
                  auto bond_itr( itr_p->First()->GetBonds().Begin()),
                  bond_itr_end( itr_p->First()->GetBonds().End());
                  bond_itr != bond_itr_end;
                  ++bond_itr
              )
              {
                if( bond_itr->GetTargetAtom().GetElementType() == chemistry::GetElementTypes().e_Hydrogen)
                {
                  proton_pos = bond_itr->GetTargetAtom().GetPosition();
                }
              }
            }
            // get the angle
            float hbond_angle( math::Angle::Degree( linal::ProjAngle( proton_pos, mol_pos, pocket_pos)));
            size_t angle_bin( hbond_angle / m_AngleSize);
            size_t n_angle_bins( 360 / m_AngleSize);
            STORAGE( ( bin_distance * n_angle_bins) + angle_bin) += 1;
          }
        }
      } // end compute products

      // output atom indices
      io::OFStream output;
      if( m_GetMolAAtomIndices.size() > 0)
      {
        io::File::MustOpenOFStream( output, m_GetMolAAtomIndices);
        output << mol_a_indices;
        io::File::CloseClearFStream( output);
      }
      if( m_GetMolBAtomIndices.size() > 0)
      {
        io::File::MustOpenOFStream( output, m_GetMolBAtomIndices);
        output << mol_b_indices;
        io::File::CloseClearFStream( output);
      }
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeInterHBondCode::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Represents hydrogen bond interactions in relative distance and angle correlation bins"
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
        "angle size",
        "size of each angle increment in degrees",
        io::Serialization::GetAgentWithRange( &m_AngleSize, 1, 180),
        "8"
      );
      parameters.AddInitializer
      (
        "steps",
        "# of steps/bins (each of size = step size) used in the function",
        io::Serialization::GetAgentWithRange( &m_NumberSteps, 1, 1000000),
        "8"
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
    bool MoleculeInterHBondCode::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      // check static initialization
      if( command::CommandState::IsInStaticInitialization())
      {
        return true;
      }

      // Set properties
      m_AtomHBD = CheminfoProperty("Atom_HbondDonors");
      m_AtomHBA = CheminfoProperty("Multiply(Atom_HbondAcceptors,Not(Atom_HbondDonors),LessEqual(lhs=BondTypeCount,rhs=Constant(2)))");

      // Ensure bin size/resolution specified
      if( m_NumberSteps == 0 || m_StepSize == 0.0)
      {
        ERR_STREAM << "m_NumberSteps equals zero - 3DA code will be empty!";
        return false;
      }

      // Ensure that selected angle discretization is perfect multiple of 360 degrees
      if( 360 % m_AngleSize != 0)
      {
        BCL_MessageStd( "Specified angle size must be an even multiple of 360 degrees");
        return false;
      }

      // Initialize our raw value bins and get smoothing coefficients
      m_DiscreteCode = linal::Vector< float>( GetNormalSizeOfFeatures(), 0.0);

      // Try reading in protein binding pocket from a filename specified in code object file
      if( !m_MolBFilename.empty())
      {
        // map filename to pocket file
        BCL_MessageVrb( "Using reference conformer in CodeObjectFile to generate pocket 3DA");

        // create a temporary inner map
        storage::Map< util::ObjectDataLabel, linal::Vector< float> > hbd_map, hba_map;

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
          m_PocketHBD = m_AtomHBD->CollectValuesOnEachElementOfObject( m_MolB);
          m_PocketHBA = m_AtomHBA->CollectValuesOnEachElementOfObject( m_MolB);

          // save everything to our map of maps
          hbd_map[ m_AtomHBD.GetLabel()] = m_PocketHBD;
          hba_map[ m_AtomHBD.GetLabel()] = m_PocketHBA;
          s_Pockets[ m_MolBFilename] = std::make_pair( m_MolB, hbd_map);
          s_Pockets[ m_MolBFilename] = std::make_pair( m_MolB, hba_map);
        }
        // if we have seen this pocket but not this property, then associate property with pocket
        else if
        (
            s_Pockets.Has( m_MolBFilename) &&
            (
                !s_Pockets.Find( m_MolBFilename)->second.Second().Has( m_AtomHBD.GetLabel()) ||
                !s_Pockets.Find( m_MolBFilename)->second.Second().Has( m_AtomHBA.GetLabel())
            )
        )
        {
          // get pocket from map
          BCL_MessageVrb( "Retrieving pocket from cache: " + m_MolBFilename);
          m_MolB = s_Pockets.Find( m_MolBFilename)->second.First();

          // compute new properties to associate with the pocket
          m_PocketHBD = m_AtomHBD->CollectValuesOnEachElementOfObject( m_MolB);
          m_PocketHBA = m_AtomHBA->CollectValuesOnEachElementOfObject( m_MolB);

          // save everything to our map of maps
          hbd_map[ m_AtomHBD.GetLabel()] = m_PocketHBD;
          hbd_map[ m_AtomHBA.GetLabel()] = m_PocketHBA;
          s_Pockets[ m_MolBFilename] = std::make_pair( m_MolB, hbd_map);
          s_Pockets[ m_MolBFilename] = std::make_pair( m_MolB, hba_map);
        }
        // If we have seen this pocket/property pair before, then get everything from the map
        else
        {
          BCL_MessageVrb( "Retrieving pocket from cache: " + m_MolBFilename);
          m_MolB = s_Pockets.Find( m_MolBFilename)->second.First();
          m_PocketHBD = s_Pockets.Find( m_MolBFilename)->second.Second()[ m_AtomHBD.GetLabel()];
          m_PocketHBA = s_Pockets.Find( m_MolBFilename)->second.Second()[ m_AtomHBA.GetLabel()];
        }

        // Require that the protein has atoms and properties
        BCL_Assert( m_MolB.GetNumberAtoms() > 0, "The protein indicated by the misc property id label contains no atoms");
        BCL_Assert( m_PocketHBD.GetSize() > 0, "The protein HBD property vector is empty");
        BCL_Assert( m_PocketHBA.GetSize() > 0, "The protein HBA property vector is empty");
      } // end trying to read MolB from the object file
      return true;
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    // occurs after the ReadInitializerSuccessHook
    void MoleculeInterHBondCode::SetObjectHook()
    {
      // Set properties
      m_AtomHBD = CheminfoProperty("Atom_HbondDonors");
      m_AtomHBA = CheminfoProperty("Multiply(Atom_HbondAcceptors,Not(Atom_HbondDonors),LessEqual(lhs=BondTypeCount,rhs=Constant(2)))");

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
          m_PocketHBD = m_AtomHBD->CollectValuesOnEachElementOfObject( m_MolB);
          m_PocketHBA = m_AtomHBA->CollectValuesOnEachElementOfObject( m_MolB);
        }
        // If the pocket is not specified in the CodeObjectFile
        // AND is specified as an MDL property on the input SDF
        else if( !pocket_name.empty())
        {
          // map filename to pocket file
          BCL_MessageVrb( "Using MiscProperty pocket ID label to generate pocket 3DA");

          // create a temporary inner map
          storage::Map< util::ObjectDataLabel, linal::Vector< float> > hbd_map, hba_map;

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
            m_PocketHBD = m_AtomHBD->CollectValuesOnEachElementOfObject( m_MolB);
            m_PocketHBA = m_AtomHBA->CollectValuesOnEachElementOfObject( m_MolB);

            // save everything to our map of maps
            hbd_map[ m_AtomHBD.GetLabel()] = m_PocketHBD;
            hbd_map[ m_AtomHBA.GetLabel()] = m_PocketHBA;
            s_Pockets[ pocket_name] = std::make_pair( m_MolB, hbd_map);
            s_Pockets[ pocket_name] = std::make_pair( m_MolB, hba_map);
          }
          // if we have seen this pocket but not this property, then associate property with pocket
          else if
          (
              s_Pockets.Has( pocket_name) &&
              (
                  !s_Pockets.Find( pocket_name)->second.Second().Has( m_AtomHBD.GetLabel()) ||
                  !s_Pockets.Find( pocket_name)->second.Second().Has( m_AtomHBA.GetLabel())
              )
          )
          {
            // get pocket from map
            BCL_MessageVrb( "Retrieving pocket from cache: " + pocket_name);
            m_MolB = s_Pockets.Find( pocket_name)->second.First();

            // compute new properties to associate with the pocket
            m_PocketHBD = m_AtomHBD->CollectValuesOnEachElementOfObject( m_MolB);
            m_PocketHBA = m_AtomHBA->CollectValuesOnEachElementOfObject( m_MolB);

            // save everything to our map of maps
            hbd_map[ m_AtomHBD.GetLabel()] = m_PocketHBD;
            hbd_map[ m_AtomHBA.GetLabel()] = m_PocketHBA;
            s_Pockets[ pocket_name] = std::make_pair( m_MolB, hbd_map);
            s_Pockets[ pocket_name] = std::make_pair( m_MolB, hba_map);
          }
          // If we have seen this pocket/property pair before, then get everything from the map
          else
          {
            BCL_MessageVrb( "Retrieving pocket from cache: " + pocket_name);
            m_MolB = s_Pockets.Find( pocket_name)->second.First();
            m_PocketHBD = s_Pockets.Find( pocket_name)->second.Second()[ m_AtomHBD.GetLabel()];
            m_PocketHBA = s_Pockets.Find( pocket_name)->second.Second()[ m_AtomHBA.GetLabel()];
          }

          // Require that the protein has atoms and properties
          BCL_Assert( m_MolB.GetNumberAtoms() > 0, "The protein indicated by the misc property id label contains no atoms");
          BCL_Assert( m_PocketHBD.GetSize() > 0, "The protein HBD property vector is empty");
          BCL_Assert( m_PocketHBA.GetSize() > 0, "The protein HBA property vector is empty");
        }
        s_Mutex.Unlock();
      }
      // Require that the protein has atoms and properties
      BCL_Assert( m_MolB.GetNumberAtoms() > 0, "The protein indicated by the misc property id label contains no atoms");
      BCL_Assert( m_PocketHBD.GetSize() > 0, "The protein HBD property vector is empty");
      BCL_Assert( m_PocketHBA.GetSize() > 0, "The protein HBA property vector is empty");
    } // end SetObjectHook
  } // namespace descriptor
} // namespace bcl
