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
#include "descriptor/bcl_descriptor_molecule_3da_pair_convolution.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_side_chain_factory.h"
#include "chemistry/bcl_chemistry_aa_fragment_complete.h"
#include "command/bcl_command_command_state.h"
#include "descriptor/bcl_descriptor_molecule_3da_smooth_sign_occlusion_code.h"
#include "io/bcl_io_ifstream.h"
#include "iterate/bcl_iterate_reflecting.h"
#include "math/bcl_math_running_min_max.h"
#include "pdb/bcl_pdb_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> Molecule3DAPairConvolution::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new Molecule3DAPairConvolution()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Molecule3DAPairConvolution::Molecule3DAPairConvolution() :
      m_NumberSteps( 20),
      m_StepSize( 0.50),
      m_WindowSize( 4),
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

    // intiialize static member data
    storage::Map< std::string, storage::Pair< chemistry::FragmentComplete, storage::Map< util::ObjectDataLabel, linal::Vector< float> > > >
        Molecule3DAPairConvolution::s_Pockets =
            storage::Map< std::string, storage::Pair< chemistry::FragmentComplete, storage::Map< util::ObjectDataLabel, linal::Vector< float> > > >();
    sched::Mutex Molecule3DAPairConvolution::s_Mutex = sched::Mutex();

    //! @brief constructor from number of steps, and mapped atom property
    Molecule3DAPairConvolution::Molecule3DAPairConvolution
    (
      const CheminfoProperty &ATOM_PROPERTY,
      const size_t NUMBER_STEPS,
      const float STEP_SIZE,
      const size_t WINDOW_SIZE
    ) :
      m_AtomProperty( ATOM_PROPERTY),
      m_NumberSteps( NUMBER_STEPS),
      m_StepSize( STEP_SIZE),
      m_WindowSize( WINDOW_SIZE),
      m_InternalDescriptorSize( 1)
    {
      BCL_Assert
      (
        ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
        "Failed to create " + GetClassIdentifier()
      );
    }

    //! @brief virtual copy constructor
    //! @return pointer to new Molecule3DAPairConvolution
    Molecule3DAPairConvolution *Molecule3DAPairConvolution::Clone() const
    {
      return new Molecule3DAPairConvolution( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &Molecule3DAPairConvolution::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &Molecule3DAPairConvolution::GetAlias() const
    {
      static const std::string s_name( "3DAPairConvolution");
      return s_name;
    }

    //! @brief get number of steps of code
    //! @return number of steps in 2da code
    size_t Molecule3DAPairConvolution::GetNumberSteps() const
    {
      return m_NumberSteps;
    }

    //! @brief get step size of code
    //! @return step size of 3DA code
    float Molecule3DAPairConvolution::GetStepSize() const
    {
      return m_StepSize;
    }

    //! @brief get atom property of code
    //! @return atom property mapped in 2da code
    const CheminfoProperty &Molecule3DAPairConvolution::GetAtomProperty() const
    {
      return m_AtomProperty;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    iterate::Generic< Base< chemistry::AtomConformationalInterface, float> > Molecule3DAPairConvolution::GetInternalDescriptors()
    {
      return
          iterate::Generic< Base< chemistry::AtomConformationalInterface, float> >
      (
        &m_3DASmoothSign,
        &m_3DASmoothSign + 1
      );
    }

    //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
    //! @return the cache preference, assuming this feature has its normal dimension setting
    CachePreference Molecule3DAPairConvolution::GetNormalCachePreference() const
    {
      return e_PreferCache;
    }

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void Molecule3DAPairConvolution::Calculate( linal::VectorReference< float> &STORAGE)
    {
      //! Initial 3DA for molecule 1 property to be used in convolution
      linal::Vector< float> mol_a_3da;
      mol_a_3da = m_3DASmoothSign.SumOverObject( *this->GetCurrentObject());
      linal::Vector< float> mol_a_3da_windowavg( mol_a_3da.GetSize(), 0.0);
      linal::Vector< float>::iterator itr_mol_a_3da_windowavg( mol_a_3da_windowavg.Begin());

      for( auto mol_a_3da_itr( mol_a_3da.Begin()); mol_a_3da_itr != mol_a_3da.End(); ++++++mol_a_3da_itr, ++itr_mol_a_3da_windowavg)
      {
        math::RunningAverage< float> average_nn, average_pp, average_pn;
        iterate::Reflecting< const float> itr_refl_mol_aa( mol_a_3da.Begin(), mol_a_3da.End(), mol_a_3da_itr);
        for( size_t win_pos( 0); win_pos < m_WindowSize; ++win_pos, ++itr_refl_mol_aa)
        {
          // Account for the signed nature of the 3DAs (--/++/+-) in the window weighting
          average_nn.AddWeightedObservation(*itr_refl_mol_aa,m_WindowWeights(win_pos));
          ++itr_refl_mol_aa;
          average_pp.AddWeightedObservation(*itr_refl_mol_aa,m_WindowWeights(win_pos));
          ++itr_refl_mol_aa;
          average_pn.AddWeightedObservation(*itr_refl_mol_aa,m_WindowWeights(win_pos));
        }

        iterate::Reflecting< const float> itr_refl_mol_ab( mol_a_3da.Begin(), mol_a_3da.End(), mol_a_3da_itr);
        --itr_refl_mol_ab;
        for( size_t win_pos( 1); win_pos < m_WindowSize; ++win_pos, --itr_refl_mol_ab)
        {
          // Account for the signed nature of the 3DAs (--/++/+-) in the window weighting
          average_pn.AddWeightedObservation(*itr_refl_mol_ab,m_WindowWeights(win_pos));
          --itr_refl_mol_ab;
          average_pp.AddWeightedObservation(*itr_refl_mol_ab,m_WindowWeights(win_pos));
          --itr_refl_mol_ab;
          average_nn.AddWeightedObservation(*itr_refl_mol_ab,m_WindowWeights(win_pos));
        }
        *itr_mol_a_3da_windowavg = average_nn.GetAverage();
        ++itr_mol_a_3da_windowavg;
        *itr_mol_a_3da_windowavg = average_pp.GetAverage();
        ++itr_mol_a_3da_windowavg;
        *itr_mol_a_3da_windowavg = average_pn.GetAverage();
      }

      // compare weighted window averages of molecule A 3DAs with each 3DA from molecule B
      // The bins correspond to --/-- , --/++ , --/-+ , ++/-- , ++/++ , ++/-+ , -+/-- , -+/++ , -+/-+
      // Indices correspond to    0       1       2       3       4       5       6       7       8
      for( size_t mol_a_3da_avg_index( 0); mol_a_3da_avg_index < mol_a_3da_windowavg.GetSize(); ++++++mol_a_3da_avg_index)
      {
        for( size_t i( 0); i < 9; ++i)
        {
          STORAGE( 3 * mol_a_3da_avg_index + i) = math::Sqrt( mol_a_3da_windowavg( mol_a_3da_avg_index + ( i / 3)) * m_MolB3DA( mol_a_3da_avg_index + ( i % 3)));
        }
      }
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer Molecule3DAPairConvolution::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Relates the autocorrelations of two independent molecules"
      );

      parameters.AddInitializer
      (
        "property",
        "property over which to calculate the function",
        io::Serialization::GetAgent( &m_AtomProperty)
      );
      parameters.AddInitializer
      (
        "steps",
        "# of steps/bins (each of size = step size) used in the function",
        io::Serialization::GetAgentWithRange( &m_NumberSteps, 1, 1000000),
        "20"
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
        "window size",
        "size of each step in multiples of step size",
        io::Serialization::GetAgent( &m_WindowSize),
        "4.0"
      );
      parameters.AddInitializer
      (
        "weighting",
        "method of generating weights for the windowing function",
        io::Serialization::GetAgent( &m_WindowWeightsCreator),
        "Hann"
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
      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool Molecule3DAPairConvolution::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      // Setup 3DAs
      m_3DASmoothSign = Molecule3DASmoothSignCode( m_AtomProperty, m_NumberSteps, m_StepSize, 100.0, false, true);

      // Ensure bin size/resolution specified
      if( m_NumberSteps == 0 || m_StepSize == 0.0)
      {
        ERR_STREAM << "m_NumberSteps equals zero - 3DA code will be empty!";
        return false;
      }

      // Check that a property is provided
      if( m_AtomProperty.IsDefined() && m_AtomProperty->GetNormalSizeOfFeatures() != 1)
      {
        ERR_STREAM
           << "Expected a property that returned 1 properties per atom, but property returns "
           << m_AtomProperty->GetNormalSizeOfFeatures()
           << " values per atom ( property was "
           << m_AtomProperty->GetString()
           << ")";
        return false;
      }

      // Check that a property is provided
      if( !m_AtomProperty.IsDefined())
      {
        return false;
      }

      // Read in protein binding pocket / molecule B
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
          ReduceProteinPocket( m_MolB);

          // compute properties
          BCL_MessageVrb( "Saving properties of base type: " + util::Format()( m_AtomProperty.GetLabel()));
          m_MolB3DA = ( m_3DASmoothSign.SumOverObject( m_MolB));

          // save everything to our map of maps
          temp_map[ m_AtomProperty.GetLabel()] = m_MolB3DA;
          s_Pockets[ m_MolBFilename] = std::make_pair( m_MolB, temp_map);
        }
        // if we have seen this pocket but not this property, then associate property with pocket
        else if
        (
            s_Pockets.Has( m_MolBFilename) &&
            !s_Pockets.Find( m_MolBFilename)->second.Second().Has( m_AtomProperty.GetLabel())
        )
        {
          // get pocket from map
          BCL_MessageVrb( "Retrieving pocket from cache: " + m_MolBFilename);
          m_MolB = s_Pockets.Find( m_MolBFilename)->second.First();

          // compute new properties to associate with the pocket
          BCL_MessageVrb( "Saving properties of base type: " + util::Format()( m_AtomProperty.GetLabel()));
          m_MolB3DA = ( m_3DASmoothSign.SumOverObject( m_MolB));

          // save everything to our map of maps
          temp_map[ m_AtomProperty.GetLabel()] = m_MolB3DA;
          s_Pockets[ m_MolBFilename] = std::make_pair( m_MolB, temp_map);
        }
        // If we have seen this pocket/property pair before, then get everything from the map
        else
        {
          BCL_MessageVrb( "Retrieving pocket from cache: " + m_MolBFilename);
          m_MolB = s_Pockets.Find( m_MolBFilename)->second.First();
          BCL_MessageVrb( "Retrieving properties of base type: " + util::Format()( m_AtomProperty.GetLabel()));
          m_MolB3DA = ( m_3DASmoothSign.SumOverObject( m_MolB));
        }

        // Require that the protein has atoms and properties
        BCL_Assert( m_MolB.GetNumberAtoms() > 0, "The protein indicated by the misc property id label contains no atoms");
        BCL_Assert( m_MolB3DA.GetSize() > 0, "The protein property vector is empty");

      } // end trying to read MolB from the object file

      // get the coefficients for the half window
      m_WindowWeights = m_WindowWeightsCreator->operator()( m_WindowSize);
      return true;
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void Molecule3DAPairConvolution::SetObjectHook()
    {
      // Get filename from MDL property
      util::SiPtr< const chemistry::ConformationInterface> current_mol( this->GetCurrentObject());
      const chemistry::ConformationInterface &molecule( *current_mol);
      std::string pocket_name( util::Strip( molecule.GetMDLProperty( m_MiscPropertyString), " \t\n\r"));

      // If the pocket is not specified in the CodeObjectFile
      // NOR is specified as an MDL property on the input SDF
      // AND we have seen 1 pocket already,
      // then assume we are using the same pocket for every calculation
      if( m_MolBFilename.empty())
      {
        s_Mutex.Lock();
        if( s_Pockets.GetSize() == 1 && pocket_name.empty())
        {
          BCL_MessageVrb( "Using pre-computed coordinates and properties");
          pocket_name = *( s_Pockets.GetKeysAsVector().Begin());
          m_MolB = s_Pockets.Find( pocket_name)->second.First();
          m_MolB3DA = ( m_3DASmoothSign.SumOverObject( m_MolB));
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
            ReduceProteinPocket( m_MolB);

            // compute properties
            BCL_MessageVrb( "Saving properties of base type: " + util::Format()( m_AtomProperty.GetLabel()));
            m_MolB3DA = ( m_3DASmoothSign.SumOverObject( m_MolB));

            // save everything to our map of maps
            temp_map[ m_AtomProperty.GetLabel()] = m_MolB3DA;
            s_Pockets[ pocket_name] = std::make_pair( m_MolB, temp_map);
          }
          // if we have seen this pocket but not this property, then associate property with pocket
          else if
          (
              s_Pockets.Has( pocket_name) &&
              !s_Pockets.Find( pocket_name)->second.Second().Has( m_AtomProperty.GetLabel())
          )
          {
            // get pocket from map
            BCL_MessageVrb( "Retrieving pocket from cache: " + pocket_name);
            m_MolB = s_Pockets.Find( pocket_name)->second.First();

            // compute new properties to associate with the pocket
            BCL_MessageVrb( "Saving properties of base type: " + util::Format()( m_AtomProperty.GetLabel()));
            m_MolB3DA = ( m_3DASmoothSign.SumOverObject( m_MolB));

            // save everything to our map of maps
            temp_map[ m_AtomProperty.GetLabel()] = m_MolB3DA;
            s_Pockets[ pocket_name] = std::make_pair( m_MolB, temp_map);
          }
          // If we have seen this pocket/property pair before, then get everything from the map
          else
          {
            BCL_MessageVrb( "Retrieving pocket from cache: " + pocket_name);
            m_MolB = s_Pockets.Find( pocket_name)->second.First();
            BCL_MessageVrb( "Retrieving properties of base type: " + util::Format()( m_AtomProperty.GetLabel()));
            m_MolB3DA = ( m_3DASmoothSign.SumOverObject( m_MolB));
          }

          // Require that the protein has atoms and properties
          BCL_Assert( m_MolB.GetNumberAtoms() > 0, "The protein indicated by the misc property id label contains no atoms");
          BCL_Assert( m_MolB3DA.GetSize() > 0, "The protein property vector is empty");

        }
        s_Mutex.Unlock();
      }
      // Require that the protein has atoms and properties
      BCL_Assert( m_MolB.GetNumberAtoms() > 0, "The protein indicated by the misc property id label contains no atoms");
      BCL_Assert( m_MolB3DA.GetSize() > 0, "The protein property vector is empty");
    }

    //! @brief reduce pocket to atoms contacting the ligand
    void Molecule3DAPairConvolution::ReduceProteinPocket( chemistry::FragmentComplete &POCKET)
    {
    }

  } // namespace descriptor
} // namespace bcl
