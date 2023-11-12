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
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "command/bcl_command_command_state.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "descriptor/bcl_descriptor_protein_ligand_correlation_dnn.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_atom_sigma_charge.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProteinLigandCorrelationDNN::ProteinLigandCorrelationDNN() :
        m_Molecule( chemistry::FragmentComplete()),
        m_MDLPocketName( "receptor_filename"),
        m_ProteinPocketFilename( std::string()),
        m_ApplicabilityWeight( false),
        m_DockANNScore( false),
        m_ReportBindingFreeEnergy( false),
        m_Temperature( float( 298.75))
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

    //! @brief specify applicability domain model constructor
    ProteinLigandCorrelationDNN::ProteinLigandCorrelationDNN
    (
      const bool &APPLICABILITY_WEIGHT,
      const bool &DOCK_ANN_SCORE
    ) :
      m_Molecule( chemistry::FragmentComplete()),
      m_MDLPocketName( "receptor_filename"),
      m_ProteinPocketFilename( std::string()),
      m_ApplicabilityWeight( APPLICABILITY_WEIGHT),
      m_DockANNScore( DOCK_ANN_SCORE),
      m_ReportBindingFreeEnergy( false),
      m_Temperature( float( 298.75))
    {
      BCL_Assert
      (
        ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
        "Failed to create " + GetClassIdentifier()
      );
    }

    //! @brief detailed constructor
    ProteinLigandCorrelationDNN::ProteinLigandCorrelationDNN
    (
      const std::string &MDL_PROPERTY,
      const std::string &POCKET_FILENAME,
      const bool &APPLICABILITY_WEIGHT,
      const bool &DOCK_ANN_SCORE
    ) :
      m_Molecule( chemistry::FragmentComplete()),
      m_MDLPocketName( MDL_PROPERTY),
      m_ProteinPocketFilename( POCKET_FILENAME),
      m_ApplicabilityWeight( APPLICABILITY_WEIGHT),
      m_DockANNScore( DOCK_ANN_SCORE),
      m_ReportBindingFreeEnergy( false),
      m_Temperature( float( 298.75))
    {
      BCL_Assert
      (
        ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
        "Failed to create " + GetClassIdentifier()
      );
    }

    //! @brief detailed constructor with binding free energy conversion
    ProteinLigandCorrelationDNN::ProteinLigandCorrelationDNN
    (
      const std::string &MDL_PROPERTY,
      const std::string &POCKET_FILENAME,
      const bool &APPLICABILITY_WEIGHT,
      const bool &REPORT_DG_BIND,
      const bool &DOCK_ANN_SCORE,
      const float &TEMPERATURE
    ) :
      m_Molecule( chemistry::FragmentComplete()),
      m_MDLPocketName( MDL_PROPERTY),
      m_ProteinPocketFilename( POCKET_FILENAME),
      m_ApplicabilityWeight( APPLICABILITY_WEIGHT),
      m_DockANNScore( DOCK_ANN_SCORE),
      m_ReportBindingFreeEnergy( REPORT_DG_BIND),
      m_Temperature( TEMPERATURE)
    {
      BCL_Assert
      (
        ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
        "Failed to create " + GetClassIdentifier()
      );
    }

    //! @brief virtual copy constructor
    //! @return pointer to new ProteinLigandCorrelationDNN
    ProteinLigandCorrelationDNN *ProteinLigandCorrelationDNN::Clone() const
    {
      return new ProteinLigandCorrelationDNN( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &ProteinLigandCorrelationDNN::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &ProteinLigandCorrelationDNN::GetAlias() const
    {
      // make some fancy names for our models
      static const std::string s_name( "AffinityNet"), s_ad_name( "AffinityNetAD"), dock_s_name( "DockANNScore");

      // pick a model to use
      if( m_DockANNScore)
      {
        // computes product of AffinityNet score and PredictionInfo for index 0 of BCL-DockANNScore
        return dock_s_name;
      }
      else if( m_ApplicabilityWeight)
      {
        // computes product of AffinityNet score and 1 - AD score (limited at [0.0, 1.0])
        return s_ad_name;
      }
      else
      {
        // computes the AffinityNet score
        return s_name;
      }
    }

    // initialize static member data
    sched::Mutex ProteinLigandCorrelationDNN::s_Mutex = sched::Mutex();

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void ProteinLigandCorrelationDNN::Calculate( linal::VectorReference< float> &STORAGE)
    {
      // Get our molecule
      util::SiPtr< const chemistry::ConformationInterface> current_mol( this->GetCurrentObject());
      m_Molecule = *current_mol;

      // if not pocket filename specified, try to get it from MDL property
      s_Mutex.Lock();
      std::string receptor_filename;
      if( !m_ProteinPocketFilename.size())
      {
        // get property
        receptor_filename = m_Molecule.GetStoredProperties().GetMDLProperty( m_MDLPocketName);

        // Set the filename as the MDL property value
        m_Molecule.GetStoredPropertiesNonConst().SetMDLProperty( m_MDLPocketName, receptor_filename);
      }
      else
      {
        m_Molecule.GetStoredPropertiesNonConst().SetMDLProperty( m_MDLPocketName, m_ProteinPocketFilename);
      }

      // setup static affinity prediction descriptors
      static util::Implementation< Base< chemistry::AtomConformationalInterface, float> > affinity_net, affinity_net_ad_score;
      affinity_net =
          (
              "PredictionMean(storage=File(prefix=model,directory=" + model::Model::AddModelPath
              ( "cheminfo/qsar/BCL-AffinityNet/models/") + "))"
          );

      affinity_net_ad_score =
          (
              "Subtract(lhs=1,rhs=Limit(PredictionMean(storage=File(prefix=RS,directory=" + model::Model::AddModelPath
              ( "cheminfo/qsar/BCL-AffinityNet/") + ")),max=1.0,min=0.0))"
          );
      s_Mutex.Unlock();

      // compute the predicted affinity
      s_Mutex.Lock();
      float model_result( affinity_net->SumOverObject( m_Molecule)( 0));
      s_Mutex.Unlock();

      // compute the local ppv of pose being < 1.0 Angstrom from native
      if( m_DockANNScore)
      {
        s_Mutex.Lock();
        // create the descriptor to compute the docking score
        static util::Implementation< Base< chemistry::AtomConformationalInterface, float> > dock_ann_score
        (
          "Partial(PredictionInfo(predictor=File(prefix=model,directory=" + model::Model::AddModelPath
          ( "cheminfo/qsar/BCL-DockANNScore/models/") + ") ,metrics(LocalPPV)),indices(0))"
        );
        float dock_local_ppv( dock_ann_score->SumOverObject( m_Molecule)( 0));
        model_result *= dock_local_ppv;
        s_Mutex.Unlock();
      }
      // compute the AffinityNet AD score
      else if( m_ApplicabilityWeight)
      {
        s_Mutex.Lock();
        float ad_result( affinity_net_ad_score->SumOverObject( m_Molecule)( 0));
        model_result *= ad_result;
        s_Mutex.Unlock();
      }

      // conversion to binding free energy
      if( m_ReportBindingFreeEnergy && !m_DockANNScore)
      {
        model_result = logf( pow( 10.0, -1.0 * model_result)) * 0.001987 * m_Temperature;
      }

      // store final prediction
      STORAGE( 0) = model_result;
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void ProteinLigandCorrelationDNN::SetObjectHook()
    {
      // If no MDL specified then set to default
      // I decided that to avoid ugliness, the MDL property used must be the same across all the molecules in the SDF
      // so this is good to keep in the SetObjectHook
      s_Mutex.Lock();
      if( m_MDLPocketName.empty())
      {
        m_MDLPocketName = "receptor_filename";
      }
      s_Mutex.Unlock();
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool ProteinLigandCorrelationDNN::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProteinLigandCorrelationDNN::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates binding affinity (in units of pKd)"
      );
      if( m_DockANNScore)
      {
        parameters.SetClassDescription
        (
          "Calculates a docking score based on the local ppv classifying a pose as within 1.0 Angstroms "
          "of the native pose weighted by the predicted affinity of the pose"
        );
      }
      else if( m_ApplicabilityWeight)
      {
        parameters.SetClassDescription
        (
          "Calculates binding affinity (in units of pKd) weighted by"
          " 1.0 minus the model applicability domain score"
        );
      }
      parameters.AddInitializer
      (
        "misc_property_id",
        "misc property name of corresponding protein binding pocket for current molecule",
        io::Serialization::GetAgent( &m_MDLPocketName),
        "receptor_filename"
      );
      parameters.AddInitializer
      (
        "receptor_filename",
        "the filename of the protein binding pocket of interest",
        io::Serialization::GetAgent( &m_ProteinPocketFilename),
        ""
      );
      parameters.AddInitializer
      (
        "report_dg",
        "reports the result as binding free energy in kcal/mol at the specified temperature",
        io::Serialization::GetAgent( &m_ReportBindingFreeEnergy),
        "0"
      );
      parameters.AddInitializer
      (
        "temperature",
        "temperature for binding free energy conversion; no effect if report_dg = 0",
        io::Serialization::GetAgent( &m_Temperature),
        "298.75"
      );
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
