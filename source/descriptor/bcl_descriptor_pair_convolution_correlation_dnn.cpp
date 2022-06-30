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
#include "descriptor/bcl_descriptor_pair_convolution_correlation_dnn.h"

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
    PairConvolutionCorrelationDNN::PairConvolutionCorrelationDNN() :
        m_Molecule( chemistry::FragmentComplete()),
        m_MDLPocketName( "bcl_pocket_id_080"),
        m_ProteinPocketFilename( std::string()),
        m_ApplicabilityWeight( false)
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
    PairConvolutionCorrelationDNN::PairConvolutionCorrelationDNN
    (
      const bool &APPLICABILITY_WEIGHT
    ) :
      m_Molecule( chemistry::FragmentComplete()),
      m_MDLPocketName( "bcl_pocket_id_080"),
      m_ProteinPocketFilename( std::string()),
      m_ApplicabilityWeight( APPLICABILITY_WEIGHT)
    {
      BCL_Assert
      (
        ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
        "Failed to create " + GetClassIdentifier()
      );
    }

    //! @brief detailed constructor
    PairConvolutionCorrelationDNN::PairConvolutionCorrelationDNN
    (
      const std::string &MDL_PROPERTY,
      const std::string &POCKET_FILENAME,
      const bool &APPLICABILITY_WEIGHT
    ) :
      m_Molecule( chemistry::FragmentComplete()),
      m_MDLPocketName( MDL_PROPERTY),
      m_ProteinPocketFilename( POCKET_FILENAME),
      m_ApplicabilityWeight( APPLICABILITY_WEIGHT)
    {
      BCL_Assert
      (
        ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger()),
        "Failed to create " + GetClassIdentifier()
      );
    }

    //! @brief virtual copy constructor
    //! @return pointer to new PairConvolutionCorrelationDNN
    PairConvolutionCorrelationDNN *PairConvolutionCorrelationDNN::Clone() const
    {
      return new PairConvolutionCorrelationDNN( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &PairConvolutionCorrelationDNN::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &PairConvolutionCorrelationDNN::GetAlias() const
    {
      static const std::string s_name( "PCC-DNN"), s_ad_name( "PCC-AD-DNN");
       return m_ApplicabilityWeight ? s_ad_name : s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void PairConvolutionCorrelationDNN::Calculate( linal::VectorReference< float> &STORAGE)
    {
      // create the descriptor to generate the raw binding affinity prediction
      util::Implementation< Base< chemistry::AtomConformationalInterface, float> > plc_descriptor
      (
        "PredictionMean(storage=File(prefix=model,directory=" + model::Model::AddModelPath
        ( "cheminfo/qsar/pcc-dnn_casf2016_sans_coreset/") + "))"
      );

      // compute the predicted affinity
      float model_result( plc_descriptor->SumOverObject( m_Molecule)( 0));

      if( m_ApplicabilityWeight)
      {
        // create the descriptor to generate the applicability domain score
        util::Implementation< Base< chemistry::AtomConformationalInterface, float> > ad_descriptor
        (
          "Subtract(lhs=1,rhs=PredictionMean(storage=File(prefix=PCC,directory=" + model::Model::AddModelPath
          ( "cheminfo/qsar/pcc-dnn_casf2016_sans_coreset/") + ")))"
        );

        float ad_result( ad_descriptor->SumOverObject( m_Molecule)( 0));
        model_result *= ad_result;
      }

      // store final prediction
      STORAGE( 0) = model_result;
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void PairConvolutionCorrelationDNN::SetObjectHook()
    {
      // Get our molecule
      util::SiPtr< const chemistry::ConformationInterface> current_mol( this->GetCurrentObject());
      m_Molecule = *current_mol;

      // If no MDL specified then set to default
      if( m_MDLPocketName.empty())
      {
        m_MDLPocketName = "bcl_pocket_id_080";
      }

      // Set the filename as the MDL property value
      m_Molecule.GetStoredPropertiesNonConst().SetMDLProperty( m_MDLPocketName, m_ProteinPocketFilename);
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool PairConvolutionCorrelationDNN::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer PairConvolutionCorrelationDNN::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates binding affinity (in units of pKd)"
      );
      if( m_ApplicabilityWeight)
      {
        parameters.SetClassDescription
        (
          "Calculates binding affinity (in units of pKd) weighted by"
          " 1.0 minus the model applicability domain score"
        );
      }
      parameters.AddInitializer
      (
        "misc property id",
        "misc property name of corresponding protein binding pocket for current molecule",
        io::Serialization::GetAgent( &m_MDLPocketName),
        "bcl_pocket_id_080"
      );
      parameters.AddInitializer
      (
        "receptor filename",
        "the filename of the protein binding pocket of interest",
        io::Serialization::GetAgent( &m_ProteinPocketFilename),
        ""
      );
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
