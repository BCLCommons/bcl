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
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "descriptor/bcl_descriptor_molecule_xpka.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_command_state.h"
#include "descriptor/bcl_descriptor_atom_sigma_charge.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor to initialize the class
    //! @param  ACID_MODEL true if acid model to be used. Initialized to true.
    MoleculeXpKa::MoleculeXpKa( const bool &ACID_MODEL) :
        m_AcidModel( ACID_MODEL)
      {
      }

    //! @brief virtual copy constructor
    //! @return pointer to new MoleculeXpKa
    MoleculeXpKa *MoleculeXpKa::Clone() const
    {
      return new MoleculeXpKa( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &MoleculeXpKa::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
//    const std::string &MoleculeXpKa::GetAlias() const
//    {
//      static const std::string s_name( "XpKa");
//      return s_name;
//    }

    const std::string &MoleculeXpKa::GetAlias() const
    {
      static const std::string s_acid_name( "XpKaAcid"), s_base_name( "XpKaBase");
      return m_AcidModel ? s_acid_name : s_base_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void MoleculeXpKa::Calculate( linal::VectorReference< float> &STORAGE)
    {
      // Get our molecule
      util::SiPtr< const chemistry::ConformationInterface> current_mol( this->GetCurrentObject());
      chemistry::FragmentComplete molecule( *current_mol);

//    // PREDICTIONMEAN - create the descriptor to generate the raw prediction values
//      static util::Implementation< Base< chemistry::AtomConformationalInterface, float> > descriptor
//      (
//        "PredictionMean(storage=File(prefix=model,directory=" + model::Model::AddModelPath
//        ( "cheminfo/qspr/ionizability/chembl_all_acid_x128_32.005_025.05.class_Sigmoid/") + "))"
//      );
//
//      linal::Vector< float> model_results( descriptor->SumOverObject( molecule));
//      STORAGE( 0) = model_results( 0);

      // PREDICTIONINFO - create the descriptor to generate the raw prediction values
//      if( m_AcidModel)
//      {
//        static util::Implementation< Base< chemistry::AtomConformationalInterface, float> > descriptor
//        (
//          "PredictionInfo(predictor=File(directory=" + model::Model::AddModelPath
//          ( "cheminfo/qspr/ionizability/chembl_all_acid_x128_32.005_025.05.class_Sigmoid/") + ", prefix=model),metrics(LocalPPV))"
//        );
//        linal::Vector< float> model_results( descriptor->SumOverObject( molecule));
//        BCL_Assert
//        (
//          model_results.GetSize() == GetNormalSizeOfFeatures(),
//          "Error! Size of features returned by model is not equal to " + util::Format()( GetNormalSizeOfFeatures())
//        );
//        for
//        (
//            size_t itr( 0);
//            itr != model_results.GetSize();
//            ++itr
//        )
//          if( model_results( itr) == model_results.Max())
//          {
//            STORAGE = itr;
//          }
//      }
//      else
//      {
//        static util::Implementation< Base< chemistry::AtomConformationalInterface, float> > descriptor
//        (
//          "PredictionInfo(predictor=File(directory=" + model::Model::AddModelPath
//          ( "cheminfo/qspr/ionizability/chembl_all_base_pred_allBins.class.realx_32.005_025/") + ", prefix=model),metrics(LocalPPV))"
//        );
//        linal::Vector< float> model_results( descriptor->SumOverObject( molecule));
//        BCL_Assert
//        (
//          model_results.GetSize() == GetNormalSizeOfFeatures(),
//          "Error! Size of features returned by model is not equal to " + util::Format()( GetNormalSizeOfFeatures())
//        );
//        for
//        (
//            size_t itr( 0);
//            itr != model_results.GetSize();
//            ++itr
//        )
//          if( model_results( itr) == model_results.Max())
//          {
//            STORAGE = itr;
//          }
//      }

      static util::Implementation< Base< chemistry::AtomConformationalInterface, float> > descriptor
      (
        "PredictionInfo(predictor=File(directory=" + model::Model::AddModelPath
        (
          m_AcidModel ?
              "cheminfo/qspr/ionizability/chembl_all_acid_x128_32.005_025.05.class_Sigmoid/" :
              "cheminfo/qspr/ionizability/chembl_all_base_pred_allBins.class.realx_32.005_025/"
        ) + ", prefix=model),metrics(LocalPPV))"
      );

      linal::Vector< float> model_results( descriptor->SumOverObject( molecule));

//      BCL_Assert
//      (
//        model_results.GetSize() == GetNormalSizeOfFeatures(),
//        "Error! Size of features returned by model is not equal to " + util::Format()( GetNormalSizeOfFeatures())
//      );
//      STORAGE = linal::VectorReference< float>( model_results.Begin(), model_results.End());
      for
      (
          size_t itr( 0);
          itr != model_results.GetSize();
          ++itr
      )
      {
        // store each localPPV:
        STORAGE( itr) = model_results( itr);
//
////        //store bin with max localPPV:
////        if( model_results( itr) == model_results.Max())
////        {
////          STORAGE = itr;
////        }
//
      }
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool MoleculeXpKa::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      if( command::CommandState::IsInStaticInitialization())
      {
        return true;
      }
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeXpKa::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates pKa using a multi-tasking deep neural network trained on "
        + std::string( m_AcidModel ? "acids" : "bases")
      );
      //set which model to use
//      parameters.AddInitializer
//      (
//        "model_selector",
//        "True for Acid, False for Base",
//        io::Serialization::GetAgent(&m_AcidModel),
//        "True"
//      );
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
