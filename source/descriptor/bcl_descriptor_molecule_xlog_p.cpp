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
#include "descriptor/bcl_descriptor_molecule_xlog_p.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_command_state.h"
#include "descriptor/bcl_descriptor_atom_sigma_charge.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief virtual copy constructor
    //! @return pointer to new MoleculeXLogP
    MoleculeXLogP *MoleculeXLogP::Clone() const
    {
      return new MoleculeXLogP( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &MoleculeXLogP::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &MoleculeXLogP::GetAlias() const
    {
      static const std::string s_name( "XLogP");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void MoleculeXLogP::Calculate( linal::VectorReference< float> &STORAGE)
    {
      // Get our molecule
      util::SiPtr< const chemistry::ConformationInterface> current_mol( this->GetCurrentObject());
      chemistry::FragmentComplete molecule( *current_mol);

      // create the descriptor to generate the raw prediction values
      static util::Implementation< Base< chemistry::AtomConformationalInterface, float> > descriptor
      (
        "PredictionMean(storage=File(prefix=model,directory=" + model::Model::AddModelPath
        ( "cheminfo/qspr/solubility/all_logp_logs_dgsolv.Scalar_Mol2D.rand.2x128-32_005_025_005/") + "))"
      );

      linal::Vector< float> model_results( descriptor->SumOverObject( molecule));
      STORAGE( 0) = model_results( 0);
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool MoleculeXLogP::ReadInitializerSuccessHook
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
    io::Serializer MoleculeXLogP::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates LogP using a multi-tasking deep neural network trained to predict LogP, LogS, and hydration free energy"
      );
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
