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
#include "assemble/bcl_assemble_collector_protein_model_conformations.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> CollectorProteinModelConformations::s_Instance
    (
      util::Enumerated< find::CollectorInterface< util::SiPtrList< const ProteinModel>, ProteinModel> >::AddInstance
      (
        new CollectorProteinModelConformations()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    CollectorProteinModelConformations::CollectorProteinModelConformations() :
      m_ConsiderCurrentConformation( false)
    {
    }

    //! @brief constructor taking parameters
    //! @param CONSIDER_CURRENT if true the current conformation will be considered with the other conformations
    CollectorProteinModelConformations::CollectorProteinModelConformations
    (
      const bool CONSIDER_CURRENT
    ) :
      m_ConsiderCurrentConformation( CONSIDER_CURRENT)
    {
    }

    //! @brief Clone function
    //! @return pointer to new CollectorProteinModelConformations
    CollectorProteinModelConformations *CollectorProteinModelConformations::Clone() const
    {
      return new CollectorProteinModelConformations( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CollectorProteinModelConformations::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &CollectorProteinModelConformations::GetAlias() const
    {
      static const std::string s_alias( "CollectorProteinModelConformations");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer CollectorProteinModelConformations::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Collects all conformations of a protein.");
      serializer.AddInitializer
      (
        "consider current",
        "consider current conformation",
        io::Serialization::GetAgent( &m_ConsiderCurrentConformation)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! Collect the t_ReturnType objects in t_ArgumentType
    //! @param PROTEIN_MODEL entity that contains a t_ReturnType
    //! @return returns Group of the collected t_ReturnType objects
    util::SiPtrList< const ProteinModel>
    CollectorProteinModelConformations::Collect( const ProteinModel &PROTEIN_MODEL) const
    {
      // get the conformations in the protein model
      const ProteinEnsemble &conformation_ensemble( PROTEIN_MODEL.GetConformationalEnsemble());

      // to hold the list of conformations from protein model
      util::SiPtrList< const ProteinModel> ensemble_list;

      // fill up the list of conformations from the protein model
      for
      (
        util::ShPtrVector< ProteinModel>::const_iterator
          model_itr( conformation_ensemble.GetEnsembleData().Begin()),
          model_itr_end( conformation_ensemble.GetEnsembleData().End());
        model_itr != model_itr_end;
        ++model_itr
      )
      {
        ensemble_list.PushBack( util::SiPtr< const ProteinModel>( *model_itr));
      }

      // true if want to add the current conformation of the protein model
      if( m_ConsiderCurrentConformation)
      {
        ensemble_list.PushBack( util::ToSiPtr( PROTEIN_MODEL));
      }

      // return the list of conformations
      return ensemble_list;
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace assemble
} // namespace bcl
