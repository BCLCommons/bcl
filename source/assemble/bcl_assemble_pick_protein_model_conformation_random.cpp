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
#include "assemble/bcl_assemble_pick_protein_model_conformation_random.h"

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
    const util::SiPtr< const util::ObjectInterface> PickProteinModelConformationRandom::s_Instance
    (
      util::Enumerated< find::PickInterface< util::SiPtr< const ProteinModel>, util::SiPtrList< const ProteinModel> > >::AddInstance( new PickProteinModelConformationRandom())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PickProteinModelConformationRandom::PickProteinModelConformationRandom()
    {
    }

    //! @brief Clone function
    //! @return pointer to new PickProteinModelConformationRandom
    PickProteinModelConformationRandom *PickProteinModelConformationRandom::Clone() const
    {
      return new PickProteinModelConformationRandom( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PickProteinModelConformationRandom::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &PickProteinModelConformationRandom::GetAlias() const
    {
      static const std::string s_alias( "PickProteinModelConformationRandom");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer PickProteinModelConformationRandom::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Picks a random conformation of a protein model.");

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Pick returns a random SSE from the domain argument
    //! @param PROTEIN_MODEL model which the locator will Pick from
    //! @return returns SiPtr to protein model conformation from protein model conformation ensemble
    util::SiPtr< const ProteinModel>
    PickProteinModelConformationRandom::Pick( const util::SiPtrList< const ProteinModel> &ENSEMBLE) const
    {
      // if the ensemble is empty return empty siptr
      if( ENSEMBLE.IsEmpty())
      {
        return util::SiPtr< const ProteinModel>();
      }

      // get random index between 0 and ensemble size
      size_t index( random::GetGlobalRandom().Random( size_t( 0), size_t( ENSEMBLE.GetSize() - 1)));

      // get SiPtr to the index of the ensemble
      util::SiPtrList< const ProteinModel>::const_iterator
        ensemble_itr( ENSEMBLE.Begin()), ensemble_itr_end( ENSEMBLE.End());

      storage::AdvanceIterator( ensemble_itr, ensemble_itr_end, index);
      const util::SiPtr< const ProteinModel> located_model( *ensemble_itr);

      return located_model;
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
