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
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "fold/bcl_fold_mutate_protein_model_switch_conformation.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateProteinModelSwitchConformation::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateProteinModelSwitchConformation())
    );
  
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelSwitchConformation::MutateProteinModelSwitchConformation() :
      m_ConformationLocator(),
      m_Scheme( GetStaticClassName< MutateProteinModelSwitchConformation>())
    {
    }

    //! @brief constructor taking parameters
    //! @param LOCATOR the method for locating an alternative conformation
    //! @param SCHEME Scheme to be used
    MutateProteinModelSwitchConformation::MutateProteinModelSwitchConformation
    (
      const util::ShPtr
      <
        find::LocatorInterface< util::SiPtr< const assemble::ProteinModel>, assemble::ProteinModel>
      > &CONFORMATION,
      const std::string &SCHEME
    ) :
      m_ConformationLocator( CONFORMATION),
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateProteinModelSwitchConformation
    MutateProteinModelSwitchConformation *MutateProteinModelSwitchConformation::Clone() const
    {
      return new MutateProteinModelSwitchConformation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateProteinModelSwitchConformation::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
    //! @param PROTEIN_MODEL protein model interest
    //! @return MutateResult with ProteinModel after the mutate
    math::MutateResult< assemble::ProteinModel>
    MutateProteinModelSwitchConformation::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // static empty model
      const math::MutateResult< assemble::ProteinModel> empty_result( util::ShPtr< assemble::ProteinModel>(), *this);

      // try to locate an alternative conformation
      util::SiPtr< const assemble::ProteinModel> located_conformation( m_ConformationLocator->Locate( PROTEIN_MODEL));

      // if no new conformation could be found
      if( !located_conformation.IsDefined())
      {
        return empty_result;
      }

      // copy the protein model
      util::ShPtr< assemble::ProteinModel> new_model( located_conformation->Clone());

      // copy the conformational ensemble from PROTEIN_MODEL
      assemble::ProteinEnsemble conformations( PROTEIN_MODEL.GetConformationalEnsemble());

      // find the index of the located conformation in the list of conformations
      const size_t index( conformations.GetEnsembleData().Find( located_conformation));

      // if the index of the located conformation could not be determined
      if( index >= conformations.GetSize())
      {
        BCL_MessageDbg( "could not determine index of located conformation");
        return empty_result;
      }

      // remove the located conformation from the list of them since it will be the conformation represented by the
      // protein model
      conformations.GetEnsembleData().RemoveElements( index);

      // clone PROTEIN_MODEL
      util::ShPtr< assemble::ProteinModel> old_model_copy( PROTEIN_MODEL.Clone());

      // copy the conformational ensemble from old_model_copy
      assemble::ProteinEnsemble old_model_copy_conformations( PROTEIN_MODEL.GetConformationalEnsemble());

      // reset the conformations in old_model_copy_conformations
      old_model_copy_conformations.Reset();

      // set the conformations in old_model_copy so that it does not have any conformational ensemble on its own
      old_model_copy->SetConformationalEnsemble( old_model_copy_conformations);

      // add PROTEIN_MODEL to the conformational ensemble
      conformations.InsertElement( old_model_copy);

      // set the conformations in the new model
      new_model->SetConformationalEnsemble( conformations);

      // return the new model
      return math::MutateResult< assemble::ProteinModel>( new_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateProteinModelSwitchConformation::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ConformationLocator, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateProteinModelSwitchConformation::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ConformationLocator, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace fold
} // namespace bcl
