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
#include "fold/bcl_fold_mutate_protein_ensemble_remove.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateProteinEnsembleRemove::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateProteinEnsembleRemove())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinEnsembleRemove::MutateProteinEnsembleRemove()
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateProteinEnsembleRemove
    MutateProteinEnsembleRemove *MutateProteinEnsembleRemove::Clone() const
    {
      return new MutateProteinEnsembleRemove( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateProteinEnsembleRemove::GetClassIdentifier() const
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
    //! @param PROTEIN_ENSEMBLE ensemble of proteins to be mutated
    //! @return MutateResult with ProteinEnsemble after the mutate
    math::MutateResult< assemble::ProteinEnsemble>
    MutateProteinEnsembleRemove::operator()( const assemble::ProteinEnsemble &PROTEIN_ENSEMBLE) const
    {
      // returned if the mutate cannot occur
      static util::ShPtr< assemble::ProteinEnsemble> s_empty;

      // true if the protein ensemble has size 0 or 1, should not remove anything
      if( PROTEIN_ENSEMBLE.GetSize() < 2)
      {
        return math::MutateResult< assemble::ProteinEnsemble>( s_empty, *this);
      }

      // make a new copy of the current ensemble
      util::ShPtr< assemble::ProteinEnsemble> new_ensemble( PROTEIN_ENSEMBLE.Clone());

      // remove a random model from the current ensemble
      new_ensemble->RemoveElement
      (
        random::GetGlobalRandom().Iterator( new_ensemble->Begin(), new_ensemble->End(), new_ensemble->GetSize())
      );

      // return mutate result with new ensemble
      return math::MutateResult< assemble::ProteinEnsemble>( new_ensemble, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateProteinEnsembleRemove::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateProteinEnsembleRemove::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace fold
  
} // namespace bcl
