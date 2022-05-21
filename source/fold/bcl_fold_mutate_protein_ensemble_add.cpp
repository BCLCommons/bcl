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
#include "fold/bcl_fold_mutate_protein_ensemble_add.h"

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
    const util::SiPtr< const util::ObjectInterface> MutateProteinEnsembleAdd::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateProteinEnsembleAdd())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinEnsembleAdd::MutateProteinEnsembleAdd() :
      m_ProteinPool()
    {
    }

    //! @brief constructor taking a pool of all possible protein models that could be swapped
    //! @param PROTEIN_POOL a pool of all possible protein models that could be swapped
    MutateProteinEnsembleAdd::MutateProteinEnsembleAdd
    (
      const util::ShPtr< assemble::ProteinEnsemble> &PROTEIN_POOL
    ) :
      m_ProteinPool( PROTEIN_POOL)
    {
      // sort the ensemble
      m_ProteinPool->Sort( std::less< util::ShPtr< assemble::ProteinModel> >());
    }

    //! @brief Clone function
    //! @return pointer to new MutateProteinEnsembleSwap
    MutateProteinEnsembleAdd *MutateProteinEnsembleAdd::Clone() const
    {
      return new MutateProteinEnsembleAdd( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MutateProteinEnsembleAdd::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
    //! @param PROTEIN_ENSEMBLE ensemble of proteins to be mutated
    //! @return MutateResult with ProteinEnsemble after the mutate
    math::MutateResult< assemble::ProteinEnsemble>
    MutateProteinEnsembleAdd::operator()( const assemble::ProteinEnsemble &PROTEIN_ENSEMBLE) const
    {
      // returned if the mutate cannot occur
      static util::ShPtr< assemble::ProteinEnsemble> s_empty;

      // true if pool is undefined
      if( !m_ProteinPool.IsDefined())
      {
        return math::MutateResult< assemble::ProteinEnsemble>( s_empty, *this);
      }

      // true if the pool is empty
      if( m_ProteinPool->IsEmpty())
      {
        return math::MutateResult< assemble::ProteinEnsemble>( s_empty, *this);
      }

      // make a new copy of the current ensemble
      util::ShPtr< assemble::ProteinEnsemble> new_ensemble( PROTEIN_ENSEMBLE.Clone());

      // true if there are no models currently in the ensemble
      if( PROTEIN_ENSEMBLE.IsEmpty())
      {
        // can just add random model from m_ProteinPool
        new_ensemble->InsertElement
        (
          *random::GetGlobalRandom().Iterator( m_ProteinPool->Begin(), m_ProteinPool->End(), m_ProteinPool->GetSize())
        );

        // return the mutated ensemble
        return math::MutateResult< assemble::ProteinEnsemble>( new_ensemble, *this);
      }

      // holds the models that are in the pool but not in "PROTEIN_ENSEMBLE"
      assemble::ProteinEnsemble swappable_models
      (
        m_ProteinPool->Difference
        (
          *new_ensemble, std::less< util::ShPtr< assemble::ProteinModel> >()
        )
      );

      // true if the there are no models eligible to be swapped (maybe PROTEIN_ENSEMBLE already contains all models)
      if( swappable_models.IsEmpty())
      {
        BCL_MessageDbg( "no models to swap");
        return math::MutateResult< assemble::ProteinEnsemble>( s_empty, *this);
      }

      // add a random model from the pool of swappable models into the new ensemble
      new_ensemble->InsertElement
      (
        *random::GetGlobalRandom().Iterator( swappable_models.Begin(), swappable_models.End(), swappable_models.GetSize())
      );

      // return the mutated ensemble
      return math::MutateResult< assemble::ProteinEnsemble>( new_ensemble, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateProteinEnsembleAdd::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ProteinPool, ISTREAM);

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateProteinEnsembleAdd::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ProteinPool, OSTREAM, INDENT);

      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
