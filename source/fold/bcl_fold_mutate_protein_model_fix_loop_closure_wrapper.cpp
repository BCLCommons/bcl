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
#include "fold/bcl_fold_mutate_protein_model_fix_loop_closure_wrapper.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_geometry_packing.h"
#include "assemble/bcl_assemble_voxel_grid_aa.h"
#include "fold/bcl_fold_mutate_protein_model_sse_pair_clash.h"
#include "fold/bcl_fold_mutate_protein_model_sse_pair_fix_loop_closure.h"
#include "io/bcl_io_directory.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "score/bcl_score_aa_pair_hi_res_clash.h"
#include "score/bcl_score_loop_closure.h"
#include "score/bcl_score_protein_model_membrane_topology.h"
#include "score/bcl_score_protein_model_sse_linear_loop_proximity.h"
// external includes - sorted alphabetically

// Define this to visualize all clash removals
//#define DEBUG_BCL_FOLD_MUTATE_CLASH_RESOLVER

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateProteinModelFixLoopClosureWrapper::MutateProteinModelFixLoopClosureWrapper
    (
      const math::MutateInterface< assemble::ProteinModel> &MUTATE_A,
      const std::string &SCHEME
    ) :
      m_PrimaryMutate
      (
        MUTATE_A,
        MutateProteinModelSSEPairClash(),
        false,
        MUTATE_A.GetScheme()
      )
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutateProteinModelFixLoopClosureWrapper
    MutateProteinModelFixLoopClosureWrapper *MutateProteinModelFixLoopClosureWrapper::Clone() const
    {
      return new MutateProteinModelFixLoopClosureWrapper( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateProteinModelFixLoopClosureWrapper::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
    //! @param PROTEIN_MODEL protein model interest
    //! @return MutateResult with ProteinModel after the mutate
    math::MutateResult< assemble::ProteinModel> MutateProteinModelFixLoopClosureWrapper::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      static util::ShPtr< assemble::ProteinModel> s_empty_model;
      // initialize empty model

      // max times to retry mutate. Since mutates are very fast relative to scoring, it is fast to run a given mutate
      // many times before giving up
      static size_t s_MaxRetries( 10);

      // loop closure score
      static score::LoopClosure s_closure( size_t( 1), 0.0, 1.0);

      // sse clash score
      static score::AAPairHiResClash s_Clash;

      // linear loop proximity score
      static score::ProteinModelSSELinearLoopProximity s_lin_loop_prox;

      // membrane topology score
      static score::ProteinModelMembraneTopology s_topology;
      static bool s_MPTopologyInitialized( s_topology.InitializeFromFlag());

      for( size_t i( 0); i < s_MaxRetries; ++i)
      {
        math::MutateResult< assemble::ProteinModel> result( m_PrimaryMutate( PROTEIN_MODEL));
        if( !result.GetArgument().IsDefined())
        {
          return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
        }
        if
        (
          s_closure( *result.GetArgument()) < 0.5
          && s_lin_loop_prox( *result.GetArgument()) < 1.5
          && s_Clash( *result.GetArgument()) < 0.5
          && ( !s_MPTopologyInitialized || s_topology( *result.GetArgument()) < 0.5)
        )
        {
          return math::MutateResult< assemble::ProteinModel>( result.GetArgument(), *this);
        }
      }

      // return
      return math::MutateResult< assemble::ProteinModel>( s_empty_model, *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateProteinModelFixLoopClosureWrapper::Read( std::istream &ISTREAM)
    {
      // return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &MutateProteinModelFixLoopClosureWrapper::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
