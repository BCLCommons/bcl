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

#ifndef BCL_FOLD_MUTATE_PROTEIN_MODEL_LOOP_DOMAIN_GROW_H_
#define BCL_FOLD_MUTATE_PROTEIN_MODEL_LOOP_DOMAIN_GROW_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "find/bcl_find_locator_interface.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "math/bcl_math_mutate_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateProteinModelLoopDomainGrow
    //! @brief is for growing the non rigid portions of a loop domain into a protein model
    //! @details It grows the non rigid portions in the N- to C- terminus direction and connect the rigid portions
    //! to the grown regions of the loop domain as necessary.
    //!
    //! @see @link example_fold_mutate_protein_model_loop_domain_grow.cpp @endlink
    //! @author alexanns
    //! @date Sep 8, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateProteinModelLoopDomainGrow :
      public math::MutateInterface< assemble::ProteinModel>
    {

    private:

    //////////
    // data //
    //////////

      //! the locator which will be used to find a Loop Domain in a protein model
      util::ShPtr< find::LocatorInterface< util::ShPtr< LoopDomain>, assemble::DomainInterface> > m_LocatorLoopDomain;

      //! the method that will be used in order to generate phi and psi angles as the loop domain is grown
      util::ShPtr< math::FunctionInterfaceSerializable< MutationResidue, storage::VectorND< 2, double> > > m_PhiPsiGenerator;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutateProteinModelLoopDomainGrow();

      //! @brief constructor taking member variable parameters
      //! @param LOOP_DOMAIN_LOCATOR the locator which will be used to find a Loop Domain in a protein model
      //! @param PHI_PSI_GENERATOR method to be used in order to generate phi and psi angles as the loop domain is grown
      MutateProteinModelLoopDomainGrow
      (
        const util::ShPtr< find::LocatorInterface< util::ShPtr< LoopDomain>, assemble::DomainInterface> > &LOOP_DOMAIN_LOCATOR,
        const util::ShPtr
        <
          math::FunctionInterfaceSerializable< MutationResidue, storage::VectorND< 2, double> >
        > &PHI_PSI_GENERATOR
      );

      //! @brief Clone function
      //! @return pointer to new MutateProteinModelLoopDomainGrow
      MutateProteinModelLoopDomainGrow *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator taking an ARGUMENT and returning a mutated object of t_ArgumentType
      //! @param PROTEIN_MODEL Argument of interest
      //! @return MutateResult that results from mutating to the argument
      math::MutateResult< assemble::ProteinModel> operator()
      (
        const assemble::ProteinModel &PROTEIN_MODEL
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief ConnectNonRigidSSE connects a non rigid c-terminal sse to an n-terminal anchor sse
      //! @param N_TERMINAL_SSE the sse which will have a c-terminal portion connected to it
      //! @param C_TERMINAL_SSE the non-rigid c-terminal sse which will be connected to "N_TERMINAL_SSE"
      util::ShPtr< assemble::SSE> ConnectNonRigidSSE
      (
        const assemble::SSE &N_TERMINAL_SSE, const assemble::SSE &C_TERMINAL_SSE
      ) const;

      //! @brief ConnectRigidSequence connects a non rigid c-terminal sse to an n-terminal anchor sse
      //! @param N_TERMINAL_SSE the sse which will have a c-terminal portion connected to it
      //! @param C_TERMINAL_SEQUENCE
      //! @param C_TERMINAL_SSES the non-rigid c-terminal sse which will be connected to "N_TERMINAL_SSE"
      storage::Pair< util::ShPtr< assemble::SSE>, storage::List< LoopSegment> > ConnectRigidSequence
      (
        const assemble::SSE &N_TERMINAL_SSE,
        const biol::AASequence &C_TERMINAL_SEQUENCE,
        const storage::List< LoopSegment> &C_TERMINAL_SSES
      ) const;

    }; // class MutateProteinModelLoopDomainGrow

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_PROTEIN_MODEL_LOOP_DOMAIN_GROW_H_
