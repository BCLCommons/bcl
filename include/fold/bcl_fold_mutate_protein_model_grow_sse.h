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

#ifndef BCL_FOLD_MUTATE_PROTEIN_MODEL_GROW_SSE_H_
#define BCL_FOLD_MUTATE_PROTEIN_MODEL_GROW_SSE_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "biol/bcl_biol.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence_flexibility.h"
#include "find/bcl_find_locator_interface.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "math/bcl_math_mutate_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateProteinModelGrowSSE
    //! @brief grows an sse one residue at a time while trying to optimize placement of each additional residue
    //! @details Uses an sse of interest which has all the residues in it (with or without coordinates) and grows it
    //!          from the first residue to the last residue. For each residue a random phi psi angle is chosen and
    //!          applied. Several phi psi trials can be done and the best one according to scores will be chosen.
    //!          This class is most likely to be used in conjunction with loop building.
    //!
    //! @see @link example_fold_mutate_protein_model_grow_sse.cpp @endlink
    //! @author alexanns
    //! @date Jan 13, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateProteinModelGrowSSE :
      public math::MutateInterface< assemble::ProteinModel>
    {

    private:

    //////////
    // data //
    //////////

      //! the locator which will be used to find the sse which is going to be grown
      util::ShPtr< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> > m_LocatorSSE;

      //! the method that will be used in order to generate phi and psi angles as the loop domain is grown
      util::ShPtr< math::FunctionInterfaceSerializable< MutationResidue, storage::VectorND< 2, double> > > m_PhiPsiGenerator;

      //! locator to find the residue to which the n-terminus of the growing sse will be anchored (connected to)
      util::ShPtr< find::LocatorInterface< util::SiPtr< const biol::AABase>, assemble::ProteinModel> > m_AnchorAA;

      //! growing direction of the sse
      biol::AASequenceFlexibility::DirectionEnum m_GrowingDirection;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutateProteinModelGrowSSE();

      //! @brief constructor taking parameters
      //! @param LOCATOR_SSE the locator which will be used to find the sse which is going to be grown
      //! @param PHI_PSI_GEN method used in order to generate phi and psi angles as the loop domain is grown
      //! @param ANCHOR_AA locator to residue to which the n-terminus of the growing sse will be anchored (connected to)
      //! @param GROWING_DIRECTION growing direction of the sse
      MutateProteinModelGrowSSE
      (
        const util::ShPtr
        <
          find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface>
        > &LOCATOR_SSE,
        const util::ShPtr
        <
          math::FunctionInterfaceSerializable< MutationResidue, storage::VectorND< 2, double> >
        > &PHI_PSI_GEN,
        const util::ShPtr
        <
          find::LocatorInterface< util::SiPtr< const biol::AABase>, assemble::ProteinModel>
        > &ANCHOR_AA,
        const biol::AASequenceFlexibility::SequenceDirection GROWING_DIRECTION
      );

      //! @brief Clone function
      //! @return pointer to new MutateProteinModelGrowSSE
      MutateProteinModelGrowSSE *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
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

      //! @brief operator taking an ARGUMENT and returning a mutated object of t_ArgumentType
      //! @param PROTEIN_MODEL Argument of interest
      //! @return MutateResult that results from mutating to the argument
      math::MutateResult< assemble::ProteinModel> GrowTowardsCTerminus
      (
        const assemble::ProteinModel &PROTEIN_MODEL
      ) const;

      //! @brief operator taking an ARGUMENT and returning a mutated object of t_ArgumentType
      //! @param PROTEIN_MODEL Argument of interest
      //! @return MutateResult that results from mutating to the argument
      math::MutateResult< assemble::ProteinModel> GrowTowardsNTerminus
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

    public:

      //! @brief static function that grows all the missing coil regions in the given ProteinModel
      //!        splits the coil and grows first part N to C and the second part C to N
      //! @param PROTEIN_MODEL ProteinModel of interest
      //! @param PHI_PSI_GEN method used in order to generate phi and psi angles as the loop domain is grown
      //! @return ShPtr to ProteinModel with the loops grown
      static util::ShPtr< assemble::ProteinModel> GrowAllMissingCoilsBidirectional
      (
        const assemble::ProteinModel &PROTEIN_MODEL,
        const util::ShPtr< math::FunctionInterfaceSerializable< MutationResidue, storage::VectorND< 2, double> > > &PHI_PSI_GEN
      );

      //! @brief grows coordinates for specified sses
      //! @param PROTEIN_MODEL ProteinModel of interest
      //! @param PHI_PSI_GEN method used in order to generate phi and psi angles as the loop domain is grown
      //! @param LOCATOR_SSES the specific sses that will be grown
      //! @return ptr to protein model that has the grown sses as specifed
      static util::ShPtr< assemble::ProteinModel> GrowSpecifiedCoilsBidirectional
      (
        const assemble::ProteinModel &START_MODEL,
        const util::ShPtr< math::FunctionInterfaceSerializable< MutationResidue, storage::VectorND< 2, double> > > &PHI_PSI_GEN,
        const util::SiPtrVector< const assemble::LocatorSSE> &LOCATOR_SSES
      );

      //! @brief grows a given sse in a protein model
      //! @param SSE the sse in the protein model that will be grown
      //! @param MODEL the protein model in which the sse will be grown
      //! @param DIRECTION the direction the sse needs to be grown
      //! @param PHI_PSI_GEN the method for generating the phi and psi angles to be assigned to the sse
      static void GrowSSE
      (
        const assemble::SSE &SSE, util::ShPtr< assemble::ProteinModel> &MODEL,
        const biol::AASequenceFlexibility::SequenceDirection &DIRECTION,
        const util::ShPtr< math::FunctionInterfaceSerializable< MutationResidue, storage::VectorND< 2, double> > > &PHI_PSI_GEN
      );

    }; // class MutateProteinModelGrowSSE

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_PROTEIN_MODEL_GROW_SSE_H_ 
