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

#ifndef BCL_FOLD_MUTATE_AA_SEQUENCE_GROW_H_
#define BCL_FOLD_MUTATE_AA_SEQUENCE_GROW_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_function_interface_serializable.h"
#include "math/bcl_math_mutate_interface.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_si_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateAASequenceGrow
    //! @brief is for extending a sequence in the c-terminal direction.
    //! @details It extends a sequence residue by residue with given phi and psi angles.
    //!
    //! @see @link example_fold_mutate_aa_sequence_grow.cpp @endlink
    //! @author alexanns
    //! @date Sep 5, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateAASequenceGrow :
      public math::MutateInterface< biol::AASequence>
    {

    private:

    //////////
    // data //
    //////////

      //! the method that will produce phi and psi angles as the sequence grows
      util::ShPtr< math::FunctionInterfaceSerializable< MutationResidue, storage::VectorND< 2, double> > > m_PhiPsiGenerator;

      //! the starting sequence which will be grown onto
      util::SiPtr< const biol::AASequence> m_AnchorSequence;

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
      MutateAASequenceGrow();

      //! @brief constructor from member variables
      //! @param PHI_PSI_GENERATOR the method that will produce phi and psi angles as the sequence grows
      //! @param N_TERMINAL_AA_SEQUENCE the starting sequence which will be grown onto
      MutateAASequenceGrow
      (
        const util::ShPtr
        <
          math::FunctionInterfaceSerializable< MutationResidue, storage::VectorND< 2, double> >
        > &PHI_PSI_GENERATOR,
        const util::SiPtr< const biol::AASequence> &N_TERMINAL_AA_SEQUENCE
      );

      //! @brief Clone function
      //! @return pointer to new MutateAASequenceInitializeCoordinates
      MutateAASequenceGrow *Clone() const;

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
      //! @param SEQUENCE_TO_GROW Argument of interest
      //! @return MutateResult that results from mutating to the argument
      math::MutateResult< biol::AASequence> operator()( const biol::AASequence &SEQUENCE_TO_GROW) const;

      //! @brief grows the given sequence by connecting the cterminus to the anchor residue then growing nterminally
      //! @param SEQUENCE_TO_GROW sequence for which coordinates will be added
      //! @return MutateResult contains the sequence with newly assigned coordinates grown onto the anchor residue
      math::MutateResult< biol::AASequence> GrowTowardsNTerminus( const biol::AASequence &SEQUENCE_TO_GROW) const;

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

    }; // class MutateAASequenceGrow

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_AA_SEQUENCE_GROW_H_ 
