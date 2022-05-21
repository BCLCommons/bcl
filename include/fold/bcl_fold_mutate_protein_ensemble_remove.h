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

#ifndef BCL_FOLD_MUTATE_PROTEIN_ENSEMBLE_REMOVE_H_
#define BCL_FOLD_MUTATE_PROTEIN_ENSEMBLE_REMOVE_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_mutate_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateProteinEnsembleRemove
    //! @brief Mutates a protein ensemble by randomly removing one protein
    //! @details Randomly removes one protein from an ensemble in order to mutate the ensemble
    //!
    //! @see @link example_fold_mutate_protein_ensemble_remove.cpp @endlink
    //! @author alexanns
    //! @date Feb 12, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateProteinEnsembleRemove :
      public math::MutateInterface< assemble::ProteinEnsemble>
    {

    private:

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      MutateProteinEnsembleRemove();

      //! @brief Clone function
      //! @return pointer to new MutateProteinEnsembleRemove
      MutateProteinEnsembleRemove *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

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
      operator()( const assemble::ProteinEnsemble &PROTEIN_ENSEMBLE) const;

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

    private:

    }; // class MutateProteinEnsembleRemove

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_PROTEIN_ENSEMBLE_REMOVE_H_ 
