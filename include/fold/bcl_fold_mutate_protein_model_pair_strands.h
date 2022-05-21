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

#ifndef BCL_FOLD_MUTATE_PROTEIN_MODEL_PAIR_STRANDS_H_
#define BCL_FOLD_MUTATE_PROTEIN_MODEL_PAIR_STRANDS_H_

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
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateProteinModelPairStrands
    //! @brief mutate for pairing unpaired strands in protein models
    //! @details This mutates locates unpaired strands in the given ProteinModel. It there are 2 or more then it pairs
    //! them with each other. If there is only one unpaired strand and one or more sheets, then it moves the unpaired
    //! strand and adds its to a randomly selected beta-sheet.
    //!
    //! @see @link example_fold_mutate_protein_model_pair_strands.cpp @endlink
    //! @author karakam
    //! @date May 4, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateProteinModelPairStrands :
      public math::MutateInterface< assemble::ProteinModel>
    {

    private:

    //////////
    // data //
    //////////

      //! scheme of this mutate
      std::string m_Scheme;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor, uses get static class name as scheme name
      MutateProteinModelPairStrands();

      //! @brief constructor from a scheme
      //! @param SCHEME scheme of this mutate
      explicit MutateProteinModelPairStrands( const std::string &SCHEME);

      //! @brief Clone function
      //! @return pointer to new MutateProteinModelPairStrands
      MutateProteinModelPairStrands *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns scheme of this mutate
      //! @return scheme of this mutate
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking an ARGUMENT and returning a mutate object of t_ArgumentType
      //! @param PROTEIN_MODEL protein model interest
      //! @return MutateResult with ProteinModel after the mutate
      math::MutateResult< assemble::ProteinModel> operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

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

    }; // class MutateProteinModelPairStrands

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_PROTEIN_MODEL_PAIR_STRANDS_H_ 
