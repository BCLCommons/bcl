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

#ifndef BCL_FOLD_MUTATE_PROTEIN_MODEL_ADD_SHEET_FROM_TEMPLATE_H_
#define BCL_FOLD_MUTATE_PROTEIN_MODEL_ADD_SHEET_FROM_TEMPLATE_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_fold_placement_domain.h"
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_fold_template.h"
#include "math/bcl_math_mutate_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateProteinModelAddSheetFromTemplate
    //! @brief This mutate was designed to be used in the de novo folding process. It selects a random number of strands
    //! from the SSE pool and fits them to a sheet template in the database before adding them to the protein model.
    //!
    //! @author fischea
    //! @date April 17, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateProteinModelAddSheetFromTemplate :
      public math::MutateInterface< assemble::ProteinModel>
    {

    //////////
    // data //
    //////////

    public:

      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    private:

      util::ShPtr< PlacementInterface< assemble::Domain, assemble::ProteinModel> > m_Placer;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief construct from placer
      //! @param PLACER object used to place the sheet in the protein model
      MutateProteinModelAddSheetFromTemplate
      (
        const PlacementInterface< assemble::Domain, assemble::ProteinModel> &PLACER = PlacementDomain()
      );

      //! @brief returns a pointer to a new MutateProteinModelAddSheetFromTemplate
      //! @return pointer to a new MutateProteinModelAddSheetFromTemplate
      MutateProteinModelAddSheetFromTemplate *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

      //! @brief returns the scheme for this mutate
      //! @return the scheme for this mutate
      const std::string &GetScheme() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief adds a sheet to the given protein model and returns the mutated model
      //! @detail a random number of strands is selected from the SSE pool of the given model and fitted to a sheet
      //! template in the library before being added to the given protein model
      //! @param MODEL protein model to add the sheet to
      //! @return the mutated protein model with additional information regarding the mutation
      math::MutateResult< assemble::ProteinModel> operator()( const assemble::ProteinModel &MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read object from input stream
      //! @param ISTREAM input stream to read object from
      //! @return input stream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write object into  output stream
      //! @param OSTREAM output stream to write object into
      //! @param INDENT number of indentations to separate members
      //! @return output stream object was written into
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class MutateProteinModelAddSheetFromTemplate

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_PROTEIN_MODEL_ADD_SHEET_FROM_TEMPLATE_H_
