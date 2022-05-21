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

#ifndef BCL_FOLD_MUTATE_SHEET_FIT_TO_TEMPLATE_H_
#define BCL_FOLD_MUTATE_SHEET_FIT_TO_TEMPLATE_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sheet_template_handler.h"
#include "math/bcl_math_mutate_interface.h"
#include "storage/bcl_storage_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateSheetFitToTemplate
    //! @brief This mutate fits a given sheet to a randomly selected sheet template
    //! @details This mutate allows fitting of sheets into sheet templates. the sheet template is thereby chosen from a
    //! database of sheets un der the condictions that the template has the same number of strands and those strands
    //! have a similar size like the strands that shall be fitted. This approach allows bending of all strands in a
    //! beta sheet by fitting the geometries of all strands to the geometries represented in a template
    //!
    //! @see @link example_fold_mutate_sheet_fit_to_template.cpp @endlink
    //! @author karakam, fischea
    //! @date Oct 21, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateSheetFitToTemplate :
      public math::MutateInterface< assemble::Domain>
    {

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief returns a pointer to a new MutateSheetFitToTemplate
      //! @return pointer to a new MutateSheetFitToTemplate
      MutateSheetFitToTemplate *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the name of this class
      //! @return the name of this class
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief fits the given sheet to a randomly selected template
      //! @param SHEET sheet to be fitted
      //! @return the fitted sheet and a boolean indicating if the application of the mutate was successful
      math::MutateResult< assemble::Domain> operator()( const assemble::Domain &SHEET) const;

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

    }; // class MutateSheetFitToTemplate

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_SHEET_FIT_TO_TEMPLATE_H_ 
