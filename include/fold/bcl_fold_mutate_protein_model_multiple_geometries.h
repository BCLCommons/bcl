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

#ifndef BCL_FOLD_MUTATE_PROTEIN_MODEL_MULTIPLE_GEOMETRIES_H_
#define BCL_FOLD_MUTATE_PROTEIN_MODEL_MULTIPLE_GEOMETRIES_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_fold_template.h"
#include "math/bcl_math_mutate_result.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MutateProteinModelMultipleGeometries
    //! @brief adds SSEs using fragment geometries read in from a fold template.
    //! @details This allows for fine placement of SSEs into a protein model during a folding run.
    //!
    //! @see @link example_fold_mutate_protein_model_multiple_geometries.cpp @endlink
    //! @author weinerbe
    //! @date Apr 26, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MutateProteinModelMultipleGeometries :
      public math::MutateInterface< assemble::ProteinModel>
    {

    private:

    //////////
    // data //
    //////////

      //! Fold template to use instead of of handler
      assemble::FoldTemplate m_FoldTemplate;

      //! probability of fitting into a small, equal, or large template, respectively
      storage::VectorND< 3, double> m_TemplateSizeProbabilities;

      //! scheme
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

      //! @brief default constructor
      MutateProteinModelMultipleGeometries();

      //! @brief constructor from a fold template and scheme
      //! @param FOLD_TEMPLATE fold template to be used
      //! @param SCHEME Scheme to be used
      MutateProteinModelMultipleGeometries
      (
        const assemble::FoldTemplate &FOLD_TEMPLATE,
        const std::string &SCHEME = GetStaticClassName< MutateProteinModelMultipleGeometries>()
      );

      //! @brief constructor from a fold template handler and scheme
      //! @param PROBABILITIES probability of fitting into a small, equal, or large template, respectively
      //! @param SCHEME Scheme to be used
      MutateProteinModelMultipleGeometries
      (
        const storage::VectorND< 3, double> &PROBABILITIES, // = storage::VectorND< 3, double>( 0.0, 1.0, 0.0),
        const std::string &SCHEME = GetStaticClassName< MutateProteinModelMultipleGeometries>()
      );

      //! @brief Clone function
      //! @return pointer to new MutateProteinModelMultipleGeometries
      MutateProteinModelMultipleGeometries *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this mutate
      //! @return the scheme for this mutate
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator taking a protein model and returning a mutate object of protein model type
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

    }; // class MutateProteinModelMultipleGeometries

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_MUTATE_PROTEIN_MODEL_MULTIPLE_GEOMETRIES_H_ 
