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

#ifndef BCL_ASSEMBLE_COLLECTOR_PROTEIN_MODEL_CONFORMATION_BY_SCORE_H_
#define BCL_ASSEMBLE_COLLECTOR_PROTEIN_MODEL_CONFORMATION_BY_SCORE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "find/bcl_find_collector_interface.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "util/bcl_util_binary_function_interface_serializable.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CollectorProteinModelConformationByScore
    //! @brief TODO: add a brief comment
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_assemble_collector_protein_model_conformation_by_score.cpp @endlink
    //! @author alexanns
    //! @date Apr 15, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CollectorProteinModelConformationByScore :
      public find::CollectorInterface< util::SiPtrList< const ProteinModel>, ProteinModel>
    {

    private:

    //////////
    // data //
    //////////

      util::Implementation< math::FunctionInterfaceSerializable< ProteinModel, double> > m_ScoreFunction;

      //! number of conformations to collect
      size_t m_NumberToCollect;

      //! the method for sorting scores to determine whether the best or worst scored models are collected
      //! use less than to get the best by score, use greater than to get the worst by score
      util::Implementation< util::BinaryFunctionInterfaceSerializable< double, double, bool> > m_Sort;

      //! if true the current conformation held by the protein model will be considered with the other conformations
      bool m_ConsiderCurrentConformation;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      CollectorProteinModelConformationByScore();

      //! @brief constructor taking parameters
      //! @param SCORE the scoring function to be used
      //! @param NUM_TO_COLLECT the number of conformations to collect
      //! @param SORT the method for sorting scores to determine whether the best or worst scored models are collected
      //! @param CONSIDER_CURRENT if true the current conformation will be considered with the other conformations
      CollectorProteinModelConformationByScore
      (
        const util::ShPtr< math::FunctionInterfaceSerializable< ProteinModel, double> > &SCORE,
        const size_t NUM_TO_COLLECT,
        const util::ShPtr< util::BinaryFunctionInterfaceSerializable< double, double, bool> > &SORT,
        const bool CONSIDER_CURRENT
      );

      //! @brief constructor taking parameters
      //! @param NUM_TO_COLLECT the number of conformations to collect
      //! @param SORT the method for sorting scores to determine whether the best or worst scored models are collected
      //! @param CONSIDER_CURRENT if true the current conformation will be considered with the other conformations
      CollectorProteinModelConformationByScore
      (
        const size_t NUM_TO_COLLECT,
        const util::ShPtr< util::BinaryFunctionInterfaceSerializable< double, double, bool> > &SORT,
        const bool CONSIDER_CURRENT
      );

      //! @brief Clone function
      //! @return pointer to new CollectorProteinModelConformationByScore
      CollectorProteinModelConformationByScore *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! Collect the t_ReturnType objects in t_ArgumentType
      //! @param PROTEIN_MODEL entity that contains a t_ReturnType
      //! @return returns Group of the collected t_ReturnType objects
      virtual util::SiPtrList< const ProteinModel> Collect( const ProteinModel &PROTEIN_MODEL) const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

    private:

    }; // class CollectorProteinModelConformationByScore

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_COLLECTOR_PROTEIN_MODEL_CONFORMATION_BY_SCORE_H_ 
