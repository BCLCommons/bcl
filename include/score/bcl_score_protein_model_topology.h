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

#ifndef BCL_SCORE_PROTEIN_MODEL_TOPOLOGY_H_
#define BCL_SCORE_PROTEIN_MODEL_TOPOLOGY_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinModelTopology
    //! @brief return the percentage of sse contacts recovered from one given model to the other one
    //! @details This class implements the sse contact recovery quality measure for comparing protein models. 
    //! It measures the percentage of sse contacts from the first given protein model (ideally a native/template model) 
    //! that were recovered in the second given model.
    //!
    //! @author heinzes1
    //! @date Sep 26, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinModelTopology :
      public math::BinaryFunctionInterfaceSerializable< assemble::ProteinModel, assemble::ProteinModel, double>
    {
    public:
      
      //! f=f-score i.e. harmonic mean of tpr and ppv, tpr=true positive rate, ppv=positive predictive value
      enum score_result { score_f, score_tpr, score_ppv};
      
    private:

    //////////
    // data //
    //////////
      
      score_result m_Result; //!< which score result to return

    public:
      
      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor allowing to define the score result type, and default constructor
      ProteinModelTopology( score_result RESULT = score_f);

      //! @brief Clone function
      //! @return pointer to new ProteinModelTopology
      ProteinModelTopology *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief calculate and return the percentage of recovered contacts
      //! @param TEMPLATE the known correct protein model
      //! @param MODEL the protein model to be evaluated against the TEMPLATE
      //! @return the percentage of recovered contacts
      double operator()( const assemble::ProteinModel &TEMPLATE, const assemble::ProteinModel &MODEL) const;

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
      
    /////////////////////
    // helper funtions //
    /////////////////////
      
    private:
      
      typedef storage::Map< util::SiPtr< const assemble::SSEGeometryInterface>, util::SiPtr< const assemble::SSEGeometryInterface> > SSEMap;
      
      //! @brief calculates a mapping template_sse->model_sse based on the Q3 score
      //! @param TEMPLATE_SSES all template SSEs (copy, no const ref, is modified internally)
      //! @param MODEL_SSES all model SSEs (copy, no const ref, is modified internally)
      //! @return the mapping
      static SSEMap CalculateSSEMapping
      (
        util::SiPtrVector< const assemble::SSEGeometryInterface> TEMPLATE_SSES,
        util::SiPtrVector< const assemble::SSEGeometryInterface> MODEL_SSES
      );       
    }; // class ProteinModelTopology

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_PROTEIN_MODEL_TOPOLOGY_H_
