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

#ifndef BCL_SCORE_PROTEIN_MODEL_FRAGMENT_TOPOLOGY_H_
#define BCL_SCORE_PROTEIN_MODEL_FRAGMENT_TOPOLOGY_H_

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
    //! @class ProteinModelFragmentTopology
    //! @brief return the fragment-based topoloy score between two protein models
    //! @details This class implements the fragment-based topology metric for comparing protein models. 
    //! It measures the similarity in SSE fragment contacts between 
    //! that were recovered in the second given model.
    //!
    //! @author fooksams, mendenjl
    //! @date Oct 10, 2016
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinModelFragmentTopology :
      public math::BinaryFunctionInterfaceSerializable< assemble::ProteinModel, assemble::ProteinModel, double>
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

      //! @brief Clone function
      //! @return pointer to new ProteinModelFragmentTopology
      ProteinModelFragmentTopology *Clone() const;
      
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
    };      
  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_PROTEIN_MODEL_FRAGMENT_TOPOLOGY_H_
